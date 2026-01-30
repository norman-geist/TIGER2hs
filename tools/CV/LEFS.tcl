# -----------------------------------------------------------------------------
# LEFS - Local Environment Fingerprint Similarity
#
# LEFS.tcl  (Tcl 8.5)  — NAMD scriptedColvar-compatible
#
# REQUIREMENT:
#   calc_LEFS must accept ONLY (coordsA coordsB).
#   All parameters are taken from GLOBAL variables (with sane defaults).
#
# COORDINATE FORMAT:
#   coordsA/coordsB are flat lists: {x1 y1 z1 x2 y2 z2 ... xN yN zN}
#
# OUTPUT:
#   Similarity in [0,1] where 1.0 = identical local environments.
#
# GLOBAL PARAMS (optional, defaults below):
#   ::LEFS_dmax        (Å)    default 12.0
#   ::LEFS_K           (int)  default 6
#   ::LEFS_rbfWidth    (Å)    default 1.0
#   ::LEFS_min_sep     (int)  default 3
#
# Notes:
# - Segment-free, alignment-free.
# - Uses a cell-list (grid) neighbor search for near-linear scaling.
# -----------------------------------------------------------------------------

namespace eval LEFS {
    proc _get_params {} {
        # Defaults
        set dmax     12.0
        set K        6
        set rbfWidth 1.0
        set min_sep  3

        # Override from globals if present
        if {[info exists ::LEFS_dmax]}     { set dmax     $::LEFS_dmax }
        if {[info exists ::LEFS_K]}        { set K        $::LEFS_K }
        if {[info exists ::LEFS_rbfWidth]} { set rbfWidth $::LEFS_rbfWidth }
        if {[info exists ::LEFS_min_sep]}  { set min_sep  $::LEFS_min_sep }

        # Basic validation (keep strict, fail fast)
        if {$dmax <= 0.0}   { error "LEFS: ::LEFS_dmax must be > 0" }
        if {$K < 1}         { error "LEFS: ::LEFS_K must be >= 1" }
        if {$rbfWidth <= 0.0} { error "LEFS: ::LEFS_rbfWidth must be > 0" }
        if {$min_sep < 0}   { error "LEFS: ::LEFS_min_sep must be >= 0" }

        return [list $dmax $K $rbfWidth $min_sep]
    }

    # Build a simple cell list: maps integer 3D cell index -> list of atom indices.
    proc _build_cell_list {coords cellSize} {
        set n [expr {[llength $coords] / 3}]
        array set cells {}

        for {set i 0} {$i < $n} {incr i} {
            set x [lindex $coords [expr {3*$i}]]
            set y [lindex $coords [expr {3*$i + 1}]]
            set z [lindex $coords [expr {3*$i + 2}]]

            set ix [expr {int(floor($x / $cellSize))}]
            set iy [expr {int(floor($y / $cellSize))}]
            set iz [expr {int(floor($z / $cellSize))}]

            set key "${ix},${iy},${iz}"
            if {[info exists cells($key)]} {
                lappend cells($key) $i
            } else {
                set cells($key) [list $i]
            }
        }
        return [array get cells] ;# dict-like list usable with "array set"
    }

    # Compute K-bin radial fingerprint for a single index i.
    proc _fingerprint_i {coords i cellsDict cellSize dmax K rbfWidth min_sep} {
        # Bin centers in (0, dmax):
        # mu_b = (b+1)*dmax/(K+1)
        set centers {}
        for {set b 0} {$b < $K} {incr b} {
            lappend centers [expr {double($b+1) * $dmax / double($K+1)}]
        }

        # Initialize fingerprint bins
        set fp {}
        for {set b 0} {$b < $K} {incr b} { lappend fp 0.0 }

        # i coordinates
        set xi [lindex $coords [expr {3*$i}]]
        set yi [lindex $coords [expr {3*$i + 1}]]
        set zi [lindex $coords [expr {3*$i + 2}]]

        # i cell
        set ix [expr {int(floor($xi / $cellSize))}]
        set iy [expr {int(floor($yi / $cellSize))}]
        set iz [expr {int(floor($zi / $cellSize))}]

        # Expand dict into array for fast lookup
        array set cells $cellsDict

        set dmax2    [expr {$dmax * $dmax}]
        set inv_rbf  [expr {1.0 / double($rbfWidth)}]

        # Iterate 27 neighboring cells
        for {set dx -1} {$dx <= 1} {incr dx} {
            for {set dy -1} {$dy <= 1} {incr dy} {
                for {set dz -1} {$dz <= 1} {incr dz} {
                    set key "[expr {$ix+$dx}],[expr {$iy+$dy}],[expr {$iz+$dz}]"
                    if {![info exists cells($key)]} { continue }

                    foreach j $cells($key) {
                        if {$j == $i} { continue }
                        if {[expr {abs($j - $i)}] < $min_sep} { continue }

                        set xj [lindex $coords [expr {3*$j}]]
                        set yj [lindex $coords [expr {3*$j + 1}]]
                        set zj [lindex $coords [expr {3*$j + 2}]]

                        set dxv [expr {$xi - $xj}]
                        set dyv [expr {$yi - $yj}]
                        set dzv [expr {$zi - $zj}]
                        set r2  [expr {$dxv*$dxv + $dyv*$dyv + $dzv*$dzv}]
                        if {$r2 > $dmax2} { continue }
                        set r [expr {sqrt($r2)}]

                        # Smooth cutoff weight: w = exp(-(r/dmax)^6)
                        set tcut [expr {$r / double($dmax)}]
                        set w [expr {exp(-pow($tcut, 6.0))}]

                        # Accumulate into radial basis bins
                        for {set b 0} {$b < $K} {incr b} {
                            set mu [lindex $centers $b]
                            set tt [expr {($r - $mu) * $inv_rbf}]
                            set add [expr {$w * exp(-$tt*$tt)}]
                            lset fp $b [expr {[lindex $fp $b] + $add}]
                        }
                    }
                }
            }
        }

        # Normalize to sum 1 (if nonzero)
        set s 0.0
        for {set b 0} {$b < $K} {incr b} {
            set s [expr {$s + [lindex $fp $b]}]
        }
        if {$s > 0.0} {
            set invs [expr {1.0 / $s}]
            for {set b 0} {$b < $K} {incr b} {
                lset fp $b [expr {[lindex $fp $b] * $invs}]
            }
        }
        return $fp
    }

    # Core computation (called by global wrapper proc calc_LEFS)
    proc _calc {coordsA coordsB} {
        if {[llength $coordsA] != [llength $coordsB]} {
            error "calc_LEFS: Coordinate lists must be the same length"
        }

        foreach {dmax K rbfWidth min_sep} [_get_params] break

        set n [expr {[llength $coordsA] / 3}]
        if {$n < 2} { return 1.0 }

        # Cell size: dmax
        set cellSize $dmax

        # Build cell lists
        set cellsA [_build_cell_list $coordsA $cellSize]
        set cellsB [_build_cell_list $coordsB $cellSize]

        # Per-index fingerprint similarity:
        # L1(fpA, fpB) in [0,2] -> sim_i = 1 - 0.5*L1 in [0,1]
        set simSum 0.0
        set count  0

        for {set i 0} {$i < $n} {incr i} {
            set fpA [_fingerprint_i $coordsA $i $cellsA $cellSize $dmax $K $rbfWidth $min_sep]
            set fpB [_fingerprint_i $coordsB $i $cellsB $cellSize $dmax $K $rbfWidth $min_sep]

            set l1 0.0
            for {set b 0} {$b < $K} {incr b} {
                set l1 [expr {$l1 + abs([lindex $fpA $b] - [lindex $fpB $b])}]
            }

            set sim_i [expr {1.0 - 0.5*$l1}]
            if {$sim_i < 0.0} { set sim_i 0.0 }
            if {$sim_i > 1.0} { set sim_i 1.0 }

            set simSum [expr {$simSum + $sim_i}]
            incr count
        }

        return [expr {$simSum / double($count)}]
    }
}

# -----------------------------------------------------------------------------
# NAMD scriptedColvar entry point:
# Must accept ONLY coordinates.
# -----------------------------------------------------------------------------
proc calc_LEFS {coordsA coordsB} {
    return [LEFS::_calc $coordsA $coordsB]
}

# -----------------------------------------------------------------------------
# Example: set globals (do this in your NAMD Tcl context before calling)
# -----------------------------------------------------------------------------
# set ::LEFS_dmax 12.0
# set ::LEFS_K 6
# set ::LEFS_rbfWidth 1.0
# set ::LEFS_min_sep 3
#
# Then NAMD scriptedColvar calls:
#   set s [calc_LEFS $coordsA $coordsB]
# -----------------------------------------------------------------------------
