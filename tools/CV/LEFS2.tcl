# -----------------------------------------------------------------------------
# LEFS2.tcl  (Tcl 8.5) — NAMD scriptedColvar-compatible
#
# LEFS2 = "identity-aware" Local Environment Fingerprint Similarity
#
#   LEFS2: per-residue fingerprint is a (bin,neighbor-index) weighted map:
#          it encodes not just "how much stuff at distance r" but also "WHO"
#          (atom index j) contributes at that distance around i.
#
# This makes LEFS2 more topology/identity specific (harder for swaps/permutations
# to score well), at the cost of more work per evaluation.
#
# REQUIREMENT:
#   calc_LEFS2 must accept ONLY (coordsA coordsB).
#   All parameters are GLOBAL variables (defaults provided).
#
# COORDINATE FORMAT:
#   coordsA/coordsB are flat lists: {x1 y1 z1 x2 y2 z2 ...}
#
# OUTPUT:
#   Similarity in [0,1] where 1.0 = identical local (bin,neighbor) fingerprints.
#
# GLOBAL PARAMS (optional):
#   ::LEFS2_dmax        (Å)    default 12.0
#   ::LEFS2_K           (int)  default 6
#   ::LEFS2_rbfWidth    (Å)    default 1.0
#   ::LEFS2_min_sep     (int)  default 3
#
# Notes:
# - Uses a cell-list neighbor search (grid) for near-linear behavior in practice.
# - LEFS2 is heavier than LEFS because fingerprints contain neighbor identities.
# -----------------------------------------------------------------------------

namespace eval LEFS2 {
    proc _get_params {} {
        set dmax     12.0
        set K        6
        set rbfWidth 1.0
        set min_sep  3

        if {[info exists ::LEFS2_dmax]}     { set dmax     $::LEFS2_dmax }
        if {[info exists ::LEFS2_K]}        { set K        $::LEFS2_K }
        if {[info exists ::LEFS2_rbfWidth]} { set rbfWidth $::LEFS2_rbfWidth }
        if {[info exists ::LEFS2_min_sep]}  { set min_sep  $::LEFS2_min_sep }

        if {$dmax <= 0.0}     { error "LEFS2: ::LEFS2_dmax must be > 0" }
        if {$K < 1}           { error "LEFS2: ::LEFS2_K must be >= 1" }
        if {$rbfWidth <= 0.0} { error "LEFS2: ::LEFS2_rbfWidth must be > 0" }
        if {$min_sep < 0}     { error "LEFS2: ::LEFS2_min_sep must be >= 0" }

        return [list $dmax $K $rbfWidth $min_sep]
    }

    # Cell list: key "ix,iy,iz" -> list of indices
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
        return [array get cells]
    }

    # Identity-aware fingerprint for residue/atom i:
    # Returns a dict-like list {key value key value ...} where key = "b:j"
    # b = radial bin index, j = neighbor atom index
    #
    # Weight per neighbor:
    #   weight = w_cut(r) * rbf_b(r)
    # with smooth cutoff w_cut(r)=exp(-(r/dmax)^6) and Gaussian RBF shells.
    #
    # Then normalize the whole map so sum(weights)=1, which makes L1 bounded.
    proc _fingerprint_i {coords i cellsDict cellSize dmax K rbfWidth min_sep} {
        # Bin centers in (0, dmax): mu_b = (b+1)*dmax/(K+1)
        set centers {}
        for {set b 0} {$b < $K} {incr b} {
            lappend centers [expr {double($b+1) * $dmax / double($K+1)}]
        }

        # i coords
        set xi [lindex $coords [expr {3*$i}]]
        set yi [lindex $coords [expr {3*$i + 1}]]
        set zi [lindex $coords [expr {3*$i + 2}]]

        # i cell
        set ix [expr {int(floor($xi / $cellSize))}]
        set iy [expr {int(floor($yi / $cellSize))}]
        set iz [expr {int(floor($zi / $cellSize))}]

        array set cells $cellsDict

        set dmax2   [expr {$dmax * $dmax}]
        set inv_rbf [expr {1.0 / double($rbfWidth)}]

        # Accumulator map
        array set acc {}
        set sumw 0.0

        # Iterate 27 neighboring cells
        for {set dx -1} {$dx <= 1} {incr dx} {
            for {set dy -1} {$dy <= 1} {incr dy} {
                for {set dz -1} {$dz <= 1} {incr dz} {
                    set keycell "[expr {$ix+$dx}],[expr {$iy+$dy}],[expr {$iz+$dz}]"
                    if {![info exists cells($keycell)]} { continue }

                    foreach j $cells($keycell) {
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

                        # Smooth cutoff weight
                        set tcut [expr {$r / double($dmax)}]
                        set wcut [expr {exp(-pow($tcut, 6.0))}]

                        # Contribute into each bin with identity j
                        for {set b 0} {$b < $K} {incr b} {
                            set mu [lindex $centers $b]
                            set tt [expr {($r - $mu) * $inv_rbf}]
                            set wbin [expr {exp(-$tt*$tt)}]
                            set w [expr {$wcut * $wbin}]
                            if {$w <= 0.0} { continue }

                            set k "b${b}:j${j}"
                            if {[info exists acc($k)]} {
                                set acc($k) [expr {$acc($k) + $w}]
                            } else {
                                set acc($k) $w
                            }
                            set sumw [expr {$sumw + $w}]
                        }
                    }
                }
            }
        }

        # Normalize so sum(weights)=1 (if any)
        if {$sumw > 0.0} {
            set invs [expr {1.0 / $sumw}]
            foreach k [array names acc] {
                set acc($k) [expr {$acc($k) * $invs}]
            }
        }

        return [array get acc]
    }

    # L1 distance between two normalized maps {k v ...} and {k v ...}
    # Returns value in [0,2] if both maps sum to 1.
    proc _l1_maps {mapA mapB} {
        array set A $mapA
        array set B $mapB

        set l1 0.0
        array set seen {}

        foreach k [array names A] {
            set a $A($k)
            if {[info exists B($k)]} {
                set b $B($k)
                set l1 [expr {$l1 + abs($a - $b)}]
            } else {
                set l1 [expr {$l1 + abs($a)}]
            }
            set seen($k) 1
        }
        foreach k [array names B] {
            if {[info exists seen($k)]} { continue }
            set l1 [expr {$l1 + abs($B($k))}]
        }
        return $l1
    }

    proc _calc {coordsA coordsB} {
        if {[llength $coordsA] != [llength $coordsB]} {
            error "calc_LEFS2: Coordinate lists must be the same length"
        }

        foreach {dmax K rbfWidth min_sep} [_get_params] break

        set n [expr {[llength $coordsA] / 3}]
        if {$n < 2} { return 1.0 }

        # Cell size: dmax
        set cellSize $dmax
        set cellsA [_build_cell_list $coordsA $cellSize]
        set cellsB [_build_cell_list $coordsB $cellSize]

        # Per-index similarity:
        #   L1 in [0,2] -> sim_i = 1 - 0.5*L1 in [0,1]
        set simSum 0.0
        set count  0

        for {set i 0} {$i < $n} {incr i} {
            set fA [_fingerprint_i $coordsA $i $cellsA $cellSize $dmax $K $rbfWidth $min_sep]
            set fB [_fingerprint_i $coordsB $i $cellsB $cellSize $dmax $K $rbfWidth $min_sep]

            set l1 [_l1_maps $fA $fB]
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
# NAMD scriptedColvar entry point: must accept ONLY coordinates.
# -----------------------------------------------------------------------------
proc calc_LEFS2 {coordsA coordsB} {
    return [LEFS2::_calc $coordsA $coordsB]
}

# -----------------------------------------------------------------------------
# Example globals (set in your NAMD Tcl context):
# -----------------------------------------------------------------------------
# set ::LEFS2_dmax     12.0
# set ::LEFS2_K        6
# set ::LEFS2_rbfWidth 1.0
# set ::LEFS2_min_sep  3
# -----------------------------------------------------------------------------
