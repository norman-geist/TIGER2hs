# -----------------------------------------------------------------------------
# SAXS CRYSOL score for NAMD scriptedColvar / VMD Tcl
# -----------------------------------------------------------------------------
#
# Required public entry point:
#
#   calc_saxs_score coords
#
# Coordinate format:
#   coords = flat list {x1 y1 z1 x2 y2 z2 ...}
#
# Global settings:
#   ::saxs_exp_dat   Experimental SAXS .dat file
#   ::crysol_bin     Path to CRYSOL executable, or "crysol" if in PATH
#   ::saxs_tmpdir    Temporary working directory
#
# Example:
#   set ::saxs_exp_dat "/media/urzstore/geistn/Beta2/SAXS/B2GP1_SAXS_digi.dat"
#   set ::crysol_bin   "/griffinrt/bin/ATSAS/ATSAS-3.2.1-1/bin/crysol"
#   set ::saxs_tmpdir  $env(TMPDIR)
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Standalone fallback for VMD
# -----------------------------------------------------------------------------
if {[llength [info commands myReplica]] == 0} {
    proc myReplica {} {
        return 0
    }
    
    # Convert a VMD atom selection to a flat coordinate list:
    # {x1 y1 z1 x2 y2 z2 ...}
    proc make_coords {seltext} {
        set sel [atomselect top $seltext]
    
        set coords {}
        foreach xyz [$sel get {x y z}] {
            lappend coords {*}$xyz
        }
        return $coords
    }
    set coords [make_coords all]
    
    set ::saxs_exp_dat "/media/urzstore/geistn/Beta2/SAXS/B2GP1_SAXS_digi.dat"
    set ::crysol_bin   "/griffinrt/bin/ATSAS/ATSAS-3.2.1-1/bin/crysol"
    set ::saxs_tmpdir  $env(TMPDIR)
}

namespace eval CRYSOL {
    proc _get_params {} {
        set exp_dat ""
        set crysol_bin "crysol"
        set tmpdir "/tmp"

        if {[info exists ::saxs_exp_dat]} {
            set exp_dat $::saxs_exp_dat
        }
        if {[info exists ::crysol_bin]} {
            set crysol_bin $::crysol_bin
        }
        if {[info exists ::saxs_tmpdir]} {
            set tmpdir $::saxs_tmpdir
        }

        if {$exp_dat eq ""} {
            error "CRYSOL: ::saxs_exp_dat is not set"
        }
        if {![file exists $exp_dat]} {
            error "CRYSOL: experimental SAXS file does not exist: $exp_dat"
        }

        if {$crysol_bin eq ""} {
            error "CRYSOL: ::crysol_bin is not set"
        }
        if {![file exists $crysol_bin] && [auto_execok $crysol_bin] eq ""} {
            error "CRYSOL: CRYSOL executable not found: $crysol_bin"
        }

        if {$tmpdir eq ""} {
            error "CRYSOL: ::saxs_tmpdir is empty"
        }
        if {![file isdirectory $tmpdir]} {
            error "CRYSOL: temporary directory does not exist: $tmpdir"
        }

        return [list $exp_dat $crysol_bin $tmpdir]
    }

    # Write dummy-atom PDB from flat coordinate list:
    # x1 y1 z1 x2 y2 z2 ...
    proc write_dummy_pdb {coords pdbfile} {
        if {[expr {[llength $coords] % 3}] != 0} {
            error "CRYSOL: coordinate list length is not divisible by 3"
        }

        set fh [open $pdbfile w]

        set natoms [expr {[llength $coords] / 3}]
        for {set i 0} {$i < $natoms} {incr i} {
            set x [lindex $coords [expr {3*$i    }]]
            set y [lindex $coords [expr {3*$i + 1}]]
            set z [lindex $coords [expr {3*$i + 2}]]

            puts $fh [format "ATOM  %5d  CA  ALA A%4d    %8.3f%8.3f%8.3f  1.00  0.00           C" \
                      [expr {$i+1}] [expr {$i+1}] $x $y $z]
        }

        puts $fh "TER"
        puts $fh "END"
        close $fh
    }

    proc parse_crysol_score {text} {
        foreach line [split $text "\n"] {
            if {[regexp -nocase {Chi-square of fit[ .]*:\s*([0-9.+-Ee]+)} $line -> val]} {
                return $val
            }
        }
        error "CRYSOL: could not parse CRYSOL Chi-square of fit"
    }

    proc _calc {coords} {
        foreach {exp_dat crysol_bin tmpdir} [_get_params] break

        set pdbfile [file join $tmpdir "saxs_[myReplica].pdb"]

        write_dummy_pdb $coords $pdbfile

        set olddir [pwd]
        cd $tmpdir

        set cmd [list $crysol_bin $pdbfile $exp_dat]

        if {[catch {exec {*}$cmd 2>@1} out]} {
            cd $olddir
            error "CRYSOL: CRYSOL failed: $out"
        }

        cd $olddir

        set chi [parse_crysol_score $out]

        return $chi
    }
}

# -----------------------------------------------------------------------------
# NAMD scriptedColvar entry point: must accept ONLY coordinates.
# -----------------------------------------------------------------------------
proc calc_CRYSOL {coords} {
    return [CRYSOL::_calc $coords]
}