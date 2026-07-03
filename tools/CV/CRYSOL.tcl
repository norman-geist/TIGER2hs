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
#   ::exp_dat   Experimental SAXS .dat file
#   ::crysol_bin     Path to CRYSOL executable, or "crysol" if in PATH
#   ::tmpdir    Temporary working directory
#   ::ref_pdb   Reference PDB whose ATOM/HETATM records provide atom/residue/element metadata
#
# Example:
#   set ::exp_dat "saxa.dat"
#   set ::crysol_bin   "/griffinrt/bin/ATSAS/ATSAS-3.2.1-1/bin/crysol"
#   set ::tmpdir  $env(TMPDIR)
#   set ::ref_pdb   "/path/to/reference_annotated.pdb"
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Standalone fallback for VMD
# -----------------------------------------------------------------------------
if {[llength [info commands myReplica]] == 0} {
    proc myReplica {} {
        return 0
    }
    
    set ::tmpdir  $env(TMPDIR)
    set ::exp_dat "/media/urzstore/geistn/vWF/ng_vwf154/1_ADAMTS13_SmallerPeak_trunc.dat"
    set ::crysol_bin   "/griffinrt/bin/ATSAS/ATSAS-3.2.1-1/bin/crysol"
    set ::ref_pdb [file join $::tmpdir "saxs_ref_[myReplica].pdb"]

    # Convert a VMD atom selection to a flat coordinate list:
    # {x1 y1 z1 x2 y2 z2 ...}
    proc make_coords {seltext} {
        if {[llength [info commands atomselect]] == 0} {
            error "CRYSOL: make_coords requires VMD atomselect"
        }
        set sel [atomselect top $seltext]
        $sel writepdb $::ref_pdb

        set coords {}
        foreach xyz [$sel get {x y z}] {
            lappend coords {*}$xyz
        }
        return $coords
    }

    set coords [make_coords "protein and noh"]
}

namespace eval CRYSOL {
    proc _get_params {} {
        set exp_dat ""
        set crysol_bin "crysol"
        set tmpdir "/tmp"

        if {[info exists ::exp_dat]} {
            set exp_dat $::exp_dat
        }
        if {[info exists ::crysol_bin]} {
            set crysol_bin $::crysol_bin
        }
        if {[info exists ::tmpdir]} {
            set tmpdir $::tmpdir
        }

        if {$exp_dat eq ""} {
            error "CRYSOL: ::exp_dat is not set"
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
            error "CRYSOL: ::tmpdir is empty"
        }
        if {![file isdirectory $tmpdir]} {
            error "CRYSOL: temporary directory does not exist: $tmpdir"
        }

        return [list $exp_dat $crysol_bin $tmpdir]
    }

    # Return ATOM/HETATM records from a reference PDB. These records define the
    # atom names, residue names, chain IDs, residue numbers, occupancies, B-factors,
    # and element columns that CRYSOL needs for chemically meaningful scattering.
    proc read_reference_atom_records {ref_pdb} {
        if {$ref_pdb eq ""} {
            error "CRYSOL: ::ref_pdb is empty"
        }
        if {![file exists $ref_pdb]} {
            error "CRYSOL: reference PDB does not exist: $ref_pdb"
        }

        set fh [open $ref_pdb r]
        set records {}
        while {[gets $fh line] >= 0} {
            set recname [string range $line 0 5]
            if {[string match "ATOM  " $recname] || [string match "HETATM" $recname]} {
                # Pad short lines so fixed-width ranges below are always safe.
                set line [format "%-80s" $line]
                lappend records $line
            }
        }
        close $fh

        if {[llength $records] == 0} {
            error "CRYSOL: no ATOM/HETATM records found in reference PDB: $ref_pdb"
        }
        return $records
    }

    # Write a chemically annotated PDB by taking atom/residue/element metadata from
    # a reference PDB and replacing only the x/y/z coordinate fields.
    # Coordinate order MUST match the ATOM/HETATM order in the reference PDB.
    proc write_reference_pdb {coords ref_pdb pdbfile} {
        if {[expr {[llength $coords] % 3}] != 0} {
            error "CRYSOL: coordinate list length is not divisible by 3"
        }

        set records [read_reference_atom_records $ref_pdb]
        set natoms [expr {[llength $coords] / 3}]
        set nref [llength $records]

        if {$natoms != $nref} {
            error "CRYSOL: coordinate/reference atom count mismatch: coords=$natoms reference=$nref. The colvar coordinate order and selection must match ::ref_pdb ATOM/HETATM order."
        }

        set fh [open $pdbfile w]
        for {set i 0} {$i < $natoms} {incr i} {
            set x [lindex $coords [expr {3*$i    }]]
            set y [lindex $coords [expr {3*$i + 1}]]
            set z [lindex $coords [expr {3*$i + 2}]]

            set line [lindex $records $i]
            set left  [string range $line 0 29]
            set right [string range $line 54 79]
            puts $fh [format "%s%8.3f%8.3f%8.3f%s" $left $x $y $z $right]
        }

        puts $fh "TER"
        puts $fh "END"
        close $fh
    }

    # Fallback retained for testing only. Do not use this for atomistic protein SAXS.
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

        if {[info exists ::ref_pdb] && $::ref_pdb ne ""} {
            write_reference_pdb $coords $::ref_pdb $pdbfile
        } else {
            # Fallback for standalone testing only. For production CRYSOL scoring,
            # set ::ref_pdb so CRYSOL receives chemically meaningful atoms.
            write_dummy_pdb $coords $pdbfile
        }

        set olddir [pwd]
        cd $tmpdir

        set cmd [list $crysol_bin $pdbfile $exp_dat]

        if {[catch {exec {*}$cmd --implicit-hydrogen=1 2>@1} out]} {
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