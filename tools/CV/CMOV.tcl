#Contact Map Overlap Vector

#Default variables
if {![info exists CMOV_rc]} {
    set CMOV_rc 6; #decay length of Gaussian contacts
}
if {![info exists CMOV_min_sep]} {
    set CMOV_min_sep 3; #skip close sequence neighbors
}
if {![info exists CMOV_dmax]} {
    set CMOV_dmax 12; #max contact distance
}

proc calc_CMOV {args} {
    global CMOV_rc CMOV_min_sep CMOV_dmax
    
    set coordsA [lindex $args 0]
    set coordsB [lindex $args 1]
    
    set rc      $CMOV_rc
    set min_sep $CMOV_min_sep
    set dmax    $CMOV_dmax
    
    set n_atoms [expr {[llength $coordsA] / 3}]
    if {[llength $coordsB] != [llength $coordsA]} {
        error "Coordinate lists must be the same length"
    }

    set sum_diff 0.0
    set n_pairs 0

    for {set i 0} {$i < $n_atoms} {incr i} {
        for {set j [expr {$i + $min_sep}]} {$j < $n_atoms} {incr j} {
            set dxA [expr {[lindex $coordsA [expr {$i*3}]] - [lindex $coordsA [expr {$j*3}]]}]
            set dyA [expr {[lindex $coordsA [expr {$i*3+1}]] - [lindex $coordsA [expr {$j*3+1}]]}]
            set dzA [expr {[lindex $coordsA [expr {$i*3+2}]] - [lindex $coordsA [expr {$j*3+2}]]}]
            set distA [expr {sqrt($dxA*$dxA + $dyA*$dyA + $dzA*$dzA)}]

            set dxB [expr {[lindex $coordsB [expr {$i*3}]] - [lindex $coordsB [expr {$j*3}]]}]
            set dyB [expr {[lindex $coordsB [expr {$i*3+1}]] - [lindex $coordsB [expr {$j*3+1}]]}]
            set dzB [expr {[lindex $coordsB [expr {$i*3+2}]] - [lindex $coordsB [expr {$j*3+2}]]}]
            set distB [expr {sqrt($dxB*$dxB + $dyB*$dyB + $dzB*$dzB)}]

            # Optional sparsity cutoff
            if {$distA > $dmax && $distB > $dmax} {
                continue
            }

            # Gaussian contact value
            set cA [expr {exp(-($distA*$distA)/($rc*$rc))}]
            set cB [expr {exp(-($distB*$distB)/($rc*$rc))}]
            set diff [expr {abs($cA - $cB)}]

            set sum_diff [expr {$sum_diff + $diff}]
            incr n_pairs
        }
    }

    if {$n_pairs == 0} {
        return 0.0
    }

    set cmov [expr {1.0 - $sum_diff / double($n_pairs)}]
    return $cmov
}
