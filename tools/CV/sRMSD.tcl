#Default variables
if {![info exists sRMSD_window]} {
    set sRMSD_window 12
}

# Utility: subtract two vectors
proc vec_sub {a b} {
    list [expr {[lindex $a 0] - [lindex $b 0]}] \
         [expr {[lindex $a 1] - [lindex $b 1]}] \
         [expr {[lindex $a 2] - [lindex $b 2]}]
}

# Compute center of geometry
proc center_of_mass {coords} {
    set n [llength $coords]
    set cx 0.0
    set cy 0.0
    set cz 0.0
    foreach p $coords {
        foreach {x y z} $p break
        set cx [expr {$cx + $x}]
        set cy [expr {$cy + $y}]
        set cz [expr {$cz + $z}]
    }
    return [list [expr {$cx/$n}] [expr {$cy/$n}] [expr {$cz/$n}]]
}

# Subtract center from coords
proc translate_coords {coords center} {
    set newcoords {}
    foreach p $coords {
        set x [expr {[lindex $p 0] - [lindex $center 0]}]
        set y [expr {[lindex $p 1] - [lindex $center 1]}]
        set z [expr {[lindex $p 2] - [lindex $center 2]}]
        lappend newcoords [list $x $y $z]
    }
    return $newcoords
}

# Matrix multiplication: 3x3 * 3x3
proc mat_mult {A B} {
    set C {}
    for {set i 0} {$i < 3} {incr i} {
        set row {}
        for {set j 0} {$j < 3} {incr j} {
            set sum 0.0
            for {set k 0} {$k < 3} {incr k} {
                set sum [expr {$sum + [lindex $A $i $k] * [lindex $B $k $j]}]
            }
            lappend row $sum
        }
        lappend C $row
    }
    return $C
}

# Matrix transpose 3xN
proc transpose_matrix {M} {
    set result {}
    for {set j 0} {$j < 3} {incr j} {
        set row {}
        foreach r $M {
            lappend row [lindex $r $j]
        }
        lappend result $row
    }
    return $result
}

# Multiply matrix with vector
proc matvec_mult {M v} {
    set result {}
    foreach row $M {
        set s 0.0
        for {set i 0} {$i < 3} {incr i} {
            set s [expr {$s + [lindex $row $i] * [lindex $v $i]}]
        }
        lappend result $s
    }
    return $result
}

# Rotate coords by 3x3 matrix
proc rotate_coords {coords R} {
    set result {}
    foreach p $coords {
        lappend result [matvec_mult $R $p]
    }
    return $result
}

# Kabsch rotation
proc kabsch_rotation {P Q} {
    set H [list {0.0 0.0 0.0} {0.0 0.0 0.0} {0.0 0.0 0.0}]
    for {set i 0} {$i < [llength $P]} {incr i} {
        foreach {px py pz} [lindex $P $i] break
        foreach {qx qy qz} [lindex $Q $i] break
        lset H 0 0 [expr {[lindex $H 0 0] + $qx * $px}]
        lset H 0 1 [expr {[lindex $H 0 1] + $qx * $py}]
        lset H 0 2 [expr {[lindex $H 0 2] + $qx * $pz}]
        lset H 1 0 [expr {[lindex $H 1 0] + $qy * $px}]
        lset H 1 1 [expr {[lindex $H 1 1] + $qy * $py}]
        lset H 1 2 [expr {[lindex $H 1 2] + $qy * $pz}]
        lset H 2 0 [expr {[lindex $H 2 0] + $qz * $px}]
        lset H 2 1 [expr {[lindex $H 2 1] + $qz * $py}]
        lset H 2 2 [expr {[lindex $H 2 2] + $qz * $pz}]
    }

    # Approximate orthonormalization (normalize rows)
    set U {}
    for {set i 0} {$i < 3} {incr i} {
        set row [lindex $H $i]
        set norm 0.0
        foreach x $row { set norm [expr {$norm + $x*$x}] }
        set norm [expr {sqrt($norm)}]
        set norm_row {}
        foreach x $row { lappend norm_row [expr {$x / $norm}] }
        lappend U $norm_row
    }

    return [transpose_matrix $U]
}

# RMSD between aligned fragments
proc rmsd {A B} {
    set sum 0.0
    set N [llength $A]
    for {set i 0} {$i < $N} {incr i} {
        foreach {x1 y1 z1} [lindex $A $i] break
        foreach {x2 y2 z2} [lindex $B $i] break
        set dx [expr {$x1 - $x2}]
        set dy [expr {$y1 - $y2}]
        set dz [expr {$z1 - $z2}]
        set sum [expr {$sum + $dx*$dx + $dy*$dy + $dz*$dz}]
    }
    return [expr {sqrt($sum / $N)}]
}

# Main function: calc_sRMSD
proc calc_sRMSD {args} {

    global sRMSD_window
    
    set window $sRMSD_window
    set coordsA [lindex $args 0]
    set coordsB [lindex $args 1]

    if {[llength $coordsA] != [llength $coordsB]} {
        error "Coordinate list lengths must match"
    }
    set natoms [expr {[llength $coordsA] / 3}]
    set vecsA {}
    set vecsB {}
    for {set i 0} {$i < $natoms} {incr i} {
        lappend vecsA [lrange $coordsA [expr {$i*3}] [expr {$i*3+2}]]
        lappend vecsB [lrange $coordsB [expr {$i*3}] [expr {$i*3+2}]]
    }

    set total 0.0
    set count 0
    for {set i 0} {$i <= ($natoms - $window)} {incr i} {
        set fragA [lrange $vecsA $i [expr {$i+$window-1}]]
        set fragB [lrange $vecsB $i [expr {$i+$window-1}]]

        set comA [center_of_mass $fragA]
        set comB [center_of_mass $fragB]

        set fragA0 [translate_coords $fragA $comA]
        set fragB0 [translate_coords $fragB $comB]

        set R [kabsch_rotation $fragA0 $fragB0]
        set fragB_aligned [rotate_coords $fragB0 $R]

        set val [rmsd $fragA0 $fragB_aligned]
        set total [expr {$total + $val}]
        incr count
    }

    return [expr {$total / $count}]
}
