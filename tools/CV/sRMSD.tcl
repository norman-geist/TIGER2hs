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

# Compute dot product of two 3D vectors
proc dot {a b} {
    expr {[lindex $a 0]*[lindex $b 0] + [lindex $a 1]*[lindex $b 1] + [lindex $a 2]*[lindex $b 2]}
}

# Compute cross product of two 3D vectors
proc cross {a b} {
    list \
        [expr {[lindex $a 1]*[lindex $b 2] - [lindex $a 2]*[lindex $b 1]}] \
        [expr {[lindex $a 2]*[lindex $b 0] - [lindex $a 0]*[lindex $b 2]}] \
        [expr {[lindex $a 0]*[lindex $b 1] - [lindex $a 1]*[lindex $b 0]}]
}

# Normalize a 3D vector
proc normalize {v} {
    set norm [expr {sqrt([dot $v $v])}]
    if {$norm == 0} {
        return {0 0 0}
    }
    list \
        [expr {[lindex $v 0] / $norm}] \
        [expr {[lindex $v 1] / $norm}] \
        [expr {[lindex $v 2] / $norm}]
}

# Build 3x3 zero matrix
proc zero_matrix {} {
    return {{0.0 0.0 0.0} {0.0 0.0 0.0} {0.0 0.0 0.0}}
}

# Add outer product of vectors a^T * b to matrix M
proc accumulate_covariance {M a b} {
    for {set i 0} {$i < 3} {incr i} {
        for {set j 0} {$j < 3} {incr j} {
            lset M $i $j [expr {[lindex $M $i $j] + [lindex $a $i] * [lindex $b $j]}]
        }
    }
    return $M
}

# Transpose a 3x3 matrix
proc transpose3x3 {M} {
    list \
        [list [lindex $M 0 0] [lindex $M 1 0] [lindex $M 2 0]] \
        [list [lindex $M 0 1] [lindex $M 1 1] [lindex $M 2 1]] \
        [list [lindex $M 0 2] [lindex $M 1 2] [lindex $M 2 2]]
}

# Determinant of 3x3 matrix via scalar triple product (more stable)
proc determinant3x3 {M} {
    set u [list [lindex $M 0 0] [lindex $M 1 0] [lindex $M 2 0]]
    set v [list [lindex $M 0 1] [lindex $M 1 1] [lindex $M 2 1]]
    set w [list [lindex $M 0 2] [lindex $M 1 2] [lindex $M 2 2]]
    return [dot $u [cross $v $w]]
}

# Approximate orthonormalization (polar decomposition using Gram-Schmidt)
proc orthonormalize {M} {
    # Columns of M
    set u0 [list [lindex $M 0 0] [lindex $M 1 0] [lindex $M 2 0]]
    set u1 [list [lindex $M 0 1] [lindex $M 1 1] [lindex $M 2 1]]
    set u2 [list [lindex $M 0 2] [lindex $M 1 2] [lindex $M 2 2]]

    set e0 [normalize $u0]

    set proj01 [expr {[dot $u1 $e0]}]
    set temp1 [list \
        [expr {[lindex $u1 0] - $proj01 * [lindex $e0 0]}] \
        [expr {[lindex $u1 1] - $proj01 * [lindex $e0 1]}] \
        [expr {[lindex $u1 2] - $proj01 * [lindex $e0 2]}]]
    set e1 [normalize $temp1]

    set e2 [normalize [cross $e0 $e1]]

    # Construct matrix with orthonormalized columns
    return [transpose3x3 [list $e0 $e1 $e2]]
}

# Main Kabsch rotation function (no reflections)
proc kabsch_rotation {P Q} {
    if {[llength $P] != [llength $Q]} {
        error "Point sets must be same length"
    }

    set H [zero_matrix]
    for {set i 0} {$i < [llength $P]} {incr i} {
        set p [lindex $P $i]
        set q [lindex $Q $i]
        set H [accumulate_covariance $H $q $p]  ;# H = sum(q_i * p_i^T)
    }
    
    set R [orthonormalize $H]

    # Reject mirror solution if det(R) < 0
    #puts "determinant: [determinant3x3 $R]"
    if {[determinant3x3 $R] < 0} {
        # Flip last column
        for {set i 0} {$i < 3} {incr i} {
            lset R $i 2 [expr {-[lindex $R $i 2]}]
        }
    }

    return $R
}


# RMSD between aligned fragments
proc rmsd {P Q} {
    set n [llength $P]
    set sum 0.0
    for {set i 0} {$i < $n} {incr i} {
        set a [lindex $P $i]
        set b [lindex $Q $i]
        set dx [expr {[lindex $a 0] - [lindex $b 0]}]
        set dy [expr {[lindex $a 1] - [lindex $b 1]}]
        set dz [expr {[lindex $a 2] - [lindex $b 2]}]
        set sum [expr {$sum + $dx*$dx + $dy*$dy + $dz*$dz}]
    }
    return [expr {sqrt($sum / $n)}]
}

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
