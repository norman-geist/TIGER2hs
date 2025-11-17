# TM-Score function in TCL
# Takes two lists of XYZ coordinates as input: coordsA and coordsB
# coordsA and coordsB are lists of XYZ coordinates in the format {x1 y1 z1} {x2 y2 z2} ...
# Returns the TM-Score for the two structuress

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
        for {set i 0} {$i < 3} {incr i} {
            lset R $i 2 [expr {-[lindex $R $i 2]}]
        }
    }

    return $R
}

proc calc_TM-score {args} {

    set coordsA [lindex $args 0]
    set coordsB [lindex $args 1]
    
    if {[llength $coordsA] != [llength $coordsB]} {
        error "Coordinate list lengths must match"
    }
    
    set num_atoms [expr {[llength $coordsA] / 3}]
    set vecsA {}
    set vecsB {}
    for {set i 0} {$i < $num_atoms} {incr i} {
        lappend vecsA [lrange $coordsA [expr {$i*3}] [expr {$i*3+2}]]
        lappend vecsB [lrange $coordsB [expr {$i*3}] [expr {$i*3+2}]]
    }
    
    #align structures
    set comA [center_of_mass $vecsA]
    set comB [center_of_mass $vecsB]  
    set vecsA0 [translate_coords $vecsA $comA]
    set vecsB0 [translate_coords $vecsB $comB]
    set R [kabsch_rotation $vecsA0 $vecsB0]
    set vecsB0A [rotate_coords $vecsB0 $R]

    # Define the reference value d_0 (based on the size of the proteins)
    set d_0 [expr {1.24 * pow($num_atoms, 1.0/3.0)}]
    
    set sum_tm_score 0.0
    set n 0

    # Loop over all atoms and compute RMSD for each pair
    for {set i 0} {$i < $num_atoms} {incr i} {
        # Get the XYZ coordinates for the current atom in both structures
        lassign [lindex $vecsA0 $i] x1 y1 z1
        lassign [lindex $vecsB0A $i] x2 y2 z2

        # Compute the squared difference in positions
        set rmsd [expr {pow($x1 - $x2, 2) + pow($y1 - $y2, 2) + pow($z1 - $z2, 2)}]

        # Compute the TM-score contribution for this atom
        set tm_contrib [expr {1.0 / (1.0 + $rmsd / ($d_0 * $d_0))}]
        
        # Sum the TM-score contributions
        set sum_tm_score [expr {$sum_tm_score + $tm_contrib}]
        incr n
    }

    # Return the normalized TM-score
    return [expr {1 - ($sum_tm_score / $n)}]
}

