#!/usr/bin/env tclsh
# Analytical, steady-state estimate of completed baseline segment lengths.
# Tcl 8.5+. No packages or simulation files required.

proc usage {} {
    puts {Estimate continuous baseline segments and milestone injection coverage.

Usage:
  tclsh estimate_shooting_segments.tcl --run-ps PS --replicas N \
      --basereps B --shootingreps S --shootingcycle C [options]

Required:
  --run-ps PS          Baseline time per TIGER run in ps (>0).
  --replicas N         Total replicas (>=1).
  --basereps B         Baseline replicas (1..N).
  --shootingreps S     Replicas injected per shooting cycle (0..N).
  --shootingcycle C    TIGER runs between shooting cycles (>=1).

Optional:
  --runs R            TIGER runs for segment-count estimates. Default: 1000.
  --threshold-ns NS   Long-segment threshold in ns (>0). Default: 1.
  --exchange-prob E   Baseline/elevated pair acceptance (0..1). Default: 0.2.
  --retention P       Measured baseline retention per run (0..1); overrides
                      the exchange estimate. Default: calculated from pairing.
  --milestones M      Path milestone count. Default: 0 (coverage disabled).
  --coverage-trials T Coverage simulations (>=1). Default: 500.
  --seed SEED         Coverage seed (1..2147483646). Default: 42.
  --help              Show this help.

Example:
  tclsh estimate_shooting_segments.tcl --run-ps 64 --replicas 16 \
      --basereps 8 --shootingreps 2 --shootingcycle 25 \
      --runs 5000 --milestones 128

Estimates:
  Segment length is the long-run mean. Counts refer to completed segments
  across all baseline replicas during --runs, starting with fresh segments.
  Segments still open at the end are not counted. The ns fraction describes
  the long-run distribution; finite-run counts also account for run duration.
  Coverage reports average shooting cycles and TIGER runs to inject each
  target fraction of milestones at least once, plus segment counts by then.
  Time is simulated baseline-run equivalent, not wall time.

Assumptions:
  Random disjoint exchange pairs; one random replica is unpaired for odd N.
  Injection replicas are drawn from the full ensemble. Each receives a
  separate milestone draw weighted by 1/(previous injections+1), with counts
  updated immediately. Coverage starts with all injection counts at zero.
  First injection follows C runs; exchange precedes injection at a boundary.
  Estimates describe sampling opportunities, not kinetic convergence.}
}

proc finite {value} {
    if {![string is double -strict $value]} {return 0}
    if {[catch {expr {abs($value) < Inf}} answer]} {return 0}
    return $answer
}

proc estimate_segments {runps n b s cycle p} {
    set q [expr {double($s)/$n}]
    set exchange [expr {1.0-$p}]
    set injection [expr {$p*$q/$cycle}]
    set hazard [expr {$exchange+$injection}]
    if {$hazard == 0.0} {
        set meanruns Inf
        set meanps Inf
    } else {
        set meanruns [expr {1.0/$hazard}]
        set meanps [expr {$runps/$hazard}]
    }
    return [list $q $exchange $injection $meanruns $meanps]
}

proc segment_tail {p q cycle minruns} {
    # Starts immediately after an injection boundary have weight (1-p+p*q);
    # other phases have weight (1-p). Survival spans minruns-1 boundaries.
    set total [expr {$cycle*(1.0-$p)+$p*$q}]
    if {$total == 0.0} {return "undefined"}
    set m [expr {$minruns-1}]
    set periods [expr {$m/$cycle}]
    set remainder [expr {$m%$cycle}]
    return [expr {pow($p,$m)*pow(1.0-$q,$periods)*
        (1.0-$remainder*(1.0-$p)*$q/$total)}]
}

proc expected_counts {run_count baseline p q cycle minruns} {
    # Each boundary closes baseline segments through exchange, or through
    # injection of survivors. Do not double-count exchange plus injection.
    set total [expr {$baseline*((1.0-$p)*$run_count +
        $p*$q*($run_count/$cycle))}]
    set long 0.0
    if {$run_count >= $minruns} {
        # Group ending boundaries by phase modulo the injection interval.
        set phases [expr {min($cycle,$run_count-$minruns+1)}]
        set exchange_survival [expr {pow($p,$minruns-1)}]
        for {set offset 0} {$offset < $phases} {incr offset} {
            set end [expr {$minruns+$offset}]
            set multiplicity [expr {1+($run_count-$end)/$cycle}]
            set prior_injections [expr {($end-1)/$cycle-($end-$minruns)/$cycle}]
            set survival [expr {$exchange_survival*pow(1.0-$q,$prior_injections)}]
            set hazard [expr {1.0-$p+($end%$cycle == 0 ? $p*$q : 0.0)}]
            set long [expr {$long+$baseline*$multiplicity*$survival*$hazard}]
        }
    }
    return [list $total $long]
}

proc coverage_simulations {milestones injections trials seed} {
    set draws $injections
    set targets {}
    foreach fraction {0.5 0.9 0.95 0.99 1.0} {
        lappend targets [expr {int(ceil($fraction*$milestones))}]
    }
    set samples [lrepeat 5 {}]
    set rng $seed
    for {set trial 0} {$trial < $trials} {incr trial} {
        set counts [lrepeat $milestones 0]
        set weights [lrepeat $milestones 1.0]
        set total [expr {double($milestones)}]
        set seen 0
        set draw 0
        set nexttarget 0
        while {$seen < $milestones} {
            # Private Park-Miller stream, independent of Tcl's global rand().
            set rng [expr {(wide($rng)*16807)%2147483647}]
            set key [expr {($rng/2147483647.0)*$total}]
            set selected [expr {$milestones-1}]
            for {set i 0} {$i < $milestones} {incr i} {
                set key [expr {$key-[lindex $weights $i]}]
                if {$key < 0.0} {set selected $i; break}
            }
            set old [lindex $counts $selected]
            if {$old == 0} {incr seen}
            set updated [expr {$old+1}]
            set weight [expr {1.0/($updated+1.0)}]
            set total [expr {$total+$weight-[lindex $weights $selected]}]
            lset counts $selected $updated
            lset weights $selected $weight
            incr draw
            while {$nexttarget < 5 && $seen >= [lindex $targets $nexttarget]} {
                set cycles [expr {($draw+$draws-1)/$draws}]
                lset samples $nexttarget [concat [lindex $samples $nexttarget] [list $cycles]]
                incr nexttarget
            }
        }
    }
    return $samples
}

proc print_coverage {milestones injections interval runps trials seed baseline p q} {
    puts "\nMILESTONE COVERAGE ESTIMATE"
    puts "Pool: $milestones milestones; weight: 1/(injection count + 1)."
    if {$injections == 0} {
        puts "No injection coverage: shootingreps is zero."
        return
    }
    puts stderr "Estimating milestone coverage..."
    set samples [coverage_simulations $milestones $injections $trials $seed]
    puts [format "%-18s %16s %16s %18s %16s" "Milestones seeded" "Shooting cycles" "TIGER runs" "ns/replica eq." "Segments"]
    set index 0
    foreach target {50 90 95 99 100} {
        set sum 0.0
        foreach value [lindex $samples $index] {set sum [expr {$sum+$value}]}
        set cycles [expr {wide(ceil($sum/$trials))}]
        set runs [expr {$cycles*$interval}]
        set segments [lindex [expected_counts $runs $baseline $p $q $interval 1] 0]
        puts [format "%-18s %16d %16d %18.3f %16.1f" "$target%" $cycles $runs [expr {$runs*$runps/1000.0}] $segments]
        incr index
    }
    puts "Average estimates from $trials simulations (seed $seed), rounded up to whole cycles."
    puts "Segments: expected completed baseline segments by the listed run count."
    puts "Seeded means injected at least once; it does not mean adequately sampled."
    puts "Starts with zero injection counts. Time excludes TIGER overhead and is not wall time."
}

proc main {arguments} {
    array set opt {--retention auto --exchange-prob 0.2 --threshold-ns 1 --milestones 0 --coverage-trials 500 --seed 42 --runs 1000}
    set allowed {--run-ps --replicas --shootingreps --shootingcycle --basereps --retention --exchange-prob --threshold-ns --milestones --coverage-trials --seed --runs}
    while {[llength $arguments]} {
        set key [lindex $arguments 0]
        if {$key eq "--help" || $key eq "-h"} {usage; return}
        if {[lsearch -exact $allowed $key] < 0} {error "Unknown option: $key (use --help)"}
        if {[llength $arguments] < 2} {error "Missing value for $key"}
        set opt($key) [lindex $arguments 1]
        set arguments [lrange $arguments 2 end]
    }
    foreach key {--run-ps --replicas --shootingreps --shootingcycle --basereps} {
        if {![info exists opt($key)]} {error "Missing required option: $key (use --help)"}
    }
    foreach key {--replicas --shootingreps --shootingcycle --basereps} {
        if {![string is integer -strict $opt($key)]} {error "$key must be an integer"}
    }
    set run_count $opt(--runs)
    if {![string is integer -strict $run_count] || $run_count < 1} {
        error "--runs must be a positive integer"
    }
    set trials $opt(--coverage-trials)
    set seed $opt(--seed)
    if {![string is integer -strict $trials] || $trials < 1} {
        error "--coverage-trials must be a positive integer"
    }
    if {![string is integer -strict $seed] || $seed < 1 || $seed > 2147483646} {
        error "--seed must be an integer from 1 to 2147483646"
    }
    set milestones $opt(--milestones)
    if {![string is integer -strict $milestones] || $milestones < 0} {
        error "--milestones must be a positive integer (0 disables coverage)"
    }
    set n $opt(--replicas)
    set b $opt(--basereps)
    set s $opt(--shootingreps)
    set c $opt(--shootingcycle)
    set dt $opt(--run-ps)
    set threshold $opt(--threshold-ns)
    if {![finite $threshold] || $threshold <= 0} {error "--threshold-ns must be finite and positive"}
    if {$n < 1 || $b < 1 || $b > $n || $s < 0 || $s > $n || $c < 1} {
        error "Require N>=1, 1<=B<=N, 0<=S<=N and C>=1"
    }
    if {![finite $dt] || $dt <= 0} {error "--run-ps must be finite and positive"}
    set exchangeprob $opt(--exchange-prob)
    if {![finite $exchangeprob] || $exchangeprob < 0 || $exchangeprob > 1} {
        error "--exchange-prob must be between 0 and 1"
    }
    if {$opt(--retention) eq "auto"} {
        if {$n == 1} {
            set p 1.0
        } elseif {$n % 2 == 0} {
            set p [expr {1.0-$exchangeprob*double($n-$b)/($n-1)}]
        } else {
            set p [expr {1.0-$exchangeprob*double($n-$b)/$n}]
        }
        set assumption "random disjoint pairing; cross-temperature acceptance $exchangeprob"
        if {$n > 1 && $n % 2} {append assumption "; one random replica unpaired"}
    } else {
        set p $opt(--retention)
        if {![finite $p] || $p < 0 || $p > 1} {error "--retention must be between 0 and 1"}
        set assumption "user-supplied retention; --exchange-prob is overridden"
    }
    lassign [estimate_segments $dt $n $b $s $c $p] q ex inj runs duration
    puts "\nBASELINE SEGMENT ESTIMATE"
    puts [format "%-36s %12d" "Total replicas:" $n]
    puts [format "%-36s %12d" "Baseline replicas:" $b]
    puts [format "%-36s %12.3f ps" "Baseline time per TIGER run:" $dt]
    puts [format "%-36s %12d runs" "Injection interval:" $c]
    puts [format "%-36s %12d" "Replicas injected per cycle:" $s]
    puts [format "%-36s %12.4f" "Baseline retention probability:" $p]
    puts [format "%-36s %12.4f" "Reseeding probability per cycle:" $q]
    puts "\nAssumption: $assumption"
    if {$duration == Inf} {
        puts "No segment termination predicted; duration is limited by total runtime."
    } else {
        puts [format "\n%-36s %12.3f runs" "Mean segment length:" $runs]
        puts [format "%-36s %12.3f ps" "Mean segment duration:" $duration]
        puts [format "%-36s %12.6f ns" "Mean segment duration:" [expr {$duration/1000.0}]]
        puts [format "%-36s %11.1f%%" "Termination attributed to exchange:" [expr {100*$ex/($ex+$inj)}]]
        puts [format "%-36s %11.1f%%" "Termination attributed to injection:" [expr {100*$inj/($ex+$inj)}]]
        if {$p < 1} {
            puts [format "%-36s %12.3f ps" "Mean without milestone injections:" [expr {$dt/(1-$p)}]]
        } else {
            puts "Without injections, no termination is predicted."
        }
    }
    set minruns [expr {max(1,wide(ceil(1000.0*$threshold/$dt)))}]
    set tail [segment_tail $p $q $c $minruns]
    puts ""
    if {$tail eq "undefined"} {
        puts "Fraction reaching >= $threshold ns: undefined (no completed segments predicted)."
    } else {
        puts [format "%-36s %11.4g%%" "Segments reaching >= $threshold ns:" [expr {100.0*$tail}]]
        puts [format "%-36s %12d runs" "Minimum length to reach threshold:" $minruns]
        puts "Fraction of segments, not time; periodic-injection model."
    }
    lassign [expected_counts $run_count $b $p $q $c $minruns] segments longsegments
    puts "\nEXPECTED SEGMENTS OVER $run_count TIGER RUNS"
    puts [format "%-36s %12.1f" "Completed segments (all baseline):" $segments]
    puts [format "%-36s %12.1f" "Completed segments >= $threshold ns:" $longsegments]
    puts "Counts start from fresh segments and exclude segments still open at the end."
    if {$milestones > 0} {
        print_coverage $milestones $s $c $dt $trials $seed $b $p $q
    }
    puts "\nModel estimate, not guaranteed achievable sampling; see --help for assumptions."
}

if {[file normalize [info script]] eq [file normalize $argv0]} {
    if {[catch {main $argv} message]} {
        puts stderr "ERROR: $message"
        exit 1
    }
}
