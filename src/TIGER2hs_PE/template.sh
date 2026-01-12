#!/bin/bash -l
################# "REMD schwarm-TIGER2hsPE"  ###################
# NOTE Requires topology and pdb for both explicit and implict solvation
# NOTE expecting ${jobname}_HmassR.top ${jobname}.pdb
#################################################################
jobname="ng_beta28"						# jobname
basedir="/scratch-emmy/usr/mvbdelcm/${jobname}"				# base directory
cd $basedir
numreplicas="8"						# number of replicas (incl. 1 energy replica)
mintemp="300"							# lower bound for temperature
maxtemp="450"							# upper bound for temperature
tigerheat="0"							# tiger2 heat steps
tigersample="4000"						# tiger2 sample steps
tigerquench="2000"						# tiger2 quenching steps
tigersolute="1-352"		                # resid selection for solute
tigersolvent="353-60549"				# resid selection for solvent
tigerignore=""							# resid selection to ignore part of the solute during shellsearch
tigeromm="remd/omm_impl_spe_charmm.py"  # OpenMM implicit solvent energy evaluation
tigerimplplatform="CPU"                 # OpenMM platform used for energy calculation [CPU, CUDA, OpenCL]
tigerimplgb="OBC2"                      # Amber implicit solvent model [HCT, OBC1, OBC2, GBn, GBn2] 
tigerimplsaltcon="0.15"                 # implicit solvent salt conc (ignored for pbc systems)
tigerimplpbc="1"                        # activate pbc for implicit solvent calculations
tigerimpltop="${jobname}_impl.psf"      # amber topology for the implicit solvent system
tigershell="0"                          # number if shell solvents
tigerspace="3"
tigerconheat="1" 
ommpre=""
ommsuff="-pf $(cat omm_charmm_params.dat)"
minruns="2"						        #perform a minimization of this number of runs before REMD; 0 for off; no exchange of course
randruns="0"                            #perform MD at maxtemp for randomization of replicas
randtemp="600"
numruns="500000"						#number of exchange
runsperrestart="10"						#number of exchange attempts before writing a set of restarts
#-
swarmseedreps="4"                       #Replicas to be seeded with current best state
swarmpoolreps="2"                       #Replicas to be seeded from structure pool
swarmpoolsize="16"                      #Size of structure pool of previous best states
swarmpooldiversity="0.05"               #Percent diversity threshold to add states to the pool
swarmcycle="1"                          #Number of TIGER runs between swarm cycles
swarmdynamic="0"                        #Increase swarmcycle based on success rate of finding better structures
swarmcolvar="min(test),max(test2)"      #Name of colvars to overwatch (min(name) or max(name)) seperated by ","
swarmlimits="test=60"                   #Hard capping of upper/lower limit in colvar colvar=limit seperated by ","
swarmpermitresets="1"                   #New best states in one CV override in other CVs too (reasonable in multi-CV)
#-
CellX="125"							    #Cell dimensions 0 for none
CellY="125"
CellZ="125"
pmeon="1"   							#0 for off 1 for on
pmegridspacing="1"
timestep="4"
fullelectfrequency="1"
stepspercycle="20"
cutoff="10"
switchdist="9"							#0 for off
pairlistdist="12"
rigidbonds="all"
dielectric="1.0"
twoawayx="0"							#0 for off 1 for on
twoawayy="0"
remdpressuregen="1"						#generate temperature dependend target pressures for langevinpiston
langevinpiston="1"						#make sure replica.namd sets langevinpistontemp
langevinpistontarget="1.01325"
langevinpistondecay="100"
langevinpistonperiod="200"
constraints="0"
constraintsfile="${jobname}_constraints.pdb"
constraintscol="B"
constraintscaling="1"
extrabonds="0"
extrabondsfile="ionbonds.dat"
wrapWater="1"							#some kind of wrapping may be required while using seref as diffusing solvent may leave pdb format
wrapAll="0"							    #when simulation is running for long time.
wrapNearest="1"
#######################################################################

cat >  ${jobname}_colvar.conf <<END
colvarsTrajFrequency 500
indexFile index.ndx

#Overwatch distance between DI and DV as minimal swarm colvar --
colvar {
  name test
  distance {
     group1 {
      indexGroup DI
    }
    group2 {
      indexGroup DV
    }
    forceNoPBC yes
  }
}
#-----------------------------------------------------------

#Overwatch distance between DII and DIV as minimal swarm colvar --
colvar {
  name test2
  distance {
     group1 {
      indexGroup DII
    }
    group2 {
      indexGroup DIV
    }
    forceNoPBC yes
  }
}
#-----------------------------------------------------------

#Keep DI intact but allow for induced fit --------
colvar {
  name DI_rmsd
  rmsd {
    atoms {
      indexGroup DI_colvar
    }
    refPositionsFile ${jobname}_DI_colvar.pdb
    refPositionsCol B
  }
}

harmonicWalls {
  name DI_rmsd_potential
  colvars DI_rmsd
  upperWallConstant 10
  upperWalls 2
}
#--------------------------------------------------

#Keep DII intact but allow for induced fit --------
colvar {
  name DII_rmsd
  rmsd {
    atoms {
      indexGroup DII_colvar
    }
    refPositionsFile ${jobname}_DII_colvar.pdb
    refPositionsCol B
  }
}

harmonicWalls {
  name DII_rmsd_potential
  colvars DII_rmsd
  upperWallConstant 10
  upperWalls 2
}
#---------------------------------------------------

#Keep DIII intact but allow for induced fit --------
colvar {
  name DIII_rmsd
  rmsd {
    atoms {
      indexGroup DIII_colvar
    }
    refPositionsFile ${jobname}_DIII_colvar.pdb
    refPositionsCol B
  }
}

harmonicWalls {
  name DIII_rmsd_potential
  colvars DIII_rmsd
  upperWallConstant 10
  upperWalls 2
}
#--------------------------------------------------

#Keep DIV intact but allow for induced fit --------
colvar {
  name DIV_rmsd
  rmsd {
    atoms {
      indexGroup DIV_colvar
    }
    refPositionsFile ${jobname}_DIV_colvar.pdb
    refPositionsCol B
  }
}

harmonicWalls {
  name DIV_rmsd_potential
  colvars DIV_rmsd
  upperWallConstant 10
  upperWalls 2
}
#--------------------------------------------------

#Keep DV intact but allow for induced fit ---------
colvar {
  name DV_rmsd
  rmsd {
    atoms {
      indexGroup DV_colvar
    }
    refPositionsFile ${jobname}_DV_colvar.pdb
    refPositionsCol B
  }
}

harmonicWalls {
  name DV_rmsd_potential
  colvars DV_rmsd
  upperWallConstant 10
  upperWalls 2
}
#--------------------------------------------------
END

cat >  ${jobname}_rand.conf <<END

END

############################# BINARY OPTIONS #########################
############################# CPU VERSION ############################
module load anaconda2
conda activate t2hs_openmm

module load openmpi-4.0.4
module load namd-2.14-cuda9
# module load namd-2.11
PRE="--mca btl ^openib --mca btl_tcp_if_include eth0"
# module swap PrgEnv-cray PrgEnv-intel
# module load craype-hugepages8M
# module load namd/namd-2.11
ADD="++ppn $SLURM_CPUS_PER_TASK +devicesperreplica 1"
######################################################################
############################# BRAIN VERSION ############################
# module load anaconda/anaconda2
# conda activate t2hs_openmm
# 
# export OPENMM_CPU_THREADS=32
# 
# module load openmpi/openmpi-3.0.0
# module load namd/namd-2.13-cuda9
# #module load namd/namd-2.14
# # module load namd-2.11
# PRE="--bind-to none taskset -c 0-127"
# # module swap PrgEnv-cray PrgEnv-intel
# # module load craype-hugepages8M
# # module load namd/namd-2.11
# ADD="++ppn $SLURM_CPUS_PER_TASK +devicesperreplica 1"
######################################################################
############################ HLRN VERSION ############################
# #export OPENMM_CPU_THREADS=32
# 
# export TMPDIR=$LOCAL_TMPDIR
# export TMP=$LOCAL_TMPDIR
# module load openmpi/openmpi-4.0.2
# module load namd/namd-2.13
# #ADD="++ppn $SLURM_CPUS_PER_TASK +setcpuaffinity"
# #ADD="+setcpuaffinity"
# PRE="--bind-to none"
# 
# #must be last to not mess up Lua Modules
# module load anaconda/anaconda2
# conda activate t2hs_openmm
######################################################################


#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
######################################################################
######## DO NOT CHANGE BELOW UNLESS YOU KNOW WHAT U R DOING ##########
######################################################################
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

############################# Prepearation hacks #####################
cd $basedir
remd/make_output_dirs.sh output $numreplicas
cd $basedir
unset SGE_ROOT
######################################################################

############################# BASE CONFIG ############################
cat > ${jobname}_base.namd <<END
####### basic input #######
structure ${jobname}_HmassR.psf
paraTypeCharmm on
parameters toppar_water_ions_namd.str
parameters toppar/par_all36_carb.prm
parameters toppar/par_all36_cgenff.prm
parameters toppar/par_all36_lipid.prm
parameters toppar/par_all36m_prot.prm
parameters toppar/par_all36_na.prm
parameters toppar/par_interface.prm
parameters toppar/toppar_all36_carb_glycolipid.str
parameters toppar/toppar_all36_carb_glycopeptide.str
parameters toppar/toppar_all36_carb_imlab.str
parameters toppar/toppar_all36_label_fluorophore.str
parameters toppar/toppar_all36_lipid_cardiolipin.str
parameters toppar/toppar_all36_lipid_detergent.str
parameters toppar/toppar_all36_lipid_ether.str
parameters toppar/toppar_all36_lipid_hmmm.str
parameters toppar/toppar_all36_lipid_inositol.str
parameters toppar/toppar_all36_lipid_lps.str
parameters toppar/toppar_all36_lipid_miscellaneous.str
parameters toppar/toppar_all36_lipid_model.str
parameters toppar/toppar_all36_lipid_prot.str
parameters toppar/toppar_all36_na_nad_ppi.str
parameters toppar/toppar_all36_na_rna_modified.str
parameters toppar/toppar_all36_prot_fluoro_alkanes.str
parameters toppar/toppar_all36_prot_heme.str
parameters toppar/toppar_all36_prot_modify_res.str
parameters toppar/toppar_all36_prot_na_combined.str
parameters toppar/toppar_all36_prot_retinol.str
parameters toppar/toppar_all36_prot_stapling.str
###########################

# pbc 
cellBasisVector1     ${CellX} 0      0
cellBasisVector2     0      ${CellY} 0
cellBasisVector3     0      0      ${CellZ}
cellOrigin           0      0      0 

# simulation options

# Basic dynamics
exclude                 scaled1-4
1-4scaling              1
dielectric              ${dielectric}

# Simulation space partitioning
if {${switchdist} > 0} {
  switching         on
  switchdist        ${switchdist}
}
cutoff                  ${cutoff}
pairlistdist        ${pairlistdist}
rigidbonds      ${rigidbonds}

# Multiple timestepping
timestep                ${timestep}
stepspercycle           ${stepspercycle}
fullElectFrequency      ${fullelectfrequency}

if {${pmeon}} {
  # PME
  pme yes
  pmegridspacing ${pmegridspacing}
}
#usecuda2 no

# temperature control
langevin on
langevinDamping 10.0

# Constant Pressure Control (variable volume)
if {${langevinpiston}} {
  LangevinPiston on
  if {!${remdpressuregen}} {
    LangevinPistonTarget ${langevinpistontarget}
  }
  LangevinPistonPeriod ${langevinpistonperiod}
  LangevinPistonDecay  ${langevinpistondecay}
  useGroupPressure yes
}

# constraints
if {${constraints}} {
  constraints     on
  consref   ${constraintsfile}
  conskfile ${constraintsfile}
  conskcol  ${constraintscol}
  constraintScaling ${constraintscaling}
}

# Extra bonds
if {${extrabonds}} {
  extrabonds     on
  extraBondsFile   ${extrabondsfile}
}

if {${wrapAll}} {
  wrapAll on
}

if {${wrapWater}} {
  wrapWater on
}

if {${wrapNearest}} {
  wrapNearest on
}

if {${twoawayx}} {
  twoawayx yes
}
if {${twoawayy}} {
  twoawayy yes
}
END
######################################################################

######################## PARALLEL TEMPERING CONF #####################
cat > ${jobname}_remd.conf <<END
set num_replicas ${numreplicas}
set min_temp ${mintemp}
set max_temp ${maxtemp}
set tigerheat ${tigerheat}
set tigersample ${tigersample}
set tigerquench ${tigerquench}
set tigersolute ${tigersolute}
set tigersolvent ${tigersolvent}
set tigerignore "${tigerignore}"
set tigerconheat ${tigerconheat}
set tigershell ${tigershell}
set tigerspace ${tigerspace}
set tigeromm ${tigeromm}
set tigerimplgb ${tigerimplgb}
set tigerimplsaltcon ${tigerimplsaltcon}
set tigerimplpbc ${tigerimplpbc}
set tigerimpltop ${tigerimpltop}
set tigerimplplatform ${tigerimplplatform}
set ommpre "${ommpre}"
set ommsuff "${ommsuff}"
set num_runs ${numruns}
set runs_per_restart ${runsperrestart}
set swarmseedreps ${swarmseedreps}
set swarmpoolreps ${swarmpoolreps}
set swarmpoolsize ${swarmpoolsize}
set swarmpooldiversity ${swarmpooldiversity}
set swarmcycle ${swarmcycle}
set swarmdynamic ${swarmdynamic}
set swarmcolvar "${swarmcolvar}"
set swarmlimits "${swarmlimits}"
set swarmpermitresets $swarmpermitresets
set namd_config_file "${jobname}_base.namd"
set colvar_config_file "${jobname}_colvar.conf"
set output_root "output/%s/${jobname}"
set remdpressuregen ${remdpressuregen}
set minruns ${minruns}
set randruns ${randruns}
set rand_temp ${randtemp}
set colvar_config_file_rand "${jobname}_rand.conf"
set initial_pdb    "${jobname}.pdb"
END
######################################################################

############################# JOBS ###################################
jobnum=0
while ls output/${jobname}.job${jobnum}.restart*.tcl > /dev/null 2>&1
do
  jobnum=$(($jobnum+1))
  restartfile=$(ls -tr output/${jobname}.job$(($jobnum-1)).restart*.tcl | tail -n 1)
done

cat > job${jobnum}.conf <<END
source ${jobname}_remd.conf

if {$jobnum < 1} {
  margin 5
} else {
  margin 5
  source [format ${restartfile} ""]
}

if { ! [catch numPes] } { source remd/replica.namd }
END
mpirun $PRE namd2 $ADD +idlepoll +replicas $numreplicas job${jobnum}.conf +stdout output/%d/job${jobnum}.%d.log >> output/job${jobnum}.log 2> output/job${jobnum}.out
######################################################################
