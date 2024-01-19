#!/bin/sh
 
# allocate 4 nodes
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#
# job name
#SBATCH --job-name=Sch_mk+GVm


# jobs always start in $HOME
# change to work directory
cd /home/woody/gwpa/gwpa005h/MorphoSim-Empirical/


# run your revbayes job
filename="Schoch_et_Sues_2013a_paleobiodb.nex"
num_states="4" 
runID="Schoch"
model="mk+GVmultistate"

bash sim-start $filename $num_states, $runID, $model