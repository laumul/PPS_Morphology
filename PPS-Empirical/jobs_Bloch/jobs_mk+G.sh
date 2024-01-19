#!/bin/sh
 
# allocate 4 nodes
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#
# job name
#SBATCH --job-name=Blo_mk+G


# jobs always start in $HOME
# change to work directory
cd /home/woody/gwpa/gwpa005h/MorphoSim-Empirical/


# run your revbayes job
filename="Bloch_etal_2001a_paleobiodb.nex"
num_states="6" 
runID="Bloch"
model="mk+G"

bash sim-start $filename $num_states, $runID, $model