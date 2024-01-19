#!/bin/sh
 
# allocate 4 nodes
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#
# job name
#SBATCH --job-name=Shos_mk+Gm


# jobs always start in $HOME
# change to work directory
cd /home/woody/gwpa/gwpa005h/MorphoSim-Empirical/


# run your revbayes job
filename="Shoshani_etal_2006a_paleobiodb.nex"
num_states="6" 
runID="Shoshani"
model="mk+Gmultistate"

bash sim-start $filename $num_states, $runID, $model