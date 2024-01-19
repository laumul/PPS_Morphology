#!/bin/sh
 
# allocate 4 nodes
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#
# job name
#SBATCH --job-name=Bou_mk+Gm


# jobs always start in $HOME
# change to work directory
cd /home/woody/gwpa/gwpa005h/MorphoSim-Empirical/


# run your revbayes job
filename="Bourdon_etal_2009a_paleobiodb.nex"
num_states="3" 
runID="Bourdon"
model="mk+Gmultistate"

bash sim-start $filename $num_states, $runID, $model