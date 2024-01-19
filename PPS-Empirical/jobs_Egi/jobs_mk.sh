#!/bin/sh
 
# allocate 4 nodes
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#
# job name
#SBATCH --job-name=Egi_mk


# jobs always start in $HOME
# change to work directory
cd /home/woody/gwpa/gwpa005h/MorphoSim-Empirical/


# run your revbayes job
filename="Egi_etal_2005a_paleobiodb.nex"
num_states="5" 
runID="Egi"
model="mk"

bash sim-start $filename $num_states, $runID, $model