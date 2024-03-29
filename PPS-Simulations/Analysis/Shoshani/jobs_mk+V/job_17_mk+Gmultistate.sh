#!/bin/sh
 
# allocate 4 nodes
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#
# job name
#SBATCH --job-name=ShV_mk+Gm_17


# jobs always start in $HOME
# change to work directory
cd /home/woody/gwpa/gwpa005h/MorphoSim

filename=17.nex
numchar=6
runID=17_mk+V
data=data_mk+V
out=output_mk+V

# run your revbayes job
bash sim-start.sh $filename $numchar $runID $data $out mk+Gmultistate




