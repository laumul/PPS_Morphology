#!/bin/sh
 
# allocate 4 nodes
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#
# job name
#SBATCH --job-name=ShVm_mk+V_1


# jobs always start in $HOME
# change to work directory
cd /home/woody/gwpa/gwpa005h/MorphoSim

filename=1.nex
numchar=6
runID=1_mk+Vmultistate
data=data_mk+Vmultistate
out=output_mk+Vmultistate

# run your revbayes job
bash sim-start.sh $filename $numchar $runID $data $out mk+V




