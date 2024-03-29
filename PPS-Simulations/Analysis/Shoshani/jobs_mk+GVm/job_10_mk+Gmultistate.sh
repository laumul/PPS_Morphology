#!/bin/sh
 
# allocate 4 nodes
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#
# job name
#SBATCH --job-name=ShGVm_mk+Gm_10


# jobs always start in $HOME
# change to work directory
cd /home/woody/gwpa/gwpa005h/MorphoSim

filename=10.nex
numchar=6
runID=10_mk+GVmultistate
data=data_mk+GVmultistate
out=output_mk+GVmultistate

# run your revbayes job
bash sim-start.sh $filename $numchar $runID $data $out mk+Gmultistate




