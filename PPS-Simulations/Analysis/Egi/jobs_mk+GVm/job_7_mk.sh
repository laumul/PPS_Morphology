#!/bin/sh
 
# allocate 4 nodes
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#
# job name
#SBATCH --job-name=EgiGVm_mk+mk_7


# jobs always start in $HOME
# change to work directory
cd /home/woody/gwpa/gwpa005h/Egi/

filename=7.nex
numchar=5
runID=7_mk+GVmultistate
data=data_mk+GVmultistate
out=output_mk+GVmultistate

# run your revbayes job
bash sim-start.sh $filename $numchar $runID $data $out mk




