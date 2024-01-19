#!/bin/sh
 
# allocate 4 nodes
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#
# job name
#SBATCH --job-name=SS_ShGV_mk+GVm_14


# jobs always start in $HOME
# change to work directory
cd /home/woody/gwpa/gwpa005h/SS-MorphoSim-Shoshani/


# run your revbayes job

filename=14.nex
model=mk+GVmultistate
data=mk+GV

bash stepping_stone.sh $filename $model $data

