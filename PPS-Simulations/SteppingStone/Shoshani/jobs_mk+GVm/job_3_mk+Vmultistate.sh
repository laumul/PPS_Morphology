#!/bin/sh
 
# allocate 4 nodes
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#
# job name
#SBATCH --job-name=SS_ShGVm_mk+Vm_3


# jobs always start in $HOME
# change to work directory
cd /home/woody/gwpa/gwpa005h/SS-MorphoSim-Shoshani/


# run your revbayes job

filename=3.nex
model=mk+Vmultistate
data=mk+GVmultistate

bash stepping_stone.sh $filename $model $data

