#!/bin/sh
 
# allocate 4 nodes
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#
# job name
#SBATCH --job-name=SS_ShV_mk+Vm_15


# jobs always start in $HOME
# change to work directory
cd /home/woody/gwpa/gwpa005h/SS-MorphoSim-Shoshani/


# run your revbayes job

filename=15.nex
model=mk+Vmultistate
data=mk+V

bash stepping_stone.sh $filename $model $data

