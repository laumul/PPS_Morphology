#!/bin/sh
 
# allocate 4 nodes
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#
# job name
#SBATCH --job-name=SS_EgiVm_mk+V_3


# jobs always start in $HOME
# change to work directory
cd /home/woody/gwpa/gwpa005h/SS-MorphoSim-Egi/


# run your revbayes job

filename=3.nex
model=mk+V
data=mk+Vmultistate

bash stepping_stone.sh $filename $model $data

