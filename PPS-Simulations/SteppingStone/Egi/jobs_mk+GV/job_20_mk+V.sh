#!/bin/sh
 
# allocate 4 nodes
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#
# job name
#SBATCH --job-name=SS_EgiGV_mk+V_20


# jobs always start in $HOME
# change to work directory
cd /home/woody/gwpa/gwpa005h/SS-MorphoSim-Egi/


# run your revbayes job

filename=20.nex
model=mk+V
data=mk+GV

bash stepping_stone.sh $filename $model $data

