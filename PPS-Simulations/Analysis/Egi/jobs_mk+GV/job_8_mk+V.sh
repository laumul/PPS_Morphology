#!/bin/sh
 
# allocate 4 nodes
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#
# job name
#SBATCH --job-name=EgiGV_mk+V_8


# jobs always start in $HOME
# change to work directory
cd /home/woody/gwpa/gwpa005h/Egi/

filename=8.nex
numchar=5
runID=8_mk+GV
data=data_mk+GV
out=output_mk+GV

# run your revbayes job
bash sim-start.sh $filename $numchar $runID $data $out mk+V




