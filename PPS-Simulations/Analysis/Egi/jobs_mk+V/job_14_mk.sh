#!/bin/sh
 
# allocate 4 nodes
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#
# job name
#SBATCH --job-name=EgiV_mk+mk_14


# jobs always start in $HOME
# change to work directory
cd /home/woody/gwpa/gwpa005h/Egi/

filename=14.nex
numchar=5
runID=14_mk+V
data=data_mk+V
out=output_mk+V

# run your revbayes job
bash sim-start.sh $filename $numchar $runID $data $out mk




