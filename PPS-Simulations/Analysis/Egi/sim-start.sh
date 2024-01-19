#!/bin/bash 

#23.08.2023

rb=/home/hpc/gwpa/gwpa005h/rb2/revbayes/projects/cmake/rb

filename=$1
num_states=$2
runID=$3
data=$4
out=$5
model=$6
iterations=10000

file="Convergence.csv"
counter=1

### pps simulation workflow
## 1. normal mcmc
while [ ! -f $out/$model/$runID/output_$model/$file ] && [ $counter -lt  7 ]
do
#
it=$iterations
rb_command="filename=\"$filename\"; model=\"$model\"; num_states=\"$num_states\"; runID=\"$runID\"; prior=\"$prior\"; data=\"$data\"; out=\"$out\";  it=\"$it\"; source(\"scripts/pps_MCMC_Simulation.Rev\");"
#
echo  $rb_command | $rb
#

Rscript CheckConvergence.r $out $model $runID

if [[ "$counter" -eq 1 ]]
then
iterations=50000
else 

let iterations=iterations+50000
fi
let counter=counter+1
done

echo -e "$it" >> $out/$model/$runID/input.txt

#


## 2. simulate nex files
##
#Rscript MorphoSim.r $filename $runID $model $data $out
##
## 3. mcmc for simulated nex file
rb_command="filename=\"$filename\"; model=\"$model\"; num_states=\"$num_states\"; runID=\"$runID\"; data=\"$data\"; out=\"$out\"; source(\"scripts/PosteriorPredictive_MCMC.Rev\");"
##
##
#echo  $rb_command | $rb
##
### 4. tree summaries
rb_command="filename=\"$filename\"; model=\"$model\"; num_states=\"$num_states\"; runID=\"$runID\"; data=\"$data\"; out=\"$out\";  source(\"scripts/PosteriorPredictive_TreeSummary.Rev\");"
##
##
#echo  $rb_command | $rb
##
### 5. P-values
rb_command="filename=\"$filename\"; model=\"$model\"; num_states=\"$num_states\"; runID=\"$runID\"; data=\"$data\"; out=\"$out\";  source(\"scripts/PosteriorPredictive_Pvalues.Rev\");"
##
#echo $rb_command | rb
#
#
#
#
##### stich multistate files back together
#Rscript MorphoSim.r $filename $runID $model $data $out
##
#
#
#export runID
#export out
#export model
#export data
#export filename
#bash combine_nex.sh
#
#
#

echo -e "$filename\n$num_states" >> $out/$model/$runID/input.txt
#



 


