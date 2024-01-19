## run all datasets using all models and compare the outputs



files=$(ls data)

# how many matrices to include
numtaxa=122




# folder name
counter=1

# line in results file to take info from
line=2



for i in $files
do

for j in {1..7}
do
model=$(cat models.csv | sed -n $j'p')



#while [ $ct -le $numtaxa ]
#do
# get the name of the nexus file
filename=$(cat dataset.csv | awk '{print $1}' | cut -f1 -d";" | sed -n $line'p')


# get the number of character states
num_states=$(cat datasets.csv | awk '{print $NF}' | sed 's/;//g' | sed -n $line'p')


# run ID name
ID=$(echo $filename | awk -F'_' '{print $1}')
runID=$(echo $counter"_"$ID)


## pps simulation workflow
# 1. normal mcmc
rb_command="filename=\"$filename\"; model=\"$model\"; num_states=\"$num_states\"; runID=\"$runID\"; prior=\"$prior\";  source(\"scripts/pps_MCMC_Simulation.Rev\");"

echo  $rb_command | rb

done

let line=line+1
let counter=counter+1 
done
