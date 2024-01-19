### 06.05.2022

### Analysis just for stepping stone stuff

# define where the revbayes executable is
rb=/home/hpc/gwpa/gwpa005h/rb2/revbayes/projects/cmake/rb


filename=$1
model=$2
data=$3

# get the number of character states
num_states=6

# run ID name
#ID=$(echo $filename | awk -F'_' '{print $1}')
runID=$(echo $filename  | tr -d .nex)


## pps simulation workflow
# 1. normal mcmc
rb_command="filename=\"$filename\"; model=\"$model\"; num_states=\"$num_states\"; runID=\"$runID\"; data=\"$data\";  source(\"scripts/Steppingstone.Rev\");"

echo $rb_command | $rb



 

 


