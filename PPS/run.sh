## How to run PPS from terminal 

# step one: run initial MCMC infernce 
# once the model has been specified and the number of states set in the scritps 

rb 
source("scripts/pps_MCMC.Rev")


# step two: Simulate data in R 

Rscript MorphoSim.r $filename $model $NumSim

#e.g. Rscript MorphoSim.r Egi_etal_2005a_paleobiodb.nex mk+GV 5

# you will get a warning output here. You can ignore that. 

# step three: Run inference on simulated data sets 

rb
source("scripts/pps_MCMC_Simulations.Rev")

# step four: summarise results in revbayes 
rb
source("scripts/pps_TreeSummary.Rev")


# step five: summary statistics in R 

Rscript SummaryStats.r $filename