## Data and scripts from the study "Assessing the Adequacy of Morphological Models used in Palaeobiology"


**1. Empirical-Inf**

This contains everything for the empirical analysis. 
Data directroy contains all the morphological data sets analysed. 
Scripts direcroy contains all rev scripts for inference and Rscripts for down stream analysis. 

`Start.sh` will run all data sets under the 7 models.

**2. PPS Simulations**

This contains the set up for the simulation study. 

###### Simulation Directory

 
The Data directory contains the output from the empirical inference of Egi and Shoshani for 
four different models. Simulated data sets are in the data_model directories.
Scripts contains all the rev scripts used for the inference. 
Data was simualted using the `Sim.r` file.To simulate data, you need to specify the model, data set, and runID. This requires eidting the Sim.r file directly, lines 41-45


###### Analysis Directory
Contains two directories, Egi and Shsoahni. Each contains the anaylsis of the data simulated based on that data set. Each directoy contains the same set up.
To start the analysis run the scripts in the jobs_"model" files.
There are four R scripts used in the analysis 

- `sim-start.sh`: this contains all the commands for each individual run
- `CheckConvergence.r`: this file ensures that the initail MCMC has reach convergence before simulating new data sets
- `Anaylsis.r`: this file calculates the tests statsitics. Need to specify on line 3 and 11 which simulation set up (which model) you are analysising.
- `Cumulative.r`: this file was used to assess the number of replicatie data sets are necessary to determine the adequacy of a model. This was onyl used once, after which we simulated 500 data sets for the rest of the models
