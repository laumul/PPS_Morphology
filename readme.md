## Data and scripts from the study "Assessing the Adequacy of Morphological Models used in Palaeobiology"


**1. Empirical-Inf**

This contains everything for the empirical analysis. 
Data directroy contains all the morphological data sets analysed. 
Scripts direcroy contains all rev scripts for inference and Rscripts for down stream analysis. 

`Start.sh` will run all data sets under the 7 models.

**2. PPS Simulations**

This contains the set up for the simulation study. 

###### Simulation 

 
The data directory contains the output from the empirical inference of Egi and Shoshani for 
four different models. Simulated data sets are in the data_model directories.
Scripts contains all the rev scripts used for the inference. 
Data was simualted using the `Sim.r` file.To simulate data, you need to specify the model, data set, and runID. This requires eidting the Sim.r file directly, lines 41-45


###### Analysis 
Contains two directories, Egi and Shsoahni. Each contains the posterior predictive simlations workflow. 
To start the analysis run the scripts in the jobs_"model" files. All rev scripts used are in the scripts folder.
There are four other scripts used in the analysis 

- `sim-start.sh:` this contains all the commands for each individual run
- `CheckConvergence.r:` this file ensures that the initail MCMC has reach convergence before simulating new data sets
- `Anaylsis.r:` this file calculates the tests statsitics. Need to specify on line 3 and 11 which simulation set up (which model) you are analysising.
- `Cumulative.r:` this file was used to assess the number of replicatie data sets are necessary to determine the adequacy of a model. This was onyl used once, after which we simulated 500 data sets for the rest of the models


###### Stepping Stone
Contains two directories, Egi and Shsoahni. Each contains the scripts for the stepping stone analysis. To start the analysis run the scripts in the jobs_"model" files. All rev scripts used are in the scripts folder.
There are two other scripts used in the analysis 

- `stepping_stone.sh:` his contains all the commands for each individual run
- `ml_ss.sh`: gets all the marginal likelihoods and puts them into a cvs file. Requires you provide the model name as an argument.


**3. PPS Empirical**
