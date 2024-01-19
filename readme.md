## Data and scripts from the study "Assessing the Adequacy of Morphological Models used in Palaeobiology"


**1. Empirical-Inf**

This contains everything for the empirical analysis. 
Data directory contains all the morphological data sets analysed. 
Scripts directory contains all revbayes scripts for inference and Rscripts for down stream analysis. 

`Start.sh` will run all data sets under the 7 models.

**2. PPS Simulations**

This contains the set up for the simulation study. 

###### Simulation 

 
The data directory contains the output from the empirical inference of Egi and Shoshani for 
four different models. Simulated data sets are in the data_model directories.
Scripts contains all the revbayes scripts used for the inference. 
Data was simulated using the `Sim.r` file. To simulate data, you need to specify the model, data set, and runID. This requires editing the Sim.r file directly, lines 41-45


###### Analysis 
Contains two directories, Egi and Shsoahni. Each contains the posterior predictive simulations’ workflow. 
To start the analysis run the scripts in the jobs_"model" files. All revbayes scripts used are in the scripts folder.
There are five other scripts used in the analysis 

- `sim-start.sh:` this contains all the commands for each individual run
- `CheckConvergence.r:` this file ensures that the initial MCMC has reach convergence before simulating new data sets
- `MorphoSim.r:` simulates data sets in R using the phangorn R package.
- `Anaylsis.r:` this file calculates the tests statistics. Need to specify on line 3 and 11 which simulation set up (which model) you are analysing.
- `Cumulative.r:` this file was used to assess the number of replicates data sets are necessary to determine the adequacy of a model. This was only used once, after which we simulated 500 data sets for the rest of the models


###### Stepping Stone
Contains two directories, Egi and Shsoahni. Each contains the scripts for the stepping stone analysis. To start the analysis run the scripts in the jobs_"model" files. All revbayes scripts used are in the scripts folder.
There are two other scripts used in the analysis 

- `stepping_stone.sh:` his contains all the commands for each individual run
- `ml_ss.sh`: gets all the marginal likelihoods and puts them into a cvs file. Requires you provide the model name as an argument.


**3. PPS Empirical**

The scripts required to carry out posterior predictive simulations on the empirical data sets. 
All the revbayes scripts are in the scripts directory. Running the files in the jobs_"datasetname" directories will start the analysis.
There are four other scripts used in the analysis

- `MorphoSim.r:` simulates data sets in R using the phangorn R package.
- `Test_Statistics.r:` calculates all test statistics
- `Fig-6.r:` calculates consistency index and retention index test stats
- `P-vals.r:` calculates the posterior predictive p-values
