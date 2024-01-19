setwd("Desktop/PhD/Posterior-Prediction/Final/Empirical-PPS/output/")


files <- list.files()
data <- list.files("../data/")
for ( k in 1:length(files) ){

models <- list.files(files[k])[1:7]
Cal_pVaules <- matrix(nrow=length(models), ncol=4)
rownames(Cal_pVaules) <- models
colnames(Cal_pVaules) <- c("Lower-1-tailed", "Upper-1-tailed", "Two-tailed", "Midpoint")
emp_data <- ape::read.nexus.data(paste0("../data/", data[k]))
phylev <- unique(unlist(unique(emp_data)))

for (j in 1:length(models)){
  sim_replicaties <- c()
empt <- ape::read.nexus(paste0(files[k],"/",models[j],"/output/MCC.tre"))

emp_data_phan <-  phangorn::phyDat(emp_data, type="USER", levels=phylev, return.index = TRUE)
empCI <- phangorn::CI(empt, emp_data_phan)


reps <- list.files(paste0(files[k],"/",models[j],"/output/PPS_Morpho_post_sims/"))

for (i in 1:length(reps)){
  sim_data <- ape::read.nexus.data(paste0(files[k],"/",models[j],"/output/PPS_Morpho_post_sims/", reps[i], "/seq.nex"))
  sim_data_phan <-  phangorn::phyDat(sim_data, type="USER", levels=phylev, return.index = TRUE)
  simCI <- phangorn::CI(empt, sim_data_phan)
  sim_replicaties <- rbind(sim_replicaties, simCI)
}

hist(sim_replicaties, main = models[j], xlim = c(0,1))
abline(v=empCI)

pvalues <- pVal(empCI, sim_replicaties)
mid_p_values <- pvalues[2] + pvalues[3]*0.5
two_tailed_p_value <- 2*min(pvalues[2]+pvalues[3], pvalues[1]+pvalues[3])
Cal_pVaules[models[j], "Lower-1-tailed"] <- pvalues[2]
Cal_pVaules[models[j], "Upper-1-tailed"] <- pvalues[1]
Cal_pVaules[models[j], "Two-tailed"] <- two_tailed_p_value
Cal_pVaules[models[j], "Midpoint"] <- mid_p_values
}

write.csv(Cal_pVaules, paste0(files[k],"/CI.csv"))


for (j in 1:length(models)){
  sim_replicaties <- c()
  empt <- ape::read.nexus(paste0(files[k],"/",models[j],"/output/MCC.tre"))
  
  emp_data_phan <-  phangorn::phyDat(emp_data, type="USER", levels=phylev, return.index = TRUE)
  empCI <- phangorn::RI(empt, emp_data_phan)
  
  
  reps <- list.files(paste0(files[k],"/",models[j],"/output/PPS_Morpho_post_sims/"))
  
  for (i in 1:length(reps)){
    sim_data <- ape::read.nexus.data(paste0(files[k],"/",models[j],"/output/PPS_Morpho_post_sims/", reps[i], "/seq.nex"))
    sim_data_phan <-  phangorn::phyDat(sim_data, type="USER", levels=phylev, return.index = TRUE)
    simCI <- phangorn::RI(empt, sim_data_phan)
    sim_replicaties <- rbind(sim_replicaties, simCI)
  }
  
  hist(sim_replicaties, main = models[j], xlim = c(0,1))
  abline(v=empCI)
  
  pvalues <- pVal(empCI, sim_replicaties)
  mid_p_values <- pvalues[2] + pvalues[3]*0.5
  two_tailed_p_value <- 2*min(pvalues[2]+pvalues[3], pvalues[1]+pvalues[3])
  Cal_pVaules[models[j], "Lower-1-tailed"] <- pvalues[2]
  Cal_pVaules[models[j], "Upper-1-tailed"] <- pvalues[1]
  Cal_pVaules[models[j], "Two-tailed"] <- two_tailed_p_value
  Cal_pVaules[models[j], "Midpoint"] <- mid_p_values
}

write.csv(Cal_pVaules, paste0(files[k],"/RI.csv"))

}

pVal <- function(e,s){
  greaterThan <- 0
  lessThan <- 0
  Equal <- 0
  for (i in 1:length(s)){
    if (e < s[i]){
      greaterThan <- greaterThan + 1
    } else if (e > s[i]){
      lessThan <- lessThan + 1
    } else if (e == s[i]){
      Equal <- Equal + 1
    }
  }
  pvals <- cbind(greaterThan/length(s), lessThan/length(s), Equal/length(s))
  return(pvals)

}




  
  
  






