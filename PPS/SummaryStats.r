main <- function(){
  args <- commandArgs(trailingOnly = TRUE)

library(gtools)


mod <- as.data.frame(list.files("output/"))
colnames(mod) <- "V1"

testz <- c("TL", "RF", "CI", "RI","ged","gc")
effect_sizes <- matrix(ncol=length(mod$V1), nrow=length(testz))
colnames(effect_sizes) <- mod$V1
rownames(effect_sizes) <- testz


k=1

pdf("Test_Statistics.pdf")
par(mfrow=c(2,2))
simname <- paste0("output/",mod[k,], "/results/simulated_inference_APW-example.csv")
sim <- read.csv(simname)

final <- matrix(nrow = length(sim$simID), ncol=length(mod$V1))
colnames(final) <- mod$V1
for (k in 1:length(mod$V1)){
  
  simname <- paste0("output/",mod[k,], "/results/simulated_inference_APW-example.csv")
  empname <- paste0("output/",mod[k,], "/results/empirical_inference_APW-example.csv")
  sim <- read.csv(simname)
  emp <- read.csv(empname)
  std <- sd(sim$mean_tl)
  
  for (i in 1:length(sim$simID)){
    numsim <- sim$mean_tl[i]
    final[i,mod[k,]] <- (emp$mean_tl - numsim) / std
  }
  
  
}

mods <- stack(as.data.frame(final))

x <- boxplot(mods$values ~ mods$ind, lty=1, main = "Tree Length")
abline(h=0, lty = 2)
for (m in 1:length(mod$V1)){
  
  effect_sizes["TL", mod[m,]] <- x$stats[3,m]
}



#### Robinson Foulds

final <- matrix(nrow = length(sim$simID), ncol=length(mod$V1))
colnames(final) <- mod$V1
for (k in 1:length(mod$V1)){
  
  simname <- paste0("output/",mod[k,], "/results/simulated_inference_APW-example.csv")
  empname <- paste0("output/",mod[k,], "/results/empirical_inference_APW-example.csv")
  sim <- read.csv(simname)
  emp <- read.csv(empname)
  std <- sd(sim$mean_rf)
  
  for (i in 1:length(sim$simID)){
    numsim <- sim$mean_rf[i]
    final[i,mod[k,]] <- (emp$mean_rf - numsim) / std
  }
  
  
}

mods <- stack(as.data.frame(final))


x <- boxplot(mods$values ~ mods$ind, lty=1, main = "Robinsons Fould")
abline(h=0, lty = 2)
for (m in 1:length(mod$V1)){
  
  effect_sizes["RF", mod[m,]] <- x$stats[3,m]
}





### Consistency Index and Retention Index
k=1
for (n in 1:2){
  
  loopname <- paste0("output/", mod[k,], "/output/PPS_Morpho_post_sims/")
  loop <- list.files(loopname)
  
  empdata <- ape::read.nexus.data(paste0("data/Egi_etal_2005a_paleobiodb.nex"))
  phylev <- unique(unlist(unique( empdata)))
  models <- matrix(nrow=length(loop), ncol=length(mod$V1))
  colnames(models) <- mod$V1
  for (k in 1:length(mod$V1)){
    
    loopname <- paste0("output/",mod[k,], "/output/PPS_Morpho_post_sims/")
    loop <- list.files(loopname)
    
    empname <- paste0("output/",mod[k,] , "/output/MCC.tre" )
    emp1 <- ape::read.nexus(empname)
    empd <-  phangorn::phyDat(empdata, type="USER", levels=phylev, return.index = TRUE)
    
    if (n == 1){
      empci <- phangorn::CI(emp1, empd)
    }
    
    else {
      empci <- phangorn::RI(emp1, empd)
    }
    
    ConInd <- matrix(nrow=length(loop), ncol = 2)
    colnames(ConInd) <- c("RI" , "dist")
    for (i in 1:length(loop)){
      file <- paste0(loopname , loop[i] , "/seq.nex")
      simdata <- ape::read.nexus.data(file)
      simd <-  phangorn::phyDat(simdata, type="USER", levels=phylev, return.index = TRUE)
      
      if (n ==1 ){
        ConInd[i,1] <- phangorn::CI(emp1, simd)
      }
      else {
        ConInd[i,1] <- phangorn::RI(emp1, simd) 
      }
      
      
    }
    
    standdev <- sd(ConInd[,1])
    
    for (j in 1:length(loop)){
      
      ConInd[j,2] <- (empci - ConInd[j,1] )/ standdev
      
    }
    
    models[,mod[k,]] <- ConInd[,2]
    
  }
  
  if (n == 1 ){
    mods <- stack(as.data.frame(models))
    x <- boxplot(mods$values ~ mods$ind, lty=1, main ="Consistency Index" , par(cex.axis=0.4))
    abline(h=0, lty = 2)
    
  } else {
    mods <- stack(as.data.frame(models))
    x <- boxplot(mods$values ~ mods$ind, lty=1, main ="Retention Index" , par(cex.axis=0.4))
    abline(h=0, lty = 2)
  }
  
  
  
  for (m in 1:length(mod$V1)){
    if (n ==1){
      effect_sizes["CI", mod[m,]] <- x$stats[3,m]
    } else {
      effect_sizes["RI", mod[m,]] <- x$stats[3,m]
    }
  }
  
}


library(Claddis)

emp <- ape::read.nexus.data(paste0("data/Egi_etal_2005a_paleobiodb.nex"))

TestStats <- c("ged", "gc")

for (l in 1:length(TestStats)){

emp_matrix <- matrix(ncol=length(emp[[1]]), nrow = length(names(emp)))
rownames(emp_matrix) <- names(emp)

for (v in 1:length(emp)){
  emp_matrix[names(emp)[v],] <- as.character(emp[[v]])
}

emp_matrix[emp_matrix =="?"] <- NA
emp_matrix <- build_cladistic_matrix(emp_matrix, ordering = rep("unord", ncol(emp_matrix)))

# Get morphological distances
distances <- Claddis::calculate_morphological_distances(cladistic_matrix = emp_matrix, distance_metric = TestStats[l])

# Show distance matrix:
emp_mean =mean(distances$distance_matrix)

k=1
# create container to hold the output
loopname <- paste0("output/", mod[k,], "/output/PPS_Morpho_post_sims/")
loop <- list.files(loopname)
container = matrix(ncol = length(mod$V1), nrow=length(loop))
colnames(container) = mod$V1
es_container = matrix(ncol = length(mod$V1), nrow=length(loop))
colnames(es_container) = mod$V1


for (i in 1:length(mod$V1)){
  
  for (j in 1:length(loop)){
    
    temp_sim <- ape::read.nexus.data(paste0("output/",mod[i,], "/output/PPS_Morpho_post_sims/", 
                                            loop[j], "/seq.nex"))
    
    # create mytrix for character data
    sim_matrix <- matrix(ncol=length(temp_sim[[1]]), nrow = length(temp_sim))
    rownames(sim_matrix) <- names(temp_sim)
    
    for (v in 1:length(temp_sim)){
      sim_matrix[names(temp_sim)[v],] <- as.character(temp_sim[[v]])
    }
    
    sim_matrix[sim_matrix =="?"] <- NA
    # build cladistic matrix
    sim_matrix <- build_cladistic_matrix(sim_matrix, ordering = rep("unord", ncol(sim_matrix)))
    
    # Get morphological distances 
    distances <- Claddis::calculate_morphological_distances(cladistic_matrix = sim_matrix, distance_metric = TestStats[l])
    
    
    
    # Show distance matrix:
    mean_sim_dm <- mean(distances$distance_matrix)
    container[j,mod[i,]] <- mean_sim_dm
  }
}

for (i in 1:length(mod$V1)){
  std <- sd(container[,mod[i,]])
  for (j in 1:length(loop)){
    es_container[j,mod[i,]] <- (emp_mean - container[j,mod[i,]])/std
    
    
  }
}

mods <- stack(as.data.frame(es_container))
boxplot(mods$values ~ mods$ind, lty=1, main =TestStats[l] , par(cex.axis=0.4))
abline(h=0, lty = 2)

for (m in 1:length(mod$V1)){
  
  effect_sizes[TestStats[l], mod[m,]] <- x$stats[3,m]
}
}

effect_sizes <- as.data.frame(effect_sizes)




dev.off()  


write.csv(effect_sizes, "Test_Statistics.csv")






}

main()