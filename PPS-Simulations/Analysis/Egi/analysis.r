#setwd("Desktop/PhD/Posterior-Prediction/Analysis/Simulations/Agnolin/")

main <- function(){
  args <- commandArgs(trailingOnly = TRUE)

library(gtools)

g=3
muds <- c("mk+V", "mk+GV", "mk+GVmultistate", "mk+Vmultistate")

mod <- as.data.frame(list.files("output_mk+GVmultistate/"))
#colnames(mod) <- "V1"
#mod <- as.data.frame(mod[-5,])
colnames(mod) <- "V1"
#for (g in 1:length(muds)){
  
  dfiles <- gtools::mixedsort(sort(list.files(paste0("data_",muds[g],"/"))))

r <- mixedsort(sort(list.files(paste0("output_",muds[g],"/mk"))))

for (s in 1:1){
s=3
  message(paste0("working on sim ", s))
  
  ######## add names and paths #########
  filename <- paste0("data_" ,muds[g], "/" ,dfiles[s])
  
  # paths should have the end of the directory. i.e. output/models/"path"
  path_sims_nexus <- paste0("/",r[s],"/output/pps_Nov_post_sims/")
  path_sims_csv <- paste0("/",r[s],"/results/simulated_inference_pps_Nov.csv")
  path_emp_csv <- paste0("/",r[s],"/results/empirical_inference_pps_Nov.csv")
  path_emp_mcc <- paste0("/",r[s],"/output/MCC.tre")
  
  
  ### matrix to hold median values 
  testz <- c("TL", "RF", "CI", "RI","mord", "gc", "ged", "red")
  effect_sizes <- matrix(ncol=length(mod$V1), nrow=length(testz))
  colnames(effect_sizes) <- mod$V1
  rownames(effect_sizes) <- testz
  
  
  ### list to hold all values generated
  Results <- list()
  
  if (!file.exists(paste0("Figures/",muds[g],"/", r[s])) ){
    
    dir.create(paste0("Figures/",muds[g],"/", r[s]))
  }
  
  ### Save all to pdf
  pdf(paste0("Figures/" ,muds[g],"/", r[s],"/Test_Statistics.pdf"))
  par(mfrow=c(2,2))
  
  
  ######## Tree Length #######################################
  
  # to get dimentions for results container 
  simname <- paste0("output_",muds[g],"/mk", path_sims_csv)
  sim <- read.csv(simname)
  final <- matrix(nrow = length(sim$simID), ncol=length(mod$V1))
  colnames(final) <- mod$V1
  
  for (k in 1:length(mod$V1)){
#message(mod[k,])    
    simname <- paste0("output_",muds[g],"/",mod[k,], path_sims_csv)
    empname <- paste0("output_",muds[g],"/",mod[k,], path_emp_csv)
    sim <- read.csv(simname)
    emp <- read.csv(empname)
    # standard deviatio across simulations 
    std <- sd(sim$mean_tl)
    
    for (i in 1:length(sim$simID)){
      numsim <- sim$mean_tl[i]
      # the number of standard deviations the simulated is from the empirical value
      final[i,mod[k,]] <- (emp$mean_tl - numsim) / std
    }
    
    
  }
  
  Results[["Tree Length"]]<- final
  
  mods <- stack(as.data.frame(final))
  
  
  x <- boxplot(mods$values ~ mods$ind, lty=1)
  abline(h=0, lty = 2)
  for (m in 1:length(mod$V1)){
    effect_sizes["TL", mod[m,]] <- x$stats[3,m]
  }
  
  
  
  
  ######## Robinson Foulds #######################################
  
  # results container 
  final <- matrix(nrow = length(sim$simID), ncol=length(mod$V1))
  colnames(final) <- mod$V1
  
  for (k in 1:length(mod$V1)){
    
    simname <- paste0("output_",muds[g],"/",mod[k,], path_sims_csv)
    empname <- paste0("output_",muds[g],"/",mod[k,], path_emp_csv)
    sim <- read.csv(simname)
    emp <- read.csv(empname)
    # standard deviatio across simulations 
    std <- sd(sim$mean_rf)
    
    for (i in 1:length(sim$simID)){
      numsim <- sim$mean_rf[i]
      # the number of standard deviations the simulated is from the empirical value
      final[i,mod[k,]] <- (emp$mean_rf - numsim) / std
    }
    
    
  }
  
  Results[["Robinson Foulds"]]<- final
  mods <- stack(as.data.frame(final))
  
  x <- boxplot(mods$values ~ mods$ind, lty=1)
  abline(h=0, lty = 2)
  for (m in 1:length(mod$V1)){
    
    effect_sizes["RF", mod[m,]] <- x$stats[3,m]
  }
  
  
  
  # ###########Consistency Index & Retention Index ################
  
  for (n in 1:2){
    
    # to get dimentions for results container 
    loopname <- paste0("output_",muds[g],"/mk",path_sims_nexus)
    loop <- list.files(loopname)
    
    empdata <- ape::read.nexus.data(filename)
    phylev <- unique(unlist(unique( empdata)))
    models <- matrix(nrow=length(loop), ncol=length(mod$V1))
    colnames(models) <- mod$V1
    
    for (k in 1:length(mod$V1)){
      
      loopname <- paste0("output_",muds[g], "/",mod[k,], path_sims_nexus)
      loop <- list.files(loopname)
      
      empname <- paste0("output_",muds[g], "/",mod[k,] ,path_emp_mcc )
      emp1 <- ape::read.nexus(empname)
      empd <-  phangorn::phyDat(empdata, type="USER", levels=phylev, return.index = TRUE)
      
      if (n == 1){
        empci <- phangorn::CI(emp1, empd)
      }
      
      else {
        empci <- phangorn::RI(emp1, empd)
      }
      
      ConInd <- matrix(nrow=length(loop), ncol = 2)
      colnames(ConInd) <- c("test-stat" , "dist")
      for (i in 1:length(loop)){
        file <- paste0(loopname , loop[i] , "/seq.nex")
        simdata <- ape::read.nexus.data(file)
        phylevs <- unique(unlist(unique( simdata)))
        simd <-  phangorn::phyDat(simdata, type="USER", levels=phylevs, return.index = TRUE)
        
        if (n ==1 ){
          ConInd[i,1] <- phangorn::CI(emp1, simd)
        }
        else {
          ConInd[i,1] <- phangorn::RI(emp1, simd) 
        }
        
        
      }
      
      standdev <- sd(ConInd[,1])
      
      
      for (j in 1:length(loop)){
        
        # the number of standard deviations the simulated is from the empirical value
        ConInd[j,2] <- (empci - ConInd[j,1] )/ standdev
        
      }
      
      models[,mod[k,]] <- ConInd[,2]
      
      if (n == 1){
        
        Results[["Consistency Index"]] <- models
      }else {
        Results[["Retention Index"]] <- models
      }
      
    }
    
    if (n == 1 ){
      mods <- stack(as.data.frame(models))
      x <- boxplot(mods$values ~ mods$ind, lty=1, main ="Consistency Index" , par(cex.axis=0.4))
      abline(h=0, lty = 2)
      
    } else {
      mods <- stack(as.data.frame(models))
      x <- boxplot(mods$values ~ mods$ind, lty=1, main ="Retention Index", par(cex.axis=0.4))
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
  
  
  
  ############## Data Statistics ##########################################
  
  # the distance metric here is set to "ged" 
  # can change to "mord" "gc", or "red"
  
  library(Claddis)
  opts <- c("mord", "gc", "ged" ,"red")
  for (p in 1:length(opts)){
    
    emp <- ape::read.nexus.data(filename)
    
    
    emp_matrix <- matrix(ncol=length(emp[[1]]), nrow = length(names(emp)))
    rownames(emp_matrix) <- names(emp)
    
    for (v in 1:length(emp)){
      emp_matrix[names(emp)[v],] <- as.character(emp[[v]])
    }
    
    emp_matrix[emp_matrix =="?"] <- NA
    emp_matrix[emp_matrix =="-"] <- NA
    
    emp_matrix <- Claddis::build_cladistic_matrix(emp_matrix, ordering = rep("unord", ncol(emp_matrix)))
    
    # Get morphological distances
    distances <- Claddis::calculate_morphological_distances(cladistic_matrix = emp_matrix, distance_metric =  opts[p],)
    
    # Show distance matrix:
    emp_mean =mean(distances$distance_matrix)
    
    
    # create container to hold the output
    loopname <- paste0("output_",muds[g],"/mk",path_sims_nexus)
    loop <- list.files(loopname)
    container = matrix(ncol = length(mod$V1), nrow=length(loop))
    colnames(container) = mod$V1
    es_container = matrix(ncol = length(mod$V1), nrow=length(loop))
    colnames(es_container) = mod$V1
    
    
    for (i in 1:length(mod$V1)){
      
      for (j in 1:length(loop)){
        
        temp_sim <- ape::read.nexus.data(paste0("output_",muds[g], "/",mod[i,],path_sims_nexus, 
                                                loop[j], "/seq.nex"))
        
        # create mytrix for character data
        sim_matrix <- matrix(ncol=length(temp_sim[[1]]), nrow = length(temp_sim))
        rownames(sim_matrix) <- names(temp_sim)
        
        for (v in 1:length(temp_sim)){
          sim_matrix[names(temp_sim)[v],] <- as.character(temp_sim[[v]])
        }
        
        sim_matrix[sim_matrix =="?"] <- NA
        sim_matrix[sim_matrix =="-"] <- NA
        # build cladistic matrix
        sim_matrix <- Claddis::build_cladistic_matrix(sim_matrix, ordering = rep("unord", ncol(sim_matrix)))
        
        # Get morphological distances 
        distances <- Claddis::calculate_morphological_distances(cladistic_matrix = sim_matrix, distance_metric =  opts[p],)
        
        
        
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
    
    Results[[opts[p]]] <- es_container
    mods <- stack(as.data.frame(es_container))
    x <- boxplot(mods$values ~ mods$ind, lty=1, main = opts[p] , par(cex.axis=0.4))
    abline(h=0, lty = 2)
    
    for (m in 1:length(mod$V1)){
      
      effect_sizes[opts[p], mod[m,]] <- x$stats[3,m]
    }
    
  }
  # sum up the effect sizes across different test. Probably not the best way to do it but currently looking at it 
  effect_sizes <- as.data.frame(effect_sizes)
  
  effect_sizes["Total",] <- 0
  
  for (m in 1:length(mod$V1)){
    effect_sizes["Total",m] <- sum(abs(effect_sizes[1:8,mod[m,]]))
  }
  
  
  
  dev.off()  
  
  
  write.csv(effect_sizes, paste0("Figures/" ,muds[g],"/", r[s],"/Test_Statistics.csv"))
  
}

#}

}

main()



