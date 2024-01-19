main <- function(){
  args <- commandArgs(trailingOnly = TRUE)
  
  library(gtools)
  
  g=2
  muds <- c("mk+V", "mk+GV", "mk+GVmultistate", "mk+Vmultistate")
  
  num_reps <- 1000
  
  mod <- as.data.frame(list.files("output_mk+GV/"))
  colnames(mod) <- "V1"
  
  #for (g in 1:length(muds)){
  
  dfiles <- gtools::mixedsort(sort(list.files(paste0("data_",muds[g],"/"))))
  
  r <- mixedsort(sort(list.files(paste0("output_",muds[g],"/mk"))))
  pdf("CumulativeMeans.pdf")
  par(mfrow=c(2,2))
  for (s in 1:20){
    message(paste0("working on sim ", s))
    
    ######## add names and paths #########
    filename <- paste0("data_" ,muds[g], "/" ,dfiles[s])
    
    # paths should have the end of the directory. i.e. output/models/"path"
    path_sims_nexus <- paste0("/",r[s],"/output/pps_Nov_post_sims/")
    path_sims_csv <- paste0("/",r[s],"/results/simulated_inference_pps_Nov.csv")
    path_emp_csv <- paste0("/",r[s],"/results/empirical_inference_pps_Nov.csv")
    path_emp_mcc <- paste0("/",r[s],"/output/MCC.tre")
    
  
    
    
  
    
    
    ######## Tree Length #######################################
    
    # to get dimentions for results container 
    simname <- paste0("output_",muds[g],"/mk", path_sims_csv)
    sim <- read.csv(simname)
    final <- matrix(nrow = num_reps, ncol=length(mod$V1))
    colnames(final) <- mod$V1
    
    
  
    
    for (k in 1:length(mod$V1)){
      
      simname <- paste0("output_",muds[g],"/",mod[k,], path_sims_csv)
      empname <- paste0("output_",muds[g],"/",mod[k,], path_emp_csv)
      sim <- read.csv(simname)
      sim <- sim[1:num_reps,]
      emp <- read.csv(empname)
      # standard deviatio across simulations 
      std <- sd(sim$mean_tl)
      
      plot(cumsum(sim$mean_tl) / seq_along(sim$mean_tl), main = paste0(r[s], "_",mod[k,], "_TreeLength" ))
      abline(v=400, lty=2)
      
      cum_avg <- cumsum(sim$mean_tl) / seq_along(sim$mean_tl)
      
      Variance <- matrix(nrow = 38, ncol=1)
      colnames(Variance) <- "Mean"
      
      
      for (i in 1:38){
        if (i == 1){
          bina <- 20
        } else {
          bina <- bina + 20
        }
        
        binb <- + 20 
        Variance[i,"Mean"]<-  mean(cum_avg[bina:binb])
        
      }
      
      
      plot(Variance, type="b", main = paste0(r[s], "_",mod[k,], "_TreeLength"))
           abline(v=14, lty=2)
           
           
    }
 
    
    
    ######## Robinson Foulds #######################################
    
  
    
    for (k in 1:length(mod$V1)){
      
      simname <- paste0("output_",muds[g],"/",mod[k,], path_sims_csv)
      empname <- paste0("output_",muds[g],"/",mod[k,], path_emp_csv)
      sim <- read.csv(simname)
      sim <- sim[1:num_reps,]
      emp <- read.csv(empname)
     
      plot(cumsum(sim$mean_rf) / seq_along(sim$mean_rf), main = paste0(r[s], "_",mod[k,], "_RobinsonsFould" ))
      abline(v=400, lty=2)
      
      cum_avg <- cumsum(sim$mean_rf) / seq_along(sim$mean_rf)
      
      Variance <- matrix(nrow = 38, ncol=1)
      colnames(Variance) <- "Mean"
      
      
      for (i in 1:38){
        if (i == 1){
          bina <- 20
        } else {
          bina <- bina + 20
        }
        
        binb <- + 20 
        Variance[i,"Mean"]<-  mean(cum_avg[bina:binb])
        
      }
      
      
      plot(Variance, type="b", main = paste0(r[s], "_",mod[k,], "_RobinsonsFould"))
      abline(v=14, lty=2)
      
      
      
     
      
    }
    
    
  
    
    
    # ###########Consistency Index & Retention Index ################
    
    for (n in 1:2){
     
      
      # to get dimentions for results container 
      loopname <- paste0("output_",muds[g],"/mk",path_sims_nexus)
      loop <- list.files(loopname)
      
      empdata <- ape::read.nexus.data(filename)
      phylev <- unique(unlist(unique( empdata)))
      models <- matrix(nrow=num_reps, ncol=length(mod$V1))
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
        
        ConInd <- matrix(nrow=num_reps, ncol = 2)
        colnames(ConInd) <- c("test-stat" , "dist")
        for (i in 1:num_reps){
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
        
        
        if ( n == 1){
        plot(cumsum(ConInd[,1]) / seq_along(ConInd[,1]), main = paste0(r[s], "_",mod[k,], "_ConsistencyIndex" ))
        abline(v=400, lty=2)
        
        cum_avg <- cumsum(ConInd[,1]) / seq_along(ConInd[,1])
        
        Variance <- matrix(nrow = 38, ncol=1)
        colnames(Variance) <- "Mean"
        
        
        for (i in 1:38){
          if (i == 1){
            bina <- 20
          } else {
            bina <- bina + 20
          }
          
          binb <- + 20 
          Variance[i,"Mean"]<-  mean(cum_avg[bina:binb])
          
        }
        
        
        plot(Variance, type="b", main = paste0(r[s], "_",mod[k,], "_ConsistencyIndex"))
        abline(v=14, lty=2)
        
        
        } else {
          
          plot(cumsum(ConInd[,1]) / seq_along(ConInd[,1]), main = paste0(r[s], "_",mod[k,], "_RetentionIndex" ))
        abline(v=400, lty=2)

          cum_avg <- cumsum(ConInd[,1]) / seq_along(ConInd[,1])
          
          Variance <- matrix(nrow = 38, ncol=1)
          colnames(Variance) <- "Mean"
          
          
          for (i in 1:38){
            if (i == 1){
              bina <- 20
            } else {
              bina <- bina + 20
            }
            
            binb <- + 20 
            Variance[i,"Mean"]<-  mean(cum_avg[bina:binb])
            
          }
          
          
          plot(Variance, type="b", main = paste0(r[s], "_",mod[k,], "_RetentionIndex"))
          abline(v=14, lty=2)
                  }
                  
                  
                  

        
      }
        
      }  
    
    
    ############## Data Statistics ##########################################
    
    # the distance metric here is set to "ged" 
    # can change to "mord" "gc", or "red"
    
    library(Claddis)
    opts <- c( "gc", "ged")
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
      container = matrix(ncol = length(mod$V1), nrow=num_reps)
      colnames(container) = mod$V1
      es_container = matrix(ncol = length(mod$V1), nrow=num_reps)
      colnames(es_container) = mod$V1
      
      
      for (i in 1:length(mod$V1)){
        
        for (j in 1:num_reps){
          
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
        
        
    
  
  
      plot(cumsum(container[,i]) / seq_along(container[,i]), main = paste0(r[s], "_",mod[i,], "_",opts[p]))
      abline(v=400, lty=2)
      
      cum_avg <- cumsum(container[,i]) / seq_along(container[,i])
      
      Variance <- matrix(nrow = 38, ncol=1)
      colnames(Variance) <- "Mean"
      
      
      for (i in 1:38){
        if (i == 1){
          bina <- 20
        } else {
          bina <- bina + 20
        }
        
        binb <- + 20 
        Variance[i,"Mean"]<-  mean(cum_avg[bina:binb])
        
      }
      
      
      plot(Variance, type="b", main = paste0(r[s], "_",mod[i,], "_",opts[p]))
      abline(v=14, lty=2)
      
      }
  
}

    
  }
  
  dev.off()

}
main()



