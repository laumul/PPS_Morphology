setwd("Desktop/PhD/Posterior-Prediction/Final/Empirical-Inf/")

setwd("../Empirical-PPS/")
datas <- list.files("data")

#pdf("Results.pdf")

mod <- as.data.frame(list.files("output/Agnolin/")[1:7])
colnames(mod) <- "V1"
datasets <-  as.data.frame(list.files("output/"))
par(mar = c(1,1,0.5,0.5))
par(mfrow=c(2,2))
#par(mfrow=c(1,1))
vecp <- c("Mk", "Mk+G", "MkV", "MkV+G"," MkP+G", "MkVP", "MkVP+G")
for (w in 1:length(datas)){
  
  
  runID <- datasets[w,]
  filename <- paste0("data/", datas[w])
  
  
  testz <- c("CI", "RI")
  effect_sizes <- matrix(ncol=length(mod$V1), nrow=length(testz))
  colnames(effect_sizes) <- mod$V1
  rownames(effect_sizes) <- testz
for (n in 1:2){
    
    loopname <- paste0("output/", runID, "/mk/output", "/PPS_Morpho_post_sims/")
    loop <- list.files(loopname)
    
    empdata <- ape::read.nexus.data(filename)
    phylev <- unique(unlist(unique( empdata)))
    models <- matrix(nrow=length(loop), ncol=length(mod$V1))
    colnames(models) <- mod$V1
    for (k in 1:length(mod$V1)){
      
      loopname <- paste0("output/",  runID, "/", mod[k,],"/output/PPS_Morpho_post_sims/")
      loop <- list.files(loopname)
      
      empname <- paste0("output/",  runID, "/", mod[k,],"/output/MCC.tre" )
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
      models <- models[,c(1,2,6,4,3,7,5)]
      mods <- stack(as.data.frame(models))
      
      x <- boxplot(mods$values ~ mods$ind, lty=1,  par(cex.axis=1), horizontal = TRUE, las=1, ylab = "",
                   xlab = "", yaxt = "n")
      axis(2, at = 1:7, label = vecp, las = 2)
      rect(par("usr")[1], par("usr")[3],
           par("usr")[2], par("usr")[4],
           col = "#f7f7f7")
 
      boxplot(mods$values ~ mods$ind, add = TRUE, lty=1,  par(cex.axis=1), yaxt = "n",
              horizontal = TRUE, las=1, ylab = "", xlab = "", col = "steelblue4")
      abline(v=0, lty = 2, lwd =2)
    } else {
      
      models <- models[,c(1,2,6,4,3,7,5)]
      mods <- stack(as.data.frame(models))
      
      x <- boxplot(mods$values ~ mods$ind, lty=1,  par(cex.axis=1), yaxt = "n",
                   horizontal = TRUE, las=1, ylab = "", xlab = ""
                  )
      rect(par("usr")[1], par("usr")[3],
           par("usr")[2], par("usr")[4],
           col = "#f7f7f7")
      
      
      
      boxplot(mods$values ~ mods$ind, add = TRUE, lty=1,  par(cex.axis=1), yaxt = "n",
              horizontal = TRUE, las=1, ylab = "", xlab = "", col = "steelblue4")
      abline(v=0, lty = 2, lwd =2)
    }
    
    
    
    for (m in 1:length(mod$V1)){
      if (n ==1){
        effect_sizes["CI", mod[m,]] <- x$stats[3,m]
      } else {
        effect_sizes["RI", mod[m,]] <- x$stats[3,m]
      }
    }
    
  }
  
  
  
  effect_sizes <- as.data.frame(effect_sizes)
  write.csv(effect_sizes, paste0("output/",datasets[w,], "/test-stats.csv"))
  
  
  
  # dev.off()  
  
  
}






