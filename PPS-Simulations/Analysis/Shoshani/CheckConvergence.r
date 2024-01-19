main <- function(){
  args <- commandArgs(trailingOnly = TRUE)
  
  
  ESS_check = function(dat, parameter) {
    coda::effectiveSize(mcmcr::as.mcmc(dat[,parameter])) >= 200
  }
  
  ### create a dataframe to hold all unconverged parameters

      
      dat <- read.table(paste0(args[1], "/",args[2], "/" ,args[3],"/output_", args[2],"/pps_Nov_posterior.log"), header =TRUE)
      
      ## burnin of first 10% of trees 
      burnin <- (length(dat$Iteration) -2) * .1
      dat <- dat[burnin:length(dat$Iteration),]
      parameters <- names(dat)[3:length(names(dat))]
      
      # create matrix to hold output 
      
      ConverTest <- matrix(nrow=1, ncol=length(parameters))
      rownames(ConverTest) <- c("ESS")
      colnames(ConverTest) <- parameters
      
      
      

      for (k in 1:length(parameters)){
        
        ConverTest["ESS",parameters[k]] <- (ESS_check(dat, parameters[k]))
      
    
      }
    
      ConverTest <- as.data.frame(ConverTest)
     
      if (length(which(ConverTest == TRUE)) == length(ConverTest)) {
        write.csv(ConverTest, paste0(args[1], "/",args[2], "/" ,args[3],"/output_", args[2],"/Convergence.csv"))
      }
  
}

main()


