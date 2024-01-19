#####
# Simulated 1000 data sets for PPS

main <- function(){
  args <- commandArgs(trailingOnly = TRUE)
  
  dir=paste0(args[2], "/", args[3], "/output/")
  dir.create(paste0(dir,"/PPS_Morpho_post_sims"))
  
  sim.morpho <- function(tree, k = 2, nchar = c(10), variable.coding = FALSE, 
                         partitions = 1, alpha = NULL){
    
    if(is.null(alpha)) rates = 1
    else rates = phangorn::discrete.gamma(alpha, 4) # TODO: NEED TO CHECK WITH BEN
    
    for(p in partitions:1){ # for now, goes backwards from p max. Otherwise we have to deal with attr(seq, "levels"); or you can use cbind(x, seq)
      if(max(partitions) == 1) states = as.character(c(0:(k-1)))
      else states = as.character(c(0:p))
      
      for(i in 1:nchar[p]){
        if(nchar[p] == 0) next
        x = phangorn::simSeq(tree, l = 1, type = "USER", levels = states, rate = sample(rates, 1))
        
        if(variable.coding){
          while( length(unique(as.character(x))) == 1 ){
            x = phangorn::simSeq(tree, l = 1, type = "USER", levels = states, rate = sample(rates, 1))
          }
        }
        
        if(i == 1 && p == max(partitions)){
          seq = x
        }
        else seq = cbind(seq, x)
      }
    }
    seq
  }
  
  
  
  
  dat <- t(data.frame(ape::read.nexus.data(paste0("data/",args[1]))))

 nstate <- max(unique(na.omit(as.numeric(dat)))) +1  
  v=c()
  for (i in 1:length(dat[1,])){
    v=c(v,max(unique(na.omit(as.numeric(dat[,i])))))
    
  }
  
    part <- as.vector(table(v))
    
     char= c()
    states <- seq(0, (nstate-1), 1)
      for (k in 1:length(states)){
        
        char = c(char,length(which(v == states[k])))
      }
      
#    char = c(1,2,3,4)
    s1 <- char[1]+char[2]
    char = c(s1,char[3:nstate])
      
if (char[1] == length(v)){ 
  char = char[1]
  states= 1
 nstate=2 
}
         
    tree <- ape::read.tree(paste0(args[2],"/",args[3], "/output/PPS_Morpho_posterior.var"))
    tree <- ape::root.multiPhylo (tree, tree[[1]]$tip.label[4], resolve.root = TRUE)

    
    
    InferenceLog <- read.table(paste0(args[2],"/",args[3], "/output/PPS_Morpho_posterior.var"), header=TRUE) 
    
    
    if (args[3] == "mk"){
      
      for (i in 1:500){
        dir.create(paste0(dir,"/PPS_Morpho_post_sims/posterior_predictive_sim_", i ))
        #alpha_morpho <- 1/InferenceLog$alpha_inv[i] 
        t <- sim.morpho(tree[[i]], k=nstate, nchar = length(v), partitions = 1, variable.coding = FALSE,
                        alpha = NULL)
        
        phangorn::write.phyDat(t, paste0(dir,"/PPS_Morpho_post_sims/posterior_predictive_sim_", i,"/seq.nex"), format = "nexus" )
           
      }
    } else if (args[3] == "mk+G"){
      
      for (i in 1:500){
        dir.create(paste0(dir,"/PPS_Morpho_post_sims/posterior_predictive_sim_", i ))
        alpha_morpho <- 1/InferenceLog$alpha_inv[i] 
        t <- sim.morpho(tree[[i]], k=nstate, nchar =length(v), partitions = 1 , variable.coding = FALSE,
                        alpha = alpha_morpho)
        
        phangorn::write.phyDat(t, paste0(dir,"/PPS_Morpho_post_sims/posterior_predictive_sim_", i,"/seq.nex"), format = "nexus" )
        
      } 
    } else if (args[3] == "mk+GV"){
      
      for (i in 1:500){
        dir.create(paste0(dir,"/PPS_Morpho_post_sims/posterior_predictive_sim_", i ))
        alpha_morpho <- 1/InferenceLog$alpha_inv[i] 
        t <- sim.morpho(tree[[i]], k=nstate, nchar = length(v), partitions = 1, variable.coding = TRUE,
                        alpha = alpha_morpho)
        
        phangorn::write.phyDat(t, paste0(dir,"/PPS_Morpho_post_sims/posterior_predictive_sim_", i,"/seq.nex"), format = "nexus" )
        
      } 
    } else if (args[3] == "mk+V"){
      
      for (i in 1:500){
        dir.create(paste0(dir,"/PPS_Morpho_post_sims/posterior_predictive_sim_", i ))
        #alpha_morpho <- 1/InferenceLog$alpha_inv[i] 
        t <- sim.morpho(tree[[i]], k=nstate, nchar = length(v), partitions = 1, variable.coding = TRUE,
                        alpha = NULL)
        
        phangorn::write.phyDat(t, paste0(dir,"/PPS_Morpho_post_sims/posterior_predictive_sim_", i,"/seq.nex"), format = "nexus" )
        
      } 
    } else if (args[3] == "mk+Vmultistate"){
      
      for (i in 1:500){
        dir.create(paste0(dir,"/PPS_Morpho_post_sims/posterior_predictive_sim_", i ))
        #alpha_morpho <- 1/InferenceLog$alpha_inv[i]
	rm(Morpho_seq) 
    numbers = 2
	m=1 
  for ( s in 1:length(char)){
      if (char[s] == 0) {
  numbers=numbers+1
  }else{
   t <- sim.morpho(tree[[i]], k = numbers, nchar = char[s], partitions = 1, variable.coding = TRUE,
                      alpha = NULL)
      phangorn::write.phyDat(t, paste0(dir,"/PPS_Morpho_post_sims/posterior_predictive_sim_", i, "/", "m_morph[",m,"].nex"), format = "nexus") 
      numbers=numbers+1
	m=m+1
       if (! exists("Morpho_seq")) {
        Morpho_seq <- t
      } else {
        Morpho_seq <- cbind(t, Morpho_seq)
      }
    }
}
     phangorn::write.phyDat(Morpho_seq, paste0(dir,"/PPS_Morpho_post_sims/posterior_predictive_sim_", i,"/seq.nex"),format ="nexus" )

      }
    } else if (args[3] == "mk+GVmultistate"){
      
      for (i in 1:500){
        dir.create(paste0(dir,"/PPS_Morpho_post_sims/posterior_predictive_sim_", i ))
        alpha_morpho <- 1/InferenceLog$alpha_inv[i] 
        #t <- sim.morpho(tree[[i]], k=nstate, nchar = char, partitions = length(char), variable.coding = TRUE,
         #               alpha = alpha_morpho)
        
       # phangorn::write.phyDat(t, paste0(dir,"/PPS_Morpho_post_sims/posterior_predictive_sim_", i,"/seq.nex"), format = "nexus" )
      rm(Morpho_seq)
	 numbers = 2 
	m=1
	for ( s in 1:length(char)){
      if (char[s] == 0) {
	numbers=numbers+1
	}else{
	 t <- sim.morpho(tree[[i]], k = numbers, nchar = char[s], partitions = 1, variable.coding = TRUE,
                      alpha = alpha_morpho)
      phangorn::write.phyDat(t, paste0(dir,"/PPS_Morpho_post_sims/posterior_predictive_sim_", i, "/", "m_morph[",m,"].nex"), format ="nexus") 
      numbers=numbers+1
	m=m+1
       if (! exists("Morpho_seq")) {
        Morpho_seq <- t
      } else {
        Morpho_seq <- cbind(t, Morpho_seq)
      }
    }
}
     phangorn::write.phyDat(Morpho_seq, paste0(dir,"/PPS_Morpho_post_sims/posterior_predictive_sim_", i,"/seq.nex"), format ="nexus" )

      } 
    } else if (args[3] == "mk+Gmultistate"){
      
      for (i in 1:500){
        dir.create(paste0(dir,"/PPS_Morpho_post_sims/posterior_predictive_sim_", i ))
        alpha_morpho <- 1/InferenceLog$alpha_inv[i]
	rm(Morpho_seq) 
    numbers = 2 
	m=1
  for ( s in 1:length(char)){
      if (char[s] == 0) {
  numbers=numbers+1
  }else{
   t <- sim.morpho(tree[[i]], k = numbers, nchar = char[s], partitions = 1, variable.coding = FALSE,
                      alpha = alpha_morpho)
      phangorn::write.phyDat(t, paste0(dir,"/PPS_Morpho_post_sims/posterior_predictive_sim_", i, "/", "m_morph[",m,"].nex"), format = "nexus") 
      numbers=numbers+1
	m=m+1
       if (! exists("Morpho_seq")) {
        Morpho_seq <- t
      } else {
        Morpho_seq <- cbind(t,Morpho_seq)
      }
    }
}
     phangorn::write.phyDat(Morpho_seq, paste0(dir,"/PPS_Morpho_post_sims/posterior_predictive_sim_", i,"/seq.nex"), format = "nexus" )

      }      
    }
}
    
  
  
  main()
  
