#setwd("Desktop/PhD/Posterior-Prediction/Final/Simulations/")


######
# partition based on empirical data 


sim.morpho = function(tree, k = 2, nchar = c(10), variable.coding = FALSE, partitions = "PPS", alpha = alpha_morpho){
  # Mk
  # + V
  # + P
  # + G
  
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

model="mk+GV"
EmpData = "Egi_etal_2005a_paleobiodb.nex"
#EmpData = "Shoshani_etal_2006a_paleobiodb.nex"

runID ="1_Egi"
#setwd("../morphosim/")
dat <- t(data.frame(ape::read.nexus.data(paste0("data/",EmpData))))
nstate <- max(unique(na.omit(as.numeric(dat)))) +1


v=c()
for (i in 1:length(dat[1,])){
  v=c(v,max(unique(na.omit(as.numeric(dat[,i])))))
  
}

char= c()
states <- seq(0, nstate-1, 1)
for (k in 1:length(states)){
  
  char = c(char,length(which(v == states[k])))
}

#    char = c(1,2,3,4)
s1 <- char[1]+char[2]
char = c(s1,char[3:(nstate)])

if (char[1] == length(v)){ 
  char = char[1]
  states= 1
  nstate=2 
}

## read in var trees
tree <- ape::read.tree(paste0(model,"/",runID,"/output_",model,"/pps_Nov_posterior.var"))
tree <- ape::root.multiPhylo (tree, tree[[1]]$tip.label[4], resolve.root = TRUE)


InferenceLog <- read.table(paste0(model,"/",runID,"/output_",model,"/pps_Nov_posterior.var"), header=TRUE) 



dir.create(paste0("data_",model))

if (model == "mk+GVmultistate") {
  for (i in 1:20) {
    alpha_morpho <- 1 / InferenceLog$alpha_inv[i]
    t <- sim.morpho(tree[[i]], k = nstate, nchar = char, partitions = length(char), variable.coding = TRUE,
                    alpha = alpha_morpho)
    phangorn::write.phyDat(t, paste0("data_", model, "/", i, ".nex"), format = "nexus")
    
  }
} else if (model == "mk+GV") {
  for (i in 1:500) {
    alpha_morpho <- 1 / InferenceLog$alpha_inv[i]
    t <- sim.morpho(tree[[i]], k = nstate, nchar = length(v), partitions = 1, variable.coding = TRUE,
                    alpha = alpha_morpho)
    phangorn::write.phyDat(t, paste0("data_", model, "/", i, ".nex"), format = "nexus")
  }
} else if (model == "mk+Vmultistate") {
  for (i in 1:20) {
    alpha_morpho <- 1 / InferenceLog$alpha_inv[i]
    t <- sim.morpho(tree[[i]], k = nstate, nchar = char, partitions = length(char), variable.coding = TRUE,
                    alpha = NULL)
    phangorn::write.phyDat(t, paste0("data_", model, "/", i, ".nex"), format = "nexus")
  }
} else if (model == "mk+V") {
  for (i in 1:20) {
    alpha_morpho <- 1 / InferenceLog$alpha_inv[i]
    t <- sim.morpho(tree[[i]], k = nstate, nchar = length(v), partitions = 1, variable.coding = TRUE,
                    alpha = NULL)
    phangorn::write.phyDat(t, paste0("data_", model, "/", i, ".nex"), format = "nexus")
  }
}
