setwd("../../ur15efac/Downloads/Morpho/")
setwd("output-complete/")


pdf("Tree Lengths.pdf")
par(mfrow=c(2,2))
for (i in 1:length(datasets)){

toplt <- matrix(ncol=3, nrow=length(models))
colnames(toplt) <- c("model", "average", "sd")
rownames(toplt) <- models
for (j in 1:length(models)){
  container <- matrix(ncol = 1, nrow=4002)
  #colnames(container) <- models
  tr <- read.table(paste0("output/", datasets[i],"/" , models[j],
                          "/all-nexus_posterior.log"), header = TRUE)
  container[,1] <- tr$tree_length
  container <- as.data.frame(container)
  toplt[models[j], "model"]  <- models[j]
  toplt[models[j], "average"] <- as.numeric(mean(container$V1))
  toplt[models[j], "sd"] <- as.numeric(sd(container$V1))
}
toplt <- as.data.frame(toplt)
vec <- c("mk","mk+G","mk+V","mk+GV", "mk+Gmultistate", "mk+Vmultistate", "mk+GVmultistate")
toplt<- toplt[match(vec, toplt$model), ]

vecp <- c("Mk", "Mk+G", "Mkv", "Mkv+G"," Mkmul+G", "Mkvmul", "MkVmul+G")
## Variables for plotting
sdev <- as.numeric(toplt$sd)
avg <- as.numeric(toplt$average)
x <- 1:length(models)
s <-max(sdev)
a <- max(avg) + s
b <- min(avg) - s
l <- sum(avg) / length(models)


plot(as.numeric(toplt$average),x, pch =c(15,16,17,18,19,8,3), main=i,
     col= col,
     xlim = c(b,a), xlab = "Mean tree Length", ylab = "", yaxt = "n", )
arrows(avg-sdev,x, avg+sdev, x, length=0.05, angle=90, code=3)
points(as.numeric(toplt$average),x, pch =c(15,16,17,18,19,8,3), main="",
       col= RColorBrewer::brewer.pal(7, "Set2"))
axis(2, at=1:7, labels=vecp, cex.axis=0.8, las=2, lwd = 0, lwd.ticks = 1)


}
dev.off()










#### unknown 3
files <- list.files("output")
files <- mixedsort(sort(files))
install.packages("ape")
install.packages("phangorn")
library(ape)
#pdf("RF-distance.pdf)
  mk <- ape::read.tree(paste0("output/",files[i],"/mk/all-nexus_posterior.trees"))
  burnin <- length(mk) * 0.25
  size <-  length(mk)
  mk <- mk[burnin:size]
  
  mk_G <- ape::read.tree(paste0("output/",files[i],"/mk+G/all-nexus_posterior.trees"))
  mk_G <- mk_G[burnin:size]
  
  mk_V <- ape::read.tree(paste0("output/",files[i],"/mk+V/all-nexus_posterior.trees"))
  mk_V <- mk_V[burnin:size]
  
  mk_GV <- ape::read.tree(paste0("output/",files[i],"/mk+GV/all-nexus_posterior.trees"))
  mk_GV <- mk_GV[burnin:size]
  
  mk_Gm <- ape::read.tree(paste0("output/",files[i],"/mk+Gmultistate/all-nexus_posterior.trees"))
  mk_Gm <- mk_Gm[burnin:size]
  
  mk_Vm <- ape::read.tree(paste0("output/",files[i],"/mk+Vmultistate/all-nexus_posterior.trees"))
  mk_Vm <- mk_Vm[burnin:size]
  
  mk_GVm <- ape::read.tree(paste0("output/",files[i],"/mk+Gvmultistate//all-nexus_posterior.trees"))
  mk_GVm <- mk_GVm[burnin:size]
#par(mfrow=c(2,2))
## read in trees
i=85
#for (i in 1:2){
  #message(paste0("Working on file ", i, " out of ", length(files)))

  ### sample from posterior
  SAMPLE <- function(data){
    data[sample(length(data),1000)]
  }
  mk.sample <- SAMPLE(mk)
  mk_G.sample <- SAMPLE(mk_G)
  mk_V.sample <- SAMPLE(mk_V)
  mk_GV.sample <- SAMPLE(mk_GV)
  mk_Gm.sample <- SAMPLE(mk_Gm)
  mk_Vm.sample <- SAMPLE(mk_Vm)
  mk_GVm.sample <- SAMPLE(mk_GVm)
  #Prepare for graphic representation
  col <- c(viridis::viridis(7))#c(viridis::viridis(5))
  pch <- c(0,1,2,3,4,5,6)#c(16,1,15,17,17)
  legend <- c("mk","mk_G","mk_V","mk_GV", "mk_Gm", "mk_Vm", "mk_GVm")
  group<-rep(c(1,2,3,4,5,6,7),each=300)
  #Robinson-Foulds ####
  #RF distances ####
  Results <- function(df1,df2,ROWNAMES,COLNAMES){
    RF <- data.frame(matrix(ncol=length(df2),nrow=length(df1)))
    for (i in 1:length(df1)){
      for (j in 1:length(df2)){
        RF[i,j] <- phangorn::RF.dist(ape::unroot(df1[[i]]),ape::unroot(df2[[j]]),
                                     normalize=TRUE, check.labels=TRUE,rooted=FALSE)
        rownames(RF)[i] <- paste(ROWNAMES,i,sep="_")
        colnames(RF)[j] <- paste(COLNAMES,j,sep="_")
      }
    }
    return(RF)
  }
  RF1 <- Results(mk.sample,mk.sample,"mk","mk")
  RF2 <- Results(mk.sample,mk_G.sample,"mk","mk_G")
  RF3 <- Results(mk.sample,mk_V.sample,"mk","mk_V")
  RF4 <- Results(mk.sample,mk_GV.sample,"mk","mk_GV")
  RF5 <- Results(mk.sample,mk_Gm.sample,"mk","mk_Gm")
  RF6 <- Results(mk.sample,mk_Vm.sample,"mk","mk_Vm")
  RF7 <- Results(mk.sample,mk_GVm.sample,"mk","mk_GVm")
  RF8 <- Results(mk_G.sample,mk_G.sample,"mk_G","mk_G")
  RF9 <- Results(mk_G.sample,mk_V.sample,"mk_G","mk_V")
  RF10 <- Results(mk_G.sample,mk_GV.sample,"mk_G","mk_GV")
  RF11 <- Results(mk_G.sample,mk_GV.sample,"mk_G","mk_Gm")
  RF12 <- Results(mk_G.sample,mk_Vm.sample,"mk_G","mk_Vm")
  RF13 <- Results(mk_G.sample,mk_GVm.sample,"mk_G","mk_GVm")
  RF14 <- Results(mk_V.sample,mk_V.sample,"mk_V","mk_V")
  RF15 <- Results(mk_V.sample,mk_GV.sample,"mk_V","mk_GV")
  RF16 <- Results(mk_V.sample,mk_Gm.sample,"mk_V","mk_Gm")
  RF17 <- Results(mk_V.sample,mk_Vm.sample,"mk_V","mk_Vm")
  RF18 <- Results(mk_V.sample,mk_GVm.sample,"mk_V","mk_GVm")
  RF19 <- Results(mk_GV.sample,mk_GV.sample,"mk_GV","mk_GV")
  RF20 <- Results(mk_GV.sample,mk_Gm.sample,"mk_GV","mk_Gm")
  RF21 <- Results(mk_GV.sample,mk_Vm.sample,"mk_GV","mk_Vm")
  RF22 <- Results(mk_GV.sample,mk_GVm.sample,"mk_GV","mk_GVm")
  RF23 <- Results(mk_Gm.sample,mk_Gm.sample,"mk_Gm","mk_Gm")
  RF24 <- Results(mk_Gm.sample,mk_Vm.sample,"mk_Gm","mk_Vm")
  RF25 <- Results(mk_Gm.sample,mk_GVm.sample,"mk_Gm","mk_GVm")
  RF26 <- Results(mk_Vm.sample,mk_Vm.sample,"mk_Vm","mk_Vm")
  RF27 <- Results(mk_Vm.sample,mk_GVm.sample,"mk_Vm","mk_GVm")
  RF28 <- Results(mk_GVm.sample,mk_GVm.sample,"mk_GVm","mk_GVm")
  
  
  
  
  
  
  


