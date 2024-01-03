
#setwd("../../Final/Figures/Figure-3/")
main <- function(){
  args <- commandArgs(trailingOnly = TRUE)

library(dplyr)
library(phangorn)
library(vegan)
library(gtools)
library(RColorBrewer)


files <- mixedsort(sort(list.files("output")))


i=as.numeric(args[1])
#i=2
print(files[i]) 
nm <- gsub(".*_", "", files[i])
#par(mfrow=c(2,2))


#message(paste0("Working on file ", i, " out of ", length(files))) 
mk <- ape::read.tree(paste0("output/",files[i],"/mk/all-nexus_posterior.trees"))
b <-(length(mk) -2 ) *0.1 
mk <- mk[b:length(mk)]

mk_G <- ape::read.tree(paste0("output/",files[i],"/mk+G/all-nexus_posterior.trees")) 
b <-(length(mk_G) -2 ) *0.1
mk_G <- mk_G[b:length(mk_G)]

mk_V <- ape::read.tree(paste0("output/",files[i],"/mk+V/all-nexus_posterior.trees"))
b <-(length(mk_V) -2 ) *0.1
mk_V <- mk_V[b:length(mk_V)]

mk_GV <- ape::read.tree(paste0("output/",files[i],"/mk+GV/all-nexus_posterior.trees"))
b <-(length(mk_GV) -2 ) *0.1
mk_GV <- mk_GV[b:length(mk_GV)]

mk_Gm <- ape::read.tree(paste0("output/",files[i],"/mk+Gmultistate/all-nexus_posterior.trees"))
b <-(length(mk_Gm) -2 ) *0.1
mk_Gm <- mk_Gm[b:length(mk_Gm)]

mk_Vm <- ape::read.tree(paste0("output/",files[i],"/mk+Vmultistate/all-nexus_posterior.trees"))
b <-(length(mk_Vm) -2 ) *0.1
mk_Vm <- mk_Vm[b:length(mk_Vm)]

mk_GVm <- ape::read.tree(paste0("output/",files[i],"/mk+GVmultistate//all-nexus_posterior.trees"))
b <-(length(mk_GVm) -2 ) *0.1
mk_GVm <- mk_GVm[b:length(mk_GVm)]

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


#Graphics helpers
col <- brewer.pal(7, "Set2")
#pch <- c(0,1,2,3,4,5,6)#c(16,1,15,17,17)
pch=c(15,4,17,18,19,8,3)
legend <- c("mk","mk_G","mk_V","mk_GV", "mk_Gm", "mk_Vm", "mk_GVm")
group<-rep(c(1,2,3,4,5,6,7),each=1000)


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




join.tot <- cbind(rbind(RF1,t(RF2),t(RF3), t(RF4), t(RF5), t(RF6), t(RF7)),
                  rbind(RF2,t(RF8),t(RF9),t(RF10),t(RF11),t(RF12), t(RF13)),
                  rbind(t(RF3),t(RF9),t(RF14),t(RF15),t(RF16),t(RF17), t(RF18)),
                  rbind(t(RF4),t(RF10),t(RF15),t(RF19),t(RF20),t(RF21), t(RF22)),
                  rbind(t(RF5),t(RF11),t(RF16),t(RF20),t(RF23),t(RF24), t(RF25)),
                  rbind(t(RF6),t(RF12),t(RF17),t(RF21),t(RF24),t(RF26), t(RF27)),
                  rbind(t(RF7),t(RF13),t(RF18),t(RF22),t(RF25),t(RF27), t(RF28)))





# An object with distance information to be converted to a "dist" object  
dist.join <- stats::as.dist(join.tot, diag=TRUE)

join.group <- as.factor(rep(c("mk","mk_G","mk_V","mk_GV", "mk_Gm", "mk_Vm", "mk_GVm"),
                            each=1000))

beta.test <- vegan::betadisper(dist.join, join.group)

stats::TukeyHSD(beta.test, which = "group", ordered = FALSE, conf.level = 0.95)
permtest <- vegan::permutest(beta.test, pairwise = TRUE)


pdf(paste0("output/",files[i],"/", nm,".pdf"))
plot(beta.test, axes = c(1,2), cex = 0.5, pch=pch, col=col,
     lty = 1, lwd = 2, hull = TRUE,  ellipse = FALSE,
     segments = FALSE, label = FALSE, main = files[i],sub="",
     xlab =paste0("PCoA 1= ",round(beta.test$eig[1]/ sum(beta.test$eig)*100,2), "%"),
     ylab= paste0("PCoA 1= ",round(beta.test$eig[2]/ sum(beta.test$eig)*100,2), "%"))
legend("bottomleft", legend=legend, col=col, pch=pch, cex =0.8)

dev.off()

saveRDS(beta.test, file= paste0("output/",files[i],"/", nm,".RData"))

write.table(permtest$pairwise$permuted, paste0("output/",files[i],"/","permut_",nm,".txt"))


}

main()



