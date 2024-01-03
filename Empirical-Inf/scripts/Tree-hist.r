setwd("Desktop/PhD/Posterior-Prediction/Final/")
install.packages("ggridges")
library(ggridges)

library(ggplot2)
Files <- list.files("output/")
Results <- matrix(nrow=length(Files), ncol= 6)
mods <- list.files("output/1_Agnolin/")
mods <- mods[c(2,3,4,5,6,7)]
colnames(Results) <- mods
for (i in 1:length(Files)){
  MkTree <-read.table(paste0("output/",Files[i],"/mk/all-nexus_posterior.log"), header = TRUE)
  MkMean <- mean(MkTree$tree_length)
  for (k in 1:length(mods)){
    CompTree <- read.table(paste0("output/",Files[i],"/", mods[k],"/all-nexus_posterior.log"), header = TRUE)
    CompTreeMean <- mean(CompTree$tree_length)
    CompDiff <- CompTreeMean - MkMean
    CompPer <- CompDiff/ MkMean
    Results[i,k] <- CompPer * 100
  }
}
Results <- as.data.frame(Results)
pdf("tree-lengths.pdf")
par(mfrow=c(1,1))
for (l in 1:6){
  hist(Results[,l], main = mods[l], xlim = c(-0.5,0.4))
}
dev.off()
library(ggplot2)
Results<- Results[, c(1, 5, 3, 2, 6,4)]
ResultsGg <- stack(as.data.frame(Results))


#white background
pdf("ridgeplots.pdf")
ggplot(ResultsGg, aes(x=values, y = ind, fill= ind)) +
  geom_density_ridges() +
  theme_ridges() +
  theme(legend.position = "none") +
  geom_vline(xintercept = 0, colour= "black", linetype="dashed", size=1)
dev.off()



##### calulcate differences 

counter=0

for (i in 1:length(Results$`mk+G`)){
if(Results$`mk+V`[i] < 0 ){
counter=counter+1
}
}



counter=0

for (i in 1:length(Results$`mk+G`)){
  if(Results$`mk+G`[i] > 0 ){
    counter=counter+1
  }
}









