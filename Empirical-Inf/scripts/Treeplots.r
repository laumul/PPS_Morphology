nstall.packages("gtools")
library(gtools)

datasets <- mixedsort(sort(list.files("output")))
models <- list.files("output/1_Agnolin/")
#models <- models[c(1,2,3,4,5,6,7)]
pch=c(15,16,17,18,19,8,3)
col <- RColorBrewer::brewer.pal(7, "Set2")

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
  
  
  plot(as.numeric(toplt$average),x, pch =c(15,16,17,18,19,8,3), main=paste0(datasets[i]),
       col= col,
       xlim = c(b,a), xlab = "Mean tree Length", ylab = "", yaxt = "n", )
  arrows(avg-sdev,x, avg+sdev, x, length=0.05, angle=90, code=3)
  points(as.numeric(toplt$average),x, pch =c(15,16,17,18,19,8,3), main="",
         col= RColorBrewer::brewer.pal(7, "Set2"))
  axis(2, at=1:7, labels=vecp, cex.axis=0.8, las=2, lwd = 0, lwd.ticks = 1)
  
  
}
dev.off()