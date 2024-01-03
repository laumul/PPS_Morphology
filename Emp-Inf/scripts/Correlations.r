#### plot the variance in tree length across models with num characters

information <- read.csv2("results.csv")[1:114,]


Variance <- matrix(ncol =5, nrow=length(information$Name))
colnames(Variance) <- c("variance", "numChar", "numTaxa", "numState", "mix")
rownames(Variance) <- datasets


for (i in 1:length(datasets)){
v <- c()
mk <- ape::read.nexus(paste0("output/",datasets[i],"/mk/MCC.tre"))
v <- rbind(v, sum(mk$edge.length))
mkG <- ape::read.nexus(paste0("output/",datasets[i],"/mk+G/MCC.tre"))
v <- rbind(v, sum(mkG$edge.length))
mkV <- ape::read.nexus(paste0("output/",datasets[i],"/mk+V/MCC.tre"))
v <- rbind(v, sum(mkV$edge.length))
mkGV <- ape::read.nexus(paste0("output/",datasets[i],"/mk+GV/MCC.tre"))
v <- rbind(v, sum(mkGV$edge.length))
mkPG <- ape::read.nexus(paste0("output/",datasets[i],"/mk+Gmultistate/MCC.tre"))
v <- rbind(v, sum(mkPG$edge.length))
mkVP <- ape::read.nexus(paste0("output/",datasets[i],"/mk+Vmultistate/MCC.tre"))
v <- rbind(v, sum(mkVP$edge.length))
mkVPG <- ape::read.nexus(paste0("output/",datasets[i],"/mk+GVmultistate/MCC.tre"))
v <- rbind(v, sum(mkVPG$edge.length))

Variance[datasets[i],"variance"] <-  max(v) - min(v)
Variance[datasets[i], "numChar"] <- information[i,"Num.Characters"]
Variance[datasets[i], "numTaxa"] <- information[i,"Num.taxa"]
Variance[datasets[i], "numState"] <- information[i,"Num.states"]
Variance[datasets[i], "mix"] <- (information[i,"Num.states"] * information[i,"Num.Characters"]) / information[i,"Num.taxa"]


}
par(mfrow=c(1,1))

  #Num char
Variance <- as.data.frame(Variance)
plot(Variance$numChar ,Variance$variance, xlab = "Number of Characters", ylab = "Variance in tree length across models")
fit <- lm(variance ~ numChar, data = Variance)
abline(fit)
summary_fit <- summary(fit)
text(75, 11, paste0("R = ",round(summary_fit$r.squared,4)),  cex =0.7)
text(75, 10, paste0("Pvalue = ",round(summary_fit$coefficients[8],4)),  cex =0.7)  

# numTaxa
plot(Variance$numTaxa ,Variance$variance, xlab = "Number of Taxa", ylab = "Variance in tree length across models")
fit <- lm(variance ~ numTaxa, data = Variance)


abline(fit)
summary_fit <- summary(fit)
text(20, 11, paste0("R = ",round(summary_fit$r.squared,4)), cex =0.7)
text(20, 10, paste0("Pvalue = ",round(summary_fit$coefficients[8],4)),  cex =0.7)     
   

## Num State
plot(Variance$numState ,Variance$variance, xlab = "Number of States", ylab = "Variance in tree length across models")
fit <- lm(variance ~ numState, data = Variance)

abline(fit)
summary_fit <- summary(fit)
text(3, 11, paste0("R = ",round(summary_fit$r.squared,4)), cex =0.7) 
text(3, 10, paste0("Pvalue = ",round(summary_fit$coefficients[8],4)),  cex =0.7) 
    

# mix

Variance <- as.data.frame(Variance)
plot(Variance$mix ,Variance$variance, xlab = "Number of Taxa", ylab = "Variance in tree length across models")
fit <- lm(variance ~ mix, data = Variance)

#plot(Variance$variance ,Variance$numChar,Variance$variance, xlab = "Number of Characters", ylab = "Variance in tree length across models")
#fit <- lm(numChar~ variance, data = Variance)
abline(fit)
summary_fit <- summary(fit)
round(summary_fit$r.squared,4)
text(20, 11, paste0("R = ",round(summary_fit$r.squared,4)))
     
     
     