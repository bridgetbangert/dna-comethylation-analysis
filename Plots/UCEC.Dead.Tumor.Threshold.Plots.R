setwd("/home/bmb191/math5376D/Final.Project/Output.Files")

# Compare the thresholds
UCEC_dt_counts_0.8 <- unlist(read.table("Correlation.Counts/UCEC.dead.tumor.high.corr.cnts.0.8.threshold.txt", header=FALSE))
UCEC_dt_counts_0.9 <- unlist(read.table("Correlation.Counts/UCEC.dead.tumor.high.corr.cnts.0.9.threshold.txt", header=FALSE))
UCEC_dt_distances_0.8 <- unlist(read.table("Correlation.Distances/UCEC.dead.tumor.high.corr.distances.threshold.0.8.txt", header=FALSE))
UCEC_dt_distances_0.9 <- unlist(read.table("Correlation.Distances/UCEC.dead.tumor.high.corr.distances.threshold.0.9.txt", header=FALSE))

UCEC_dn_counts_0.8 <- unlist(read.table("Correlation.Counts/UCEC.dead.normal.high.corr.cnts.0.8.threshold.txt", header=FALSE))
UCEC_dn_counts_0.9 <- unlist(read.table("Correlation.Counts/UCEC.dead.normal.high.corr.cnts.0.9.threshold.txt", header=FALSE))
UCEC_dn_distances_0.8 <- unlist(read.table("Correlation.Distances/UCEC.dead.normal.high.corr.distances.threshold.0.8.txt", header=FALSE))
UCEC_dn_distances_0.9 <- unlist(read.table("Correlation.Distances/UCEC.dead.normal.high.corr.distances.threshold.0.9.txt", header=FALSE))

pdf("Plots/UCEC.dead.threshold.comparison.counts.pdf",16,16)
par(mfrow=c(1,2))
plot(
  density(UCEC_dn_counts_0.8^(1/4)),
  col="blue", 
  xlim=c(1,7),
  ylim=c(0,1.1),
  ylab=NA,
  xlab="0.9 high correlation threshold",
  main=NA
)
legend("topright",legend=c("Normal","Tumor"),col=c("blue","red"),lty=1,lwd=2,cex=0.75,bty="n")
lines(density(UCEC_dt_counts_0.8^(1/4)), col="red")

plot(
  density(UCEC_dn_counts_0.9^(1/4)),
  col="blue", 
  xlim=c(1,7),
  ylim=c(0,1.1),
  ylab=NA,
  xlab="0.8 high correlation threshold",
  main=NA
)
legend("topright",legend=c("Normal","Tumor"),col=c("blue","red"),lty=1,lwd=2,cex=0.75,bty="n")
lines(density(UCEC_dt_counts_0.9^(1/4)), col="red")

title("UCEC Dead Data High Correlation Threshold Counts Comparison", outer=TRUE, line=-2)
dev.off()

pdf("Plots/UCEC.dead.threshold.comparison.distances.pdf",16,16)
par(mfrow=c(1,2))
plot(
  density(UCEC_dn_distances_0.8),
  col="blue", 
  #xlim=c(1,7),
  ylim=c(0,2.25*10^(-8)),
  ylab=NA,
  xlab="0.9 high correlation threshold",
  main=NA
)
legend("topright",legend=c("Normal","Tumor"),col=c("blue","red"),lty=1,lwd=2,cex=0.75,bty="n")
lines(density(UCEC_dt_distances_0.8), col="red")

plot(
  density(UCEC_dn_distances_0.9),
  col="blue", 
  #xlim=c(1,7),
  ylim=c(0,2.25*10^(-8)),
  ylab=NA,
  xlab="0.8 high correlation threshold",
  main=NA
)
legend("topright",legend=c("Normal","Tumor"),col=c("blue","red"),lty=1,lwd=2,cex=0.75,bty="n")
lines(density(UCEC_dt_distances_0.9), col="red")

title("UCEC Dead Data High Correlation Threshold Distances Comparison", outer=TRUE, line=-2)
dev.off()

summary(UCEC_dn_counts_0.9)
summary(UCEC_dn_counts_0.8)
summary(UCEC_dt_counts_0.9)
summary(UCEC_dt_counts_0.8)

summary(UCEC_dn_distances_0.9)
summary(UCEC_dn_distances_0.8)
summary(UCEC_dt_distances_0.9)
summary(UCEC_dt_distances_0.8)
