# ------------------------------------------------------------
# This R file will generate the plots of the 
# the start points in the X chromosome
#
# arg[1] - chrX data filename
# arg[2] - output file name (will be stored in Output.Files/
# Distance.Matrices/ directory)
# ------------------------------------------------------------
setwd("/home/bmb191/math5376D/Final.Project/Output.Files")

UCEC_dt_counts_0.8 <- unlist(read.table("Correlation.Counts/UCEC.dead.tumor.high.corr.cnts.0.8.threshold.txt", header=FALSE))
UCEC_dt_counts_0.9 <- unlist(read.table("Correlation.Counts/UCEC.dead.tumor.high.corr.cnts.0.9.threshold.txt", header=FALSE))
UCEC_dt_distances_0.8 <- unlist(read.table("Correlation.Distances/UCEC.dead.tumor.high.corr.distances.threshold.0.8.txt", header=FALSE))
UCEC_dt_distances_0.9 <- unlist(read.table("Correlation.Distances/UCEC.dead.tumor.high.corr.distances.threshold.0.9.txt", header=FALSE))

pdf("Plots/UCEC.dead.tumor.threshold.comparison.counts.pdf",16,16)
plot(density(UCEC_dt_counts_0.9),col="blue", main="UCEC high correlation count densities - comparing high correlation thresholds")
legend("topright",legend=c("0.9 threshold","0.8 threshold"),col=c("blue","red"),lty=1,lwd=2)
lines(density(UCEC_dt_counts_0.8), col="red")
dev.off()

pdf("Plots/UCEC.dead.tumor.threshold.comparison.distances.pdf",16,16)
plot(density(UCEC_dt_distances_0.9),col="blue", main="UCEC high correlation distance densities - comparing high correlation thresholds")
legend("topright",legend=c("0.9 threshold","0.8 threshold"),col=c("blue","red"),lty=1,lwd=2)
lines(density(UCEC_dt_distances_0.8), col="red")
dev.off()

summary(UCEC_dt_counts_0.9)
summary(UCEC_dt_counts_0.8)

summary(UCEC_dt_distances_0.9)
summary(UCEC_dt_distances_0.8)
