# ------------------------------------------------------------
# This R code will generate a very high-level summary of
# the distribution of the correlation count matrices for
# abs, negative, and positive.
# args[1] - correlation matrix filepath
# ------------------------------------------------------------
args <- commandArgs(TRUE)
setwd("/home/bmb191/math5376D/Final.Project/Output.Files/Plots")

# ------------------------------------------------------------
# Load in correlation count data and set diagonals to 0.
# ------------------------------------------------------------
# Absolute value
corr_cnts <- read.table(
	paste0("../Correlation.Counts/", args[1]), 
	header=FALSE)[,1]
diag(corr_cnts) <- 0

# Positive
pos_corr_cnts <- read.table(
	paste0("../Correlation.Counts/positive.", args[1]), 
	header=FALSE)[,1]
diag(pos_corr_cnts) <- 0

# Negative
neg_corr_cnts <- read.table(
	paste0("../Correlation.Counts/negative.", args[1]), 
	header=FALSE)[,1]
diag(neg_corr_cnts) <- 0

# ------------------------------------------------------------
# Output summary statistics for each count distribution
# ------------------------------------------------------------
summary(corr_cnts)
summary(pos_corr_cnts)
summary(neg_corr_cnts)

# ------------------------------------------------------------
# Output how many of the CG sites posses each low value
# ------------------------------------------------------------
for (i in 0:5) {
	print(paste("Count of correlation counts equal to",i,":",sum(corr_cnts==i)))
	print(paste("Count of positive correlation counts equal to",i,":",sum(pos_corr_cnts==i)))
	print(paste("Count of negative correlation counts equal to",i,":",sum(neg_corr_cnts==i)))
}

# ------------------------------------------------------------
# Plot the density of the highly correlated pair counts
# ------------------------------------------------------------
pdf(paste0(args[2],'.high.corr.cnts.density.pdf'))
plot(
	density((corr_cnts)^(1/4)),
	main="High Corr Counts Density"
)
dev.off()

# ------------------------------------------------------------
# Plot the density of the positive highly correlated pair 
# counts
# ------------------------------------------------------------
pdf(paste0(args[2],'.positive.high.corr.cnts.density.pdf'))
plot(
	density((pos_corr_cnts)^(1/4)),
	main="Positive High Corr Counts Density",
	xlim=c()
)
dev.off()

# ------------------------------------------------------------
# Plot the density of the negative highly correlated pair 
# counts
# ------------------------------------------------------------
pdf(paste0(args[2],'.negative.high.corr.cnts.density.pdf'))
plot(
	density((neg_corr_cnts)^(1/4)),
	main="Negative High Corr Counts Density"
)
dev.off()
