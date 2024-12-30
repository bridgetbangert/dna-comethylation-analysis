# ------------------------------------------------------------
# This R code will generate the scatterplots comparing
# the original to the data with noise.
# ------------------------------------------------------------
# args[1] - Original data filename
# args[2] - Noise data filename
# args[3] - Output filename
# args[4] - Cancer type
# args[5] - Alive/Dead
# args[6] - Normal/Tumor
# ------------------------------------------------------------
args <- commandArgs(TRUE)
setwd("/home/bmb191/math5376D/Final.Project/Output.Files")

orig_data <- as.numeric(unlist(read.table(paste0("Correlation.Counts/",args[1]),header=FALSE)))
noise_data <- as.numeric(unlist(read.table(paste0("Correlation.Counts/",args[2]),header=FALSE)))

pdf(paste0("Plots/",args[3],".pdf"))
plot(
	orig_data,
	noise_data,
	xlab="Original Data",
	ylab="Data with Noise",
	main=paste(args[4],args[5],args[6])
)
dev.off()