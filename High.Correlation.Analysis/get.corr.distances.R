# ------------------------------------------------------------
# This R file will grab all of the distinct distances of the
# highly correlated CG sites as a vector
#
# arg[1] - correlation matrix file
# arg[2] - distance matrix file
# arg[3] - output file name (will be save in the Output.Files/
# directory)
# ------------------------------------------------------------
args <- commandArgs(TRUE)
setwd("/home/bmb191/math5376D/Final.Project/")

# load in data
cg_corr <- read.table(
	paste0("Output.Files/Correlation.Matrices/",args[1]), 
	header=FALSE
)
cd_dist <- read.table(
	paste0("Output.Files/Distance.Matrices/",args[2]), 
	header=FALSE
)

# ------------------------------------------------------------
# Get the indexes of the locations with highly correlated
# CG sites.
# Create a vector to store the array of the distances of
# highly correlated CG pairs
# ------------------------------------------------------------
high_corr_indices <- which(abs(cg_corr) >= 0.8, arr.ind = TRUE)
high_corr_indices <- high_corr_indices[high_corr_indices[, 1] >= high_corr_indices[, 2], ]
high_corr_distances <- cd_dist[high_corr_indices]

# ------------------------------------------------------------
# Write result to a file in the leap server
# ------------------------------------------------------------
write.table(
	high_corr_distances,
	paste0("Output.Files/Correlation.Distances/",args[3]), 
	quote=FALSE, 
	sep="\t", 
	row.names=FALSE, 
	col.names=FALSE
)
