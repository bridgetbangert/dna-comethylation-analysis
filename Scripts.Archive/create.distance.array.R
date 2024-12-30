# ------------------------------------------------------------
# This R code will generate the vector of distinct distances
# of all of the CG sites.
#
# args[1] - input distance matrix filename
# args[2] - output distance vector filename
# ------------------------------------------------------------
args <- commandArgs(TRUE)
setwd("/home/bmb191/math5376D/Final.Project/")

# ------------------------------------------------------------
# Get alive subjects' data
# ------------------------------------------------------------
dist_matrix <- read.table(
	paste0("Output.Files/Distance.Matrices/",args[1]),
	header=FALSE
)

write.table(
	which(dist_matrix > 0, arr.ind = TRUE),
	paste0("Output.Files/Correlation.Distances/",args[2]), 
	quote=FALSE, 
	sep="\t", 
	row.names=FALSE, 
	col.names=FALSE
)
