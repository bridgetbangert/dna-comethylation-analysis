# ------------------------------------------------------------
# This R file will compute the distance matrix of all of
# the start points in the X chromosome
#
# arg[1] - chrX data filename
# arg[2] - output file name (will be stored in Output.Files/
# Distance.Matrices/ directory)
# ------------------------------------------------------------
args <- commandArgs(TRUE)
setwd("/home/bmb191/math5376D/Final.Project/")
cg_sites <- read.table(paste0("Data/",args[1]), header=TRUE)

# ------------------------------------------------------------
# Make start column numerical in original data
# ------------------------------------------------------------
cg_sites <- cbind(
	cg_sites[1:2], apply(cg_sites[3:4],2, as.numeric), 
	cg_sites[5:length(colnames(cg_sites))]
)

# ------------------------------------------------------------
# Initialize matrix to grab distances
# ------------------------------------------------------------
start_site_list <- cg_sites[[3]]
distance_matrix <- matrix(0,nrow=nrow(cg_sites),ncol=nrow(cg_sites))

# ------------------------------------------------------------
# Put distances of each CG pair into a matrix
# ------------------------------------------------------------
for (i in 1:nrow(cg_sites)) { 
  site_start <- cg_sites[[3]][i]
  for (j in 1:nrow(cg_sites)) { 
	  distance_matrix[i,j] <- abs(site_start - start_site_list[j])
  } 
} 

# ------------------------------------------------------------
# Write the matrix to a file in the server
# ------------------------------------------------------------
write.table(
	distance_matrix,
	paste0("Output.Files/Distance.Matrices/",args[2]), 
	quote=FALSE, 
	sep="\t", 
	row.names=FALSE, 
	col.names=FALSE
)