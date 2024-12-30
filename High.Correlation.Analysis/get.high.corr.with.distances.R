# ------------------------------------------------------------
# This R file will grab all of the distinct distances of the
# highly correlated CG sites as a vector
#
# args[1] - correlation matrix file
# args[2] - distance matrix file
# args[3] - output index file name (will be saved in the 
# Output.Files/Distance.Indexes directory)
# args[4] - output distance file name (will be saved in the 
# Output.Files/Correlation.Distances directory)
# args[5] - high correlation threshold
# ------------------------------------------------------------
args <- commandArgs(TRUE)
setwd("/home/bmb191/math5376D/Final.Project/Output.Files")

# ------------------------------------------------------------
# This function will get the indexes of records that have
# high correlations in the matrix
# ------------------------------------------------------------
get_indexes <- function(cmat) {
	diag(cmat) <- 0
	return(which(abs(cmat)>=args[5], arr.ind=TRUE))
}
get_indexes_pos <- function(cmat) {
	diag(cmat) <- 0
	return(which(cmat>=args[5], arr.ind=TRUE))
}
get_indexes_neg <- function(cmat) {
	diag(cmat) <- 0
	threshold <- as.numeric(args[5])
	return(which(cmat<=(-1*threshold), arr.ind=TRUE))
}

# ------------------------------------------------------------
# This function will get the distances of the records that 
# have high correlations in the matrix
# ------------------------------------------------------------
get_dist <- function(cmat, indexes) {
	return(cmat[indexes])
}

# ------------------------------------------------------------
# Read correlation matrix
# ------------------------------------------------------------
corr_matrix <- read.table(paste0("Correlation.Matrices/",args[1]), header=FALSE)
dist_matrix <- read.table(paste0("Distance.Matrices/",args[2]), header=FALSE)

index_list <- get_indexes(corr_matrix)
write.table(
	index_list, 
	paste0("High.Corr.Distance.Indeces/",args[3]), 
	row.names=FALSE, 
	col.names=FALSE
)
index_list_pos <- get_indexes_pos(corr_matrix)
write.table(
	index_list_pos, 
	paste0("High.Corr.Distance.Indeces/positive.",args[3]), 
	row.names=FALSE, 
	col.names=FALSE
)
index_list_neg <- get_indexes_neg(corr_matrix)
write.table(
	index_list_neg, 
	paste0("High.Corr.Distance.Indeces/negative.",args[3]), 
	row.names=FALSE, 
	col.names=FALSE
)

dist_list <- get_dist(dist_matrix, index_list)
write.table(
	dist_list, 
	paste0("Correlation.Distances/",args[4]), 
	row.names=FALSE, 
	col.names=FALSE
)
dist_list_pos <- get_dist(dist_matrix, index_list_pos)
write.table(
	dist_list_pos, 
	paste0("Correlation.Distances/positive.",args[4]), 
	row.names=FALSE, 
	col.names=FALSE
)
dist_list_neg <- get_dist(dist_matrix, index_list_neg)
write.table(
	dist_list_neg, 
	paste0("Correlation.Distances/negative.",args[4]), 
	row.names=FALSE, 
	col.names=FALSE
)
