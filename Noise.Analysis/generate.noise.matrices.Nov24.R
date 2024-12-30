# ------------------------------------------------------------
# This script is to generate its input matrix with noise
# with the seed of the user's choice and output its 
# correlation matrix.
#
# args[1] - input filepath
# args[2] - starting column
# args[3] - ending column
# args[4] - output filepath
# args[5] - seed
# ------------------------------------------------------------
setwd("/home/bmb191/math5376D/Final.Project")
args <- commandArgs(TRUE)
set.seed(args[5])

generate_noise_matrix <- function(input.mat, noise.min, noise.max) {
	n.row <- nrow(input.mat)
	n.col <- as.numeric(args[3]) - as.numeric(args[2]) + 1
	return(matrix(
		runif(
			n.row*n.col, 
			min=noise.min, 
			max=noise.max
		), 
		nrow=n.row, 
		ncol=n.col)
	)
}

get_data_with_noise <- function(input.mat, noise.min, noise.max) {
	return(
		generate_noise_matrix(input.mat, noise.min, noise.max) 
		+ input.mat[ ,as.numeric(args[2]):as.numeric(args[3])]
	)
}

input_data <- read.table(paste0("Data/", args[1]), header = TRUE)
noise_matrix <- get_data_with_noise(input_data, -0.01, 0.01)
matrix_numerics <- apply(noise_matrix,1,as.numeric)
matrix_corr <- cor(matrix_numerics, use ="everything", method = "spearman")

write.table(
	round(matrix_numerics,4), 
	paste0("Data/",args[1],".with.noise.txt"),
	quote=FALSE, 
	sep="\t", 
	row.names=FALSE, 
	col.names=FALSE
)

write.table(
	round(matrix_corr,4), 
	paste0("Output.Files/Correlation.Matrices/",args[4]),
	quote=FALSE, 
	sep="\t", 
	row.names=FALSE, 
	col.names=FALSE
)

