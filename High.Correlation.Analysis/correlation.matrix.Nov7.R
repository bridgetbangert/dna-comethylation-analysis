# ------------------------------------------------------------
# This is an updated version of the original correlation
# matrix computation file. This one will use the
# "everything" param to compute the matrix, so that it matches
# with Madison's results.
# ------------------------------------------------------------
setwd("/home/bmb191/math5376D/Final.Project")

# ------------------------------------------------------------
# Get alive subjects' data
# ------------------------------------------------------------
alive_normal <- read.table(
	"Data/chrX9498cg.BRCA.53Alive.Normal.380355cg.75col.w.Header.txt",
	header=T
)
alive_tumor <- read.table(
	"Data/chrX9498cg.BRCA.53Alive.Tumor.380355cg.75col.w.Header.txt",
	header=T
)

# ------------------------------------------------------------
# Get dead subjects' data
# ------------------------------------------------------------
dead_normal <- read.table(
	"Data/chrX9498cg.BRCA.32Dead.Normal.380355cg.54col.w.Header.txt", 
	header=T
)
dead_tumor <- read.table(
	"Data/chrX9498cg.BRCA.32Dead.Tumor.380355cg.54col.w.Header.txt", 
	header=T
)

# ------------------------------------------------------------
# Make the methylation signals numeric types and calculate
# the correlation matrix using the Spearman method, then
# write to file in the leap server.
# ------------------------------------------------------------

# ------------------------------------------------------------
# Alive normal data
# ------------------------------------------------------------
alive_normal_numerics <- apply(alive_normal[,5:57],1,as.numeric)
alive_normal_corr <- cor(alive_normal_numerics, use ="everything", method = "spearman")
write.table(
	alive_normal_corr, 
	"Output.Files/Correlation.Matrices/BRCA.alive.normal.corr.matrix.txt", 
	quote=FALSE, 
	sep="\t", 
	row.names=FALSE, 
	col.names=FALSE
)

# ------------------------------------------------------------
# Alive tumor data
# ------------------------------------------------------------
alive_tumor_numerics <- apply(alive_tumor[,5:57],1,as.numeric)
alive_tumor_corr <- cor(alive_tumor_numerics, use ="everything", method = "spearman")
write.table(
	alive_tumor_corr, 
	"Output.Files/Correlation.Matrices/BRCA.alive.tumor.corr.matrix.txt", 
	quote=FALSE, 
	sep="\t", 
	row.names=FALSE, 
	col.names=FALSE
)

# ------------------------------------------------------------
# Dead normal data
# ------------------------------------------------------------
dead_normal_numerics <- apply(dead_normal[,5:36],1,as.numeric)
dead_normal_corr <- cor(dead_normal_numerics, use ="everything", method = "spearman")
write.table(
	dead_normal_corr, 
	"Output.Files/Correlation.Matrices/BRCA.dead.normal.corr.matrix.txt", 
	quote=FALSE, 
	sep="\t", 
	row.names=FALSE, 
	col.names=FALSE
)

# ------------------------------------------------------------
# Dead tumor data
# ------------------------------------------------------------
dead_tumor_numerics <- apply(dead_tumor[,5:36],1,as.numeric)
dead_tumor_corr <- cor(dead_tumor_numerics, use ="everything", method = "spearman")
write.table(
	dead_tumor_corr, 
	"Output.Files/Correlation.Matrices/BRCA.dead.tumor.corr.matrix.txt", 
	quote=FALSE, 
	sep="\t", 
	row.names=FALSE, 
	col.names=FALSE
)
