# ------------------------------------------------------------
# Original correlation matrix code. This did not match
# Madison's results, but this is because she used the
# "everything" use parameter.
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
# alive normal correlation matrix
# ------------------------------------------------------------
alive_normal_numerics <- apply(alive_normal[,5:57],1,as.numeric)
alive_normal_corr_matrix <- cor(alive_normal_numerics, method="spearman")
write.table(
	alive_normal_corr_matrix, 
	"Output.Files/alive.normal.corr.matrix.txt", 
	row.names=FALSE, 
	col.names=FALSE
)

# ------------------------------------------------------------
# alive tumor correlation matrix
# ------------------------------------------------------------
alive_tumor_numerics <- apply(alive_tumor[,5:57],1,as.numeric)
alive_tumor_corr_matrix <- cor(alive_tumor_numerics, method="spearman")
write.table(
	alive_tumor_corr_matrix, 
	"Output.Files/alive.tumor.corr.matrix.txt", 
	row.names=FALSE, 
	col.names=FALSE
)

###############################################

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
# dead normal correlation matrix
# ------------------------------------------------------------
dead_normal_numerics <- apply(dead_normal[,5:36],1,as.numeric)
dead_normal_corr_matrix <- cor(dead_normal_numerics, method="spearman")
write.table(
  dead_normal_corr_matrix, 
  "Output.Files/dead.normal.corr.matrix.txt", 
  row.names=FALSE, 
  col.names=FALSE
)

# ------------------------------------------------------------
# dead tumor correlation matrix
# ------------------------------------------------------------
dead_tumor_numerics <- apply(dead_tumor[,5:36],1,as.numeric)
dead_tumor_corr_matrix <- cor(dead_tumor_numerics, method="spearman")
write.table(
	dead_tumor_corr_matrix, 
	"Output.Files/dead.tumor.corr.matrix.txt", 
	row.names=FALSE, 
	col.names=FALSE
)

###############################################

# ------------------------------------------------------------
# get counts for each row alive matrix with |corr|>.8
# ------------------------------------------------------------
high_corr_counts <- function(corr_matrix) {
	apply(corr_matrix, 1, function(x) sum(abs(x)>.8))
}

# ------------------------------------------------------------
# call the high_corr_counts function for the alive data
# ------------------------------------------------------------
alive_tumor_high_corr_cnts <- high_corr_counts(alive_tumor_corr_matrix)
alive_normal_high_corr_cnts <- high_corr_counts(alive_normal_corr_matrix)

# ------------------------------------------------------------
# boxplots for the alive data
# ------------------------------------------------------------
boxplot(
	alive_tumor_high_corr_cnts,
	alive_normal_high_corr_cnts,
	main="Alive high correlation counts",
	names=c("Tumor","Normal")
)
boxplot(
  alive_tumor_high_corr_cnts,
  alive_normal_high_corr_cnts,
  main="Alive high correlation counts (removed outliers)",
  names=c("Tumor","Normal"),
  outline=FALSE
)

# ------------------------------------------------------------
# histograms for the alive data
# ------------------------------------------------------------
hist(
  alive_normal_high_corr_cnts,
  main="Alive high correlation frequencies (non-restricted range)",
  xlab="High Correlation Counts",
  ylab="",
  col=rgb(0, 0, 1, alpha=0.5),
  breaks=100
)
hist(
  alive_tumor_high_corr_cnts,
  col=rgb(1, 0, 0, alpha=0.5),
  breaks=100,
  add=TRUE
)
legend(
  "topright", 
  legend=c("Normal","Tumor"), 
  fill=c(rgb(0, 0, 1, alpha=0.5), rgb(1, 0, 0, alpha=0.5))
)
hist(
  alive_normal_high_corr_cnts,
  main="Alive high correlation frequencies (Restricted range)",
  xlab="High Correlation Counts",
  ylab="",
  ylim=c(0,100),
  col=rgb(0, 0, 1, alpha=0.5),
  breaks=100
)
hist(
  alive_tumor_high_corr_cnts,
  col=rgb(1, 0, 0, alpha=0.5),
  breaks=100,
  add=TRUE
)
legend(
  "topright", 
  legend=c("Normal","Tumor"), 
  fill=c(rgb(0, 0, 1, alpha=0.5), rgb(1, 0, 0, alpha=0.5))
)

# ------------------------------------------------------------
# call the high_corr_counts function for the dead data
# ------------------------------------------------------------
dead_tumor_high_corr_cnts <- high_corr_counts(dead_tumor_corr_matrix)
dead_normal_high_corr_cnts <- high_corr_counts(dead_normal_corr_matrix)

# ------------------------------------------------------------
# boxplots for the dead data
# ------------------------------------------------------------
boxplot(
  dead_tumor_high_corr_cnts,
  dead_normal_high_corr_cnts,
  main="Dead high correlation counts",
  names=c("Tumor","Normal")
)
boxplot(
  dead_tumor_high_corr_cnts,
  dead_normal_high_corr_cnts,
  main="Dead high correlation counts (removed outliers)",
  names=c("Tumor","Normal"),
  outline=FALSE
)
hist(
  dead_normal_high_corr_cnts,
  main="Dead high correlation frequencies (non-restricted range)",
  xlab="High Correlation Counts",
  ylab="",
  col=rgb(0, 0, 1, alpha=0.5),
  breaks=100
)
hist(
  dead_tumor_high_corr_cnts,
  col=rgb(1, 0, 0, alpha=0.5),
  breaks=100,
  add=TRUE
)
legend(
  "topright", 
  legend=c("Normal","Tumor"), 
  fill=c(rgb(0, 0, 1, alpha=0.5), rgb(1, 0, 0, alpha=0.5))
)

# ------------------------------------------------------------
# histograms for the dead data
# ------------------------------------------------------------
hist(
  dead_normal_high_corr_cnts,
  main="Dead high correlation frequencies (Restricted range)",
  xlab="High Correlation Counts",
  ylab="",
  ylim=c(0,100),
  col=rgb(0, 0, 1, alpha=0.5),
  breaks=100
)
hist(
  dead_tumor_high_corr_cnts,
  col=rgb(1, 0, 0, alpha=0.5),
  breaks=100,
  add=TRUE
)
legend(
  "topright", 
  legend=c("Normal","Tumor"), 
  fill=c(rgb(0, 0, 1, alpha=0.5), rgb(1, 0, 0, alpha=0.5))
)
