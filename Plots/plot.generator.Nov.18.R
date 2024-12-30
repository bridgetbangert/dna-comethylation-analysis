setwd("/home/bmb191/math5376D/Final.Project/Output.Files")

setClass("CountsAnalysis",
         slots = list(
           group="character",
           normal_cnts="list",
           tumor_cnts="list",
           normal_cnts_pos="list",
           tumor_cnts_pos="list",
           normal_cnts_neg="list",
           tumor_cnts_neg="list"
))
setClass("DistancesAnalysis",
         slots = list(
           group="character",
           normal_dist="list",
           tumor_dist="list",
           normal_dist_pos="list",
           tumor_dist_pos="list",
           normal_dist_neg="list",
           tumor_dist_neg="list"
))

load_data <- function(analysis_subset, sign, cell_health, measure) {
  return(read.table(
    get_directory(analysis_subset, sign, cell_health, measure),
    header=FALSE
  ))
}

get_directory <- function(analysis_subset, sign, cell_health, measure) {
  if (sign != "") { sign = paste0(sign,".") }
  if (measure == "counts") { 
    directory="Correlation.Counts/"
    filename=".high.corr.cnts.txt"
  }
  else if (measure == "distances") { 
    directory="Correlation.Distances/"
    filename=".high.corr.distances.txt" 
  }
  else { print("Measure not found. Please try again"); return() }
  return(paste0(
    directory,
    sign,
    analysis_subset,
    ".",
    cell_health,
    filename
  ))
}

get_group <- function(analysis_subset) {
  grp_str <- strsplit(analysis_subset,"\\.")[[1]]
  disease <- grp_str[1]
  status <- grp_str[2]
  return(paste(
    grp_str[1], 
    paste0(
      toupper(substring(status, 1, 1)), 
      substring(status, 2)
    )
  ))
}

CountsAnalysis <- function(analysis_subset) {
  new("CountsAnalysis", 
    group=get_group(analysis_subset),
    normal_cnts=load_data(analysis_subset,"","normal","counts"),
    tumor_cnts=load_data(analysis_subset,"","tumor","counts"),
    normal_cnts_pos=load_data(analysis_subset,"positive","normal","counts"),
    tumor_cnts_pos=load_data(analysis_subset,"positive","tumor","counts"),
    normal_cnts_neg=load_data(analysis_subset,"negative","normal","counts"),
    tumor_cnts_neg=load_data(analysis_subset,"negative","tumor","counts")
  )
}

DistancesAnalysis <- function(analysis_subset) {
  new("DistancesAnalysis", 
      group=get_group(analysis_subset),
      normal_dist=load_data(analysis_subset,"","normal","distances"),
      tumor_dist=load_data(analysis_subset,"","tumor","distances"),
      normal_dist_pos=load_data(analysis_subset,"positive","normal","distances"),
      tumor_dist_pos=load_data(analysis_subset,"positive","tumor","distances"),
      normal_dist_neg=load_data(analysis_subset,"negative","normal","distances"),
      tumor_dist_neg=load_data(analysis_subset,"negative","tumor","distances")
  )
}

setGeneric("boxplot_comparison_cnts", function(object, scale) standardGeneric("boxplot_comparison_cnts"))
setMethod("boxplot_comparison_cnts", "CountsAnalysis", function(object, scale) {
  boxplot(unlist(object@normal_cnts)^(scale), unlist(object@tumor_cnts)^(scale), 
          names = c("Normal Counts", "Tumor Counts"),
          main = paste(object@group,"Normal vs. Tumor Count Comparison"))
})

setGeneric("pos_boxplot_comparison_cnts", function(object, scale) standardGeneric("pos_boxplot_comparison_cnts"))
setMethod("pos_boxplot_comparison_cnts", "CountsAnalysis", function(object, scale) {
  boxplot(unlist(object@normal_cnts_pos)^(scale), unlist(object@tumor_cnts_pos)^(scale), 
          names = c("Normal Counts", "Tumor Counts"),
          main = paste(object@group,"Normal vs. Tumor Positive Highly Correlated Counts Comparison"))
})

setGeneric("neg_boxplot_comparison_cnts", function(object, scale) standardGeneric("neg_boxplot_comparison_cnts"))
setMethod("neg_boxplot_comparison_cnts", "CountsAnalysis", function(object, scale) {
  boxplot(unlist(object@normal_cnts_neg)^(scale), unlist(object@tumor_cnts_neg)^(scale), 
          names = c("Normal Counts", "Tumor Counts"),
          main = paste(object@group,"Normal vs. Tumor Negative Highly Correlated Counts Comparison"))
})

setGeneric("density_comparison_cnts", function(object, scale) standardGeneric("density_comparison_cnts"))
setMethod("density_comparison_cnts", "CountsAnalysis", function(object, scale) {
  plot(density(
      unlist(object@normal_cnts)^(scale)), 
      col="blue",
      main =paste(object@group,"Normal vs. Tumor Count Comparison"))
  legend("topright",legend=c("Normal","Tumor"),col=c("blue","red"),lty=1,lwd=2)
  lines(
    density(unlist(object@tumor_cnts)^(scale)),
    col="red")
})

setGeneric("pos_density_comparison_cnts", function(object, scale) standardGeneric("pos_density_comparison_cnts"))
setMethod("pos_density_comparison_cnts", "CountsAnalysis", function(object, scale) {
  plot(density(
      unlist(object@normal_cnts_pos)^(scale)), 
      col="blue",
      main =paste(object@group,"Normal vs. Tumor Positive Highly Correlated Count Comparison"))
  legend("topright",legend=c("Normal","Tumor"),col=c("blue","red"),lty=1,lwd=2)
  lines(
    density(unlist(object@tumor_cnts_pos)^(scale)),
    col="red")
})

setGeneric("neg_density_comparison_cnts", function(object, scale) standardGeneric("neg_density_comparison_cnts"))
setMethod("neg_density_comparison_cnts", "CountsAnalysis", function(object, scale) {
  plot(density(
      unlist(object@normal_cnts_neg)^(scale)), 
      col="blue",
      main =paste(object@group,"Normal vs. Tumor Negative Highly Correlated Count Comparison"))
  legend("topright",legend=c("Normal","Tumor"),col=c("blue","red"),lty=1,lwd=2)
  lines(
    density(unlist(object@tumor_cnts_neg)^(scale)),
    col="red")
})

setGeneric("boxplot_comparison_dist", function(object, scale) standardGeneric("boxplot_comparison_dist"))
setMethod("boxplot_comparison_dist", "DistancesAnalysis", function(object, scale) {
  boxplot(unlist(object@normal_dist)^(scale), unlist(object@tumor_dist)^(scale), 
          names = c("Normal Counts", "Tumor Counts"),
          main = paste(object@group,"Normal vs. Tumor Distances Comparison"))
})

setGeneric("pos_boxplot_comparison_dist", function(object, scale) standardGeneric("pos_boxplot_comparison_dist"))
setMethod("pos_boxplot_comparison_dist", "DistancesAnalysis", function(object, scale) {
  boxplot(unlist(object@normal_dist_pos)^(scale), unlist(object@tumor_dist_pos)^(scale), 
          names = c("Normal Counts", "Tumor Counts"),
          main = paste(object@group,"Normal vs. Tumor Positive Highly Correlated Distances Comparison"))
})

setGeneric("neg_boxplot_comparison_dist", function(object, scale) standardGeneric("neg_boxplot_comparison_dist"))
setMethod("neg_boxplot_comparison_dist", "DistancesAnalysis", function(object, scale) {
  boxplot(unlist(object@normal_dist_neg)^(scale), unlist(object@tumor_dist_neg)^(scale), 
          names = c("Normal Counts", "Tumor Counts"),
          main = paste(object@group,"Normal vs. Tumor Negative Highly Correlated Distances Comparison"))
})

setGeneric("density_comparison_dist", function(object, scale) standardGeneric("density_comparison_dist"))
setMethod("density_comparison_dist", "DistancesAnalysis", function(object, scale) {
  plot(density(
    unlist(object@normal_dist)^(scale)), 
    col="blue",
    main =paste(object@group,"Normal vs. Tumor Distances Comparison"))
  legend("topright",legend=c("Normal","Tumor"),col=c("blue","red"),lty=1,lwd=2)
  lines(
    density(unlist(object@tumor_dist)^(scale)),
    col="red")
})

setGeneric("pos_density_comparison_dist", function(object, scale) standardGeneric("pos_density_comparison_dist"))
setMethod("pos_density_comparison_dist", "DistancesAnalysis", function(object, scale) {
  plot(density(
    unlist(object@normal_dist_pos)^(scale)), 
    col="blue",
    main =paste(object@group,"Normal vs. Tumor Positive Highly Correlated Distances Comparison"))
  legend("topright",legend=c("Normal","Tumor"),col=c("blue","red"),lty=1,lwd=2)
  lines(
    density(unlist(object@tumor_dist_pos)^(scale)),
    col="red")
})

setGeneric("neg_density_comparison_dist", function(object, scale) standardGeneric("neg_density_comparison_dist"))
setMethod("neg_density_comparison_dist", "DistancesAnalysis", function(object, scale) {
  plot(density(
    unlist(object@normal_dist_neg)^(scale)), 
    col="blue",
    main =paste(object@group,"Normal vs. Tumor Negative Highly Correlated Distances Comparison"))
  legend("topright",legend=c("Normal","Tumor"),col=c("blue","red"),lty=1,lwd=2)
  lines(
    density(unlist(object@tumor_dist_neg)^(scale)),
    col="red")
})

setGeneric("print_summaries_cnts", function(object) standardGeneric("print_summaries_cnts"))
setMethod("print_summaries_cnts", "CountsAnalysis", function(object) {
  print(paste("Distribution of",object@group,"alive counts:"))
  print(summary(unlist(object@normal_cnts)))
  print(paste("Distribution of",object@group,"tumor counts:"))
  print(summary(unlist(object@tumor_cnts)))
})

setGeneric("pos_print_summaries_cnts", function(object) standardGeneric("pos_print_summaries_cnts"))
setMethod("pos_print_summaries_cnts", "CountsAnalysis", function(object) {
  print(paste("Distribution of",object@group,"positive highly correlated alive counts:"))
  print(summary(unlist(object@normal_cnts_pos)))
  print(paste("Distribution of",object@group,"positive highly correlated tumor counts:"))
  print(summary(unlist(object@tumor_cnts_pos)))
})

setGeneric("neg_print_summaries_cnts", function(object) standardGeneric("neg_print_summaries_cnts"))
setMethod("neg_print_summaries_cnts", "CountsAnalysis", function(object) {
  print(paste("Distribution of",object@group,"negative highly correlated alive counts:"))
  print(summary(unlist(object@normal_cnts_neg)))
  print(paste("Distribution of",object@group,"negative highly correlated tumor counts:"))
  print(summary(unlist(object@tumor_cnts_neg)))
})

setGeneric("print_summaries_dist", function(object) standardGeneric("print_summaries_dist"))
setMethod("print_summaries_dist", "DistancesAnalysis", function(object) {
  print(paste("Distribution of",object@group,"alive distances:"))
  print(summary(unlist(object@normal_dist)))
  print(paste("Distribution of",object@group,"tumor distances:"))
  print(summary(unlist(object@tumor_dist)))
})

setGeneric("pos_print_summaries_dist", function(object) standardGeneric("pos_print_summaries_dist"))
setMethod("pos_print_summaries_dist", "DistancesAnalysis", function(object) {
  print(paste("Distribution of",object@group,"positive highly correlated alive distances:"))
  print(summary(unlist(object@normal_dist_pos)))
  print(paste("Distribution of",object@group,"positive highly correlated tumor distances:"))
  print(summary(unlist(object@tumor_dist_pos)))
})

setGeneric("neg_print_summaries_dist", function(object) standardGeneric("neg_print_summaries_dist"))
setMethod("neg_print_summaries_dist", "DistancesAnalysis", function(object) {
  print(paste("Distribution of",object@group,"negative highly correlated alive distances:"))
  print(summary(unlist(object@normal_dist_neg)))
  print(paste("Distribution of",object@group,"negative highly correlated tumor distances:"))
  print(summary(unlist(object@tumor_dist_neg)))
})

# counts plots
BRCA_alive_cnts <- CountsAnalysis("BRCA.alive")
BRCA_dead_cnts <- CountsAnalysis("BRCA.dead")
UCEC_alive_cnts <- CountsAnalysis("UCEC.alive")
UCEC_dead_cnts <- CountsAnalysis("UCEC.dead")

pdf("Plots/counts.boxplots.pdf",16,16)
par(mfrow=c(2,2))
boxplot_comparison_cnts(BRCA_alive_cnts, 1/4)
boxplot_comparison_cnts(BRCA_dead_cnts, 1/4)
boxplot_comparison_cnts(UCEC_alive_cnts, 1/4)
boxplot_comparison_cnts(UCEC_dead_cnts, 1/4)
par(mfrow=c(1,1))
dev.off()

pdf("Plots/positive.counts.boxplots.pdf",16,16)
par(mfrow=c(2,2))
pos_boxplot_comparison_cnts(BRCA_alive_cnts, 1/4)
pos_boxplot_comparison_cnts(BRCA_dead_cnts, 1/4)
pos_boxplot_comparison_cnts(UCEC_alive_cnts, 1/4)
pos_boxplot_comparison_cnts(UCEC_dead_cnts, 1/4)
par(mfrow=c(1,1))
dev.off()

pdf("Plots/negative.counts.boxplots.pdf",16,16)
par(mfrow=c(2,2))
neg_boxplot_comparison_cnts(BRCA_alive_cnts, 1/4)
neg_boxplot_comparison_cnts(BRCA_dead_cnts, 1/4)
neg_boxplot_comparison_cnts(UCEC_alive_cnts, 1/4)
neg_boxplot_comparison_cnts(UCEC_dead_cnts, 1/4)
par(mfrow=c(1,1))
dev.off()

pdf("Plots/counts.density.plots.pdf",16,16)
par(mfrow=c(2,2))
density_comparison_cnts(BRCA_alive_cnts, 1/4)
density_comparison_cnts(BRCA_dead_cnts, 1/4)
density_comparison_cnts(UCEC_alive_cnts, 1/4)
density_comparison_cnts(UCEC_dead_cnts, 1/4)
par(mfrow=c(1,1))
dev.off()

pdf("Plots/positive.counts.density.plots.pdf",16,16)
par(mfrow=c(2,2))
pos_density_comparison_cnts(BRCA_alive_cnts, 1/4)
pos_density_comparison_cnts(BRCA_dead_cnts, 1/4)
pos_density_comparison_cnts(UCEC_alive_cnts, 1/4)
pos_density_comparison_cnts(UCEC_dead_cnts, 1/4)
par(mfrow=c(1,1))
dev.off()

pdf("Plots/negative.counts.density.plots.pdf",16,16)
par(mfrow=c(2,2))
neg_density_comparison_cnts(BRCA_alive_cnts, 1/4)
neg_density_comparison_cnts(BRCA_dead_cnts, 1/4)
neg_density_comparison_cnts(UCEC_alive_cnts, 1/4)
neg_density_comparison_cnts(UCEC_dead_cnts, 1/4)
par(mfrow=c(1,1))
dev.off()

# summaries
print_summaries_cnts(BRCA_alive_cnts)
print_summaries_cnts(BRCA_dead_cnts)
print_summaries_cnts(UCEC_alive_cnts)
print_summaries_cnts(UCEC_dead_cnts)

#  positive summaries
pos_print_summaries_cnts(BRCA_alive_cnts)
pos_print_summaries_cnts(BRCA_dead_cnts)
pos_print_summaries_cnts(UCEC_alive_cnts)
pos_print_summaries_cnts(UCEC_dead_cnts)

# negative summaries
neg_print_summaries_cnts(BRCA_alive_cnts)
neg_print_summaries_cnts(BRCA_dead_cnts)
neg_print_summaries_cnts(UCEC_alive_cnts)
neg_print_summaries_cnts(UCEC_dead_cnts)

# distances plots
BRCA_alive_dist <- DistancesAnalysis("BRCA.alive")
BRCA_dead_dist <- DistancesAnalysis("BRCA.dead")
UCEC_alive_dist <- DistancesAnalysis("UCEC.alive")
UCEC_dead_dist <- DistancesAnalysis("UCEC.dead")

pdf("Plots/distances.boxplots.pdf",16,16)
par(mfrow=c(2,2))
boxplot_comparison_dist(BRCA_alive_dist, 1/4)
boxplot_comparison_dist(BRCA_dead_dist, 1/4)
boxplot_comparison_dist(UCEC_alive_dist, 1/4)
boxplot_comparison_dist(UCEC_dead_dist, 1/4)
par(mfrow=c(1,1))
dev.off()

pdf("Plots/positive.distances.boxplots.pdf",16,16)
par(mfrow=c(2,2))
pos_boxplot_comparison_dist(BRCA_alive_dist, 1/4)
pos_boxplot_comparison_dist(BRCA_dead_dist, 1/4)
pos_boxplot_comparison_dist(UCEC_alive_dist, 1/4)
pos_boxplot_comparison_dist(UCEC_dead_dist, 1/4)
par(mfrow=c(1,1))
dev.off()

pdf("Plots/negative.distances.boxplots.pdf",16,16)
par(mfrow=c(2,2))
neg_boxplot_comparison_dist(BRCA_alive_dist, 1/4)
neg_boxplot_comparison_dist(BRCA_dead_dist, 1/4)
neg_boxplot_comparison_dist(UCEC_alive_dist, 1/4)
neg_boxplot_comparison_dist(UCEC_dead_dist, 1/4)
par(mfrow=c(1,1))
dev.off()

pdf("Plots/distances.density.plots.pdf",16,16)
par(mfrow=c(2,2))
density_comparison_dist(BRCA_alive_dist, 1/4)
density_comparison_dist(BRCA_dead_dist, 1/4)
density_comparison_dist(UCEC_alive_dist, 1/4)
density_comparison_dist(UCEC_dead_dist, 1/4)
par(mfrow=c(1,1))
dev.off()

pdf("Plots/positive.distances.density.plots.pdf",16,16)
par(mfrow=c(2,2))
pos_density_comparison_dist(BRCA_alive_dist, 1/4)
pos_density_comparison_dist(BRCA_dead_dist, 1/4)
pos_density_comparison_dist(UCEC_alive_dist, 1/4)
pos_density_comparison_dist(UCEC_dead_dist, 1/4)
par(mfrow=c(1,1))
dev.off()

pdf("Plots/negative.distances.density.plots.pdf",16,16)
par(mfrow=c(2,2))
neg_density_comparison_dist(BRCA_alive_dist, 1/4)
neg_density_comparison_dist(BRCA_dead_dist, 1/4)
neg_density_comparison_dist(UCEC_alive_dist, 1/4)
neg_density_comparison_dist(UCEC_dead_dist, 1/4)
par(mfrow=c(1,1))
dev.off()

# summaries
print_summaries_dist(BRCA_alive_dist)
print_summaries_dist(BRCA_dead_dist)
print_summaries_dist(UCEC_alive_dist)
print_summaries_dist(UCEC_dead_dist)

# positive summaries
pos_print_summaries_dist(BRCA_alive_dist)
pos_print_summaries_dist(BRCA_dead_dist)
pos_print_summaries_dist(UCEC_alive_dist)
pos_print_summaries_dist(UCEC_dead_dist)

# negative summaries
neg_print_summaries_dist(BRCA_alive_dist)
neg_print_summaries_dist(BRCA_dead_dist)
neg_print_summaries_dist(UCEC_alive_dist)
neg_print_summaries_dist(UCEC_dead_dist)

