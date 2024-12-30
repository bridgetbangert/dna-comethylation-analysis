setwd("/home/bmb191/math5376D/Final.Project/Output.Files")

create_distance_plots <- function(corr_type, file_suffix) {
  if (corr_type == "Absolute") { output_path_corr <- "" }
  else {output_path_corr <- paste0(tolower(corr_type),".") }
  pdf(paste0("Plots/BRCA.",output_path_corr,"highly.correlated.distances.boxplot",file_suffix,".pdf"))
  layout(matrix(c(1,2,3,1,4,5), ncol=2, nrow=3), heights=c(1,1,14), widths=c(16,16))
  par(mai=rep(0,4))
  plot.new()
  text(0.5,0.5,paste("Fourth Root Distances of BRCA",corr_type,"Highly Correlated Data"),cex=1.5,font=1)
  plot.new()
  text(0.5,0.5,"BRCA Alive",cex=1,font=1)
  par(mai=c(0.3,0.3,0,0.15))
  boxplot(
    BRCA_an^(1/4), 
    BRCA_at^(1/4),
    names=c("Normal","Tumor")
  )
  par(mai=rep(0,4))
  plot.new()
  text(0.5,0.5,"BRCA Dead",cex=1,font=1)
  par(mai=c(0.3,0.15,0,0.3))
  boxplot(
    BRCA_dn^(1/4), 
    BRCA_dt^(1/4),
    names=c("Normal","Tumor")
  )
  dev.off()
  
  pdf(paste0("Plots/BRCA.",output_path_corr,"highly.correlated.distances.density.plot",file_suffix,".pdf"))
  layout(matrix(c(1,2,3,4,1,5,6,4), ncol=2, nrow=4), heights=c(1,1,14,1), widths=c(16,16))
  par(mai=rep(0,4))
  plot.new()
  text(0.5,0.5,paste("Fourth Root Distances of BRCA",corr_type,"Highly Correlated Data"),cex=1.5,font=1)
  plot.new()
  text(0.5,0.5,"BRCA Alive",cex=1,font=1)
  par(mai=c(0.3,0.3,0,0.15))
  plot(
    density((BRCA_an)^(1/4)),
    col="blue",
    ylim=c(0,.03),
    main=NA
  )
  lines(density((BRCA_at)^(1/4)), col="red")
  par(mai=rep(0,4))
  plot.new()
  legend(x="center", ncol=2,legend=c("Normal","Tumor"),
         fill=c("blue","red"), bty="n")
  plot.new()
  text(0.5,0.5,"BRCA Dead",cex=1,font=1)
  par(mai=c(0.3,0.15,0,0.3))
  plot(
    density((BRCA_dn)^(1/4)),
    col="blue",
    ylim=c(0,.03),
    main=NA
  )
  lines(density((BRCA_dt)^(1/4)), col="red")
  dev.off()
  
  pdf(paste0("Plots/UCEC.",output_path_corr,"highly.correlated.distances.boxplot",file_suffix,".pdf"))
  layout(matrix(c(1,2,3,1,4,5), ncol=2, nrow=3), heights=c(1,1,14), widths=c(16,16))
  par(mai=rep(0,4))
  plot.new()
  text(0.5,0.5,"Fourth Root Distances of UCEC Highly Correlated Data",cex=1.5,font=1)
  plot.new()
  text(0.5,0.5,"UCEC Alive",cex=1,font=1)
  par(mai=c(0.3,0.3,0,0.15))
  boxplot(
    UCEC_an^(1/4), 
    UCEC_at^(1/4),
    names=c("Normal","Tumor")
  )
  par(mai=rep(0,4))
  plot.new()
  text(0.5,0.5,"UCEC Dead",cex=1,font=1)
  par(mai=c(0.3,0.15,0,0.3))
  boxplot(
    UCEC_dn^(1/4), 
    UCEC_dt^(1/4),
    names=c("Normal","Tumor")
  )
  dev.off()
  
  pdf(paste0("Plots/UCEC.",output_path_corr,"highly.correlated.distances.density.plot",file_suffix,".pdf"))
  layout(matrix(c(1,2,3,4,1,5,6,4), ncol=2, nrow=4), heights=c(1,1,14,1), widths=c(16,16))
  par(mai=rep(0,4))
  plot.new()
  text(0.5,0.5,"Fourth Root Distances of UCEC Highly Correlated Data",cex=1.5,font=1)
  plot.new()
  text(0.5,0.5,"UCEC Alive",cex=1,font=1)
  par(mai=c(0.3,0.3,0,0.15))
  plot(
    density((UCEC_an)^(1/4)),
    col="blue",
    ylim=c(0,.03),
    main=NA
  )
  lines(density((UCEC_at)^(1/4)), col="red")
  par(mai=rep(0,4))
  plot.new()
  legend(x="center", ncol=2,legend=c("Normal","Tumor"),
         fill=c("blue","red"), bty="n")
  plot.new()
  text(0.5,0.5,"UCEC Dead",cex=1,font=1)
  par(mai=c(0.3,0.15,0,0.3))
  plot(
    density((UCEC_dn)^(1/4)),
    col="blue",
    ylim=c(0,.03),
    main=NA
  )
  lines(density((UCEC_dt)^(1/4)), col="red")
  dev.off()
}

create_count_plots <- function(corr_type, file_suffix) {
  output_path_corr <- paste0(tolower(corr_type),".")
  if file_suffix != "" {
  	noise_ind <- "(Noise)"
  }
  else { noise_ind <- "" }
  pdf(paste0("Plots/BRCA.",output_path_corr,"highly.correlated.counts.boxplot",file_suffix,".pdf"))
  layout(matrix(c(1,2,3,1,4,5), ncol=2, nrow=3), heights=c(1,1,14), widths=c(16,16))
  par(mai=rep(0,4))
  plot.new()
  text(0.5,0.5,paste("Fourth Root Counts of BRCA",corr_type,"Highly Correlated Data"),cex=1.5,font=1)
  plot.new()
  text(0.5,0.5,"BRCA Alive",cex=1,font=1)
  par(mai=c(0.3,0.3,0,0.15))
  boxplot(
    BRCA_an^(1/4), 
    BRCA_at^(1/4),
    names=c("Normal","Tumor"),
    ylim=c(0,5)
  )
  par(mai=rep(0,4))
  plot.new()
  text(0.5,0.5,"BRCA Dead",cex=1,font=1)
  par(mai=c(0.3,0.15,0,0.3))
  boxplot(
    BRCA_dn^(1/4), 
    BRCA_dt^(1/4),
    names=c("Normal","Tumor"),
    ylim=c(0,5)
  )
  dev.off()
  
  pdf(paste0("Plots/BRCA.",output_path_corr,"highly.correlated.counts.density.plot",file_suffix,".pdf"))
  layout(matrix(c(1,2,3,4,1,5,6,4), ncol=2, nrow=4), heights=c(1,1,14,1), widths=c(16,16))
  par(mai=rep(0,4))
  plot.new()
  text(0.5,0.5,paste("Fourth Root Counts of BRCA",corr_type,"Highly Correlated Data",noise_ind),cex=1.5,font=1)
  plot.new()
  text(0.5,0.5,"BRCA Alive",cex=1,font=1)
  par(mai=c(0.3,0.3,0,0.15))
  plot(
    density((BRCA_an)^(1/4)),
    col="blue",
    main=NA,
    ylim=c(0,3)
  )
  lines(density((BRCA_at)^(1/4)), col="red")
  par(mai=rep(0,4))
  plot.new()
  legend(x="center", ncol=2,legend=c("Normal","Tumor"),
         fill=c("blue","red"), bty="n")
  plot.new()
  text(0.5,0.5,"BRCA Dead",cex=1,font=1)
  par(mai=c(0.3,0.15,0,0.3))
  plot(
    density((BRCA_dn)^(1/4)),
    col="blue",
    main=NA,
    ylim=c(0,3)
  )
  lines(density((BRCA_dt)^(1/4)), col="red")
  dev.off()
  
  
  pdf(paste0("Plots/UCEC.",output_path_corr,"highly.correlated.counts.boxplot",file_suffix,".pdf"))
  layout(matrix(c(1,2,3,1,4,5), ncol=2, nrow=3), heights=c(1,1,14), widths=c(16,16))
  par(mai=rep(0,4))
  plot.new()
  text(0.5,0.5,paste("Fourth Root Counts of",corr_type,"UCEC Highly Correlated Data",noise_ind),cex=1.5,font=1)
  plot.new()
  text(0.5,0.5,"UCEC Alive",cex=1,font=1)
  par(mai=c(0.3,0.3,0,0.15))
  boxplot(
    UCEC_an^(1/4), 
    UCEC_at^(1/4),
    names=c("Normal","Tumor")
  )
  par(mai=rep(0,4))
  plot.new()
  text(0.5,0.5,"UCEC Dead",cex=1,font=1)
  par(mai=c(0.3,0.15,0,0.3))
  boxplot(
    UCEC_dn^(1/4), 
    UCEC_dt^(1/4),
    names=c("Normal","Tumor")
  )
  dev.off()
  
  
  pdf(paste0("Plots/UCEC.",output_path_corr,"highly.correlated.counts.density.plot",file_suffix,".pdf"))
  layout(matrix(c(1,2,3,4,1,5,6,4), ncol=2, nrow=4), heights=c(1,1,14,1), widths=c(16,16))
  par(mai=rep(0,4))
  plot.new()
  text(0.5,0.5,paste("Fourth Root Counts of UCEC",corr_type,"Highly Correlated Data",noise_ind),cex=1.5,font=1)
  plot.new()
  text(0.5,0.5,"UCEC Alive",cex=1,font=1)
  par(mai=c(0.3,0.3,0,0.15))
  plot(
    density((UCEC_an)^(1/4)),
    col="blue",
    ylim=c(0,5),
    main=NA
  )
  lines(density((UCEC_at)^(1/4)), col="red")
  par(mai=rep(0,4))
  plot.new()
  legend(x="center", ncol=2,legend=c("Normal","Tumor"),
         fill=c("blue","red"), bty="n")
  plot.new()
  text(0.5,0.5,"UCEC Dead",cex=1,font=1)
  par(mai=c(0.3,0.15,0,0.3))
  plot(
    density((UCEC_dn)^(1/4)),
    col="blue",
    ylim=c(0,5),
    main=NA
  )
  lines(density((UCEC_dt)^(1/4)), col="red")
  dev.off()
}

get_summaries <- function() {
  print(summary(BRCA_an))
  print(summary(BRCA_at))
  print(summary(BRCA_dn))
  print(summary(BRCA_dt))
  print(summary(UCEC_an))
  print(summary(UCEC_at))
  print(summary(UCEC_dn))
  print(summary(UCEC_dt))
}

BRCA_an <- as.numeric(unlist(read.table("Correlation.Distances/BRCA.alive.normal.high.corr.distances.txt", header=FALSE)))
BRCA_at <- as.numeric(unlist(read.table("Correlation.Distances/BRCA.alive.tumor.high.corr.distances.txt", header=FALSE)))
BRCA_dn <- as.numeric(unlist(read.table("Correlation.Distances/BRCA.dead.normal.high.corr.distances.txt", header=FALSE)))
BRCA_dt <- as.numeric(unlist(read.table("Correlation.Distances/BRCA.dead.tumor.high.corr.distances.txt", header=FALSE)))

UCEC_an <- as.numeric(unlist(read.table("Correlation.Distances/UCEC.alive.normal.high.corr.distances.txt", header=FALSE)))
UCEC_at <- as.numeric(unlist(read.table("Correlation.Distances/UCEC.alive.tumor.high.corr.distances.txt", header=FALSE)))
UCEC_dn <- as.numeric(unlist(read.table("Correlation.Distances/UCEC.dead.normal.high.corr.distances.txt", header=FALSE)))
UCEC_dt <- as.numeric(unlist(read.table("Correlation.Distances/UCEC.dead.tumor.high.corr.distances.txt", header=FALSE)))

create_distance_plots("Absolute","")
get_summaries()

BRCA_an <- as.numeric(unlist(read.table("Correlation.Distances/positive.BRCA.alive.normal.high.corr.distances.txt", header=FALSE)))
BRCA_at <- as.numeric(unlist(read.table("Correlation.Distances/positive.BRCA.alive.tumor.high.corr.distances.txt", header=FALSE)))
BRCA_dn <- as.numeric(unlist(read.table("Correlation.Distances/positive.BRCA.dead.normal.high.corr.distances.txt", header=FALSE)))
BRCA_dt <- as.numeric(unlist(read.table("Correlation.Distances/positive.BRCA.dead.tumor.high.corr.distances.txt", header=FALSE)))

UCEC_an <- as.numeric(unlist(read.table("Correlation.Distances/positive.UCEC.alive.normal.high.corr.distances.txt", header=FALSE)))
UCEC_at <- as.numeric(unlist(read.table("Correlation.Distances/positive.UCEC.alive.tumor.high.corr.distances.txt", header=FALSE)))
UCEC_dn <- as.numeric(unlist(read.table("Correlation.Distances/positive.UCEC.dead.normal.high.corr.distances.txt", header=FALSE)))
UCEC_dt <- as.numeric(unlist(read.table("Correlation.Distances/positive.UCEC.dead.tumor.high.corr.distances.txt", header=FALSE)))

create_distance_plots("Positive","")
get_summaries()

BRCA_an <- as.numeric(unlist(read.table("Correlation.Distances/negative.BRCA.alive.normal.high.corr.distances.txt", header=FALSE)))
BRCA_at <- as.numeric(unlist(read.table("Correlation.Distances/negative.BRCA.alive.tumor.high.corr.distances.txt", header=FALSE)))
BRCA_dn <- as.numeric(unlist(read.table("Correlation.Distances/negative.BRCA.dead.normal.high.corr.distances.txt", header=FALSE)))
BRCA_dt <- as.numeric(unlist(read.table("Correlation.Distances/negative.BRCA.dead.tumor.high.corr.distances.txt", header=FALSE)))

UCEC_an <- as.numeric(unlist(read.table("Correlation.Distances/negative.UCEC.alive.normal.high.corr.distances.txt", header=FALSE)))
UCEC_at <- as.numeric(unlist(read.table("Correlation.Distances/negative.UCEC.alive.tumor.high.corr.distances.txt", header=FALSE)))
UCEC_dn <- as.numeric(unlist(read.table("Correlation.Distances/negative.UCEC.dead.normal.high.corr.distances.txt", header=FALSE)))
UCEC_dt <- as.numeric(unlist(read.table("Correlation.Distances/negative.UCEC.dead.tumor.high.corr.distances.txt", header=FALSE)))

create_distance_plots("Negative","")
get_summaries()

BRCA_an <- as.numeric(unlist(read.table("Correlation.Counts/BRCA.alive.normal.high.corr.cnts.txt", header=FALSE)))
BRCA_at <- as.numeric(unlist(read.table("Correlation.Counts/BRCA.alive.tumor.high.corr.cnts.txt", header=FALSE)))
BRCA_dn <- as.numeric(unlist(read.table("Correlation.Counts/BRCA.dead.normal.high.corr.cnts.txt", header=FALSE)))
BRCA_dt <- as.numeric(unlist(read.table("Correlation.Counts/BRCA.dead.tumor.high.corr.cnts.txt", header=FALSE)))

UCEC_an <- as.numeric(unlist(read.table("Correlation.Counts/UCEC.alive.normal.high.corr.cnts.txt", header=FALSE)))
UCEC_at <- as.numeric(unlist(read.table("Correlation.Counts/UCEC.alive.tumor.high.corr.cnts.txt", header=FALSE)))
UCEC_dn <- as.numeric(unlist(read.table("Correlation.Counts/UCEC.dead.normal.high.corr.cnts.txt", header=FALSE)))
UCEC_dt <- as.numeric(unlist(read.table("Correlation.Counts/UCEC.dead.tumor.high.corr.cnts.txt", header=FALSE)))

create_count_plots("Absolute","")
get_summaries()

BRCA_an <- as.numeric(unlist(read.table("Correlation.Counts/negative.BRCA.alive.normal.high.corr.cnts.txt", header=FALSE)))
BRCA_at <- as.numeric(unlist(read.table("Correlation.Counts/negative.BRCA.alive.tumor.high.corr.cnts.txt", header=FALSE)))
BRCA_dn <- as.numeric(unlist(read.table("Correlation.Counts/negative.BRCA.dead.normal.high.corr.cnts.txt", header=FALSE)))
BRCA_dt <- as.numeric(unlist(read.table("Correlation.Counts/negative.BRCA.dead.tumor.high.corr.cnts.txt", header=FALSE)))

UCEC_an <- as.numeric(unlist(read.table("Correlation.Counts/negative.UCEC.alive.normal.high.corr.cnts.txt", header=FALSE)))
UCEC_at <- as.numeric(unlist(read.table("Correlation.Counts/negative.UCEC.alive.tumor.high.corr.cnts.txt", header=FALSE)))
UCEC_dn <- as.numeric(unlist(read.table("Correlation.Counts/negative.UCEC.dead.normal.high.corr.cnts.txt", header=FALSE)))
UCEC_dt <- as.numeric(unlist(read.table("Correlation.Counts/negative.UCEC.dead.tumor.high.corr.cnts.txt", header=FALSE)))

create_count_plots("Positive","")
get_summaries()

BRCA_an <- as.numeric(unlist(read.table("Correlation.Counts/negative.BRCA.alive.normal.high.corr.cnts.txt", header=FALSE)))
BRCA_at <- as.numeric(unlist(read.table("Correlation.Counts/negative.BRCA.alive.tumor.high.corr.cnts.txt", header=FALSE)))
BRCA_dn <- as.numeric(unlist(read.table("Correlation.Counts/negative.BRCA.dead.normal.high.corr.cnts.txt", header=FALSE)))
BRCA_dt <- as.numeric(unlist(read.table("Correlation.Counts/negative.BRCA.dead.tumor.high.corr.cnts.txt", header=FALSE)))

UCEC_an <- as.numeric(unlist(read.table("Correlation.Counts/negative.UCEC.alive.normal.high.corr.cnts.txt", header=FALSE)))
UCEC_at <- as.numeric(unlist(read.table("Correlation.Counts/negative.UCEC.alive.tumor.high.corr.cnts.txt", header=FALSE)))
UCEC_dn <- as.numeric(unlist(read.table("Correlation.Counts/negative.UCEC.dead.normal.high.corr.cnts.txt", header=FALSE)))
UCEC_dt <- as.numeric(unlist(read.table("Correlation.Counts/negative.UCEC.dead.tumor.high.corr.cnts.txt", header=FALSE)))

create_count_plots("Negative","")
get_summaries()

BRCA_an <- as.numeric(unlist(read.table("Correlation.Counts/negative.BRCA.alive.normal.high.corr.cnts.txt", header=FALSE)))
BRCA_at <- as.numeric(unlist(read.table("Correlation.Counts/negative.BRCA.alive.tumor.high.corr.cnts.txt", header=FALSE)))
BRCA_dn <- as.numeric(unlist(read.table("Correlation.Counts/negative.BRCA.dead.normal.high.corr.cnts.txt", header=FALSE)))
BRCA_dt <- as.numeric(unlist(read.table("Correlation.Counts/negative.BRCA.dead.tumor.high.corr.cnts.txt", header=FALSE)))

UCEC_an <- as.numeric(unlist(read.table("Correlation.Counts/negative.UCEC.alive.normal.high.corr.cnts.txt", header=FALSE)))
UCEC_at <- as.numeric(unlist(read.table("Correlation.Counts/negative.UCEC.alive.tumor.high.corr.cnts.txt", header=FALSE)))
UCEC_dn <- as.numeric(unlist(read.table("Correlation.Counts/negative.UCEC.dead.normal.high.corr.cnts.txt", header=FALSE)))
UCEC_dt <- as.numeric(unlist(read.table("Correlation.Counts/negative.UCEC.dead.tumor.high.corr.cnts.txt", header=FALSE)))

create_count_plots("Absolute","")
get_summaries()

BRCA_an <- as.numeric(unlist(read.table("Correlation.Counts/BRCA.alive.normal.high.corr.with.noise.cnts.txt", header=FALSE)))
BRCA_at <- as.numeric(unlist(read.table("Correlation.Counts/BRCA.alive.tumor.high.corr.with.noise.cnts.txt", header=FALSE)))
BRCA_dn <- as.numeric(unlist(read.table("Correlation.Counts/BRCA.dead.normal.high.corr.with.noise.cnts.txt", header=FALSE)))
BRCA_dt <- as.numeric(unlist(read.table("Correlation.Counts/BRCA.dead.tumor.high.corr.with.noise.cnts.txt", header=FALSE)))

#UCEC_an <- as.numeric(unlist(read.table("Correlation.Counts/UCEC.alive.normal.high.corr.with.noise.cnts.txt", header=FALSE)))
#UCEC_at <- as.numeric(unlist(read.table("Correlation.Counts/UCEC.alive.tumor.high.corr.with.noise.cnts.txt", header=FALSE)))
#UCEC_dn <- as.numeric(unlist(read.table("Correlation.Counts/UCEC.dead.normal.high.corr.with.noise.cnts.txt", header=FALSE)))
#UCEC_dt <- as.numeric(unlist(read.table("Correlation.Counts/UCEC.dead.tumor.high.corr.with.noise.cnts.txt", header=FALSE)))

create_count_plots("Absolute",".with.noise")
get_summaries()

