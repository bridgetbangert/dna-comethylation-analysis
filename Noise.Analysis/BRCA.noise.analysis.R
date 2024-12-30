setwd("/home/bmb191/math5376D/Final.Project/Output.Files")

# --------------------------------------------------------
# Creating the plots comparing the high correlation counts
# --------------------------------------------------------
# Reading in the data
# --------------------------------------------------------

BRCA_an <- as.numeric(unlist(read.table("Correlation.Counts/BRCA.alive.normal.high.corr.cnts.txt",header=FALSE)))
BRCA_an_noise <- as.numeric(unlist(read.table("Correlation.Counts/BRCA.alive.normal.high.corr.with.noise.cnts.txt",header=FALSE)))
BRCA_at <- as.numeric(unlist(read.table("Correlation.Counts/BRCA.alive.tumor.high.corr.cnts.txt",header=FALSE)))
BRCA_at_noise <- as.numeric(unlist(read.table("Correlation.Counts/BRCA.alive.tumor.high.corr.with.noise.cnts.txt",header=FALSE)))
BRCA_dn <- as.numeric(unlist(read.table("Correlation.Counts/BRCA.dead.normal.high.corr.cnts.txt",header=FALSE)))
BRCA_dn_noise <- as.numeric(unlist(read.table("Correlation.Counts/BRCA.dead.normal.high.corr.with.noise.cnts.txt",header=FALSE)))
BRCA_dt <- as.numeric(unlist(read.table("Correlation.Counts/BRCA.dead.tumor.high.corr.cnts.txt",header=FALSE)))
BRCA_dt_noise <- as.numeric(unlist(read.table("Correlation.Counts/BRCA.dead.tumor.high.corr.with.noise.cnts.txt",header=FALSE)))

# --------------------------------------------------------
# Creating the pdf file
# --------------------------------------------------------
pdf("Plots/BRCA.data.high.corr.cnts.comparison.with.added.noise.pdf",16,16)
par(mfrow=c(2,2))

# --------------------------------------------------------
# Alive normal data
# --------------------------------------------------------
plot(
  density((BRCA_an)^(1/4)),
  col="red",
  xlim=c(0,4),
  ylim=c(0,3),
  main=" "
)
title("Alive Normal", adj=0.5, line=1, cex.main=0.75, outer=FALSE)
legend(
  "topright",
  legend=c("Without noise","With noise"),
  col=c("red","blue"),
  cex=0.7,
  lty=1,
  lwd=2
)
lines(density((BRCA_an_noise)^(1/4)),col="blue")
# --------------------------------------------------------
# Alive tumor data
# --------------------------------------------------------
plot(density(BRCA_at^(1/4)),
     col="red",
     main=" ",
     xlim=c(0,4),
     ylim=c(0,3)
)
title("Alive Tumor", adj=0.5, line=1, cex.main=0.75, outer=FALSE)
legend(
  "topright",
  legend=c("Without noise","With noise"),
  col=c("red","blue"),
  cex=0.7,
  lty=1,
  lwd=2
)
lines(density(BRCA_at_noise^(1/4)),col="blue")
# --------------------------------------------------------
# Dead normal data
# --------------------------------------------------------
plot(density(BRCA_dn_noise^(1/4)),
     col="blue",
     main=" ",
     xlim=c(0,4),
     ylim=c(0,3)
)
title("Dead Normal", adj=0.5, line=1, cex.main=0.75, outer=FALSE)
legend(
  "topright",
  legend=c("Without noise","With noise"),
  col=c("red","blue"),
  cex=0.5,
  lty=1,
  lwd=2
)
lines(density(BRCA_dn^(1/4)),col="red")
# --------------------------------------------------------
# Dead tumor data
# --------------------------------------------------------
plot(density(BRCA_dt_noise^(1/4)),
     col="blue",
     main=" ",
     xlim=c(0,4),
     ylim=c(0,3)
)
title("Dead Tumor", adj=0.5, line=1, cex.main=0.75, outer=FALSE)
legend(
  "topright",
  legend=c("Without noise","With noise"),
  col=c("red","blue"),
  cex=0.5,
  lty=1,
  lwd=2
)
lines(density(BRCA_dt^(1/4)),col="red")
title("BRCA High Correlation Counts With and Without Noise (scaled by fourth root)", line = -2, outer = TRUE)
dev.off()

# --------------------------------------------------------
# Print summaries of the counts
# --------------------------------------------------------
summary(BRCA_an)
summary(BRCA_an_noise)
summary(BRCA_at)
summary(BRCA_at_noise)
summary(BRCA_dn)
summary(BRCA_dn_noise)
summary(BRCA_dt)
summary(BRCA_dt_noise)

# --------------------------------------------------------
# Creating the plots comparing the distances of sites with
# high correlation counts
# --------------------------------------------------------
# Reading in the data
# --------------------------------------------------------
BRCA_an <- as.numeric(unlist(read.table("Correlation.Distances/BRCA.alive.normal.high.corr.distances.txt",header=FALSE)))
BRCA_an_noise <- as.numeric(unlist(read.table("Correlation.Distances/BRCA.alive.normal.high.corr.with.noise.distances.txt",header=FALSE)))
BRCA_at <- as.numeric(unlist(read.table("Correlation.Distances/BRCA.alive.tumor.high.corr.distances.txt",header=FALSE)))
BRCA_at_noise <- as.numeric(unlist(read.table("Correlation.Distances/BRCA.alive.tumor.high.corr.with.noise.distances.txt",header=FALSE)))
BRCA_dn <- as.numeric(unlist(read.table("Correlation.Distances/BRCA.dead.normal.high.corr.distances.txt",header=FALSE)))
BRCA_dn_noise <- as.numeric(unlist(read.table("Correlation.Distances/BRCA.dead.normal.high.corr.with.noise.distances.txt",header=FALSE)))
BRCA_dt <- as.numeric(unlist(read.table("Correlation.Distances/BRCA.dead.tumor.high.corr.distances.txt",header=FALSE)))
BRCA_dt_noise <- as.numeric(unlist(read.table("Correlation.Distances/BRCA.dead.tumor.high.corr.with.noise.distances.txt",header=FALSE)))

max_x <- max(BRCA_an, BRCA_an_noise, BRCA_at, BRCA_at_noise, BRCA_dn, BRCA_dn_noise, BRCA_dt, BRCA_dt_noise)
max_y <- 2*max(density(c(BRCA_an, BRCA_an_noise, BRCA_at, BRCA_at_noise, BRCA_dn, BRCA_dn_noise, BRCA_dt, BRCA_dt_noise))$y)
# --------------------------------------------------------
# Creating the pdf file
# --------------------------------------------------------
pdf("Plots/BRCA.data.high.corr.distances.comparison.with.added.noise.pdf",16,16)
par(mfrow=c(2,2))

# --------------------------------------------------------
# Alive normal
# --------------------------------------------------------
plot(
  density(BRCA_an_noise),
  col="blue",
  xlim=c(0,max_x),
  ylim=c(0,max_y),
  main=" "
)
title("Alive Normal", adj=0.5, line=1, cex.main=0.75, outer=FALSE)
legend(
  "topright",
  legend=c("Without noise","With noise"),
  col=c("red","blue"),
  cex=0.7,
  lty=1,
  lwd=2
)
lines(density(BRCA_an),col="red")
# --------------------------------------------------------
# Alive tumor
# --------------------------------------------------------
plot(density(BRCA_at),
     col="red",
     main=" ",
     xlim=c(0,max_x),
     ylim=c(0,max_y)
)
title("Alive Tumor", adj=0.5, line=1, cex.main=0.75, outer=FALSE)
legend(
  "topright",
  legend=c("Without noise","With noise"),
  col=c("red","blue"),
  cex=0.7,
  lty=1,
  lwd=2
)
lines(density(BRCA_at_noise),col="blue")
# --------------------------------------------------------
# Dead normal
# --------------------------------------------------------
plot(density(BRCA_dn_noise),
     col="blue",
     main=" ",
     xlim=c(0,max_x),
     ylim=c(0,max_y)
)
title("Dead Normal", adj=0.5, line=1, cex.main=0.75, outer=FALSE)
legend(
  "topright",
  legend=c("Without noise","With noise"),
  col=c("red","blue"),
  cex=0.5,
  lty=1,
  lwd=2
)
lines(density(BRCA_dn),col="red")
# --------------------------------------------------------
# Dead tumor
# --------------------------------------------------------
plot(density(BRCA_dt_noise),
     col="blue",
     main=" ",
     xlim=c(0,max_x),
     ylim=c(0,max_y)
)
title("Dead Tumor", adj=0.5, line=1, cex.main=0.75, outer=FALSE)
legend(
  "topright",
  legend=c("Without noise","With noise"),
  col=c("red","blue"),
  cex=0.5,
  lty=1,
  lwd=2
)
lines(density(BRCA_dt),col="red")

title("BRCA High Correlation Distances With and Without Noise", line = -2, outer = TRUE)
dev.off()

summary(BRCA_an)
summary(BRCA_an_noise)
summary(BRCA_at)
summary(BRCA_at_noise)
summary(BRCA_dn)
summary(BRCA_dn_noise)
summary(BRCA_dt)
summary(BRCA_dt_noise)

# --------------------------------------------------------
# Creating the plots comparing the distances of sites with
# high correlation counts scaled by 4th root
# --------------------------------------------------------
max_x <- max(BRCA_an^(1/4), BRCA_an_noise^(1/4), BRCA_at^(1/4), BRCA_at_noise^(1/4), BRCA_dn^(1/4), BRCA_dn_noise^(1/4), BRCA_dt^(1/4), BRCA_dt_noise^(1/4))
max_y <- max(density(c(BRCA_an^(1/4), BRCA_an_noise^(1/4), BRCA_at^(1/4), BRCA_at_noise^(1/4), BRCA_dn^(1/4), BRCA_dn_noise^(1/4), BRCA_dt^(1/4), BRCA_dt_noise^(1/4)))$y)

# --------------------------------------------------------
# Creating the pdf file
# --------------------------------------------------------
pdf("Plots/BRCA.data.high.corr.distances.comparison.with.added.noise.scaled.pdf",16,16)
par(mfrow=c(2,2))

# --------------------------------------------------------
# Alive normal
# --------------------------------------------------------
lines(
  density(BRCA_an_noise^1/2),
  col="blue",
  #xlim=c(0,max_x),
  #ylim=c(0,max_y),
  main=" "
)
title("Alive Normal", adj=0.5, line=1, cex.main=0.75, outer=FALSE)
legend(
  "topright",
  legend=c("Without noise","With noise"),
  col=c("red","blue"),
  cex=0.7,
  lty=1,
  lwd=2
)
plot(density(BRCA_an^(1/2)),col="red")
# --------------------------------------------------------
# Alive tumor
# --------------------------------------------------------
plot(density(BRCA_at^(1/2)),
     col="red",
     main=" ",
     #xlim=c(0,max_x),
     #ylim=c(0,max_y)
)
title("Alive Tumor", adj=0.5, line=1, cex.main=0.75, outer=FALSE)
legend(
  "topright",
  legend=c("Without noise","With noise"),
  col=c("red","blue"),
  cex=0.7,
  lty=1,
  lwd=2
)
lines(density(BRCA_at_noise^(1/2)),col="blue")
# --------------------------------------------------------
# Dead normal
# --------------------------------------------------------
plot(density(BRCA_dn_noise^(1/4)),
     col="blue",
     main=" ",
     #xlim=c(0,max_x),
     #ylim=c(0,max_y)
)
title("Dead Normal", adj=0.5, line=1, cex.main=0.75, outer=FALSE)
legend(
  "topright",
  legend=c("Without noise","With noise"),
  col=c("red","blue"),
  cex=0.5,
  lty=1,
  lwd=2
)
lines(density(BRCA_dn^(1/4)),col="red")
# --------------------------------------------------------
# Dead tumor
# --------------------------------------------------------
plot(density(BRCA_dt_noise^(1/4)),
     col="blue",
     main=" ",
     #xlim=c(0,max_x),
     #ylim=c(0,max_y)
)
title("Dead Tumor", adj=0.5, line=1, cex.main=0.75, outer=FALSE)
legend(
  "topright",
  legend=c("Without noise","With noise"),
  col=c("red","blue"),
  cex=0.5,
  lty=1,
  lwd=2
)
lines(density(BRCA_dt^(1/4)),col="red")

title("BRCA High Correlation Distances With and Without Noise (scaled by 4th root)", line = -2, outer = TRUE)
dev.off()
