setwd("/home/bmb191/math5376D/Final.Project/Output.Files")

# ----------------------------------------------------

brca_a_t <- read.table("Correlation.Distances/BRCA.alive.tumor.high.corr.distances.txt",header=FALSE)
brca_a_n <- read.table("Correlation.Distances/BRCA.alive.normal.high.corr.distances.txt",header=FALSE)
brca_d_t <- read.table("Correlation.Distances/BRCA.dead.tumor.high.corr.distances.txt",header=FALSE)
brca_d_n <- read.table("Correlation.Distances/BRCA.dead.normal.high.corr.distances.txt",header=FALSE)

ucec_a_t <- read.table("Correlation.Distances/UCEC.alive.tumor.high.corr.distances.txt",header=FALSE)
ucec_a_n <- read.table("Correlation.Distances/UCEC.alive.normal.high.corr.distances.txt",header=FALSE)
ucec_d_t <- read.table("Correlation.Distances/UCEC.dead.tumor.high.corr.distances.threshold.0.9.txt",header=FALSE)
ucec_d_n <- read.table("Correlation.Distances/UCEC.dead.normal.high.corr.distances.threshold.0.8.txt",header=FALSE)

# ----------------------------------------------------

brca_corr_cnts_at <- read.table("Correlation.Counts/BRCA.alive.tumor.high.corr.cnts.txt",header=FALSE)
brca_corr_cnts_an <- read.table("Correlation.Counts/BRCA.alive.normal.high.corr.cnts.txt",header=FALSE)
brca_corr_cnts_dt <- read.table("Correlation.Counts/BRCA.dead.tumor.high.corr.cnts.txt",header=FALSE)
brca_corr_cnts_dn <- read.table("Correlation.Counts/BRCA.dead.normal.high.corr.cnts.txt",header=FALSE)

ucec_corr_cnts_at <- read.table("Correlation.Counts/UCEC.alive.tumor.high.corr.cnts.txt",header=FALSE)
ucec_corr_cnts_an <- read.table("Correlation.Counts/UCEC.alive.normal.high.corr.cnts.txt",header=FALSE)
ucec_corr_cnts_dt <- read.table("Correlation.Counts/UCEC.dead.tumor.high.corr.cnts.0.9.threshold.txt",header=FALSE)
ucec_corr_cnts_dn <- read.table("Correlation.Counts/UCEC.dead.normal.high.corr.cnts.0.8.threshold.txt",header=FALSE)

# ----------------------------------------------------

brca_a_t <- apply(brca_a_t, 2, as.numeric)
brca_a_n <- apply(brca_a_n, 2, as.numeric)
brca_d_t <- apply(brca_d_t, 2, as.numeric)
brca_d_n <- apply(brca_d_n, 2, as.numeric)

ucec_a_t <- apply(ucec_a_t, 2, as.numeric)
ucec_a_n <- apply(ucec_a_n, 2, as.numeric)
ucec_d_t <- apply(ucec_d_t, 2, as.numeric)
ucec_d_n <- apply(ucec_d_n, 2, as.numeric)

# ----------------------------------------------------

brca_corr_cnts_at <- apply(brca_corr_cnts_at, 2, as.numeric)
brca_corr_cnts_an <- apply(brca_corr_cnts_an, 2, as.numeric)
brca_corr_cnts_dt <- apply(brca_corr_cnts_dt, 2, as.numeric)
brca_corr_cnts_dn <- apply(brca_corr_cnts_dn, 2, as.numeric)

ucec_corr_cnts_at <- apply(ucec_corr_cnts_at, 2, as.numeric)
ucec_corr_cnts_an <- apply(ucec_corr_cnts_an, 2, as.numeric)
ucec_corr_cnts_dt <- apply(ucec_corr_cnts_dt, 2, as.numeric)
ucec_corr_cnts_dn <- apply(ucec_corr_cnts_dn, 2, as.numeric)

# ----------------------------------------------------

sum(brca_corr_cnts_an)
sum(brca_corr_cnts_at)
sum(brca_corr_cnts_dn)
sum(brca_corr_cnts_dt)

sum(ucec_corr_cnts_an)
sum(ucec_corr_cnts_at)
sum(ucec_corr_cnts_dn)
sum(ucec_corr_cnts_dt)

# ----------------------------------------------------

summary(brca_a_n)
summary(brca_a_t)
summary(brca_d_n)
summary(brca_d_t)

summary(ucec_a_n)
summary(ucec_a_t)
summary(ucec_d_n)
summary(ucec_d_t)

# ----------------------------------------------------

pdf("Plots/BRCA.alive.high.corr.distance.density.plot.pdf")
plot(density(brca_a_t^(1/4)),
    xlim=c(0,120),
    ylim=c(0,0.035),
    col="red",
    main="Distances between highly correlated CG sites (BRCA Alive)"
)
legend("topright", legend=c("Normal","Tumor"), fill=c("blue","red"))
lines(density(brca_a_n^(1/4)),
  xlim=c(0,2.5),
  ylim=c(0,3),
  col="blue"
)
dev.off()

pdf("Plots/UCEC.alive.high.corr.distance.density.plot.pdf")
plot(density(ucec_a_t^(1/4)),
    xlim=c(0,120),
    ylim=c(0,0.035),
    col="red",
    main="Distances between highly correlated CG sites (UCEC Alive)"
)
legend("topright", legend=c("Normal","Tumor"), fill=c("blue","red"))
lines(density(ucec_a_n^(1/4)),
  xlim=c(0,2.5),
  ylim=c(0,3),
  col="blue"
)
dev.off()

# ----------------------------------------------------

pdf("Plots/BRCA.dead.high.corr.distance.density.plot.pdf")
plot(density(brca_d_t^(1/4)),
  xlim=c(0,120),
  ylim=c(0,0.035),
  col="red",
  main="Fourth Root Distances between highly correlated CG sites (Dead)"
)
legend("topright", legend=c("Normal","Tumor"), fill=c("blue","red"))
lines(density(brca_d_n^(1/4)),
  xlim=c(0,2.5),
  ylim=c(0,3),
  col="blue"
)
dev.off()

pdf("Plots/UCEC.dead.high.corr.distance.density.plot.pdf")
plot(density(ucec_d_t^(1/4)),
  xlim=c(0,120),
  ylim=c(0,0.035),
  col="red",
  main="Fourth Root Distances between highly correlated CG sites (Dead)"
)
legend("topright", legend=c("Normal","Tumor"), fill=c("blue","red"))
lines(density(ucec_d_n^(1/4)),
  xlim=c(0,2.5),
  ylim=c(0,3),
  col="blue"
)
dev.off()

# ----------------------------------------------------

brca_a_t <- read.table("Correlation.Distances/BRCA.alive.tumor.high.corr.distances.txt",header=FALSE)
brca_a_n <- read.table("Correlation.Distances/BRCA.alive.normal.high.corr.distances.txt",header=FALSE)
brca_d_t <- read.table("Correlation.Distances/BRCA.dead.normal.high.corr.distances.txt",header=FALSE)
brca_d_n <- read.table("Correlation.Distances/BRCA.dead.tumor.high.corr.distances.txt",header=FALSE)

ucec_a_t <- read.table("Correlation.Distances/UCEC.alive.tumor.high.corr.distances.txt",header=FALSE)
ucec_a_n <- read.table("Correlation.Distances/UCEC.alive.normal.high.corr.distances.txt",header=FALSE)
ucec_d_t <- read.table("Correlation.Distances/UCEC.dead.normal.high.corr.distances.txt",header=FALSE)
ucec_d_n <- read.table("Correlation.Distances/UCEC.dead.tumor.high.corr.distances.txt",header=FALSE)

# ----------------------------------------------------

brca_a_n <- cbind(as.numeric(as.character(brca_a_n$V1)), "Normal")
brca_a_t <- cbind(as.numeric(as.character(brca_a_t$V1)), "Tumor")

ucec_a_n <- cbind(as.numeric(as.character(ucec_a_n$V1)), "Normal")
ucec_a_t <- cbind(as.numeric(as.character(ucec_a_t$V1)), "Tumor")

brca_a_dist_df <- data.frame(rbind(brca_a_n, brca_a_t))
colnames(brca_a_dist_df) <- c("Distance", "Group")

ucec_a_dist_df <- data.frame(rbind(ucec_a_n, ucec_a_t))
colnames(ucec_a_dist_df) <- c("Distance", "Group")

brca_a_dist_df$Distance <- as.numeric(as.character(brca_a_dist_df$Distance))
brca_a_dist_df$Group <- as.factor(brca_a_dist_df$Group)

ucec_a_dist_df$Distance <- as.numeric(as.character(ucec_a_dist_df$Distance))
ucec_a_dist_df$Group <- as.factor(ucec_a_dist_df$Group)

pdf("Plots/BRCA.alive.high.corr.distance.boxplot.pdf")
boxplot(
  Distance ~ Group,
  data = brca_a_dist_df,
  main = "Distances between highly correlated CG sites (BRCA Alive)"
)
dev.off()

pdf("Plots/UCEC.alive.high.corr.distance.boxplot.pdf")
boxplot(
  Distance ~ Group,
  data = ucec_a_dist_df,
  main = "Distances between highly correlated CG sites (UCEC Alive)"
)
dev.off()

# ----------------------------------------------------

brca_d_n <- cbind(as.numeric(as.character(brca_d_n$V1)), "Normal")
brca_d_t <- cbind(as.numeric(as.character(brca_d_t$V1)), "Tumor")

brca_d_dist_df <- data.frame(rbind(brca_d_n, brca_d_t))
colnames(brca_d_dist_df) <- c("Distance", "Group")

brca_d_dist_df$Distance <- as.numeric(as.character(brca_d_dist_df$Distance))
brca_d_dist_df$Group <- as.factor(brca_d_dist_df$Group)

pdf("Plots/BRCA.dead.high.corr.distance.boxplot.pdf")
boxplot(
  Distance ~ Group,
  data = brca_d_dist_df,
  main = "Distances between highly correlated CG sites (BRCA Dead)"
)
dev.off()

ucec_d_n <- cbind(as.numeric(as.character(ucec_d_n$V1)), "Normal")
ucec_d_t <- cbind(as.numeric(as.character(ucec_d_t$V1)), "Tumor")

ucec_d_dist_df <- data.frame(rbind(ucec_d_n, ucec_d_t))
colnames(ucec_d_dist_df) <- c("Distance", "Group")

ucec_d_dist_df$Distance <- as.numeric(as.character(ucec_d_dist_df$Distance))
ucec_d_dist_df$Group <- as.factor(ucec_d_dist_df$Group)

pdf("Plots/UCEC.dead.high.corr.distance.boxplot.pdf")
boxplot(
  Distance ~ Group,
  data = ucec_d_dist_df,
  main = "Distances between highly correlated CG sites (UCEC Dead)"
)
dev.off()

# ----------------------------------------------------

pdf("Plots/BRCA.alive.high.corr.counts.density.plot.pdf")
plot(density(brca_corr_cnts_an^(1/4)),
  xlim=c(-.25,3),
  ylim=c(0,3),
  col="red",
  main="Fourth Root Counts of highly correlated CG sites (BRCA Alive)"
)
legend("topright", legend=c("Normal","Tumor"), fill=c("red","blue"))
lines(density(brca_corr_cnts_at^(1/4)), col="blue")
dev.off()

pdf("Plots/UCEC.alive.high.corr.counts.density.plot.pdf")
plot(density(ucec_corr_cnts_an^(1/4)),
  xlim=c(-.25,3),
  ylim=c(0,3),
  col="red",
  main="Fourth Root Counts of highly correlated CG sites (UCEC Alive)"
)
legend("topright", legend=c("Normal","Tumor"), fill=c("red","blue"))
lines(density(ucec_corr_cnts_at^(1/4)), col="blue")
dev.off()

# ----------------------------------------------------

pdf("Plots/BRCA.dead.high.corr.counts.density.plot.pdf")
plot(
  density(brca_corr_cnts_dn^(1/4)),
  xlim=c(-.25,3),
  ylim=c(0,3),
  col="red",
  main="Fourth root counts of highly correlated CG sites (BRCA Dead)"
)
legend("topright", legend=c("Normal","Tumor"), fill=c("red","blue"))
lines(density(brca_corr_cnts_dt^(1/4)),col="blue")
dev.off()

pdf("Plots/UCEC.dead.high.corr.counts.density.plot.pdf")
plot(
  density(ucec_corr_cnts_dn^(1/4)),
  xlim=c(-.25,3),
  ylim=c(0,3),
  col="red",
  main="Fourth root counts of highly correlated CG sites (Dead)"
)
legend("topright", legend=c("Normal","Tumor"), fill=c("red","blue"))
lines(density(ucec_corr_cnts_dt^(1/4)),col="blue")
dev.off()

# ----------------------------------------------------

brca_noise_corr_an <- read.table("Correlation.Distances/BRCA.alive.normal.high.corr.with.noise.distances.txt",header=FALSE)
brca_noise_dist_at <- read.table("Correlation.Distances/BRCA.alive.tumor.high.corr.with.noise.distances.txt",header=FALSE)
brca_noise_dist_dt <- read.table("Correlation.Distances/BRCA.dead.tumor.high.corr.with.noise.distances.txt",header=FALSE)
brca_noise_dist_dn <- read.table("Correlation.Distances/BRCA.dead.normal.high.corr.with.noise.distances.txt",header=FALSE)

brca_noise_corr_dist_an <- read.table("Correlation.Distances/BRCA.alive.normal.high.corr.with.noise.distances.txt",header=FALSE)
brca_noise_corr_dist_at <- read.table("Correlation.Distances/BRCA.alive.tumor.high.corr.with.noise.distances.txt",header=FALSE)
brca_noise_corr_dist_dt <- read.table("Correlation.Distances/BRCA.dead.tumor.high.corr.with.noise.distances.txt",header=FALSE)
brca_noise_corr_dist_dn <- read.table("Correlation.Distances/BRCA.dead.normal.high.corr.with.noise.distances.txt",header=FALSE)

# ----------------------------------------------------

brca_noise_an <- cbind(as.numeric(as.character(brca_noise_corr_dist_an$V1)), "Normal")
brca_noise_at <- cbind(as.numeric(as.character(brca_noise_corr_dist_at$V1)), "Tumor")

brca_noise_a_dist_df <- data.frame(rbind(brca_noise_an, brca_noise_at))
colnames(brca_noise_a_dist_df) <- c("Distance", "Group")

brca_noise_a_dist_df$Distance <- as.numeric(as.character(brca_noise_a_dist_df$Distance))
brca_noise_a_dist_df$Group <- as.factor(brca_noise_a_dist_df$Group)

pdf("Plots/BRCA.cg.site.distances.boxplot.with.noise.pdf")
boxplot(
  Distance ~ Group,
  data = brca_noise_a_dist_df,
  main = "Distances between all highly correlated CG sites with added noise (BRCA Alive)"
)

brca_noise_dn <- cbind(as.numeric(as.character(brca_noise_corr_dist_dn$V1)), "Normal")
brca_noise_dt <- cbind(as.numeric(as.character(brca_noise_corr_dist_dt$V1)), "Tumor")

brca_noise_d_dist_df <- data.frame(rbind(brca_noise_dn, brca_noise_dt))
colnames(brca_noise_d_dist_df) <- c("Distance", "Group")

brca_noise_d_dist_df$Distance <- as.numeric(as.character(brca_noise_d_dist_df$Distance))
brca_noise_d_dist_df$Group <- as.factor(brca_noise_d_dist_df$Group)

boxplot(
  Distance ~ Group,
  data = brca_noise_d_dist_df,
  main = "Distances between all highly correlated CG sites with added noise (BRCA dead)"
)
dev.off()

# ----------------------------------------------------

summary(brca_noise_corr_dist_an)
summary(brca_noise_corr_dist_at)
summary(brca_noise_corr_dist_dt)
summary(brca_noise_corr_dist_dn)

# ----------------------------------------------------

pdf("Plots/BRCA.alive.cg.site.distances.density.plots.with.noise.pdf")
plot(
  density(brca_noise_an),
  col="red",
  main="Distribution of all highly correlated CG site distances with added noise (BRCA alive)"
)
legend("topright", legend=c("Normal","Tumor"), fill=c("red","blue"))
lines(density(brca_noise_at),col="blue")
dev.off()

pdf("Plots/BRCA.dead.cg.site.distances.density.plots.with.noise.pdf")
plot(
  density(brca_noise_dn),
  col="red",
  main="Distribution of all highly correlated CG site distances with added noise (BRCA dead)"
)
legend("topright", legend=c("Normal","Tumor"), fill=c("red","blue"))
lines(density(brca_noise_dt),col="blue")
dev.off()




