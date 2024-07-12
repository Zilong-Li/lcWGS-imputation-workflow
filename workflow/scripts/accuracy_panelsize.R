## saveRDS(snakemake, snakemake@output[["rds"]])
## q()

## curdir <- "/maps/projects/alab/people/rlk420/quilt2/human/UKBB_GEL_CEU/default5/"
## setwd(curdir)
## snakemake <- readRDS("results/summary/all.accuracy.down0.1x.chr20.rds")

snakemake@source("common.R")

groups <- as.numeric(snakemake@config[["refsize"]]) * 2
refsize0 <- 2 * as.integer(system(paste("bcftools query -l", snakemake@params$vcf, "|", "wc", "-l"), intern = TRUE))
groups[groups == 0] <- refsize0
nd <- length(groups)


truth <- fread(snakemake@input[["truth"]], data.table = F)
ds.truth <- sapply(seq(1, dim(truth)[2] - 1, 2), function(i) {
  rowSums(truth[, (i + 1):(i + 2)])
}) # dosage matrix: nsnps x nsamples
rownames(ds.truth) <- truth[,1]
truth <- truth[,-1] ## remove first id column
rownames(truth) <- rownames(ds.truth)

d.af <- fread(snakemake@input[["af"]], data.table = F)
af <- as.numeric(d.af[, 2])
names(af) <- d.af[, 1]
rm(d.af)
af <- af[!is.na(af)]

dl.quilt1 <- lapply(snakemake@input[["regular"]], parse.imputed.gts2)
dl.quilt2 <- lapply(snakemake@input[["zilong"]], parse.imputed.gts2)
dl.glimpse1 <- lapply(snakemake@input[["glimpse1"]], parse.imputed.gts2)
dl.glimpse2 <- lapply(snakemake@input[["glimpse2"]], parse.imputed.gts2)

bins <- sort(unique(c(
  c(0, 0.01, 0.02 , 0.05 ) / 1e2,
  c(0, 0.01, 0.02 , 0.05 ) / 1e1,
  c(0, 0.01, 0.02 , 0.05 ) / 1e0,
  seq(0.1, 0.5, length.out = 5)
)))

if(refsize0 %/% 1e3 > 500) {
  bins <- sort(unique(c(
    c(0, 0.01, 0.02 , 0.05 ) / 1e4,
    c(0, 0.01, 0.02 , 0.05 ) / 1e3,
    c(0, 0.01, 0.02 , 0.05 ) / 1e2,
    c(0, 0.01, 0.02 , 0.05 ) / 1e1,
    c(0, 0.01, 0.02 , 0.05 ) / 1e0,
    seq(0.1, 0.5, length.out = 5)
  )))
}

if(refsize0 %/% 1e3 < 10) {
  bins <- sort(unique(c(
    c(0, 0.01, 0.02 , 0.05 ) / 1e1,
    c(0, 0.01, 0.02 , 0.05 ) / 1e0,
    seq(0.1, 0.5, length.out = 5)
  )))
}

phasing_errors <- lapply(seq(length(groups)), function(i) {
  n <- ncol(dl.quilt1[[i]])
  quilt2 <- dl.quilt2[[i]][,-seq(3, n, by = 3 )] # get phased genotypes
  quilt1 <- dl.quilt1[[i]][,-seq(3, n, by = 3 )] # get phased genotypes
  glimpse2 <- dl.glimpse2[[i]][,-seq(3, n, by = 3 )] # get phased genotypes
  glimpse1 <- dl.glimpse1[[i]][,-seq(3, n, by = 3 )] # get phased genotypes
  ll <- acc_phasing(truth, quilt2 , glimpse2, quilt1, glimpse1)
  ## names(ll) <- c("QUILT2", "GLIMPSE2", "QUILT1", "GLIMPSE1")
  pse <- lapply(ll, function(ls) {
    lapply(ls, "[[", "pse")
  })
  sites <- lapply(ll, function(ls) {
    lapply(ls, "[[", "sites")
  })
  list(pse = pse, sites = sites)
})


r2_dosage_by_af <- lapply(seq(length(groups)), function(i) {
  n <- ncol(dl.quilt1[[i]])
  quilt2 <- dl.quilt2[[i]][,seq(3, n, by = 3 )] # get dosage
  quilt1 <- dl.quilt1[[i]][,seq(3, n, by = 3 )] # get dosage
  glimpse2 <- dl.glimpse2[[i]][,seq(3, n, by = 3 )] # get dosage
  glimpse1 <- dl.glimpse1[[i]][,seq(3, n, by = 3 )] # get dosage
  d <- acc_r2_by_af(ds.truth, quilt2 , glimpse2, quilt1, glimpse1, af, bins)
  colnames(d) <- c("bin","nsnps","QUILT2", "GLIMPSE2", "QUILT1", "GLIMPSE1")
  d
})

names(r2_dosage_by_af) <- paste0("refsize", as.character(groups))

saveRDS(list(r2_dosage_by_af, phasing_errors), snakemake@output[["rds"]])

## discard sites and re-assign
phasing_errors <- lapply(phasing_errors, function(out) {
  sapply(out[["pse"]], as.numeric)
})


png(paste0(snakemake@output[["rds"]], ".png"), w = 12, h = 6, units='in', res=300)

a1 <- r2_dosage_by_af[[1]]
x <- a1$bin[!sapply(a1[, 2], is.na)]
x <- log10(as.numeric(x))
labels <- 100 * bins[-1]
labels <- labels[!sapply(a1[, 2], is.na)]
ymin <- min(sapply(r2_dosage_by_af, function(d) {
  m <- as.matrix(apply(d[, -1], 2, unlist))
  min(m, na.rm = T)
}))

par(mfrow = c(1, 2))

plot(1, col = "transparent", axes = F, xlim = c(min(x), max(x)), ylim = c(0.9 * ymin, 1.0), ylab = "Aggregated R2 within each MAF bin", xlab = "Minor Allele Frequency")
nd <- length(groups)
for (i in 1:nd) {
  d <- r2_dosage_by_af[[i]]
  y <- rmna(d$QUILT2)
  lines(x, y, type = "l", lty = nd - i + 1, pch = 1, col = mycols[1])
  y <- rmna(d$GLIMPSE2)
  lines(x, y, type = "l", lty = nd - i + 1, pch = 1, col = mycols[2])
  y <- rmna(d$QUILT1)
  lines(x, y, type = "l", lty = nd - i + 1, pch = 1, col = mycols[3])
  y <- rmna(d$GLIMPSE1)
  lines(x, y, type = "l", lty = nd - i + 1, pch = 1, col = mycols[4])
}
axis(side = 1, at = x, labels = labels)
axis(side = 2, at = seq(0, 1, 0.2))
legend("bottomright", legend = paste0("N=", groups), lty = nd:1, bty = "n")
legend("topleft", legend = c("QUILT2", "GLIMPSE2", "QUILT1", "GLIMPSE1"), col = mycols, pch = 1, lwd = 1.5, cex = 1.0, xjust = 0, yjust = 1, bty = "n")


## boxplot(phasing_errors[[1]], col = mycols[1:4], ylab = "PSE %", main = paste("Ref Panel Size: N=", groups[1]))

boxplot(phasing_errors[[1]], ylab = "PSE %", main = paste("Ref Panel Size: N=", groups[1]))
nsamples <- nrow(phasing_errors[[1]])
for(i in 1:4) {
  vals <- phasing_errors[[1]][,i]
  j <- jitter(rep(i, nsamples), amount=1/4)
  points(j,  vals,  col = mycols[i],pch = 20)
}

dev.off()
