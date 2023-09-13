snakemake@source("common.R")

library(vcfppR)
library(parallel)

depths <- as.numeric(snakemake@config[["downsample"]])
refsize0 <- 2 * as.integer(system(paste("bcftools query -l", snakemake@params$vcf, "|", "wc", "-l"), intern = TRUE))

samples <- snakemake@params[["samples"]]
chrom <- snakemake@wildcards[["chrom"]]

truth <- tableGT(snakemake@params[["truth"]], chrom, samples)
truth$id <- paste(truth$chr, truth$pos, truth$ref, truth$alt, sep=":")
## make sure the truth$samples order is same as samples
gt.truth <- do.call(rbind, truth$gt)
n <- ncol(gt.truth)
truth$gt <- gt.truth[, seq(1, n, 2)] + gt.truth[, seq(2, n, 2)]

getGT <- function(vcffile, region, samples) {
  quilt2 <- tableGT(vcffile, region, samples)
  quilt2$id <- paste(quilt2$chr, quilt2$pos, quilt2$ref, quilt2$alt, sep=":")
  gt.quilt2 <- do.call(rbind, quilt2$gt)
  n <- ncol(gt.quilt2)
  quilt2$gt <- gt.quilt2[, seq(1, n, 2)] + gt.quilt2[, seq(2, n, 2)]
  quilt2
}

dl.quilt2 <- mclapply(snakemake@input[["zilong"]], getGT, mc.cores = 4,region=chrom, samples=samples)
dl.glimpse2 <- mclapply(snakemake@input[["glimpse2"]], getGT, mc.cores = 4,region=chrom, samples=samples)

f1.quilt2 <- lapply(dl.quilt2, function(test) acc_f1(truth, test))
names(f1.quilt2) <- depths
f1.glimpse2 <- lapply(dl.glimpse2, function(test) acc_f1(truth, test))
names(f1.glimpse2) <- depths

res <- list(f1.quilt2, f1.glimpse2)
names(res) <- c("QUILT2", "GLIMPSE2")
saveRDS(res, snakemake@output[["rds"]])

nsamples <- ncol(truth$gt)

out <- do.call(rbind, lapply(names(res), function(name) {
  v <- unlist(res[[name]])
  d <- data.frame(v = v, m=rep(name, length(v)), d=rep(depths, each=nsamples))
  d
}))

mycols <- c(QUILT2 = "#e69f00", GLIMPSE2 = "#d55e00", QUILT1 = "#56b4e9", GLIMPSE1 = "#cc79a7")
mycols <- c(QUILT2 = "#e69f00", GLIMPSE2 = "#d55e00")

out <- out[order(out$d, out$m),]

pdf(paste0(snakemake@output[["rds"]], ".pdf"), h=7,w=7)
a <- c(seq(2), seq(4,5), seq(7,8), seq(10,11))
boxplot(v ~ m + d, data=out, col=mycols[unique(out$m)],
        at = a, cex.lab=2.0, cex.axis = 1.5,
        xaxt="n", xlab = "Coverage", ylab = "F1-score")
b <- c(1.5, 4.5, 7.5, 10.5)
axis(side = 1, at = b, labels = unique(out$d), cex.axis = 2.0)
legend("bottomright", legend = names(mycols), fill = mycols, bty = "n", cex=2)
dev.off()


q()

groups <- as.numeric(snakemake@config[["downsample"]])
refsize0 <- 2 * as.integer(system(paste("bcftools query -l", snakemake@params$vcf, "|", "wc", "-l"), intern = TRUE))


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

dl.quilt2 <- lapply(snakemake@input[["zilong"]], parse.imputed.gts2)
dl.glimpse2 <- lapply(snakemake@input[["glimpse2"]], parse.imputed.gts2)

bins <- sort(unique(c(
  c(0, 0.01, 0.02 , 0.05 ) / 1e2,
  c(0, 0.01, 0.02 , 0.05 ) / 1e1,
  c(0, 0.01, 0.02 , 0.05 ) / 1e0,
  seq(0.1, 0.5, length.out = 5)
)))


# d0:truth, d1: quilt2, d2:glimpse2
acc_r2_by_af <- function(d0, d1, d2, af, bins, flip = TRUE, per_snp = FALSE) {
  id <- (intersect(intersect(rownames(d0), rownames(d1)), rownames(d2)))
  id <- intersect(id, names(af))
  res1 <- r2_by_freq(breaks = bins, af, truthG = d0, testDS = d1, which_snps = id, flip = flip, per_snp = per_snp)
  res2 <- r2_by_freq(breaks = bins, af, truthG = d0, testDS = d2, which_snps = id, flip = flip, per_snp = per_snp)
  as.data.frame(cbind(bin = bins[-1], nsnps = res1[, "n"], quilt2 = res1[, "simple"], glimpse2 = res2[, "simple"]))
}

r2_dosage_by_af <- lapply(seq(length(groups)), function(i) {
  n <- ncol(dl.quilt2[[i]])
  quilt2 <- dl.quilt2[[i]][,seq(3, n, by = 3 )] # get dosage
  glimpse2 <- dl.glimpse2[[i]][,seq(3, n, by = 3 )] # get dosage
  d <- acc_r2_by_af(ds.truth, quilt2 , glimpse2, af, bins, flip = TRUE)
  colnames(d) <- c("bin","nsnps","QUILT2", "GLIMPSE2")
  d
})
names(r2_dosage_by_af) <- paste0(as.character(groups), "x")

saveRDS(r2_dosage_by_af, snakemake@output[["rds"]])

pdf(paste0(snakemake@output[["rds"]], ".pdf"), w = 6, h = 6)

a1 <- r2_dosage_by_af[[1]]
x <- a1$bin[!sapply(a1[, 2], is.na)] # remove AF bin with NULL results
x <- log10(as.numeric(x))
labels <- 100 * bins[-1]
labels <- labels[!sapply(a1[, 2], is.na)]
ymin <- min(sapply(r2_dosage_by_af, function(d) {
  m <- as.matrix(apply(d[, -1], 2, unlist))
  min(m, na.rm = T)
}))

plot(1, col = "transparent", axes = F, xlim = c(min(x), max(x)), ylim = c(0.9 * ymin, 1.0), ylab = "Aggregated R2 within each MAF bin", xlab = "Minor Allele Frequency %")
nd <- length(groups)

for (i in 1:nd) {
  d <- r2_dosage_by_af[[i]]
  y <- rmna(d$QUILT2)
  lines(x, y, type = "l", lwd = i / nd * 2.5, pch = 1, col = mycols["QUILT2"])
  y <- rmna(d$GLIMPSE2)
  lines(x, y, type = "l", lwd = i / nd * 2.5, pch = 1, col = mycols["GLIMPSE2"])
}

axis(side = 1, at = x, labels = labels)
axis(side = 2, at = seq(0, 1, 0.2))
legend("bottomright", legend = paste0(groups, "x"), lwd = (1:nd) * 2.5 / nd, bty = "n")
legend("right", legend = c("QUILT2", "GLIMPSE2"), col = mycols, pch = 1, lwd = 1.5, cex = 1.0, xjust = 0, yjust = 1, bty = "n")

dev.off()
