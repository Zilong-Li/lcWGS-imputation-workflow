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

dl.quilt1 <- mclapply(snakemake@input[["regular"]], getGT, mc.cores = 4, region=chrom, samples=samples)
dl.quilt2 <- mclapply(snakemake@input[["zilong"]], getGT, mc.cores = 4,region=chrom, samples=samples)
dl.glimpse1 <- mclapply(snakemake@input[["glimpse1"]], getGT,mc.cores = 4, region=chrom, samples=samples)
dl.glimpse2 <- mclapply(snakemake@input[["glimpse2"]], getGT, mc.cores = 4,region=chrom, samples=samples)

f1.quilt1 <- lapply(dl.quilt1, function(test) acc_f1(truth, test))
names(f1.quilt1) <- depths
f1.quilt2 <- lapply(dl.quilt2, function(test) acc_f1(truth, test))
names(f1.quilt2) <- depths
f1.glimpse1 <- lapply(dl.glimpse1, function(test) acc_f1(truth, test))
names(f1.glimpse1) <- depths
f1.glimpse2 <- lapply(dl.glimpse2, function(test) acc_f1(truth, test))
names(f1.glimpse2) <- depths

res <- list(f1.quilt1, f1.quilt2, f1.glimpse1, f1.glimpse2)
names(res) <- c("QUILT1", "QUILT2", "GLIMPSE1", "GLIMPSE2")

saveRDS(res, snakemake@output[["rds"]])

## res <- readRDS("/maps/projects/alab/people/rlk420/quilt2/human/UKBB_GEL_ONT/default5/results/summary/all.accuracy.f1.panelsize0.chr20.rds")

## depths <- c("0.1x", "0.5x", "1.0x", "2.0x")

nsamples <- ncol(truth$gt)

out <- do.call(rbind, lapply(names(res), function(name) {
  v <- unlist(res[[name]])
  d <- data.frame(v = v, m=rep(name, length(v)), d=rep(depths, each=nsamples))
  d
}))

mycols <- c(QUILT2 = "#e69f00", GLIMPSE2 = "#d55e00", QUILT1 = "#56b4e9", GLIMPSE1 = "#cc79a7")

out <- out[order(out$d, out$m),]

pdf(paste0(snakemake@output[["rds"]], ".pdf"), h=7,w=7)
a <- c(seq(4), seq(6,9), seq(11,14), seq(16,19))
boxplot(v ~ m + d, data=out, col=mycols[unique(out$m)],
        at = a, cex.lab=2.0, cex.axis = 1.5,
        xaxt="n", xlab = "Coverage", ylab = "F1-score")
b <- c(3, 8, 13, 18)
axis(side = 1, at = b, labels = unique(out$d), cex.axis = 2.0)
legend("bottomright", legend = names(mycols), fill = mycols, bty = "n", cex=2)
dev.off()

q()

## Sys.setenv(DISPLAY="localhost:13.0")



