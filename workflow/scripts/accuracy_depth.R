
## saveRDS(snakemake, snakemake@output[["rds"]])
## q()

## snakemake <- readRDS("/maps/projects/alab/people/rlk420/quilt2/human/HRC_CEU/quilt-rare-common/results/summary/all.accuracy.panelsize0.chr20.rds")
## setwd("/maps/projects/alab/people/rlk420/quilt2/human/HRC_CEU/quilt-rare-common/")

snakemake@source("common.R")

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


names(r2_dosage_by_af) <- paste0(as.character(groups), "x")

saveRDS(r2_dosage_by_af, snakemake@output[["rds"]])

pdf(paste0(snakemake@output[["rds"]], ".pdf"), w = 12, h = 6)

a1 <- r2_dosage_by_af[[1]]
x <- a1$bin[!sapply(a1[, 2], is.na)] # remove AF bin with NULL results
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
  # https://stackoverflow.com/questions/33004238/r-removing-null-elements-from-a-list
  y <- rmna(d$QUILT2)
  lines(x, y, type = "l", lwd = i / nd * 2.5, pch = 1, col = mycols["QUILT2"])
  y <- rmna(d$GLIMPSE2)
  lines(x, y, type = "l", lwd = i / nd * 2.5, pch = 1, col = mycols["GLIMPSE2"])
  y <- rmna(d$QUILT1)
  lines(x, y, type = "l", lwd = i / nd * 2.5, pch = 1, col = mycols["QUILT1"])
  y <- rmna(d$GLIMPSE1)
  lines(x, y, type = "l", lwd = i / nd * 2.5, pch = 1, col = mycols["GLIMPSE1"])
}

axis(side = 1, at = x, labels = labels)
axis(side = 2, at = seq(0, 1, 0.2))
legend("bottomright", legend = paste0(groups, "x"), lwd = (1:nd) * 2.5 / nd, bty = "n")

plot(1, col = "transparent", axes = F, xlim = c(min(x), max(x)), ylim = c(0.90, 1.0), ylab = "Aggregated R2 within each MAF bin", xlab = "Minor Allele Frequency")
for (i in 1:nd) {
  d <- r2_dosage_by_af[[i]]
  y <- rmna(d$QUILT2)
  lines(x, y, type = "l", lwd = i / nd * 2.5, pch = 1, col = mycols["QUILT2"])
  y <- rmna(d$GLIMPSE2)
  lines(x, y, type = "l", lwd = i / nd * 2.5, pch = 1, col = mycols["GLIMPSE2"])
  y <- rmna(d$QUILT1)
  lines(x, y, type = "l", lwd = i / nd * 2.5, pch = 1, col = mycols["QUILT1"])
  y <- rmna(d$GLIMPSE1)
  lines(x, y, type = "l", lwd = i / nd * 2.5, pch = 1, col = mycols["GLIMPSE1"])
}
axis(side = 1, at = x, labels = labels)
axis(side = 2)
legend("bottomleft", legend = c("QUILT2", "GLIMPSE2", "QUILT1", "GLIMPSE1"), col = mycols, pch = 1, lwd = 1.5, cex = 1.0, xjust = 0, yjust = 1, bty = "n")

dev.off()
q()

## chunkfile <- "/maps/projects/alab/people/rlk420/quilt2/human/HRC_CEU/quilt-rare-common/results/refpanels/chr20.glimpse.chunks"
chunkfile <- snakemake@params[["chunks"]]
chunk.names <- read.table(chunkfile)[,4]
chunk <- lapply(strsplit(gsub(".*:","",chunk.names),"-"), as.integer)
pos <- as.integer(sapply(strsplit(names(af),":"),"[[",2))
chunk_af <- lapply(chunk, function(c) {
  af[which(pos > c[1] & pos < c[2])]
})
names(chunk_af) <- chunk.names

r2_dosage_by_af_chunk <- lapply(chunk_af, function(af) {
  all <- lapply(seq(length(groups)), function(i) {
    n <- ncol(dl.quilt1[[i]])
    quilt2 <- dl.quilt2[[i]][,seq(3, n, by = 3 )] # get dosage
    quilt1 <- dl.quilt1[[i]][,seq(3, n, by = 3 )] # get dosage
    glimpse2 <- dl.glimpse2[[i]][,seq(3, n, by = 3 )] # get dosage
    glimpse1 <- dl.glimpse1[[i]][,seq(3, n, by = 3 )] # get dosage
    d <- acc_r2_by_af(ds.truth, quilt2 , glimpse2, quilt1, glimpse1, af, bins)
    colnames(d) <- c("bin","nsnps","QUILT2", "GLIMPSE2", "QUILT1", "GLIMPSE1")
    d
  })
  names(all) <- paste0(as.character(groups), "x")
  all
})

for(c in 1:length(chunk.names)) {
  if(c %% 2 == 1) par(mfrow = c(1, 2))
  title <- paste(names(chunk_af)[c], "#", length(chunk_af[[c]]))
  acc_chunk <- r2_dosage_by_af_chunk[[c]]
  a1 <- acc_chunk[[1]]
  x <- a1$bin[!sapply(a1[, 2], is.na)] # remove AF bin with NULL results
  x <- log10(as.numeric(x))
  labels <- 100 * bins[-1]
  labels <- labels[!sapply(a1[, 2], is.na)]
  plot(1, col = "transparent", axes = F, xlim = c(min(x), max(x)), ylim = c(0, 1.0), ylab = "Aggregated R2 within each MAF bin", xlab = "Minor Allele Frequency",main = title)
  for (i in 1:nd) {
    d <- acc_chunk[[i]]
    y <- rmna(d$QUILT2)
    lines(x, y, type = "l", lwd = i / nd * 2.5, pch = 1, col = mycols["QUILT2"])
    y <- rmna(d$GLIMPSE2)
    lines(x, y, type = "l", lwd = i / nd * 2.5, pch = 1, col = mycols["GLIMPSE2"])
    y <- rmna(d$QUILT1)
    lines(x, y, type = "l", lwd = i / nd * 2.5, pch = 1, col = mycols["QUILT1"])
    y <- rmna(d$GLIMPSE1)
    lines(x, y, type = "l", lwd = i / nd * 2.5, pch = 1, col = mycols["GLIMPSE1"])
  }
  axis(side = 1, at = x, labels = labels)
  axis(side = 2, at = seq(0, 1, 0.2))
  legend("bottomright", legend = paste0(groups, "x"), lwd = (1:nd) * 2.5 / nd, bty = "n")
}

dev.off()
