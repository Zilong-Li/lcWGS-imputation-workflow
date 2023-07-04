
snakemake@source("common.R")

acc_r2_all <- function(d0, d1, d2, d3) {
  id <- intersect(intersect(intersect(rownames(d0), rownames(d1)), rownames(d2)), rownames(d3))
  y1 <- cor(as.vector(d0[id,]), as.vector(d1[id,]), use = "pairwise.complete")**2
  y2 <- cor(as.vector(d0[id,]), as.vector(d2[id,]), use = "pairwise.complete")**2
  y3 <- cor(as.vector(d0[id,]), as.vector(d3[id,]), use = "pairwise.complete")**2
  c(y1, y2, y3)
}

local_r2_by_af <- function(d0, d1, d2, d3, af, bins) {
  id <- intersect(intersect(intersect(rownames(d0), rownames(d1)), rownames(d2)), rownames(d3))
  id <- intersect(id, names(af))
  res1 <- r2_by_freq(breaks = bins, af, truthG = d0, testDS = d1, which_snps = id)
  res2 <- r2_by_freq(breaks = bins, af, truthG = d0, testDS = d2, which_snps = id)
  res3 <- r2_by_freq(breaks = bins, af, truthG = d0, testDS = d3, which_snps = id)
  as.data.frame(cbind(bin = bins[-1], regular = res1[, "simple"], mspbwt = res2[, "simple"], zilong = res3[, "simple"]))
}

groups <- as.numeric(snakemake@config[["downsample"]])

df.truth <- read.table(snakemake@input[["truth"]])
df.truth <- sapply(seq(1, dim(df.truth)[2] - 1, 2), function(i) {
  rowSums(df.truth[, (i + 1):(i + 2)])
}) # matrix: nsnps x nsamples
rownames(df.truth) <- read.table(snakemake@input[["truth"]])[,1]
af <- as.numeric(read.table(snakemake@input[["af"]])[, 2])
names(af) <- read.table(snakemake@input[["af"]])[, 1]

groups <- as.numeric(snakemake@config[["downsample"]])

dl.regular <- lapply(snakemake@input[["regular"]], parse.quilt.gts)
dl.mspbwt <- lapply(snakemake@input[["mspbwt"]], parse.quilt.gts)
dl.zilong <- lapply(snakemake@input[["zilong"]], parse.quilt.gts)

bins <- sort(unique(c(
  c(0, 0.01 / 100, 0.02 / 100, 0.05 / 100),
  c(0, 0.01 / 10, 0.02 / 10, 0.05 / 10),
  c(0, 0.01 / 1, 0.02 / 1, 0.05 / 1),
  seq(0.1, 0.5, length.out = 5)
)))

accuracy <- matrix(sapply(1:length(groups), function(i) {
  acc_r2_all(df.truth, dl.regular[[i]], dl.mspbwt[[i]], dl.zilong[[i]])
}), ncol = length(groups))


accuracy_by_af <- lapply(1:length(groups), function(i) {
  d <- local_r2_by_af(df.truth, dl.regular[[i]], dl.mspbwt[[i]], dl.zilong[[i]], af, bins)
  colnames(d) <- c("bin", "regular", "mspbwt", "zilong" )
  d
})
names(accuracy_by_af) <- paste0(as.character(groups), "x")
saveRDS(accuracy_by_af, snakemake@output[["rds"]])

## accuracy_by_af <- readRDS("/maps/projects/alab/people/rlk420/quilt2/human/HRC_CEU/quilt-rare-common/results/summary/quilt.accuracy.panelsize0.chr20.rds" )

wong <- c("#e69f00", "#d55e00", "#56b4e9", "#cc79a7", "#009e73", "#0072b2", "#f0e442")

pdf(paste0(snakemake@output[["rds"]], ".pdf"), w = 12, h = 6)

par(mfrow = c(1, 2))

plot(groups, accuracy[1, ], type = "b", lwd = 1.0, pch = 1, col = wong[1], ylab = "Aggregated R2 for the chromosome", xlab = "Samples sequencing depth", ylim = c(0.9 * min(accuracy), 1.0))
lines(groups, accuracy[2, ], type = "b", lwd = 1.0, pch = 1, col = wong[2])
lines(groups, accuracy[3, ], type = "b", lwd = 1.0, pch = 1, col = wong[3])

a1 <- accuracy_by_af[[1]]
x <- a1$bin[!sapply(a1[, 2], is.na)] # remove AF bin with NULL results
x <- log10(as.numeric(x))
labels <- 100 * bins[-1]
labels <- labels[!sapply(a1[, 2], is.na)]
ymin <- min(sapply(accuracy_by_af, function(d) {
  m <- as.matrix(apply(d[, -1], 2, unlist))
  min(m, na.rm = T)
}))

plot(1, col = "transparent", axes = F, xlim = c(min(x), max(x)), ylim = c(0, 1.0), ylab = "Aggregated R2 within each MAF bin", xlab = "Minor Allele Frequency")

nd <- length(groups)
for (i in 1:nd) {
  d <- accuracy_by_af[[i]]
  y <- rmna(d$regular)
  lines(x, y, type = "l", lty = nd - i + 1, pch = 1, col = wong[1])
  y <- rmna(d$mspbwt)
  lines(x, y, type = "l", lty = nd - i + 1, pch = 1, col = wong[2])
  y <- rmna(d$zilong)
  lines(x, y, type = "l", lty = nd - i + 1, pch = 1, col = wong[3])
}
axis(side = 1, at = x, labels = labels)
axis(side = 2)

legend("topleft", legend = c("QUILT-regular", "QUILT-mspbwt", "QUILT-zilong"), col = mycols, pch = 1, lwd = 1.5, cex = 1.1, xjust = 0, yjust = 1, bty = "n")
legend("bottomright", legend = paste0(groups, "x"), lwd = (1:nd) * 2.5 / nd, bty = "n")

chunkfile <- snakemake@params[["chunks"]]
chunk.names <- read.table(chunkfile)[,4]
chunk <- lapply(strsplit(gsub(".*:","",chunk.names),"-"), as.integer)
pos <- as.integer(sapply(strsplit(names(af),":"),"[[",2))
chunk_af <- lapply(chunk, function(c) {
  af[which(pos > c[1] & pos < c[2])]
})
names(chunk_af) <- chunk.names

accuracy_by_af_chunk <- lapply(chunk_af, function(af) {
  all <- lapply(seq(length(groups)), function(i) {
    d <- local_r2_by_af(df.truth, dl.regular[[i]], dl.mspbwt[[i]], dl.zilong[[i]], af, bins)
    colnames(d) <- c("bin", "regular", "mspbwt", "zilong" )
    d
  })
  names(all) <- paste0(as.character(groups), "x")
  all
})

for(c in 1:length(chunk.names)) {
  if(c %% 2 == 1) par(mfrow = c(1, 2))
  title <- paste(names(chunk_af)[c], "#", length(chunk_af[[c]]))
  acc_chunk <- accuracy_by_af_chunk[[c]]
  plot(1, col = "transparent", axes = F, xlim = c(min(x), max(x)), ylim = c(0, 1.0), ylab = "Aggregated R2 within each AF bin", xlab = "Allele Frequency",main = title)
  for (i in 1:nd) {
    d <- acc_chunk[[i]]
    y <- rmna(d$regular)
    lines(x, y, type = "l", lty = nd - i + 1, pch = 1, col = wong[1])
    y <- rmna(d$mspbwt)
    lines(x, y, type = "l", lty = nd - i + 1, pch = 1, col = wong[2])
    y <- rmna(d$zilong)
    lines(x, y, type = "l", lty = nd - i + 1, pch = 1, col = wong[3])
  }
  axis(side = 1, at = x, labels = labels)
  axis(side = 2, at = seq(0, 1, 0.2))
  legend("bottomright", legend = paste0(groups, "x"), lwd = (1:nd) * 2.5 / nd, bty = "n")
}

dev.off()
