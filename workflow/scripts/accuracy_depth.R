
snakemake@source("common.R")


groups <- as.numeric(snakemake@config[["downsample"]])

df.truth <- read.table(snakemake@input[["truth"]])
df.truth <- sapply(seq(1, dim(df.truth)[2] - 1, 2), function(i) {
  rowSums(df.truth[, (i + 1):(i + 2)])
}) # matrix: nsnps x nsamples
rownames(df.truth) <- read.table(snakemake@input[["truth"]])[,1]
af <- as.numeric(read.table(snakemake@input[["af"]])[, 2])
names(af) <- read.table(snakemake@input[["af"]])[, 1]

dl.quilt1 <- lapply(snakemake@input[["regular"]], parse.quilt.gts)
dl.quilt2 <- lapply(snakemake@input[["zilong"]], parse.quilt.gts)
dl.glimpse1 <- lapply(snakemake@input[["glimpse1"]], parse.quilt.gts)
dl.glimpse2 <- lapply(snakemake@input[["glimpse2"]], parse.quilt.gts)

bins <- sort(unique(c(
  c(0, 0.01 / 100, 0.02 / 100, 0.05 / 100),
  c(0, 0.01 / 10, 0.02 / 10, 0.05 / 10),
  c(0, 0.01 / 1, 0.02 / 1, 0.05 / 1),
  seq(0.1, 0.5, length.out = 5)
)))

accuracy_by_af <- lapply(seq(length(groups)), function(i) {
  d <- acc_r2_by_af(df.truth, dl.quilt2[[i]], dl.glimpse2[[i]], dl.quilt1[[i]], dl.glimpse1[[i]], af, bins)
  colnames(d) <- c("bin","QUILT2", "GLIMPSE2", "QUILT1", "GLIMPSE1")
  d
})
names(accuracy_by_af) <- paste0(as.character(groups), "x")

saveRDS(accuracy_by_af, snakemake@output[["rds"]])

## (rds <- readRDS("/maps/projects/alab/people/rlk420/quilt2/human/UKBB_GEL_CEU/bench-speed/results/summary/all.accuracy.panelsize0.chr20.rds"))

pdf(paste0(snakemake@output[["rds"]], ".pdf"), w = 12, h = 6)

a1 <- accuracy_by_af[[1]]
x <- a1$bin[!sapply(a1[, 2], is.na)] # remove AF bin with NULL results
x <- log10(as.numeric(x))
labels <- 100 * bins[-1]
labels <- labels[!sapply(a1[, 2], is.na)]
ymin <- min(sapply(accuracy_by_af, function(d) {
  m <- as.matrix(apply(d[, -1], 2, unlist))
  min(m, na.rm = T)
}))


par(mfrow = c(1, 2))
plot(1, col = "transparent", axes = F, xlim = c(min(x), max(x)), ylim = c(0.9 * ymin, 1.0), ylab = "Aggregated R2 within each AF bin", xlab = "Allele Frequency")
nd <- length(groups)

for (i in 1:nd) {
  d <- accuracy_by_af[[i]]
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

plot(1, col = "transparent", axes = F, xlim = c(min(x), max(x)), ylim = c(0.90, 1.0), ylab = "Aggregated R2 within each AF bin", xlab = "Allele Frequency")
for (i in 1:nd) {
  d <- accuracy_by_af[[i]]
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


## chunkfile <- "/maps/projects/alab/people/rlk420/quilt2/human/HRC_CEU/quilt-rare-common/results/refpanels/chr20.glimpse.chunks"
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
    d <- acc_r2_by_af(df.truth, dl.quilt2[[i]], dl.glimpse2[[i]], dl.quilt1[[i]], dl.glimpse1[[i]], af, bins)
    colnames(d) <- c("bin","QUILT2", "GLIMPSE2", "QUILT1", "GLIMPSE1")
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
