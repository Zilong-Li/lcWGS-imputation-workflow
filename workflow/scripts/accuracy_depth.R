
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
saveRDS(accuracy_by_af, snakemake@output[["rds"]])

wong <- c("#e69f00", "#d55e00", "#56b4e9", "#cc79a7", "#009e73", "#0072b2", "#f0e442")
mycols <- wong

pdf(snakemake@output[["pdf"]], w = 12, h = 6)

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
plot(1, col = "transparent", axes = F, xlim = c(min(x), max(x)), ylim = c(0.9 * ymin, 1.0), ylab = "Aggregated R2 within each MAF bin", xlab = "Minor Allele Frequency")
nd <- length(groups)

for (i in 1:nd) {
  d <- accuracy_by_af[[i]]
  # https://stackoverflow.com/questions/33004238/r-removing-null-elements-from-a-list
  y <- rmna(d$QUILT2)
  lines(x, y, type = "l", lwd = i / nd * 2.5, pch = 1, col = mycols[1])
  y <- rmna(d$GLIMPSE2)
  lines(x, y, type = "l", lwd = i / nd * 2.5, pch = 1, col = mycols[2])
  y <- rmna(d$QUILT1)
  lines(x, y, type = "l", lwd = i / nd * 2.5, pch = 1, col = mycols[3])
  y <- rmna(d$GLIMPSE1)
  lines(x, y, type = "l", lwd = i / nd * 2.5, pch = 1, col = mycols[4])
}
axis(side = 1, at = x, labels = labels)
axis(side = 2, at = seq(0, 1, 0.2))
legend("bottomright", legend = paste0(groups, "x"), lwd = (1:nd) * 2.5 / nd, bty = "n")

plot(1, col = "transparent", axes = F, xlim = c(min(x), max(x)), ylim = c(0.90, 1.0), ylab = "Aggregated R2 within each MAF bin", xlab = "Minor Allele Frequency")
for (i in 1:nd) {
  d <- accuracy_by_af[[i]]
  y <- rmna(d$QUILT2)
  lines(x, y, type = "l", lwd = i / nd * 2.5, pch = 1, col = mycols[1])
  y <- rmna(d$GLIMPSE2)
  lines(x, y, type = "l", lwd = i / nd * 2.5, pch = 1, col = mycols[2])
  y <- rmna(d$QUILT1)
  lines(x, y, type = "l", lwd = i / nd * 2.5, pch = 1, col = mycols[3])
  y <- rmna(d$GLIMPSE1)
  lines(x, y, type = "l", lwd = i / nd * 2.5, pch = 1, col = mycols[4])
}
axis(side = 1, at = x, labels = labels)
axis(side = 2)
legend("bottomleft", legend = c("QUILT2", "GLIMPSE2", "QUILT1", "GLIMPSE1"), col = mycols, pch = 1, lwd = 1.5, cex = 1.0, xjust = 0, yjust = 1, bty = "n")
dev.off()
