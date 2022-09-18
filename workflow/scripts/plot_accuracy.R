acc_r2 <- function(d) {
    x <- rowSums(d[, c("h1","h2")])
    y1 <- cor(x, d[, c("ds.qr")], use = 'pairwise.complete') ** 2
    y2 <- cor(x, d[, c("ds.qs")], use = 'pairwise.complete') ** 2
    y3 <- cor(x, d[, c("ds.qz")], use = 'pairwise.complete') ** 2
    c(y1,y2,y3)
}

deps <- as.numeric(snakemake@config[["downsample"]])

df.truth <- read.table(snakemake@input[["truth"]][[1]], col.names = c("snpid", "h1", "h2"))
dl.regular <- lapply(snakemake@input[["regular"]], read.table, col.names = c("snpid", "h1.qr", "h2.qr", "ds.qr"))
dl.mspbwt <- lapply(snakemake@input[["mspbwt"]], read.table, col.names = c("snpid", "h1.qm", "h2.qm", "ds.qm"))
dl.zilong <- lapply(snakemake@input[["zilong"]], read.table, col.names = c("snpid", "h1.qz", "h2.qz", "ds.qz"))

dl.merge <- lapply(1:length(deps), function(i) {
    merge(merge(merge(df.truth, dl.regular[[i]]), dl.mspbwt[[i]]), dl.zilong[[i]])
})

accuracy <- matrix(sapply(dl.merge, acc_r2), ncol = length(deps))

pdf(snakemake@output[[1]], w=6, h=6)
par(mfrow=c(1,2))
plot(deps, accuracy[1,], type = "b", lwd=1.5, pch = 1, ylab = "Aggregated R2 for the chromosome", xlab = "Samples sequencing depth", ylim=c(0, 1))
lines(deps, accuracy[2,], type = "b", lwd=1.5, pch = 2, col = "orange")
lines(deps, accuracy[3,], type = "b", lwd=1.5, pch = 5, col = "red")
legend("topleft", legend=c("QUILT-regular", "QUILT-mspbwt", "QUILT-zilong"), col=c("black", "orange", "red"), pch = c(1,2,5), lwd = 1.5, cex = 1.1, xjust = 0, yjust = 1, bty = "n")

dev.off()
