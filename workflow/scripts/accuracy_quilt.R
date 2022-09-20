acc_r2_all <- function(d0, d1, d2, d3) {
    truthGT <- as.vector(sapply(seq(1, dim(d0)[2] - 1, 2), function(i){rowSums(d0[,(i+1):(i+2)])}))
    idx <- seq(1, dim(d1)[2], 3) # dosage indicies
    y1 <- cor(truthGT, as.vector(as.matrix(d1[, idx[-1]])), use = 'pairwise.complete') ** 2
    idx <- seq(1, dim(d2)[2], 3) # dosage indicies
    y2 <- cor(truthGT, as.vector(as.matrix(d2[, idx[-1]])), use = 'pairwise.complete') ** 2
    idx <- seq(1, dim(d3)[2], 3) # dosage indicies
    y3 <- cor(truthGT, as.vector(as.matrix(d3[, idx[-1]])), use = 'pairwise.complete') ** 2
    c(y1,y2,y3)
}

acc_r2_by_af <- function(d0, d1, d2, d3, af, breaks) {
    truthGT <- sapply(seq(1, dim(d0)[2] - 1, 2), function(i){rowSums(d0[,(i+1):(i+2)])})  # matrix: nsnps x nsamples
    idx <- seq(1, dim(d1)[2], 3) # dosage indicies
    d1 <- as.matrix(d1[, idx[-1]])
    x <- cut(af, breaks = breaks)
    d1_cor_af <- tapply(1:length(x), x, function(w) { c(n = length(w),
                                                        nA = sum(truthGT[w,], na.rm = TRUE),
                                                        simple = cor(as.vector(truthGT[w,]), as.vector(d1[w,]), use = 'pairwise.complete') ** 2
                                                        )})
    idx <- seq(1, dim(d2)[2], 3) # dosage indicies
    d2 <- as.matrix(d2[, idx[-1]])
    d2_cor_af <- tapply(1:length(x), x, function(w) { c(n = length(w),
                                                        nA = sum(truthGT[w,], na.rm = TRUE),
                                                        simple = cor(as.vector(truthGT[w,]), as.vector(d2[w,]), use = 'pairwise.complete') ** 2
                                                        )})
    idx <- seq(1, dim(d3)[2], 3) # dosage indicies
    d3 <- as.matrix(d3[, idx[-1]])
    d3_cor_af <- tapply(1:length(x), x, function(w) { c(n = length(w),
                                                        nA = sum(truthGT[w,], na.rm = TRUE),
                                                        simple = cor(as.vector(truthGT[w,]), as.vector(d3[w,]), use = 'pairwise.complete') ** 2
                                                        )})
    as.data.frame(cbind(bin = breaks[-1], regular = sapply(d1_cor_af, "[[", "simple"),  mspbwt = sapply(d2_cor_af, "[[", "simple"), zilong = sapply(d3_cor_af, "[[", "simple")))
}

groups <- as.numeric(snakemake@config[["downsample"]])

df.truth <- read.table(snakemake@input[["truth"]])
af <- as.numeric(read.table(snakemake@input[["af"]])[,1])
af <- ifelse(af>0.5, 1-af, af)

dl.regular <- lapply(snakemake@input[["regular"]], read.table)
dl.mspbwt <- lapply(snakemake@input[["mspbwt"]], read.table)
dl.zilong <- lapply(snakemake@input[["zilong"]], read.table)

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
    acc_r2_by_af(df.truth, dl.regular[[i]], dl.mspbwt[[i]], dl.zilong[[i]], af, bins)
})

mycols <- c("black", "orange", "red", "blue" )

pdf(snakemake@output[[1]], w=12, h=6)
par(mfrow = c(1, 2))
plot(groups, accuracy[1,], type = "b", lwd=1.5, pch = 1, col = mycols[1], ylab = "Aggregated R2 for the chromosome", xlab = "Samples sequencing depth", ylim=c(0.7, 1.0))
lines(groups, accuracy[2,], type = "b", lwd=1.5, pch = 1, col = mycols[2])
lines(groups, accuracy[3,], type = "b", lwd=1.5, pch = 1, col = mycols[3])
legend("bottomright", legend=c("QUILT-regular", "QUILT-mspbwt", "QUILT-zilong"), col=c("black", "orange", "red"), pch = c(1,1,1), lwd = 1.5, cex = 1.1, xjust = 0, yjust = 1, bty = "n")

plot(1, col = "transparent", xlim = c(0, 0.5), ylim = c(0, 1.0), ylab = "Aggregated R2 within each MAF bin",  xlab = "Minor Allele Frequency")
nd <- length(groups)
for(i in 1:nd) {
    d <- accuracy_by_af[[i]]
    x <- log10(as.numeric(d$bin))
    y <- unlist(ifelse(sapply(d$regular, is.null), NA, d$regular))
    lines(x, y, type = "b", lwd = i/nd * 2.5, pch = 1, col = mycols[1])
    y <- unlist(ifelse(sapply(d$mspbwt, is.null), NA, d$mspbwt))
    lines(x, y, type = "b", lwd = i/nd * 2.5, pch = 1, col = mycols[2])
    y <- unlist(ifelse(sapply(d$zilong, is.null), NA, d$zilong))
    lines(x, y, type = "b", lwd = i/nd * 2.5, pch = 1, col = mycols[3])
}
legend("bottomright", legend=paste0(groups, "x"), lwd = (1:nd) * 2.5 / nd, bty = "n")

dev.off()
