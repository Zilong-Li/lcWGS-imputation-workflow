
acc_r2_by_af <- function(d0, d1, d2, d3, d4, af, breaks) {
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
    idx <- seq(1, dim(d4)[2], 3) # dosage indicies
    d4 <- as.matrix(d4[, idx[-1]])
    d4_cor_af <- tapply(1:length(x), x, function(w) { c(n = length(w),
                                                        nA = sum(truthGT[w,], na.rm = TRUE),
                                                        simple = cor(as.vector(truthGT[w,]), as.vector(d4[w,]), use = 'pairwise.complete') ** 2
                                                        )})
    as.data.frame(cbind(bin = breaks[-1], regular = sapply(d1_cor_af, "[[", "simple"),  mspbwt = sapply(d2_cor_af, "[[", "simple"), zilong = sapply(d3_cor_af, "[[", "simple"), glimpse = sapply(d4_cor_af, "[[", "simple")))
}

groups <- as.numeric(snakemake@config[["refsize"]])

df.truth <- read.table(snakemake@input[["truth"]])
af <- as.numeric(read.table(snakemake@input[["af"]])[,1])
af <- ifelse(af>0.5, 1-af, af)

dl.regular <- lapply(snakemake@input[["regular"]], read.table)
dl.mspbwt <- lapply(snakemake@input[["mspbwt"]], read.table)
dl.zilong <- lapply(snakemake@input[["zilong"]], read.table)
dl.glimpse <- lapply(snakemake@input[["glimpse"]], read.table)

bins <- sort(unique(c(
    c(0, 0.01 / 100, 0.02 / 100, 0.05 / 100),
    c(0, 0.01 / 10, 0.02 / 10, 0.05 / 10),
    c(0, 0.01 / 1, 0.02 / 1, 0.05 / 1),
    seq(0.1, 0.5, length.out = 5)
)))

accuracy_by_af <- lapply(1:length(groups), function(i) {
    acc_r2_by_af(df.truth, dl.regular[[i]], dl.mspbwt[[i]], dl.zilong[[i]], dl.glimpse[[i]], af, bins)
})

mycols <- c("black", "orange", "red", "blue" )

pdf(snakemake@output[[1]], w=12, h=6)
par(mfrow = c(1, 2))
x <- log10(as.numeric(accuracy_by_af[[1]]$bin))

plot(1, col = "transparent", axes = F, xlim = c(min(x), max(x)), ylim = c(0, 1.0), ylab = "Aggregated R2 within each MAF bin",  xlab = "Minor Allele Frequency")
nd <- length(groups)
for(i in 1:nd) {
    d <- accuracy_by_af[[i]]
    y <- unlist(ifelse(sapply(d$regular, is.null), NA, d$regular))
    lines(x, y, type = "b", lty = i, pch = 1, col = mycols[1])
    y <- unlist(ifelse(sapply(d$mspbwt, is.null), NA, d$mspbwt))
    lines(x, y, type = "b", lty = i, pch = 1, col = mycols[2])
    y <- unlist(ifelse(sapply(d$zilong, is.null), NA, d$zilong))
    lines(x, y, type = "b", lty = i, pch = 1, col = mycols[3])
    y <- unlist(ifelse(sapply(d$glimpse, is.null), NA, d$glimpse))
    lines(x, y, type = "b", lty = i, pch = 1, col = mycols[4])
}
axis(side = 1, at = x, labels=100*bins[-1])
axis(side = 2, at = seq(0, 1, 0.2))
legend("bottomright", legend=paste0("N=",groups), lty = 1:nd, bty = "n")

plot(1, col = "transparent", axes = F, xlim = c(min(x), max(x)), ylim = c(0.90, 1.0), ylab = "Aggregated R2 within each MAF bin",  xlab = "Minor Allele Frequency")
for(i in 1:nd) {
    d <- accuracy_by_af[[i]]
    y <- unlist(ifelse(sapply(d$regular, is.null), NA, d$regular))
    lines(x, y, type = "b", lty = i, pch = 1, col = mycols[1])
    y <- unlist(ifelse(sapply(d$mspbwt, is.null), NA, d$mspbwt))
    lines(x, y, type = "b", lty = i, pch = 1, col = mycols[2])
    y <- unlist(ifelse(sapply(d$zilong, is.null), NA, d$zilong))
    lines(x, y, type = "b", lty = i, pch = 1, col = mycols[3])
    y <- unlist(ifelse(sapply(d$glimpse, is.null), NA, d$glimpse))
    lines(x, y, type = "b", lty = i, pch = 1, col = mycols[4])
}
axis(side = 1, at = x, labels=100*bins[-1])
axis(side = 2)
legend("bottomright", legend=c("QUILT-regular", "QUILT-mspbwt", "QUILT-zilong", "GLIMPSE"), col=c("black", "orange", "red", "blue"), pch = 1, lwd = 1.5, cex = 1.1, xjust = 0, yjust = 1, bty = "n")
dev.off()
