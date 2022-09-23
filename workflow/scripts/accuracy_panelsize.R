library(data.table)

r2_by_freq <- function(breaks, af, truthG, testDS, which_snps = NULL, flip = FALSE) {
    if (flip) {
        w <- af > 0.5
        af[w] <- 1 - af[w]
        truthG[w] <- 2 - truthG[w]
        testDS[w] <- 2 - testDS[w]
    }
    if (!is.null(which_snps)) {
        af <- af[which_snps]
        truthG <- truthG[which_snps]
        testDS <- testDS[which_snps]
    }
    x <- cut(af, breaks = breaks)
    cors_per_af <- tapply(1:length(x), x, function(w) {
        c(
            n = length(w),
            nA = sum(truthG[w], na.rm = TRUE),
            simple = cor(truthG[w], testDS[w], use = 'pairwise.complete') ** 2,
            norm = cor(truthG[w] - 2 * af[w], testDS[w] - 2 * af[w], use = 'pairwise.complete') ** 2
        )
    })
    cors_per_af <- t(sapply(cors_per_af[!sapply(cors_per_af, is.null)], I))
    return(cors_per_af)
}

acc_r2_by_af <- function(d0, d1, d2, d3, d4, af, bins) {
    truthGT <- sapply(seq(1, dim(d0)[2] - 1, 2), function(i){rowSums(d0[,(i+1):(i+2)])})  # matrix: nsnps x nsamples
    d1 <- as.matrix(d1[, -1])
    res1 <- r2_by_freq(breaks = bins, af, truthG = truthGT, testDS = d1)
    d2 <- as.matrix(d2[, -1])
    res2 <- r2_by_freq(breaks = bins, af, truthG = truthGT, testDS = d2)
    d3 <- as.matrix(d3[, -1])
    res3 <- r2_by_freq(breaks = bins, af, truthG = truthGT, testDS = d3)
    d4 <- as.matrix(d4[, -1])
    res4 <- r2_by_freq(breaks = bins, af, truthG = truthGT, testDS = d4)
    as.data.frame(cbind(bin = bins[-1], regular = res1[,"simple"], mspbwt = res2[,"simple"], zilong = res3[,"simple"], glimpse = res4[,"simple"]))
}

# https://stackoverflow.com/questions/33004238/r-removing-null-elements-from-a-list
rmnull <- function(l) {
   l <- l[!sapply(l, is.null)]
   unlist(l)
}

groups <- as.numeric(snakemake@config[["refsize"]])
nd <- length(groups)

df.truth <- read.table(snakemake@input[["truth"]])
af <- as.numeric(read.table(snakemake@input[["af"]])[,1])

dl.regular <- lapply(snakemake@input[["regular"]], function(fn) {
    fread(cmd = paste("awk '{for(i=1;i<=NF;i=i+3) printf $i\" \"; print \"\"}'", fn), data.table = F)
})
dl.mspbwt <- lapply(snakemake@input[["mspbwt"]], function(fn) {
    fread(cmd = paste("awk '{for(i=1;i<=NF;i=i+3) printf $i\" \"; print \"\"}'", fn), data.table = F)
})
dl.zilong <- lapply(snakemake@input[["zilong"]], function(fn) {
    fread(cmd = paste("awk '{for(i=1;i<=NF;i=i+3) printf $i\" \"; print \"\"}'", fn), data.table = F)
})
dl.glimpse <- lapply(snakemake@input[["glimpse"]], function(fn) {
    fread(cmd = paste("awk '{for(i=1;i<=NF;i=i+3) printf $i\" \"; print \"\"}'", fn), data.table = F)
})

bins <- sort(unique(c(
    c(0, 0.01 / 100, 0.02 / 100, 0.05 / 100),
    c(0, 0.01 / 10, 0.02 / 10, 0.05 / 10),
    c(0, 0.01 / 1, 0.02 / 1, 0.05 / 1),
    seq(0.1, 0.5, length.out = 5)
)))

accuracy_by_af <- lapply(1:length(groups), function(i) {
    acc_r2_by_af(df.truth, dl.regular[[i]], dl.mspbwt[[i]], dl.zilong[[i]], dl.glimpse[[i]], af, bins)
})

saveRDS(accuracy_by_af, snakemake@output[["rds"]])

wong <- c("#e69f00", "#d55e00", "#56b4e9", "#cc79a7", "#009e73", "#0072b2", "#f0e442")
mycols <- wong[1:4]

pdf(snakemake@output[["pdf"]], w=8, h=12)
a1 <- accuracy_by_af[[1]]
x <- a1$bin[!sapply(a1[,2], is.null)]
x <- log10(as.numeric(x))
labels <- 100 * bins[-1]
labels <- labels[!sapply(a1[,2], is.null)]
ymin <- min(sapply(accuracy_by_af, function(d){
    m = as.matrix(apply(d[,-1], 2, unlist))
    min(m, na.rm = T)
}))

par(mfrow = c(2,1))
plot(1, col = "transparent", axes = F, xlim = c(min(x), max(x)), ylim = c(0.90, 1.0), ylab = "Aggregated R2 within each MAF bin",  xlab = "Minor Allele Frequency")
for(i in 1:nd) {
    d <- accuracy_by_af[[i]]
    y <- rmnull(d$regular)
    lines(x, y, type = "l", lty = nd-i+1, pch = 1, col = mycols[1])
    y <- rmnull(d$mspbwt)
    lines(x, y, type = "l", lty = nd-i+1, pch = 1, col = mycols[2])
    y <- rmnull(d$zilong)
    lines(x, y, type = "l", lty = nd-i+1, pch = 1, col = mycols[3])
    y <- rmnull(d$glimpse)
    lines(x, y, type = "l", lty = nd-i+1, pch = 1, col = mycols[4])
}
axis(side = 1, at = x, labels=labels)
axis(side = 2)
legend("bottomleft", legend=c("QUILT-regular", "QUILT-mspbwt", "QUILT-zilong", "GLIMPSE"), col=mycols, pch = 1, lwd = 1.5, cex = 1.0, xjust = 0, yjust = 1, bty = "n")

plot(1, col = "transparent", axes = F, xlim = c(min(x), max(x)), ylim = c(0.9*ymin, 1.0), ylab = "Aggregated R2 within each MAF bin",  xlab = "Minor Allele Frequency")
nd <- length(groups)
for(i in 1:nd) {
    d <- accuracy_by_af[[i]]
    y <- rmnull(d$regular)
    lines(x, y, type = "l", lty = nd-i+1, pch = 1, col = mycols[1])
    y <- rmnull(d$mspbwt)
    lines(x, y, type = "l", lty = nd-i+1, pch = 1, col = mycols[2])
    y <- rmnull(d$zilong)
    lines(x, y, type = "l", lty = nd-i+1, pch = 1, col = mycols[3])
    y <- rmnull(d$glimpse)
    lines(x, y, type = "l", lty = nd-i+1, pch = 1, col = mycols[4])
}
axis(side = 1, at = x, labels=labels)
axis(side = 2, at = seq(0, 1, 0.2))
legend("bottomright", legend=paste0("N=",groups), lty = nd:1, bty = "n")

dev.off()
