
snakemake@source("utils.R")

acc_r2_by_af <- function(d0, d1, d2, d3, d4, af, bins) {
    res1 <- r2_by_freq(breaks = bins, af, truthG = d0, testDS = d1)
    res2 <- r2_by_freq(breaks = bins, af, truthG = d0, testDS = d2)
    res3 <- r2_by_freq(breaks = bins, af, truthG = d0, testDS = d3)
    res4 <- r2_by_freq(breaks = bins, af, truthG = d0, testDS = d4)
    as.data.frame(cbind(bin = bins[-1], regular = res1[,"simple"], mspbwt = res2[,"simple"], zilong = res3[,"simple"], glimpse = res4[,"simple"]))
}


groups <- as.numeric(snakemake@config[["downsample"]])

df.truth <- read.table(snakemake@input[["truth"]])
df.truth <- sapply(seq(1, dim(df.truth)[2] - 1, 2), function(i){rowSums(df.truth[,(i+1):(i+2)])})  # matrix: nsnps x nsamples
af <- as.numeric(read.table(snakemake@input[["af"]])[,1])

dl.regular <- lapply(snakemake@input[["regular"]], function(fn) {
    d1 <- fread(cmd = paste("awk '{for(i=1;i<=NF;i=i+3) printf $i\" \"; print \"\"}'", fn), data.table = F)
    d1 <- as.matrix(sapply(d1[,-1], as.numeric))
})
dl.mspbwt <- lapply(snakemake@input[["mspbwt"]], function(fn) {
    d1 <- fread(cmd = paste("awk '{for(i=1;i<=NF;i=i+3) printf $i\" \"; print \"\"}'", fn), data.table = F)
    d1 <- as.matrix(sapply(d1[,-1], as.numeric))
})
dl.zilong <- lapply(snakemake@input[["zilong"]], function(fn) {
    d1 <- fread(cmd = paste("awk '{for(i=1;i<=NF;i=i+3) printf $i\" \"; print \"\"}'", fn), data.table = F)
    d1 <- as.matrix(sapply(d1[,-1], as.numeric))
})
dl.glimpse <- lapply(snakemake@input[["glimpse"]], function(fn) {
    d1 <- fread(cmd = paste("awk '{for(i=1;i<=NF;i=i+3) printf $i\" \"; print \"\"}'", fn), data.table = F)
    d1 <- as.matrix(sapply(d1[,-1], as.numeric))
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

pdf(snakemake@output[["pdf"]], w=12, h=6)

a1 <- accuracy_by_af[[1]]
x <- a1$bin[!sapply(a1[,2], is.na)] # remove AF bin with NULL results
x <- log10(as.numeric(x))
labels <- 100 * bins[-1]
labels <- labels[!sapply(a1[,2], is.na)]
ymin <- min(sapply(accuracy_by_af, function(d){
    m = as.matrix(apply(d[,-1], 2, unlist))
    min(m, na.rm = T)
}))


par(mfrow = c(1, 2))
plot(1, col = "transparent", axes = F, xlim = c(min(x), max(x)), ylim = c(0.9*ymin, 1.0), ylab = "Aggregated R2 within each MAF bin",  xlab = "Minor Allele Frequency")
nd <- length(groups)

for(i in 1:nd) {
    d <- accuracy_by_af[[i]]
    # https://stackoverflow.com/questions/33004238/r-removing-null-elements-from-a-list
    y <- rmnull(d$regular)
    lines(x, y, type = "l", lwd = i/nd * 2.5, pch = 1, col = mycols[1])
    y <- rmnull(d$mspbwt)
    lines(x, y, type = "l", lwd = i/nd * 2.5, pch = 1, col = mycols[2])
    y <- rmnull(d$zilong)
    lines(x, y, type = "l", lwd = i/nd * 2.5, pch = 1, col = mycols[3])
    y <- rmnull(d$glimpse)
    lines(x, y, type = "l", lwd = i/nd * 2.5, pch = 1, col = mycols[4])
}
axis(side = 1, at = x, labels=labels)
axis(side = 2, at = seq(0, 1, 0.2))
legend("bottomright", legend=paste0(groups, "x"), lwd = (1:nd) * 2.5 / nd, bty = "n")

plot(1, col = "transparent", axes = F, xlim = c(min(x), max(x)), ylim = c(0.90, 1.0), ylab = "Aggregated R2 within each MAF bin",  xlab = "Minor Allele Frequency")
for(i in 1:nd) {
    d <- accuracy_by_af[[i]]
    y <- rmnull(d$regular)
    lines(x, y, type = "l", lwd = i/nd * 2.5, pch = 1, col = mycols[1])
    y <- rmnull(d$mspbwt)
    lines(x, y, type = "l", lwd = i/nd * 2.5, pch = 1, col = mycols[2])
    y <- rmnull(d$zilong)
    lines(x, y, type = "l", lwd = i/nd * 2.5, pch = 1, col = mycols[3])
    y <- rmnull(d$glimpse)
    lines(x, y, type = "l", lwd = i/nd * 2.5, pch = 1, col = mycols[4])
}
axis(side = 1, at = x, labels=labels)
axis(side = 2)
legend("bottomleft", legend=c("QUILT-regular", "QUILT-mspbwt", "QUILT-zilong", "GLIMPSE"), col=mycols, pch = 1, lwd = 1.5, cex = 1.0, xjust = 0, yjust = 1, bty = "n")
dev.off()
