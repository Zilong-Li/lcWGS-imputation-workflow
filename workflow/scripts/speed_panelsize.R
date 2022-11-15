gettimes <- function(ss) {
  sapply(strsplit(ss, ":"), function(s) {
    s <- as.numeric(s)
    n <- length(s)
    sum(sapply(1:n, function(i) {
      s[i] * 60^(n - i)
    }))
  })
}

gnutime <- function(dl) {
  sapply(dl, function(d) {
    sum(gettimes(d[, 1]))
  })
}

gunram <- function(dl) {
  sapply(dl, function(d) {
    max(d[, 2]) / 1024 # MB units
  })
}

groups <- as.numeric(snakemake@config[["refsize"]])
nd <- length(groups)

dl.regular <- lapply(snakemake@input[["regular"]], read.table)
dl.mspbwt <- lapply(snakemake@input[["mspbwt"]], read.table)
dl.glimpse <- lapply(snakemake@input[["glimpse"]], read.table)
saveRDS(list(regular = dl.regular, mspbwt = dl.mspbwt, glimpse = dl.glimpse), snakemake@output[["rds"]])

## times <- cbind(gnutime(dl.regular), gnutime(dl.mspbwt), gnutime(dl.zilong), gnutime(dl.glimpse))
times <- cbind(gnutime(dl.glimpse), gnutime(dl.regular), gnutime(dl.mspbwt))
rams <- cbind(gunram(dl.glimpse), gunram(dl.regular), gunram(dl.mspbwt))


wong <- c("#e69f00", "#d55e00", "#56b4e9", "#cc79a7", "#009e73", "#0072b2", "#f0e442")
mycols <- wong[1:3]


pdf(snakemake@output[["pdf"]], w = 12, h = 6)
par(mfrow = c(1, 2))
plot(groups, times[, 1], type = "b", lwd = 1.0, pch = 1, col = mycols[1], ylab = "Total Time in seconds for the chromosome", xlab = "Reference panel size", ylim = c(min(times) * 0.9, max(times) * 1.1))
lines(groups, times[, 2], type = "b", lwd = 1.0, pch = 1, col = mycols[2])
lines(groups, times[, 3], type = "b", lwd = 1.0, pch = 1, col = mycols[3])
legend("topleft", legend = c("GLIMPSE", "QUILT-regular", "QUILT-mspbwt"), col = mycols, pch = 1, lwd = 1.5, cex = 1.1, xjust = 0, yjust = 1, bty = "n")

plot(groups, rams[, 1], type = "b", lwd = 1.0, pch = 1, col = mycols[1], ylab = "Maximum RAM in MBs for the chromosome", xlab = "Reference panel size", ylim = c(min(rams) * 0.9, max(rams) * 1.1))
lines(groups, rams[, 2], type = "b", lwd = 1.0, pch = 1, col = mycols[2])
lines(groups, rams[, 3], type = "b", lwd = 1.0, pch = 1, col = mycols[3])
dev.off()
