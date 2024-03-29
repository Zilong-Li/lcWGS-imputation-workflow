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

## saveRDS(snakemake,  snakemake@output[["rds"]])
## print(refsize0)
## q()

refsize0 <- as.integer(system(paste("bcftools query -l", snakemake@params$vcf, "|", "wc", "-l"), intern = TRUE))

groups <- as.numeric(snakemake@config[["refsize"]])
groups[groups == 0] <- refsize0
groups <- groups * 2
nd <- length(groups)

dl.regular <- lapply(snakemake@input[["regular"]], read.table)
dl.zilong <- lapply(snakemake@input[["zilong"]], read.table)
dl.glimpse1 <- lapply(snakemake@input[["glimpse1"]], read.table)
dl.glimpse2 <- lapply(snakemake@input[["glimpse2"]], read.table)
rds <- list(QUILT2 = dl.zilong, GLIMPSE2 = dl.glimpse2, QUILT1 = dl.regular, GLIMPSE1 = dl.glimpse1)
rds <- lapply(rds,function(l) {names(l) <- paste0("size=",groups); l} )

saveRDS(rds, snakemake@output[["rds"]])
rds <- readRDS(snakemake@output[["rds"]])


times <- cbind(gnutime(rds$QUILT2), gnutime(rds$GLIMPSE2), gnutime(rds$QUILT1), gnutime(rds$GLIMPSE1)) / 60
rams <- cbind(gunram(rds$QUILT2), gunram(rds$GLIMPSE2), gunram(rds$QUILT1), gunram(rds$GLIMPSE1)) / 1024


wong <- c("#e69f00", "#d55e00", "#56b4e9", "#cc79a7", "#009e73", "#0072b2", "#f0e442")
mycols <- wong


pdf(paste0(snakemake@output[["rds"]], ".pdf"), w = 12, h = 6)
par(mfrow = c(1, 2))
plot(groups, times[, 1], type = "b", lwd = 1.0, pch = 1, col = mycols[1], ylab = "Runtime in Minutes", xlab = "Reference panel size", ylim = c(min(times) * 0.9, max(times) * 1.1), log = 'y')
lines(groups, times[, 2], type = "b", lwd = 1.0, pch = 1, col = mycols[2])
lines(groups, times[, 3], type = "b", lwd = 1.0, pch = 1, col = mycols[3])
lines(groups, times[, 4], type = "b", lwd = 1.0, pch = 1, col = mycols[4])
legend("topleft", legend = names(rds), col = mycols, pch = 1, lwd = 1.5, cex = 1.1, xjust = 0, yjust = 1, bty = "n")

plot(groups, rams[, 1], type = "b", lwd = 1.0, pch = 1, col = mycols[1], ylab = "Maximum RAM in GBs", xlab = "Reference panel size", ylim = c(min(rams) * 0.9, max(rams) * 1.1))
lines(groups, rams[, 2], type = "b", lwd = 1.0, pch = 1, col = mycols[2])
lines(groups, rams[, 3], type = "b", lwd = 1.0, pch = 1, col = mycols[3])
lines(groups, rams[, 4], type = "b", lwd = 1.0, pch = 1, col = mycols[4])
dev.off()
