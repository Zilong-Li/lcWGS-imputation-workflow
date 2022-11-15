
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

groups <- as.numeric(snakemake@config[["downsample"]])
nd <- length(groups)

dl.regular <- lapply(snakemake@input, read.table)
times <- data.frame(gnutime(dl.regular))
rams <- data.frame(gunram(dl.regular))

saveRDS(list(time = times, ram = rams), snakemake@output[["rds"]])


wong <- c("#e69f00", "#d55e00", "#56b4e9", "#cc79a7", "#009e73", "#0072b2", "#f0e442")
mycols <- wong[1:3]


pdf(snakemake@output[["pdf"]], w = 12, h = 6)
par(mfrow = c(1, 2))
plot(groups, times[, 1], type = "b", lwd = 1.0, pch = 1, col = mycols[1], ylab = "Total Time in seconds for the chromosome", xlab = "Sequencing depth", ylim = c(min(times) * 0.9, max(times) * 1.1))
legend("topleft", legend = c(snakemake@params[["N"]]), col = mycols, pch = 1, lwd = 1.5, cex = 1.1, xjust = 0, yjust = 1, bty = "n")

plot(groups, rams[, 1], type = "b", lwd = 1.0, pch = 1, col = mycols[1], ylab = "Maximum RAM in MBs for the chromosome", xlab = "Sequencing depth", ylim = c(min(rams) * 0.9, max(rams) * 1.1))
dev.off()
