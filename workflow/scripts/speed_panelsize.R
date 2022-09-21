gettimes <- function(ss) {
   sapply(strsplit(ss, ':'), function(s){s=as.numeric(s);n=length(s); sum(sapply(1:n, function(i){s[i]*60^(n-i)}))})
}

gnutime <- function(dl) {
    sapply(dl, function(d) {
        sum(gettimes(d[,1]))
    })
}

gunram <- function(dl) {
    sapply(dl, function(d) {
        max(d[,2]) / 1024 # MB units
    })
}

groups <- as.numeric(snakemake@config[["refsize"]])
nd <- length(groups)

dl.regular <- lapply(snakemake@input[["regular"]], read.table)
dl.mspbwt <- lapply(snakemake@input[["mspbwt"]], read.table)
dl.zilong <- lapply(snakemake@input[["zilong"]], read.table)
dl.glimpse <- lapply(snakemake@input[["glimpse"]], read.table)
saveRDS(list(dl.regular, dl.glimpse), snakemake@output[["rds"]])

times <- cbind(gnutime(dl.regular), gnutime(dl.mspbwt), gnutime(dl.zilong), gnutime(dl.glimpse))
rams <- cbind(gunram(dl.regular), gunram(dl.mspbwt), gunram(dl.zilong), gunram(dl.glimpse))


wong <- c("#e69f00", "#d55e00", "#56b4e9", "#cc79a7", "#009e73", "#0072b2", "#f0e442")
mycols <- wong[1:4]


pdf(snakemake@output[["pdf"]], w=12, h=6)
par(mfrow=c(1,2))
plot(groups, times[,1], type = "b", lwd=1.0, pch = 1, col = mycols[1], ylab = "Total Time in seconds for the chromosome", xlab = "Reference panel size", ylim=c(min(times)*0.9, max(times)*1.1))
lines(groups, times[,2], type = "b", lwd=1.0, pch = 1, col = mycols[2])
lines(groups, times[,3], type = "b", lwd=1.0, pch = 1, col = mycols[3])
lines(groups, times[,4], type = "b", lwd=1.0, pch = 1, col = mycols[4])
legend("topleft", legend=c("QUILT-regular", "QUILT-mspbwt", "QUILT-zilong", "GLIMPSE"), col=mycols, pch = 1, lwd = 1.5, cex = 1.1, xjust = 0, yjust = 1, bty = "n")

plot(groups, rams[,1], type = "b", lwd=1.0, pch = 1, col = mycols[1], ylab = "Maximum RAM in MBs for the chromosome", xlab = "Reference panel size", ylim=c(min(rams)*0.9, max(rams)*1.1))
lines(groups, rams[,2], type = "b", lwd=1.0, pch = 1, col = mycols[2])
lines(groups, rams[,3], type = "b", lwd=1.0, pch = 1, col = mycols[3])
lines(groups, rams[,4], type = "b", lwd=1.0, pch = 1, col = mycols[4])
dev.off()
