
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

deps <- as.numeric(snakemake@config[["downsample"]])

dl.regular <- lapply(snakemake@input[["regular"]], read.table)
dl.mspbwt <- lapply(snakemake@input[["mspbwt"]], read.table)
dl.zilong <- lapply(snakemake@input[["zilong"]], read.table)
dl.glimpse <- lapply(snakemake@input[["glimpse"]], read.table)

times <- cbind(gnutime(dl.regular), gnutime(dl.mspbwt), gnutime(dl.zilong), gnutime(dl.glimpse))
rams <- cbind(gunram(dl.regular), gunram(dl.mspbwt), gunram(dl.zilong), gunram(dl.glimpse))


mycols <- c("black", "orange", "red", "blue" )

pdf(snakemake@output[[1]], w=12, h=6)
par(mfrow=c(1,2))
plot(deps, times[,1], type = "b", lwd=1.5, pch = 1, col = mycols[1], ylab = "Total Time in seconds for the chromosome", xlab = "Samples sequencing depth", ylim=c(min(times)*0.9, max(times)*1.1))
lines(deps, times[,2], type = "b", lwd=1.5, pch = 1, col = mycols[2])
lines(deps, times[,3], type = "b", lwd=1.5, pch = 1, col = mycols[3])
lines(deps, times[,4], type = "b", lwd=1.5, pch = 1, col = mycols[4])
legend("topleft", legend=c("QUILT-regular", "QUILT-mspbwt", "QUILT-zilong", "GLIMPSE"), col=mycols, pch = 1, lwd = 1.5, cex = 1.1, xjust = 0, yjust = 1, bty = "n")

plot(deps, rams[,1], type = "b", lwd=1.5, pch = 1, col = mycols[1], ylab = "Maximum RAM in MBs for the chromosome", xlab = "Samples sequencing depth", ylim=c(min(rams)*0.9, max(rams)*1.1))
lines(deps, rams[,2], type = "b", lwd=1.5, pch = 1, col = mycols[2])
lines(deps, rams[,3], type = "b", lwd=1.5, pch = 1, col = mycols[3])
lines(deps, rams[,4], type = "b", lwd=1.5, pch = 1, col = mycols[4])
dev.off()