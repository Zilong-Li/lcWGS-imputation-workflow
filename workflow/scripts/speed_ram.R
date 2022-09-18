
gettimes <- function(ss) {
   sapply(strsplit(ss, ':'), function(s){s=as.numeric(s);n=length(s); sum(sapply(1:n, function(i){s[i]*60^(n-i)}))})
}

quilttime <- function(dl) {
    sapply(dl, function(d) {
        sum(gettimes(d[,1]))
    })
}

quiltram <- function(dl) {
    sapply(dl, function(d) {
        max(d[,2])
    })
}

deps <- as.numeric(snakemake@config[["downsample"]])

dl.regular <- lapply(snakemake@input[["regular"]], read.table)
dl.mspbwt <- lapply(snakemake@input[["mspbwt"]], read.table)
dl.zilong <- lapply(snakemake@input[["zilong"]], read.table)

times <- cbind(quilttime(dl.regular), quilttime(dl.mspbwt), quilttime(dl.zilong))
rams <- cbind(quiltram(dl.regular), quiltram(dl.mspbwt), quiltram(dl.zilong))


# ss = c("2:07.71","1:2:07.71")
# lapply(strsplit(ss, ':'), function(s){s=as.numeric(s);n=length(s); sum(sapply(1:n, function(i){s[i]*60^(n-i)}))})

pdf(snakemake@output[[1]], w=12, h=6)
par(mfrow=c(1,2))
plot(deps, times[,1], type = "b", lwd=1.5, pch = 1, ylab = "Total Time in seconds for chr21", xlab = "Samples sequencing depth", ylim=c(min(times)*0.9, max(times)*1.1))
lines(deps, times[,2], type = "b", lwd=1.5, pch = 1, col = "orange")
lines(deps, times[,3], type = "b", lwd=1.5, pch = 1, col = "red")
legend("topleft", legend=c("QUILT-regular", "QUILT-mspbwt", "QUILT-zilong"), col=c("black", "orange", "red"), pch = c(1,1,1), lwd = 1.5, cex = 1.1, xjust = 0, yjust = 1, bty = "n")

plot(deps, rams[,1], type = "b", lwd=1.5, pch = 1, ylab = "Maximum RAM in KBs for chr21", xlab = "Samples sequencing depth", ylim=c(min(rams)*0.9, max(rams)*1.1))
lines(deps, rams[,2], type = "b", lwd=1.5, pch = 1, col = "orange")
lines(deps, rams[,3], type = "b", lwd=1.5, pch = 1, col = "red")
dev.off()
