
deps <- as.numeric(snakemake@config[["downsample"]])

dl <- lapply(snakemake@input, read.table)

times <- sapply(dl, function(d) {
    d[,1]
})

rams <- sapply(dl, function(d) {
    d[,2]
})

pdf(snakemake@output[[1]], w=12, h=6)
par(mfrow=c(1,2))
plot(deps,times[1,], type = "b", lwd=1.5, pch = 1, ylab = "Average Time in seconds", xlab = "Samples sequencing depth", ylim=c(min(times)*0.9, max(times)*1.1))
lines(deps,times[2,], type = "b", lwd=1.5, pch = 2, col = "orange")
lines(deps,times[3,], type = "b", lwd=1.5, pch = 5, col = "red")
plot(deps,rams[1,], type = "b", lwd=1.5, pch = 1, ylab = "Average RAM in KBs", xlab = "Samples sequencing depth", ylim=c(min(rams)*0.9, max(rams)*1.1))
lines(deps,rams[2,], type = "b", lwd=1.5, pch = 2, col = "orange")
lines(deps,rams[3,], type = "b", lwd=1.5, pch = 5, col = "red")
legend("bottomright", legend=c("QUILT-regular", "QUILT-mspbwt", "QUILT-zilong"), col=c("black", "orange", "red"), pch = c(1,2,5), lwd = 1.5, cex = 1.1, xjust = 0, yjust = 1, bty = "n")
dev.off()
