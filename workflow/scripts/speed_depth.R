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

dl.regular <- lapply(snakemake@input[["regular"]], read.table)
dl.zilong <- lapply(snakemake@input[["zilong"]], read.table)
dl.glimpse1 <- lapply(snakemake@input[["glimpse1"]], read.table)
dl.glimpse2 <- lapply(snakemake@input[["glimpse2"]], read.table)
rds <- list(QUILT2 = dl.zilong, GLIMPSE2 = dl.glimpse2, QUILT1 = dl.regular, GLIMPSE1 = dl.glimpse1)
rds <- lapply(rds,function(l) {names(l) <- paste0("depth=",groups,"x"); l} )

saveRDS(rds, snakemake@output[["rds"]])
rds <- readRDS(snakemake@output[["rds"]])


times <- data.frame(QUILT2 = gnutime(rds$QUILT2), GLIMPSE2 = gnutime(rds$GLIMPSE2), QUILT1 = gnutime(rds$QUILT1), GLIMPSE1 = gnutime(rds$GLIMPSE1))
rownames(times) <- groups
rams <- data.frame(QUILT2 = gunram(rds$QUILT2), GLIMPSE2 = gunram(rds$GLIMPSE2), QUILT1 = gunram(rds$QUILT1), GLIMPSE1 = gunram(rds$GLIMPSE1))
rownames(rams) <- groups

mycols <- c("#e69f00", "#d55e00", "#56b4e9", "#cc79a7", "#009e73", "#0072b2", "#f0e442")
palette(mycols)

pdf(paste0(snakemake@output[["rds"]], ".pdf"), w = 12, h = 6)
par(mfrow = c(1, 2))
barplot(t(times), beside = T, col = 1:4, ylab = "Runtime in Seconds", xlab = "Sequencing depth")
legend("topleft", legend = colnames(times), fill =  1:4)
barplot(t(rams) / 1024, beside = T, col = 1:4, ylab = "Maximun RAM in GBs", xlab = "Sequencing depth")
dev.off()
