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
rds <- lapply(rds,function(l) {names(l) <- paste0(groups,"x"); l} )

saveRDS(rds, snakemake@output[["rds"]])
rds <- readRDS(snakemake@output[["rds"]])

times <- lapply(1:length(groups), function(i) {
  lapply(rds, function(l) {
    gettimes(l[[i]][,1]) / 60 ## minutes
  })
})
names(times) <- groups

rams <- lapply(1:length(groups), function(i) {
  lapply(rds, function(l) {
    l[[i]][,2] / 1024**2
  })
})
names(rams) <- groups


mycols <- c("#e69f00", "#d55e00", "#56b4e9", "#cc79a7", "#009e73", "#0072b2", "#f0e442")
palette(mycols)

png(paste0(snakemake@output[["rds"]], ".png"), w = 12, h = 6, units="in", res=300)
for(i in 1:length(groups)) {
  par(mfrow = c(1, 2))
  boxplot(times[[i]], col = 1:4, ylab = "Runtime in Minutes", main = paste("Sequencing depth ", groups[i], "x"))
  legend("topleft", legend = names(times[[i]]), fill =  1:4)
  boxplot(rams[[i]], col = 1:4, ylab = "RAM in GBs", main = paste("Sequencing depth ", groups[i], "x"))
}
dev.off()
