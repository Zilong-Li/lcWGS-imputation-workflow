
snakemake@source("common.R")

acc_r2_all <- function(d0, d1) {
  id <- intersect(rownames(d0), rownames(d1))
  y1 <- cor(as.vector(d0[id,]), as.vector(d1[id,]), use = "pairwise.complete")**2
}

acc_r2_by_af <- function(d0, d1, af, bins) {
  id <- intersect(rownames(d0), rownames(d1))
  res <- r2_by_freq(breaks = bins, af[id], truthG = d0[id,], testDS = d1[id,], flip = TRUE)
  as.data.frame(cbind(bin = bins[-1], single = res[, "simple"], orphan = res[, "simple"]))
}

groups <- as.numeric(snakemake@config[["downsample"]])

truth <- fread(snakemake@input[["truth"]], data.table = F)
df.truth <- sapply(seq(1, dim(truth)[2] - 1, 2), function(i) {
  rowSums(truth[, (i + 1):(i + 2)])
}) # dosage matrix: nsnps x nsamples
rownames(df.truth) <- truth[,1]
truth <- truth[,-1] ## remove first id column
rownames(truth) <- rownames(df.truth)

d.af <- fread(snakemake@input[["af"]], data.table = F)
af <- as.numeric(d.af[, 2])
names(af) <- d.af[, 1]
rm(d.af)
af <- af[!is.na(af)]

## SNPs with (1-af) > 0.0005 & (1-af) < 0.001 are all imputed hom ALT and truth hom ALT. but those are stupidly easy to impute and don’t tell you anything
## af <- ifelse(af>0.5, 1-af, af)

dl.single <- lapply(snakemake@input[["single"]], parse.imputed.gts2)

bins <- sort(unique(c(
  c(0, 0.01 / 100, 0.02 / 100, 0.05 / 100),
  c(0, 0.01 / 10, 0.02 / 10, 0.05 / 10),
  c(0, 0.01 / 1, 0.02 / 1, 0.05 / 1),
  seq(0.1, 0.5, length.out = 5)
)))

phasing_errors <- lapply(seq(length(groups)), function(i) {
  id <- intersect(rownames(truth), rownames(dl.single[[i]]))
  single <- dl.single[[i]][,-seq(3, n, by = 3 )] # get phased genotypes
  ll <- acc_phasing_single_matrix(truth, single, id)
  ll
})

accuracy_by_af <- lapply(1:length(groups), function(i) {
  dl.single[[i]] <- dl.single[[i]][,seq(3, n, by = 3 )] # get dosages
  acc_r2_by_af(df.truth, dl.single[[i]], af, bins)
})

saveRDS(list(accuracy_by_af,phasing_errors), snakemake@output[["rds"]])

wong <- c("#e69f00", "#d55e00", "#56b4e9", "#cc79a7", "#009e73", "#0072b2", "#f0e442")
mycols <- wong[1:4]

pdf(snakemake@output[["pdf"]], w = 12, h = 6)
par(mfrow = c(1, 2))

phasing_errors <- matrix(sapply(phasing_errors, function(ls) {
  as.numeric(sapply(ls, "[[", "pse"))
}))

boxplot(phasing_errors, ylab = "PSE %", main = paste("Ref Panel Size: N=", snakemake@params[["N"]]))

a1 <- accuracy_by_af[[1]]
x <- a1$bin[!sapply(a1[, 2], is.na)] # remove AF bin with NA results
x <- log10(as.numeric(x))
labels <- 100 * bins[-1]
labels <- labels[!sapply(a1[, 2], is.na)]
ymin <- min(sapply(accuracy_by_af, function(d) {
  m <- as.matrix(apply(d[, -1], 2, unlist))
  min(m, na.rm = T)
}))

plot(1, col = "transparent", axes = F, xlim = c(min(x), max(x)), ylim = c(0, 1.0), ylab = "Aggregated R2 within each MAF bin", xlab = "Minor Allele Frequency")
nd <- length(groups)
for (i in 1:nd) {
  d <- accuracy_by_af[[i]]
  y <- rmna(d$single)
  lines(x, y, type = "l", lty = nd - i + 1, pch = 1, col = mycols[1])
}
axis(side = 1, at = x, labels = labels)
axis(side = 2)
legend("bottomright", legend = paste0(groups, "x"), lty = nd:1, bty = "n")

dev.off()
