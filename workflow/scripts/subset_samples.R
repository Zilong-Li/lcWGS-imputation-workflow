
ql <- paste("query", "-l", snakemake@params[["vcf"]])
size <- as.integer(snakemake@wildcards[["size"]])
allsamples <- as.character(system2("bcftools", ql, stdout = TRUE))
targetsamples <- snakemake@params[["samples"]]

exclude <- snakemake@params[["exclude"]]

if(file.exists(exclude) & exclude != "false") {
  samples <- as.character(read.table(exclude)[,1])
  targetsamples <- c(samples, targetsamples)
}

# remove target sample from the panel
allsamples <- allsamples[!allsamples %in% targetsamples]
if (size == 0) {
  subsets <- allsamples
} else {
  # random sample N pairs haplotypes
  subsets <- allsamples[sort(sample(1:length(allsamples), size))]
}
cat(subsets, file = snakemake@output[[1]], sep = "\n")

