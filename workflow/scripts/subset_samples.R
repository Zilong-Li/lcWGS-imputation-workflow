
ql <- paste("query", "-l", snakemake@params[["vcf"]] )
size <- int(snakemake@wildcards[["size"]] )
allsamples <- system2("bcftools", ql, stdout = TRUE)
targesamples <- snakemake@params[["samples"]]
# remove target sample from the panel
allsamples <- allsamples[! allsamples %in% targesamples ]
if(size == 0) {
    subsets <- allsamples
} else {
    # random sample N pairs haplotypes
    subsets <- sample(1:length(allsamples), size)
}
cat(subsets, file = snakemake@output[[1]], sep = "\n")
