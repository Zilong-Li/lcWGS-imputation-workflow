import os
import random
import pandas as pd

# dict : {"NA12878": {"bam":ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239334/NA12878.final.cram, "depth": 99999}, ...}
SAMPLES = (
    pd.read_csv(config["samples"], sep="\t", dtype=str)
    .set_index("sampleid")
    .to_dict(orient="index")
)

# dict : {"chr1": {"vcf": chr1.vcf.gz,"start":1, "end": 99999}, ...}
REFPANEL = (
    pd.read_csv(config["genome"]["refpanel"], sep="\t", dtype=str)
    .set_index("chr")
    .to_dict(orient="index")
)


def get_regions_list_per_chrom(chrom, chunksize):
    """split chr into chunks given chunksize; return a list of '[start,end]' pairs"""
    regions = []
    s, e = int(REFPANEL[chrom]["start"]), int(REFPANEL[chrom]["end"])
    n = int((e - s) / chunksize) + 1
    if (n - 1) * chunksize == e - s:
        n = n - 1
    for i in range(n):
        ps = chunksize * i + s
        pe = chunksize * (i + 1) + s - 1
        pe = e if pe > e else pe
        regions.append([ps, pe])
    return regions


def get_samples_list_comma(wildcards):
    samples_all = (
        os.popen(f"bcftools query -l { REFPANEL[wildcards.chrom]['vcf'] }")
        .read()
        .split("\n")
    )
    samples_target = SAMPLES.keys()
    if wildcards.size == 0:
        return "^" + ",".join(samples_target)
    else:
        [samples_all.remove(i) for i in samples_all]
        samples_subset = random.sample(samples_all, int(wildcards.size))
        return ",".join(samples_subset)


def get_quilt_output_regular(wildcards):
    regions = get_regions_list_per_chrom(wildcards.chrom, config["quilt"]["chunksize"])
    return [
        f"results/quilt/panelsize{wildcards.size}/{wildcards.chrom}/quilt.{wildcards.depth}x.regular.{wildcards.chrom}.{start}.{end}.vcf.gz"
        for start, end in regions
    ]


def get_quilt_output_mspbwt(wildcards):
    regions = get_regions_list_per_chrom(wildcards.chrom, config["quilt"]["chunksize"])
    return [
        f"results/quilt/panelsize{wildcards.size}/{wildcards.chrom}/quilt.{wildcards.depth}x.mspbwt.{wildcards.chrom}.{start}.{end}.vcf.gz"
        for start, end in regions
    ]


def get_quilt_output_zilong(wildcards):
    regions = get_regions_list_per_chrom(wildcards.chrom, config["quilt"]["chunksize"])
    return [
        f"results/quilt/panelsize{wildcards.size}/{wildcards.chrom}/quilt.{wildcards.depth}x.zilong.{wildcards.chrom}.{start}.{end}.vcf.gz"
        for start, end in regions
    ]
