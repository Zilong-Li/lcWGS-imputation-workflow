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
    starts, ends = [], []
    for i in range(n):
        ps = chunksize * i + s
        pe = chunksize * (i + 1) + s - 1
        pe = e if pe > e else pe
        starts.append(ps)
        ends.append(pe)
    return starts, ends


def get_samples_list_comma(wildcards):
    samples_target = SAMPLES.keys()
    size = int(wildcards.size)
    if size == 0:
        return "^" + ",".join(samples_target)
    else:
        samples_all = (
            os.popen(f"bcftools query -l { REFPANEL[wildcards.chrom]['vcf'] }")
            .read()
            .split("\n")
        )
        [samples_all.remove(i) for i in samples_target]
        samples_subset = random.sample(samples_all, size)
        return ",".join(samples_subset)


def get_quilt_output_regular(wildcards):
    starts, ends = get_regions_list_per_chrom(
        wildcards.chrom, config["quilt"]["chunksize"]
    )
    return expand(
        rules.quilt_run_regular.output, zip, start=starts, end=ends, allow_missing=True
    )


def get_quilt_output_mspbwt(wildcards):
    starts, ends = get_regions_list_per_chrom(
        wildcards.chrom, config["quilt"]["chunksize"]
    )
    return expand(
        rules.quilt_run_mspbwt.output, zip, start=starts, end=ends, allow_missing=True
    )


def get_quilt_output_zilong(wildcards):
    starts, ends = get_regions_list_per_chrom(
        wildcards.chrom, config["quilt"]["chunksize"]
    )
    return expand(
        rules.quilt_run_zilong.output, zip, start=starts, end=ends, allow_missing=True
    )


def get_glimpse_chunks(wildcards):
    chunks = (
        os.popen(
            f"GLIMPSE_chunk --input {REFPANEL[wildcards.chrom]['vcf']} --reference {REFPANEL[wildcards.chrom]['vcf']} --region {wildcards.chrom} --window-size {config['glimpse']['chunksize']} --buffer-size {config['glimpse']['buffer']} --output {REFPANEL[wildcards.chrom]['vcf']}.chunks && cat {REFPANEL[wildcards.chrom]['vcf']}.chunks"
        )
        .read()
        .split("\n")
    )
    d = dict()
    for chunk in chunks:
        tmp = chunk.split("\t")
        d[tmp[0]] = {"irg": tmp[2], "org": tmp[3]}
    return d


def get_glimpse_outputs(wildcards):
    d = get_glimpse_chunks(wildcards)
    return expand(rules.glimpse_phase.output, chunkid=d.keys(), allow_missing=True)


def collect_quilt_log_regular(wildcards):
    starts, ends = get_regions_list_per_chrom(
        wildcards.chrom, config["quilt"]["chunksize"]
    )
    return expand(
        rules.quilt_run_regular.log, zip, start=starts, end=ends, allow_missing=True
    )


def collect_quilt_log_mspbwt(wildcards):
    starts, ends = get_regions_list_per_chrom(
        wildcards.chrom, config["quilt"]["chunksize"]
    )
    return expand(
        rules.quilt_run_mspbwt.log, zip, start=starts, end=ends, allow_missing=True
    )


def collect_quilt_log_zilong(wildcards):
    starts, ends = get_regions_list_per_chrom(
        wildcards.chrom, config["quilt"]["chunksize"]
    )
    return expand(
        rules.quilt_run_zilong.log, zip, start=starts, end=ends, allow_missing=True
    )
