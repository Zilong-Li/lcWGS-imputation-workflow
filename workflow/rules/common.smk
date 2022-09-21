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

RUN = config["scenario"]

OUTDIR = "results"
OUTDIR_DOWNSAMPLE = os.path.join(OUTDIR, "downsample", "")
OUTDIR_PANEL = os.path.join(OUTDIR, "subrefs", "")
OUTDIR_QUILT = os.path.join(OUTDIR, "quilt", "")
OUTDIR_GLIMPSE = os.path.join(OUTDIR, "glimpse", "")
OUTDIR_SUMMARY = os.path.join(OUTDIR, "summary", "")
OUTDIR_REPORT = os.path.join(OUTDIR, "report", "")


def get_speed_plots():
    return expand(
        rules.plot_speed.output, chrom=config["chroms"], depth=config["downsample"]
    )


def get_accuracy_panelsize_plots():
    return expand(
        rules.plot_accuracy_panelsize.output,
        chrom=config["chroms"],
        depth=config["downsample"],
    )


def get_accuracy_depth_plots():
    return expand(
        rules.plot_accuracy_depth.output, chrom=config["chroms"], size=config["refsize"]
    )


def get_quilt_regular_results():
    return expand(
        rules.quilt_ligate_regular.output,
        chrom=config["chroms"],
        size=config["refsize"],
        depth=config["downsample"],
    )


def get_quilt_mspbwt_results():
    return expand(
        rules.quilt_ligate_mspbwt.output,
        chrom=config["chroms"],
        size=config["refsize"],
        depth=config["downsample"],
    )


def get_quilt_zilong_results():
    return expand(
        rules.quilt_ligate_zilong.output,
        chrom=config["chroms"],
        size=config["refsize"],
        depth=config["downsample"],
    )


def get_glimpse_results():
    return expand(
        rules.glimpse_ligate.output,
        chrom=config["chroms"],
        size=config["refsize"],
        depth=config["downsample"],
    )


def if_use_af_in_refpanel(wildcards):
    if REFPANEL[wildcards.chrom].get("af"):
        return REFPANEL[wildcards.chrom]["af"]
    else:
        return ""


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
    """ugly but it's good to use GLIMPSE_chunk split the chromosome first"""
    if not os.path.exists(OUTDIR_SUMMARY):
        os.makedirs(OUTDIR_SUMMARY)
    fn = os.path.join(OUTDIR_SUMMARY, f"{wildcards.chrom}.glimpse.chunks")
    if not os.path.isfile(fn):
        os.system(
            f"GLIMPSE_chunk --input {REFPANEL[wildcards.chrom]['vcf']} --region {wildcards.chrom} --window-size {config['glimpse']['chunksize']} --buffer-size {config['glimpse']['buffer']} --output {fn} "
        )
    d = dict()
    with open(fn) as f:
        for row in f:
            tmp = row.split("\t")
            d[tmp[0]] = {"irg": tmp[2], "org": tmp[3]}
    return d


def get_glimpse_chunki_irg(wildcards):
    d = get_glimpse_chunks(wildcards)
    return d[wildcards.chunkid]["irg"]


def get_glimpse_chunki_org(wildcards):
    d = get_glimpse_chunks(wildcards)
    return d[wildcards.chunkid]["org"]


def get_glimpse_phase_outputs(wildcards):
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


def collect_glimpse_log(wildcards):
    d = get_glimpse_chunks(wildcards)
    return expand(rules.glimpse_phase.log, chunkid=d.keys(), allow_missing=True)
