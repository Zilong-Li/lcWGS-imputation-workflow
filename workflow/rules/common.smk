import os
import pandas as pd

# dict : {"NA12878": {"bam":ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239334/NA12878.final.cram, "depth": 99999}, ...}

SAMPLES = (
    pd.read_csv(os.path.abspath(config["samples"]), sep="\t", dtype=str)
    .set_index("sampleid")
    .to_dict(orient="index")
)

# dict : {"chr1": {"vcf": chr1.vcf.gz,"start":1, "end": 99999}, ...}
REFPANEL = (
    pd.read_csv(os.path.abspath(config["genome"]["refpanel"]), sep="\t", dtype=str)
    .set_index("chr")
    .to_dict(orient="index")
)

OUTDIR = "results"
OUTDIR_DOWNSAMPLE = os.path.join(OUTDIR, "downsample", "")
OUTDIR_PANEL = os.path.join(OUTDIR, "refpanels", "")
OUTDIR_TRUTH = os.path.join(OUTDIR, "truth", "")
OUTDIR_QUILT1 = os.path.join(OUTDIR, "quilt1", "")
OUTDIR_QUILT2 = os.path.join(OUTDIR, "quilt2", "")
OUTDIR_GLIMPSE = os.path.join(OUTDIR, "glimpse1", "")
OUTDIR_GLIMPSE2 = os.path.join(OUTDIR, "glimpse2", "")
OUTDIR_SUMMARY = os.path.join(OUTDIR, "summary", "")
OUTDIR_REPORT = os.path.join(OUTDIR, "report", "")

###### global programs


def get_all_results():
    RUN = config["scenario"]
    if RUN == "all":
        return (
            get_speed_all_plots(),
            get_accuracy_panelsize_plots(),
            get_accuracy_depth_plots(),
        )
    elif RUN == "refpanel":
        return get_refpanel_sites()
    elif RUN == "accuracy":
        return get_accuracy_panelsize_plots(), get_accuracy_depth_plots()
    elif RUN == "F1":
        return get_accuracy_f1_plots()
    elif RUN == "V2":
        return get_speed_glimpse2_plots(), get_speed_quilt_mspbwt_plots()
    elif RUN == "speed":
        return get_speed_all_plots()
    elif RUN == "test":
        return get_quilt_accuracy()
    elif RUN == "quilt":
        return (
            get_quilt_regular_accuracy(),
            get_speed_quilt_regular_plots(),
            get_quilt_mspbwt_accuracy(),
            get_speed_quilt_mspbwt_plots(),
        )
    elif RUN == "fast":
        return (
            get_quilt_regular_accuracy(),
            get_speed_quilt_regular_plots(),
            get_quilt_mspbwt_accuracy(),
            get_speed_quilt_mspbwt_plots(),
            get_glimpse2_accuracy(),
            get_speed_glimpse2_plots(),
        )
    elif RUN == "quilt1":
        return get_quilt_regular_accuracy(), get_speed_quilt_regular_plots()
    elif RUN == "quilt2":
        return get_quilt_mspbwt_accuracy(), get_speed_quilt_mspbwt_plots()
    elif RUN == "glimpse1":
        return get_glimpse_accuracy(), get_speed_glimpse_plots()
    elif RUN == "glimpse2":
        return get_glimpse2_accuracy(), get_speed_glimpse2_plots()
    else:
        raise RuntimeError("this is an invalid scenario!")


def get_refpanel_sites():
    return expand(
        rules.concat_refpanel_sites_by_chunks.output,
        size=config["refsize"],
        chrom=config["chroms"],
    )


def get_speed_all_plots():
    return expand(
        rules.plot_speed_by_panelsize.output,
        chrom=config["chroms"],
        depth=config["downsample"],
    ), expand(
        rules.plot_speed_by_depth.output,
        chrom=config["chroms"],
        size=config["refsize"],
    )


def get_speed_quilt_zilong_plots():
    return expand(
        rules.plot_speed_quilt_zilong.output,
        chrom=config["chroms"],
        size=config["refsize"],
    )


def get_speed_quilt_regular_plots():
    return expand(
        rules.plot_speed_quilt_regular.output,
        chrom=config["chroms"],
        size=config["refsize"],
    )


def get_speed_quilt_mspbwt_plots():
    return expand(
        rules.plot_speed_quilt_mspbwt.output,
        chrom=config["chroms"],
        size=config["refsize"],
    )


def get_speed_glimpse2_plots():
    return expand(
        rules.plot_speed_glimpse2.output,
        chrom=config["chroms"],
        size=config["refsize"],
    )


def get_speed_glimpse_plots():
    return expand(
        rules.plot_speed_glimpse.output,
        chrom=config["chroms"],
        size=config["refsize"],
    )


def get_accuracy_panelsize_plots():
    return expand(
        rules.plot_accuracy_panelsize.output,
        chrom=config["chroms"],
        depth=config["downsample"],
    )


def get_accuracy_f1_plots():
    return expand(
        rules.plot_accuracy_f1.output, chrom=config["chroms"], size=config["refsize"]
    )


def get_accuracy_v2_plots():
    return expand(
        rules.plot_accuracy_v2.output, chrom=config["chroms"], size=config["refsize"]
    )


def get_accuracy_depth_plots():
    return expand(
        rules.plot_accuracy_depth.output, chrom=config["chroms"], size=config["refsize"]
    )


def get_quilt_regular_accuracy():
    return expand(
        rules.plot_quilt_regular.output,
        chrom=config["chroms"],
        size=config["refsize"],
    )


def get_quilt_mspbwt_accuracy():
    return expand(
        rules.plot_quilt_mspbwt.output,
        chrom=config["chroms"],
        size=config["refsize"],
    )


def get_quilt_accuracy():
    return expand(
        rules.plot_quilt_accuracy.output,
        chrom=config["chroms"],
        size=config["refsize"],
    )


def get_glimpse2_accuracy():
    return expand(
        rules.plot_glimpse2_accuracy.output,
        chrom=config["chroms"],
        size=config["refsize"],
    )


def get_glimpse_accuracy():
    return expand(
        rules.plot_glimpse_accuracy.output,
        chrom=config["chroms"],
        size=config["refsize"],
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


def get_glimpse_results():
    return expand(
        rules.glimpse_ligate.output,
        chrom=config["chroms"],
        size=config["refsize"],
        depth=config["downsample"],
    )


def if_exclude_samples_in_refpanel(wildcards):
    if REFPANEL[wildcards.chrom].get("exclude_samples"):
        return REFPANEL[wildcards.chrom]["exclude_samples"]
    else:
        return "false"


def if_use_af_in_refpanel(wildcards):
    if REFPANEL[wildcards.chrom].get("af"):
        # varify af file. 5 columns: chr pos ref alt af
        return REFPANEL[wildcards.chrom]["af"]
    else:
        return "false"


def if_use_quilt_map_in_refpanel(wildcards):
    if REFPANEL[wildcards.chrom].get("quilt_map"):
        # varify quilt genetic map file.
        return REFPANEL[wildcards.chrom]["quilt_map"]
    else:
        return "false"


def if_use_glimpse_map_in_refpanel(wildcards):
    if REFPANEL[wildcards.chrom].get("glimpse_map"):
        # varify glimpse genetic map file.
        return REFPANEL[wildcards.chrom]["glimpse_map"]
    else:
        return "false"


def get_chunks_by_chrom(chrom):
    """split chr into chunks given chunksize; return a dict of {id => chr:start-end}"""
    d = dict()
    chunksize = config["chunksize"]
    if REFPANEL[chrom].get("region"):
        d[0] = REFPANEL[chrom]["region"]
    else:
        s, e = int(REFPANEL[chrom]["start"]), int(REFPANEL[chrom]["end"])
        n = int((e - s) / chunksize) + 1
        if (n - 1) * chunksize == e - s:
            n = n - 1
        for i in range(n):
            ps = chunksize * i + s
            pe = chunksize * (i + 1) + s - 1
            pe = e if pe > e else pe
            d[i] = f"{chrom}:{ps}-{pe}"
    return d


def get_glimpse_chunks_out(fn):
    d = dict()
    with open(fn) as f:
        for row in f:
            """0       chr20   chr20:82590-6074391     chr20:82590-5574162     5491573 1893"""
            tmp = row.split("\t")
            d[int(tmp[0])] = tmp[3]
    return d


def get_glimpse_chunks_in(fn):
    d = dict()
    with open(fn) as f:
        for row in f:
            """0       chr20   chr20:82590-6074391     chr20:82590-5574162     5491573 1893"""
            tmp = row.split("\t")
            d[int(tmp[0])] = tmp[2]
    return d


def get_quilt_chunks(chrom):
    """use GLIMPSE_chunk if defined, otherwise split the chromosome by chunks"""
    d = dict()
    if REFPANEL[chrom].get("glimpse_chunk"):
        d = get_glimpse_chunks_out(REFPANEL[chrom]["glimpse_chunk"])
    else:
        d = get_chunks_by_chrom(chrom)
    return d


def get_quilt_chunk_region(chrom, chunkid):
    d = get_quilt_chunks(chrom)
    tmp2 = d[int(chunkid)].split(":")
    rg = tmp2[1].split("-")
    ps = max(int(rg[0]) - int(config["extra_buffer_in_quilt"]), 1)
    pe = int(rg[1]) + int(config["extra_buffer_in_quilt"])
    return f"{chrom}:{ps}-{pe}"


def get_quilt_chunk_region_start(wildcards):
    tmp = get_quilt_chunk_region(wildcards.chrom, wildcards.chunkid)
    tmp2 = tmp.split(":")
    rg = tmp2[1].split("-")
    return rg[0]


def get_quilt_chunk_region_end(wildcards):
    tmp = get_quilt_chunk_region(wildcards.chrom, wildcards.chunkid)
    tmp2 = tmp.split(":")
    rg = tmp2[1].split("-")
    return rg[1]


def get_quilt_regular_outputs(wildcards):
    d = get_quilt_chunks(wildcards.chrom)
    ids = list(map(str, d.keys()))
    return expand(rules.quilt_run_regular.output, chunkid=ids, allow_missing=True)


def get_quilt_mspbwt_outputs(wildcards):
    d = get_quilt_chunks(wildcards.chrom)
    ids = list(map(str, d.keys()))
    return expand(rules.quilt_run_mspbwt.output, chunkid=ids, allow_missing=True)


def collect_quilt_regular_logs(wildcards):
    d = get_quilt_chunks(wildcards.chrom)
    ids = list(map(str, d.keys()))
    return expand(rules.quilt_run_regular.log, chunkid=ids, allow_missing=True)


def collect_quilt_mspbwt_logs(wildcards):
    d = get_quilt_chunks(wildcards.chrom)
    ids = list(map(str, d.keys()))
    return expand(rules.quilt_run_mspbwt.log, chunkid=ids, allow_missing=True)


def get_refpanel_chunks(chrom):
    """use GLIMPSE_chunk if defined, otherwise split the chromosome by chunks"""
    extra_buffer_in_panel = 1000000
    d = dict()
    if REFPANEL[chrom].get("glimpse_chunk"):
        d = get_glimpse_chunks_in(REFPANEL[chrom]["glimpse_chunk"])
    else:
        d = get_chunks_by_chrom(chrom)
    return d


def get_refpanel_chunk_region(chrom, chunkid):
    """use GLIMPSE_chunk if defined, otherwise split the chromosome by chunks"""
    extra_buffer_in_panel = 1000000
    d = get_refpanel_chunks(chrom)
    tmp2 = d[int(chunkid)].split(":")
    rg = tmp2[1].split("-")
    ps = max(int(rg[0]) - extra_buffer_in_panel, 1)
    pe = int(rg[1]) + extra_buffer_in_panel
    return f"{chrom}:{ps}-{pe}"


def get_quilt_output_mspbwt_region1(wildcards):
    starts, ends = get_regions_list_per_chrom(
        wildcards.chrom, config["quilt2"]["chunksize"]
    )
    return expand(
        rules.quilt_run_mspbwt.output, zip, start=starts, end=ends, allow_missing=True
    )


def get_glimpse_chunks_in_and_out(chrom):
    """ugly but it's fair to use GLIMPSE_chunk split the chromosome first"""
    d = dict()
    if REFPANEL[chrom].get("region"):
        irg = f"{REFPANEL[chrom]['region']}"
        org = f"{REFPANEL[chrom]['region']}"
        d[0] = {"irg": irg, "org": org}
    else:
        if not os.path.exists(OUTDIR_PANEL):
            os.makedirs(OUTDIR_PANEL)
        fn = os.path.join(OUTDIR_PANEL, f"{chrom}.glimpse.chunks")
        if REFPANEL[chrom].get("glimpse_chunk"):
            fn = REFPANEL[chrom].get("glimpse_chunk")
        elif not os.path.isfile(fn):
            os.system(
                f"GLIMPSE_chunk --input {REFPANEL[chrom]['vcf']} --region {chrom} --window-size {config['chunksize']} --buffer-size {config['glimpse1']['buffer']} --output {fn} "
            )
        with open(fn) as f:
            for row in f:
                """0       chr20   chr20:82590-6074391     chr20:82590-5574162     5491573 1893"""
                tmp = row.split("\t")
                d[int(tmp[0])] = {"irg": tmp[2], "org": tmp[3]}
    return d


def get_glimpse_chunki_irg(wildcards):
    d = get_glimpse_chunks_in_and_out(wildcards.chrom)
    k = int(wildcards.chunkid)
    return d[k]["irg"]


def get_glimpse_chunki_org(wildcards):
    d = get_glimpse_chunks_in_and_out(wildcards.chrom)
    k = int(wildcards.chunkid)
    return d[k]["org"]


def get_glimpse_phase_outputs(wildcards):
    d = get_glimpse_chunks_in_and_out(wildcards.chrom)
    ids = list(map(str, d.keys()))
    return expand(rules.glimpse_phase.output, chunkid=ids, allow_missing=True)


def get_glimpse2_phase_outputs(wildcards):
    d = get_glimpse_chunks_in_and_out(wildcards.chrom)
    ids = list(map(str, d.keys()))
    return expand(rules.glimpse2_phase.output, chunkid=ids, allow_missing=True)


def collect_glimpse_log(wildcards):
    d = get_glimpse_chunks_in_and_out(wildcards.chrom)
    ids = list(map(str, d.keys()))
    return expand(rules.glimpse_phase.log, chunkid=ids, allow_missing=True)


def collect_glimpse2_log(wildcards):
    d = get_glimpse_chunks_in_and_out(wildcards.chrom)
    ids = list(map(str, d.keys()))
    return expand(rules.glimpse2_phase.log, chunkid=ids, allow_missing=True)


def collect_quilt_log_regular_region2(wildcards):
    starts, ends = get_regions_list_from_glimpse_chunk(wildcards.chrom, True)
    return expand(
        rules.quilt_run_regular.log, zip, start=starts, end=ends, allow_missing=True
    )


def collect_quilt_log_mspbwt_region1(wildcards):
    starts, ends = get_regions_list_per_chrom(
        wildcards.chrom, config["quilt2"]["chunksize"]
    )
    return expand(
        rules.quilt_run_mspbwt.log, zip, start=starts, end=ends, allow_missing=True
    )


def collect_quilt_log_mspbwt_region2(wildcards):
    starts, ends = get_regions_list_from_glimpse_chunk(wildcards.chrom, True)
    return expand(
        rules.quilt_run_mspbwt.log, zip, start=starts, end=ends, allow_missing=True
    )


def collect_quilt_log_zilong_region1(wildcards):
    starts, ends = get_regions_list_per_chrom(
        wildcards.chrom, config["quilt2"]["chunksize"]
    )
    return expand(
        rules.quilt_run_zilong.log, zip, start=starts, end=ends, allow_missing=True
    )


def collect_quilt_log_zilong_region2(wildcards):
    starts, ends = get_regions_list_from_glimpse_chunk(wildcards.chrom, True)
    return expand(
        rules.quilt_run_zilong.log, zip, start=starts, end=ends, allow_missing=True
    )
