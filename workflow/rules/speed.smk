
rule collect_speed_log:
    input:
        regular=collect_quilt_log_regular,
        mspbwt=collect_quilt_log_mspbwt,
        zilong=collect_quilt_log_zilong,
        glimpse=collect_glimpse_log,
    output:
        glimpse=os.path.join(
            OUTDIR_SUMMARY,
            "speed.glimpse.regular.panelsize{size}.down{depth}x.{chrom}.txt",
        ),
        regular=os.path.join(
            OUTDIR_SUMMARY,
            "speed.quilt.regular.panelsize{size}.down{depth}x.{chrom}.txt",
        ),
        mspbwt=os.path.join(
            OUTDIR_SUMMARY,
            "speed.quilt.mspbwt.panelsize{size}.down{depth}x.{chrom}.txt",
        ),
        zilong=os.path.join(
            OUTDIR_SUMMARY,
            "speed.quilt.zilong.panelsize{size}.down{depth}x.{chrom}.txt",
        ),
    log:
        os.path.join(
            OUTDIR_SUMMARY, "speed.quilt.panelsize{size}.down{depth}x.{chrom}.txt.log"
        ),
    params:
        N="collect_speed_log",
    conda:
        "../envs/quilt.yaml"
    shell:
        """
        echo {input.regular} | tr ' ' '\\n' | xargs grep -E 'Elaps|Maximum' | awk '{{print $NF}}' | sed 'N;s/\\n/ /' > {output.regular}
        echo {input.mspbwt} | tr ' ' '\\n' | xargs grep -E 'Elaps|Maximum' | awk '{{print $NF}}' | sed 'N;s/\\n/ /' > {output.mspbwt}
        echo {input.zilong} | tr ' ' '\\n' | xargs grep -E 'Elaps|Maximum' | awk '{{print $NF}}' | sed 'N;s/\\n/ /' > {output.zilong}
        echo {input.glimpse} | tr ' ' '\\n' | xargs grep -E 'Elaps|Maximum' | awk '{{print $NF}}' | sed 'N;s/\\n/ /' > {output.glimpse}
        """


rule plot_speed:
    input:
        glimpse=expand(
            rules.collect_speed_log.output.glimpse,
            depth=config["downsample"],
            allow_missing=True,
        ),
        regular=expand(
            rules.collect_speed_log.output.regular,
            depth=config["downsample"],
            allow_missing=True,
        ),
        mspbwt=expand(
            rules.collect_speed_log.output.mspbwt,
            depth=config["downsample"],
            allow_missing=True,
        ),
        zilong=expand(
            rules.collect_speed_log.output.zilong,
            depth=config["downsample"],
            allow_missing=True,
        ),
    output:
        os.path.join(OUTDIR_SUMMARY, "speed.panelsize{size}.{chrom}.pdf"),
    log:
        os.path.join(OUTDIR_SUMMARY, "speed.panelsize{size}.{chrom}.pdf.llog"),
    params:
        N="plot_speed",
    conda:
        "../envs/quilt.yaml"
    script:
        "../scripts/speed_ram.R"
