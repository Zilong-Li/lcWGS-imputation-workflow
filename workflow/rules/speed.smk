
rule summary_speed:
    input:
        regular=collect_quilt_log_regular,
        mspbwt=collect_quilt_log_mspbwt,
        zilong=collect_quilt_log_zilong,
    output:
        regular=os.path.join(
            OUTDIR_SUMMARY, "quilt.regular.panelsize{size}.down{depth}x.{chrom}.txt"
        ),
        mspbwt=os.path.join(
            OUTDIR_SUMMARY, "quilt.mspbwt.panelsize{size}.down{depth}x.{chrom}.txt"
        ),
        zilong=os.path.join(
            OUTDIR_SUMMARY, "quilt.zilong.panelsize{size}.down{depth}x.{chrom}.txt"
        ),
    log:
        os.path.join(
            OUTDIR_SUMMARY, "quilt.panelsize{size}.down{depth}x.{chrom}.txt.log"
        ),
    conda:
        "../envs/quilt.yaml"
    shell:
        """
        echo {input.regular} | tr ' ' '\\n' | xargs grep -E 'Elaps|Maximum' | awk '{{print $NF}}' | sed 'N;s/\\n/ /' > {output.regular}
        echo {input.mspbwt} | tr ' ' '\\n' | xargs grep -E 'Elaps|Maximum' | awk '{{print $NF}}' | sed 'N;s/\\n/ /' > {output.mspbwt}
        echo {input.zilong} | tr ' ' '\\n' | xargs grep -E 'Elaps|Maximum' | awk '{{print $NF}}' | sed 'N;s/\\n/ /' > {output.zilong}
        """


rule plot_speed:
    input:
        regular=expand(
            rules.summary_speed.output.regular,
            depth=config["downsample"],
            allow_missing=True,
        ),
        mspbwt=expand(
            rules.summary_speed.output.mspbwt,
            depth=config["downsample"],
            allow_missing=True,
        ),
        zilong=expand(
            rules.summary_speed.output.zilong,
            depth=config["downsample"],
            allow_missing=True,
        ),
    output:
        os.path.join(OUTDIR_SUMMARY, "quilt.speed.panelsize{size}.{chrom}.pdf"),
    log:
        os.path.join(OUTDIR_SUMMARY, "quilt.speed.panelsize{size}.{chrom}.pdf.llog"),
    conda:
        "../envs/quilt.yaml"
    script:
        "../scripts/speed_ram.R"
