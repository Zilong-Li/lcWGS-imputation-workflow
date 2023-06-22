
rule collect_quilt_regular_speed_log:
    input:
        collect_quilt_log_regular_region2,
    output:
        os.path.join(
            OUTDIR_SUMMARY,
            "speed.quilt.regular.panelsize{size}.down{depth}x.{chrom}.txt",
        ),
    log:
        os.path.join(
            OUTDIR_SUMMARY,
            "speed.quilt.regular.panelsize{size}.down{depth}x.{chrom}.txt.log",
        ),
    params:
        N="collect_quilt_regular_speed_log",
    conda:
        "../envs/quilt.yaml"
    shell:
        """
        echo {input} | tr ' ' '\\n' | xargs grep -E 'Elaps|Maximum' | awk '{{print $NF}}' | sed 'N;s/\\n/ /' > {output}
        """


rule collect_quilt_mspbwt_speed_log:
    input:
        collect_quilt_log_mspbwt_region2,
    output:
        os.path.join(
            OUTDIR_SUMMARY,
            "speed.quilt.mspbwt.panelsize{size}.down{depth}x.{chrom}.txt",
        ),
    log:
        os.path.join(
            OUTDIR_SUMMARY,
            "speed.quilt.mspbwt.panelsize{size}.down{depth}x.{chrom}.txt.log",
        ),
    params:
        N="collect_quilt_mspbwt_speed_log",
    conda:
        "../envs/quilt.yaml"
    shell:
        """
        echo {input} | tr ' ' '\\n' | xargs grep -E 'Elaps|Maximum' | awk '{{print $NF}}' | sed 'N;s/\\n/ /' > {output}
        """


rule collect_quilt_zilong_speed_log:
    input:
        collect_quilt_log_zilong_region2,
    output:
        os.path.join(
            OUTDIR_SUMMARY,
            "speed.quilt.zilong.panelsize{size}.down{depth}x.{chrom}.txt",
        ),
    log:
        os.path.join(
            OUTDIR_SUMMARY,
            "speed.quilt.zilong.panelsize{size}.down{depth}x.{chrom}.txt.log",
        ),
    params:
        N="collect_quilt_zilong_speed_log",
    conda:
        "../envs/quilt.yaml"
    shell:
        """
        echo {input} | tr ' ' '\\n' | xargs grep -E 'Elaps|Maximum' | awk '{{print $NF}}' | sed 'N;s/\\n/ /' > {output}
        """


rule collect_glimpse2_speed_log:
    input:
        collect_glimpse2_log,
    output:
        os.path.join(
            OUTDIR_SUMMARY,
            "speed.glimpse2.panelsize{size}.down{depth}x.{chrom}.txt",
        ),
    log:
        os.path.join(
            OUTDIR_SUMMARY,
            "speed.glimpse2.panelsize{size}.down{depth}x.{chrom}.txt.log",
        ),
    params:
        N="collect_glimpse2_speed_log",
    conda:
        "../envs/quilt.yaml"
    shell:
        """
        echo {input} | tr ' ' '\\n' | xargs grep -E 'Elaps|Maximum' | awk '{{print $NF}}' | sed 'N;s/\\n/ /' > {output}
        """


rule collect_glimpse_speed_log:
    input:
        collect_glimpse_log,
    output:
        os.path.join(
            OUTDIR_SUMMARY,
            "speed.glimpse.panelsize{size}.down{depth}x.{chrom}.txt",
        ),
    log:
        os.path.join(
            OUTDIR_SUMMARY,
            "speed.glimpse.panelsize{size}.down{depth}x.{chrom}.txt.log",
        ),
    params:
        N="collect_glimpse_speed_log",
    conda:
        "../envs/quilt.yaml"
    shell:
        """
        echo {input} | tr ' ' '\\n' | xargs grep -E 'Elaps|Maximum' | awk '{{print $NF}}' | sed 'N;s/\\n/ /' > {output}
        """


rule plot_speed_quilt_regular:
    input:
        expand(
            rules.collect_quilt_regular_speed_log.output,
            depth=config["downsample"],
            allow_missing=True,
        ),
    output:
        pdf=os.path.join(
            OUTDIR_SUMMARY, "quilt.speed.regular.panelsize{size}.{chrom}.pdf"
        ),
        rds=os.path.join(
            OUTDIR_SUMMARY, "quilt.speed.regular.panelsize{size}.{chrom}.rds"
        ),
    log:
        os.path.join(OUTDIR_SUMMARY, "quilt.speed.regular.panelsize{size}.{chrom}.llog"),
    params:
        N="plot_speed_qulit_regular",
    conda:
        "../envs/quilt.yaml"
    script:
        "../scripts/speed_single.R"


rule plot_speed_quilt_zilong:
    input:
        expand(
            rules.collect_quilt_zilong_speed_log.output,
            depth=config["downsample"],
            allow_missing=True,
        ),
    output:
        pdf=os.path.join(
            OUTDIR_SUMMARY, "quilt.speed.zilong.panelsize{size}.{chrom}.pdf"
        ),
        rds=os.path.join(
            OUTDIR_SUMMARY, "quilt.speed.zilong.panelsize{size}.{chrom}.rds"
        ),
    log:
        os.path.join(OUTDIR_SUMMARY, "quilt.speed.zilong.panelsize{size}.{chrom}.llog"),
    params:
        N="plot_speed_qulit_zilong",
    conda:
        "../envs/quilt.yaml"
    script:
        "../scripts/speed_single.R"


rule plot_speed_quilt_mspbwt:
    input:
        expand(
            rules.collect_quilt_mspbwt_speed_log.output,
            depth=config["downsample"],
            allow_missing=True,
        ),
    output:
        pdf=os.path.join(
            OUTDIR_SUMMARY, "quilt.speed.mspbwt.panelsize{size}.{chrom}.pdf"
        ),
        rds=os.path.join(
            OUTDIR_SUMMARY, "quilt.speed.mspbwt.panelsize{size}.{chrom}.rds"
        ),
    log:
        os.path.join(OUTDIR_SUMMARY, "quilt.speed.mspbwt.panelsize{size}.{chrom}.llog"),
    params:
        N="plot_speed_qulit_mspbwt",
    conda:
        "../envs/quilt.yaml"
    script:
        "../scripts/speed_single.R"


rule plot_speed_glimpse2:
    input:
        expand(
            rules.collect_glimpse2_speed_log.output,
            depth=config["downsample"],
            allow_missing=True,
        ),
    output:
        pdf=os.path.join(OUTDIR_SUMMARY, "glimpse2.speed.panelsize{size}.{chrom}.pdf"),
        rds=os.path.join(OUTDIR_SUMMARY, "glimpse2.speed.panelsize{size}.{chrom}.rds"),
    log:
        os.path.join(OUTDIR_SUMMARY, "glimpse2.speed.panelsize{size}.{chrom}.llog"),
    params:
        N="plot_speed_glimpse2",
    conda:
        "../envs/quilt.yaml"
    script:
        "../scripts/speed_single.R"


rule plot_speed_glimpse:
    input:
        expand(
            rules.collect_glimpse_speed_log.output,
            depth=config["downsample"],
            allow_missing=True,
        ),
    output:
        pdf=os.path.join(OUTDIR_SUMMARY, "glimpse.speed.panelsize{size}.{chrom}.pdf"),
        rds=os.path.join(OUTDIR_SUMMARY, "glimpse.speed.panelsize{size}.{chrom}.rds"),
    log:
        os.path.join(OUTDIR_SUMMARY, "glimpse.speed.panelsize{size}.{chrom}.llog"),
    params:
        N="plot_speed_glimpse",
    conda:
        "../envs/quilt.yaml"
    script:
        "../scripts/speed_single.R"


rule plot_speed_all:
    input:
        glimpse1=expand(
            rules.collect_glimpse_speed_log.output,
            size=config["refsize"],
            allow_missing=True,
        ),
        glimpse2=expand(
            rules.collect_glimpse2_speed_log.output,
            size=config["refsize"],
            allow_missing=True,
        ),
        regular=expand(
            rules.collect_quilt_regular_speed_log.output,
            size=config["refsize"],
            allow_missing=True,
        ),
        zilong=expand(
            rules.collect_quilt_zilong_speed_log.output,
            size=config["refsize"],
            allow_missing=True,
        ),
    output:
        pdf=os.path.join(OUTDIR_SUMMARY, "all.speed.down{depth}x.{chrom}.pdf"),
        rds=os.path.join(OUTDIR_SUMMARY, "all.speed.down{depth}x.{chrom}.rds"),
    log:
        os.path.join(OUTDIR_SUMMARY, "all.speed.down{depth}x.{chrom}.pdf.llog"),
    params:
        N="plot_speed_all",
    conda:
        "../envs/quilt.yaml"
    script:
        "../scripts/speed_panelsize.R"
