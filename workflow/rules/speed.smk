
rule collect_quilt_regular_speed_log:
    input:
        collect_quilt_log_regular_region2,
    output:
        os.path.join(
            OUTDIR_SUMMARY,
            "speed.quilt.regular.refsize{size}.down{depth}x.{chrom}.txt",
        ),
    log:
        os.path.join(
            OUTDIR_SUMMARY,
            "speed.quilt.regular.refsize{size}.down{depth}x.{chrom}.txt.log",
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
            "speed.quilt.mspbwt.refsize{size}.down{depth}x.{chrom}.txt",
        ),
    log:
        os.path.join(
            OUTDIR_SUMMARY,
            "speed.quilt.mspbwt.refsize{size}.down{depth}x.{chrom}.txt.log",
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
            "speed.quilt.zilong.refsize{size}.down{depth}x.{chrom}.txt",
        ),
    log:
        os.path.join(
            OUTDIR_SUMMARY,
            "speed.quilt.zilong.refsize{size}.down{depth}x.{chrom}.txt.log",
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
            "speed.glimpse2.refsize{size}.down{depth}x.{chrom}.txt",
        ),
    log:
        os.path.join(
            OUTDIR_SUMMARY,
            "speed.glimpse2.refsize{size}.down{depth}x.{chrom}.txt.log",
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
            "speed.glimpse.refsize{size}.down{depth}x.{chrom}.txt",
        ),
    log:
        os.path.join(
            OUTDIR_SUMMARY,
            "speed.glimpse.refsize{size}.down{depth}x.{chrom}.txt.log",
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
            OUTDIR_SUMMARY, "quilt.speed.regular.refsize{size}.{chrom}.rds.pdf"
        ),
        rds=os.path.join(
            OUTDIR_SUMMARY, "quilt.speed.regular.refsize{size}.{chrom}.rds"
        ),
    log:
        os.path.join(
            OUTDIR_SUMMARY, "quilt.speed.regular.refsize{size}.{chrom}.rds.llog"
        ),
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
            OUTDIR_SUMMARY, "quilt.speed.zilong.refsize{size}.{chrom}.rds.pdf"
        ),
        rds=os.path.join(OUTDIR_SUMMARY, "quilt.speed.zilong.refsize{size}.{chrom}.rds"),
    log:
        os.path.join(OUTDIR_SUMMARY, "quilt.speed.zilong.refsize{size}.{chrom}.llog"),
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
            OUTDIR_SUMMARY, "quilt.speed.mspbwt.refsize{size}.{chrom}.rds.pdf"
        ),
        rds=os.path.join(OUTDIR_SUMMARY, "quilt.speed.mspbwt.refsize{size}.{chrom}.rds"),
    log:
        os.path.join(OUTDIR_SUMMARY, "quilt.speed.mspbwt.refsize{size}.{chrom}.llog"),
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
        pdf=os.path.join(OUTDIR_SUMMARY, "glimpse2.speed.refsize{size}.{chrom}.rds.pdf"),
        rds=os.path.join(OUTDIR_SUMMARY, "glimpse2.speed.refsize{size}.{chrom}.rds"),
    log:
        os.path.join(OUTDIR_SUMMARY, "glimpse2.speed.refsize{size}.{chrom}.llog"),
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
        pdf=os.path.join(OUTDIR_SUMMARY, "glimpse.speed.refsize{size}.{chrom}.rds.pdf"),
        rds=os.path.join(OUTDIR_SUMMARY, "glimpse.speed.refsize{size}.{chrom}.rds"),
    log:
        os.path.join(OUTDIR_SUMMARY, "glimpse.speed.refsize{size}.{chrom}.llog"),
    params:
        N="plot_speed_glimpse",
    conda:
        "../envs/quilt.yaml"
    script:
        "../scripts/speed_single.R"


rule plot_speed_by_panelsize:
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
            rules.collect_quilt_mspbwt_speed_log.output,
            size=config["refsize"],
            allow_missing=True,
        ),
    output:
        rds=os.path.join(OUTDIR_SUMMARY, "all.speed.down{depth}x.{chrom}.rds"),
    log:
        os.path.join(OUTDIR_SUMMARY, "all.speed.down{depth}x.{chrom}.rds.llog"),
    params:
        N="plot_speed_by_panelsize",
        vcf=lambda wildcards: REFPANEL[wildcards.chrom]["vcf"],
    conda:
        "../envs/quilt.yaml"
    script:
        "../scripts/speed_panelsize.R"


rule plot_speed_by_depth:
    input:
        glimpse1=expand(
            rules.collect_glimpse_speed_log.output,
            depth=config["downsample"],
            allow_missing=True,
        ),
        glimpse2=expand(
            rules.collect_glimpse2_speed_log.output,
            depth=config["downsample"],
            allow_missing=True,
        ),
        regular=expand(
            rules.collect_quilt_regular_speed_log.output,
            depth=config["downsample"],
            allow_missing=True,
        ),
        zilong=expand(
            rules.collect_quilt_mspbwt_speed_log.output,
            depth=config["downsample"],
            allow_missing=True,
        ),
    output:
        rds=os.path.join(OUTDIR_SUMMARY, "all.speed.refsize{size}.{chrom}.rds"),
    log:
        os.path.join(OUTDIR_SUMMARY, "all.speed.refsize{size}.{chrom}.rds.llog"),
    params:
        N="plot_speed_by_depth",
    conda:
        "../envs/quilt.yaml"
    script:
        "../scripts/speed_depth.R"
