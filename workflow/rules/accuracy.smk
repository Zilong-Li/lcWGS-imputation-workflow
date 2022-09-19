
rule collect_truth_gts:
    """would be better to use sites in subrefs"""
    input:
        sites=lambda wildcards: expand(
            rules.subset_refpanel.output.sites,
            size=config["refsize"],
            allow_missing=True,
        ),
    output:
        gt=os.path.join(OUTDIR_SUMMARY, "truth.gts.{chrom}.txt"),
        af=os.path.join(OUTDIR_SUMMARY, "af.input.panel.{chrom}.txt"),
    log:
        os.path.join(OUTDIR_SUMMARY, "truth.gts.{chrom}.log"),
    params:
        N="collect_truth_gts",
        samples=",".join(SAMPLES.keys()),
        truth=lambda wildcards: REFPANEL[wildcards.chrom]["truth"],
        ref=lambda wildcards: REFPANEL[wildcards.chrom]["vcf"],
        ql0="%CHROM:%POS:%REF:%ALT\\n",
        ql1="%CHROM:%POS:%REF:%ALT\\t%AF\\n",
        ql2="%CHROM:%POS:%REF:%ALT\\t[\\t%GT]\\n",
        awk="NR==FNR{a[$1]=1;} NR!=FNR{if(a[$1]) print $2;}",
    conda:
        "../envs/quilt.yaml"
    shell:
        """
        bcftools +fill-tags {params.ref} -- -t AF |  bcftools query -f '{params.ql1}' > {output.af}.tmp
        awk '{params.awk}' <(bcftools query -f '{params.ql0}' {input.sites[0]}) {output.af}.tmp >{output.af}
        bcftools view -s {params.samples} -T {input.sites[0]} {params.truth} | bcftools query -f '{params.ql2}' | sed -E 's/\/|\|/\\t/g' > {output.gt}
        """


rule collect_quilt_imputed_gts:
    input:
        regular=rules.quilt_ligate_regular.output.vcf,
        mspbwt=rules.quilt_ligate_mspbwt.output.vcf,
        zilong=rules.quilt_ligate_zilong.output.vcf,
    output:
        regular=os.path.join(
            OUTDIR_SUMMARY,
            "quilt.gts.regular.panelsize{size}.down{depth}x.{chrom}.txt",
        ),
        mspbwt=os.path.join(
            OUTDIR_SUMMARY, "quilt.gts.mspbwt.panelsize{size}.down{depth}x.{chrom}.txt"
        ),
        zilong=os.path.join(
            OUTDIR_SUMMARY, "quilt.gts.zilong.panelsize{size}.down{depth}x.{chrom}.txt"
        ),
    log:
        os.path.join(
            OUTDIR_SUMMARY, "quilt.gts.panelsize{size}.down{depth}x.{chrom}.llog"
        ),
    params:
        N="collect_quilt_imputed_gts",
        samples=",".join(SAMPLES.keys()),
        ql2="%CHROM:%POS:%REF:%ALT\\t[\\t%GT\\t%DS]\\n",
    conda:
        "../envs/quilt.yaml"
    shell:
        """
        bcftools query -f '{params.ql2}' -s {params.samples} {input.regular} | sed -E 's/\/|\|/\\t/g' > {output.regular}
        bcftools query -f '{params.ql2}' -s {params.samples} {input.mspbwt} | sed -E 's/\/|\|/\\t/g' > {output.mspbwt}
        bcftools query -f '{params.ql2}' -s {params.samples} {input.zilong} | sed -E 's/\/|\|/\\t/g' > {output.zilong}
        """


rule collect_glimpse_imputed_gts:
    input:
        rules.glimpse_ligate.output.vcf,
    output:
        os.path.join(
            OUTDIR_SUMMARY,
            "glimpse.gts.panelsize{size}.down{depth}x.{chrom}.txt",
        ),
    log:
        os.path.join(
            OUTDIR_SUMMARY, "glimpse.gts.panelsize{size}.down{depth}x.{chrom}.llog"
        ),
    params:
        N="collect_glimpse_imputed_gts",
        samples=",".join(SAMPLES.keys()),
        ql2="%CHROM:%POS:%REF:%ALT\\t[\\t%GT\\t%DS]\\n",
    conda:
        "../envs/quilt.yaml"
    shell:
        """
        bcftools query -f '{params.ql2}' -s {params.samples} {input} | sed -E 's/\/|\|/\\t/g' > {output}
        """


rule plot_quilt_accuracy:
    input:
        truth=rules.collect_truth_gts.output.gt,
        af=rules.collect_truth_gts.output.af,
        regular=expand(
            rules.collect_quilt_imputed_gts.output.regular,
            depth=config["downsample"],
            allow_missing=True,
        ),
        mspbwt=expand(
            rules.collect_quilt_imputed_gts.output.mspbwt,
            depth=config["downsample"],
            allow_missing=True,
        ),
        zilong=expand(
            rules.collect_quilt_imputed_gts.output.zilong,
            depth=config["downsample"],
            allow_missing=True,
        ),
    params:
        N="plot_quilt_accuracy",
    output:
        os.path.join(OUTDIR_SUMMARY, "quilt.accuracy.panelsize{size}.{chrom}.pdf"),
    log:
        os.path.join(OUTDIR_SUMMARY, "quilt.accuracy.panelsize{size}.{chrom}.pdf.llog"),
    conda:
        "../envs/quilt.yaml"
    script:
        "../scripts/accuracy.R"


rule plot_all_accuracy:
    input:
        truth=rules.collect_truth_gts.output.gt,
        af=rules.collect_truth_gts.output.af,
        glimpse=expand(
            rules.collect_glimpse_imputed_gts.output,
            depth=config["downsample"],
            allow_missing=True,
        ),
        regular=expand(
            rules.collect_quilt_imputed_gts.output.regular,
            depth=config["downsample"],
            allow_missing=True,
        ),
        mspbwt=expand(
            rules.collect_quilt_imputed_gts.output.mspbwt,
            depth=config["downsample"],
            allow_missing=True,
        ),
        zilong=expand(
            rules.collect_quilt_imputed_gts.output.zilong,
            depth=config["downsample"],
            allow_missing=True,
        ),
    output:
        report(
            os.path.join(OUTDIR_SUMMARY, "all.accuracy.panelsize{size}.{chrom}.pdf"),
            caption="report/accuracy.rst",
        ),
    log:
        os.path.join(OUTDIR_SUMMARY, "all.accuracy.panelsize{size}.{chrom}.pdf.llog"),
    params:
        N="plot_all_accuracy",
    conda:
        "../envs/quilt.yaml"
    script:
        "../scripts/accuracy2.R"
