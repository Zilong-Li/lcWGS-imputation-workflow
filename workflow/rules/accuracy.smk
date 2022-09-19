
rule collect_imputed_gts:
    input:
        regular=rules.quilt_ligate_regular.output.vcf,
        mspbwt=rules.quilt_ligate_mspbwt.output.vcf,
        zilong=rules.quilt_ligate_zilong.output.vcf,
    output:
        truth=os.path.join(
            OUTDIR_SUMMARY,
            "truth.gts.{chrom}.panelsize{size}.down{depth}x.{chrom}.txt",
        ),
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
            OUTDIR_SUMMARY, "accuracy.panelsize{size}.down{depth}x.{chrom}.llog"
        ),
    params:
        samples=",".join(SAMPLES.keys()),
        truth=lambda wildcards: REFPANEL[wildcards.chrom]["truth"],
        ql1="%CHROM:%POS:%REF:%ALT\\t%AF\\t[\\t%GT]\\n",
        ql2="%CHROM:%POS:%REF:%ALT\\t[\\t%GT\\t%DS]\\n",
    conda:
        "../envs/quilt.yaml"
    shell:
        """
        bcftools +fill-tags {params.truth} -- -t AF |  bcftools view -s {params.samples} | bcftools query -f '{params.ql1}' | sed -E 's/\/|\|/\\t/g' > {output.truth}
        bcftools query -f '{params.ql2}' -s {params.samples} {input.regular} | sed -E 's/\/|\|/\\t/g' > {output.regular}
        bcftools query -f '{params.ql2}' -s {params.samples} {input.mspbwt} | sed -E 's/\/|\|/\\t/g' > {output.mspbwt}
        bcftools query -f '{params.ql2}' -s {params.samples} {input.zilong} | sed -E 's/\/|\|/\\t/g' > {output.zilong}
        """


rule plot_accuracy:
    input:
        truth=expand(
            rules.collect_imputed_gts.output.truth,
            depth=config["downsample"],
            allow_missing=True,
        ),
        regular=expand(
            rules.collect_imputed_gts.output.regular,
            depth=config["downsample"],
            allow_missing=True,
        ),
        mspbwt=expand(
            rules.collect_imputed_gts.output.mspbwt,
            depth=config["downsample"],
            allow_missing=True,
        ),
        zilong=expand(
            rules.collect_imputed_gts.output.zilong,
            depth=config["downsample"],
            allow_missing=True,
        ),
    output:
        os.path.join(OUTDIR_SUMMARY, "quilt.accuracy.panelsize{size}.{chrom}.pdf"),
    log:
        os.path.join(OUTDIR_SUMMARY, "quilt.accuracy.panelsize{size}.{chrom}.pdf.llog"),
    conda:
        "../envs/quilt.yaml"
    script:
        "../scripts/accuracy.R"
