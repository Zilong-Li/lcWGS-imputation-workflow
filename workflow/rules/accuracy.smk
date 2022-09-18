
OUTDIR_REPORT = "results/reports"


rule collect_imputed_gts:
    input:
        regular=rules.quilt_ligate_regular.output,
        mspbwt=rules.quilt_ligate_mspbwt.output,
        zilong=rules.quilt_ligate_zilong.output,
    output:
        truth=os.path.join(
            OUTDIR_REPORT, "truth.gts.{chrom}.panelsize{size}.down{depth}x.{chrom}.txt"
        ),
        regular=os.path.join(
            OUTDIR_REPORT, "quilt.gts.regular.panelsize{size}.down{depth}x.{chrom}.txt"
        ),
        mspbwt=os.path.join(
            OUTDIR_REPORT, "quilt.gts.mspbwt.panelsize{size}.down{depth}x.{chrom}.txt"
        ),
        zilong=os.path.join(
            OUTDIR_REPORT, "quilt.gts.zilong.panelsize{size}.down{depth}x.{chrom}.txt"
        ),
    log:
        os.path.join(
            OUTDIR_REPORT, "accuracy.panelsize{size}.down{depth}x.{chrom}.llog"
        ),
    params:
        truth=lambda wildcards: REFPANEL[wildcards.chrom]['truth'],
        samples=SAMPLES.keys(),
    conda:
        "../envs/quilt.yaml"
    shell:
        """
        bcftools query -f '%CHROM:%POS:%REF:%ALT\t[\t%GT]\n' -s {params.samples} {params.truth} | sed -E 's/\/|\|/\\t/g' > {output.truth}
        bcftools query -f '%CHROM:%POS:%REF:%ALT\t[\t%GT\t%DS]\n' -s {params.samples} {input.regular} | sed -E 's/\/|\|/\\t/g' > {output.regular}
        bcftools query -f '%CHROM:%POS:%REF:%ALT\t[\t%GT\t%DS]\n' -s {params.samples} {input.mspbwt} | sed -E 's/\/|\|/\\t/g' > {output.mspbwt}
        bcftools query -f '%CHROM:%POS:%REF:%ALT\t[\t%GT\t%DS]\n' -s {params.samples} {input.zilong} | sed -E 's/\/|\|/\\t/g' > {output.zilong}
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
        os.path.join(OUTDIR_REPORT, "quilt.accuracy.panelsize{size}.{chrom}.pdf"),
    log:
        os.path.join(OUTDIR_REPORT, "quilt.accuracy.panelsize{size}.{chrom}.pdf.llog"),
    conda:
        "../envs/quilt.yaml"
    script:
        "../scripts/plot_accuracy.R"
