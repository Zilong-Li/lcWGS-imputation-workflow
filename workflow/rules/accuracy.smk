
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
        tmp=temp(os.path.join(OUTDIR_SUMMARY, "af.input.panel.{chrom}.txt.tmp")),
        tmp2=temp(os.path.join(OUTDIR_SUMMARY, "truth.gts.{chrom}.txt.tmp")),
    log:
        os.path.join(OUTDIR_SUMMARY, "truth.gts.{chrom}.log"),
    params:
        N="collect_truth_gts",
        samples=",".join(SAMPLES.keys()),
        truth=lambda wildcards: REFPANEL[wildcards.chrom]["truth"],
        ref=lambda wildcards: REFPANEL[wildcards.chrom]["vcf"],
        af=if_use_af_in_refpanel,
        ql0="%CHROM:%POS:%REF:%ALT\\n",
        ql1="%CHROM:%POS:%REF:%ALT\\t%AF\\n",
        ql2="%CHROM:%POS:%REF:%ALT\\t[\\t%GT]\\n",
        awk="NR==FNR{a[$1]=1;} NR!=FNR{if(a[$1]){print $2;}}",
        awk2="NR==FNR{a[$1]=1;} NR!=FNR{if(a[$1]){print $0;}}",
    conda:
        "../envs/quilt.yaml"
    shell:
        """
        (
        if [[ {params.af} == "" ]];then \
            bcftools +fill-tags {params.ref} -- -t AF |  bcftools query -f '{params.ql1}' > {output.tmp} \
        else \
            perl -lane 'print join(":",@F[0..3])."\t$F[4]"' {params.af} > {output.tmp} \
        fi \
        awk '{params.awk}' <(bcftools query -f '{params.ql0}' {input.sites[0]}) {output.tmp} >{output.af} \
        bcftools view -s {params.samples} {params.truth} | bcftools query -f '{params.ql2}' | sed -E 's/\/|\|/\\t/g' > {output.tmp2} \
        awk '{params.awk2}' <(bcftools query -f '{params.ql0}' {input.sites[0]}) {output.tmp2} >{output.gt} \
        ) &> {log}
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
        "../scripts/accuracy_quilt.R"


rule plot_accuracy_panelsize:
    input:
        truth=rules.collect_truth_gts.output.gt,
        af=rules.collect_truth_gts.output.af,
        glimpse=expand(
            rules.collect_glimpse_imputed_gts.output,
            size=config["refsize"],
            allow_missing=True,
        ),
        regular=expand(
            rules.collect_quilt_imputed_gts.output.regular,
            size=config["refsize"],
            allow_missing=True,
        ),
        mspbwt=expand(
            rules.collect_quilt_imputed_gts.output.mspbwt,
            size=config["refsize"],
            allow_missing=True,
        ),
        zilong=expand(
            rules.collect_quilt_imputed_gts.output.zilong,
            size=config["refsize"],
            allow_missing=True,
        ),
    output:
        pdf=os.path.join(OUTDIR_SUMMARY, "all.accuracy.down{depth}x.{chrom}.pdf"),
        rds=os.path.join(OUTDIR_SUMMARY, "all.accuracy.down{depth}x.{chrom}.rds"),
    log:
        os.path.join(OUTDIR_SUMMARY, "all.accuracy.down{depth}x.{chrom}.llog"),
    params:
        N="plot_accuracy_panelsize",
    conda:
        "../envs/quilt.yaml"
    script:
        "../scripts/accuracy_panelsize.R"


rule plot_accuracy_depth:
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
        pdf=os.path.join(OUTDIR_SUMMARY, "all.accuracy.panelsize{size}.{chrom}.pdf"),
        rds=os.path.join(OUTDIR_SUMMARY, "all.accuracy.panelsize{size}.{chrom}.rds"),
    log:
        os.path.join(OUTDIR_SUMMARY, "all.accuracy.panelsize{size}.{chrom}.pdf.llog"),
    params:
        N="plot_accuracy_depth",
    conda:
        "../envs/quilt.yaml"
    script:
        "../scripts/accuracy_depth.R"
