

rule glimpse_prepare_glvcf:
    input:
        bams=rules.bamlist.output,
        sites=lambda wildcards: expand(
            rules.subset_refpanel.output.sites,
            size=config["refsize"],
            allow_missing=True,
        ),
        tsv=lambda wildcards: expand(
            rules.subset_refpanel.output.tsv,
            size=config["refsize"],
            allow_missing=True,
        ),
    output:
        vcf=os.path.join(OUTDIR_GLIMPSE, "glvcf", "{chrom}", "down{depth}x.{chrom}.bcf"),
        csi=os.path.join(
            OUTDIR_GLIMPSE, "glvcf", "{chrom}", "down{depth}x.{chrom}.bcf.csi"
        ),
    log:
        os.path.join(OUTDIR_GLIMPSE, "glvcf/{chrom}/down{depth}.{chrom}.bcf.llog"),
    params:
        N="glimpse_prepare_glvcf",
        fasta=config["genome"]["fasta"],
        bq=config["glimpse"]["bq"],
        mq=config["glimpse"]["mq"],
    shell:
        """
        (
        bcftools mpileup -q {params.bq} -Q {params.mq} -f {params.fasta} -I -E -A -a 'FORMAT/DP' -r {wildcards.chrom} -T {input.sites[0]} -b {input.bams} -Ou \
        | bcftools call -Aim -C alleles -T {input.tsv[0]} -Ob -o {output.vcf} && bcftools index -f {output.vcf} \
        ) &> {log}
        """


rule glimpse_phase:
    input:
        refvcf=rules.subset_refpanel.output.vcf,
        glvcf=rules.glimpse_prepare_glvcf.output.vcf,
    output:
        vcf=os.path.join(
            OUTDIR_GLIMPSE,
            "panelsize{size}",
            "{chrom}",
            "down{depth}x.{chrom}.chunks{chunkid}.bcf",
        ),
        csi=os.path.join(
            OUTDIR_GLIMPSE,
            "panelsize{size}",
            "{chrom}",
            "down{depth}x.{chrom}.chunks{chunkid}.bcf.csi",
        ),
    log:
        os.path.join(
            OUTDIR_GLIMPSE,
            "panelsize{size}",
            "{chrom}",
            "down{depth}x.{chrom}.chunks{chunkid}.bcf.llog",
        ),
    params:
        N="glimpse_phase",
        irg=get_glimpse_chunki_irg,
        org=get_glimpse_chunki_org,
    shell:
        """
        (
            /usr/bin/time -v GLIMPSE_phase \
            --input {input.glvcf} \
            --reference {input.refvcf} \
            --input-region {params.irg} \
            --output-region {params.org} \
            --output {output.vcf} && \
            bcftools index -f {output.vcf} \
        ) &> {log}
        """


rule glimpse_ligate:
    input:
        get_glimpse_phase_outputs,
    output:
        vcf=os.path.join(
            OUTDIR_GLIMPSE, "panelsize{size}", "{chrom}", "down{depth}x.{chrom}.bcf"
        ),
        lst=temp(
            os.path.join(
                OUTDIR_GLIMPSE,
                "panelsize{size}",
                "{chrom}",
                "down{depth}x.{chrom}.vcf.list",
            ),
        ),
    params:
        N="glimpse_ligate",
    shell:
        """
        echo {input} > {output.lst}
        GLIMPSE_ligate --input {output.lst} --output {output.vcf} && bcftools index -f {output.vcf}
        """
