

rule glimpse_prepare_glvcf:
    input:
        bams=rules.bamlist.output,
        sites=rules.subset_refpanel.output.sites,
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
        bcftools mpileup -Ou -q {params.bq} -Q {params.mq} -f {params.fasta} -I -E -a 'FORMAT/DP' -T {input.sites} -b {input.bams} | \
        bcftools call -Aim -C alleles -T {input.sites} -Oz -o {output.vcf} && bcftools index -f {output.vcf} &> {log}
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
    params:
        get_glimpse_chunks,
        N="glimpse_phase",
    shell:
        """
        (
            /usr/bin/time -v GLIMPSE_phase \
            --input {input.glvcf} \
            --reference {input.refvcf} \
            --input-region {params[{wildcards.chunkid}]['irg']} \
            --output-region {params[{wildcards.chunkid}]['org']} \
            --output {output.vcf}
            bcftools index -f {output.csi}
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
