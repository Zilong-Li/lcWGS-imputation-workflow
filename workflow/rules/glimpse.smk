OUTDIR_GLIMPSE = "results/glimpse/"


rule glimpse_prepare_glvcf:
    input:
        bams=rules.bamlist.output,
        sites=rules.subset_refpanel.output.sites,
    output:
        vcf=OUTDIR_GLIMPSE + "glvcf/{chrom}/down{depth}x.{chrom}.bcf",
        csi=OUTDIR_GLIMPSE + "glvcf/{chrom}/down{depth}x.{chrom}.bcf.csi",
    log:
        OUTDIR_GLIMPSE + "glvcf/{chrom}/down{depth}.{chrom}.bcf.llog",
    params:
        fasta=config["genome"]["fasta"],
        bq=config["glimpse"]["bq"],
        mq=config["glimpse"]["mq"],
    shell:
        """
        bcftools mpileup -Ou -q {params.bq} -Q {params.mq} -f {params.fasta} -I -E -a 'FORMAT/DP' -T {input.sites} -b {input.bams} | \
        bcftools call -Aim -C alleles -T {input.sites} -Oz -o {output.vcf} && bcftools index -f {output.vcf} &> {log}
        """


# rule glimpse_chunk:
#     output:
#         txt=OUTDIR_GLIMPSE + "{chrom}.chunks.txt",
#     log:
#         OUTDIR_GLIMPSE + "{chrom}.chunks.llog",
#     params:
#         vcf=lambda wildcards: REFPANEL[wildcards.chrom]['vcf'],
#         chunksize=config["glimpse"]["chunksize"],
#         buffer=config["glimpse"]["buffer"],
#     shell:
#         """
#         /usr/bin/time -v GLIMPSE_chunk \
#           --input {params.vcf} \
#           --region {wildcards.chrom} \
#           --window-size {params.chunksize} \
#           --buffer-size {params.buffer} \
#           --output {output} &> {log}
#         """


rule glimpse_phase:
    input:
        refvcf=rules.subset_refpanel.output.vcf,
        glvcf=rules.glimpse_prepare_glvcf.output.vcf,
    output:
        vcf=OUTDIR_GLIMPSE
        + "panelsize{size}/{chrom}/down{depth}x.{chrom}.chunks{chunkid}.bcf",
        csi=OUTDIR_GLIMPSE
        + "panelsize{size}/{chrom}/down{depth}x.{chrom}.chunks{chunkid}.bcf.csi",
    params:
        get_glimpse_chunks,
    shell:
        """
        (
            /usr/bin/time -v GLIMPSE_phase \
            --input {input.glvcf} \
            --reference {input.refvcf} \
            --input-region {params[0]['irg']} \
            --output-region {params[0]['org']} \
            --output {output.vcf}
            bcftools index -f {output.csi}
        ) &> {log}
        """


rule glimpse_ligate:
    input:
        get_glimpse_outputs,
    output:
        vcf=OUTDIR_GLIMPSE + "panelsize{size}/{chrom}/down{depth}x.{chrom}.bcf",
        lst=temp(
            OUTDIR_GLIMPSE + "panelsize{size}/{chrom}/down{depth}x.{chrom}.vcf.list"
        ),
    shell:
        """
        echo {input} > {output.lst}
        GLIMPSE_ligate --input {output.lst} --output {output.vcf} && bcftools index -f {output.vcf}
        """
