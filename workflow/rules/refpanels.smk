
rule subset_sample_list:
    """put sample to be used in a file in case of arguments too long"""
    output:
        samples=temp(os.path.join(OUTDIR_PANEL, "{chrom}.size{size}.kept.samples")),
    params:
        N="subset_sample_list",
        samples=SAMPLES.keys(),
        vcf=lambda wildcards: REFPANEL[wildcards.chrom]["vcf"],
    log:
        os.path.join(OUTDIR_PANEL, "{chrom}.size{size}.kept.samples.llog"),
    conda:
        "../envs/quilt.yaml"
    script:
        "../scripts/subset_samples.R"


rule subset_refpanel:
    input:
        rules.subset_sample_list.output.samples,
    output:
        vcf=os.path.join(OUTDIR_PANEL, "{chrom}.size{size}.bcf"),
        hap=os.path.join(OUTDIR_PANEL, "{chrom}.size{size}.hap.gz"),
        leg=os.path.join(OUTDIR_PANEL, "{chrom}.size{size}.legend.gz"),
        sites=os.path.join(OUTDIR_PANEL, "{chrom}.size{size}.sites.vcf.gz"),
        tsv=os.path.join(OUTDIR_PANEL, "{chrom}.size{size}.sites.tsv.gz"),
    params:
        N="subset_refpanel",
        prefix=lambda wildcards, output: os.path.splitext(output[0])[0],
        vcf=lambda wildcards: REFPANEL[wildcards.chrom]["vcf"],
    log:
        os.path.join(OUTDIR_PANEL, "{chrom}.subrefs{size}.llog"),
    conda:
        "../envs/pandas.yaml"
    shell:
        """
        ( \
            bcftools view -v snps -m2 -M2 --samples-file {input} --threads 4 {params.vcf}| bcftools norm - -d snps -Ob -o {output.vcf} --threads 4 && bcftools index -f {output.vcf} && \
            bcftools convert --haplegendsample {params.prefix} {output.vcf} && \
            bcftools view -G {output.vcf} -Oz -o {output.sites} --threads 4 && tabix -f {output.sites} && \
            bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' {output.sites} | bgzip -c > {output.tsv} && tabix -s1 -b2 -e2 {output.tsv} \
        )  &> {log}
        """
