OUTDIR_PANEL = "results/subrefs"


rule subset_refpanel:
    output:
        vcf=os.path.join(OUTDIR_PANEL, "{chrom}.size{size}.bcf"),
        hap=os.path.join(OUTDIR_PANEL, "{chrom}.size{size}.hap.gz"),
        leg=os.path.join(OUTDIR_PANEL, "{chrom}.size{size}.legend.gz"),
        sites=os.path.join(OUTDIR_PANEL, "{chrom}.size{size}.sites.vcf.gz"),
    params:
        prefix=lambda wildcards, output: os.path.splitext(output[0])[0],
        vcf=lambda wildcards: REFPANEL[wildcards.chrom]["vcf"],
        samples=get_samples_list_comma,
    log:
        os.path.join(OUTDIR_PANEL, "{chrom}.subrefs{size}.llog"),
    conda:
        "../envs/pandas.yaml"
    shell:
        """
        ( \
            bcftools view -v snps -m2 -M2 -s {params.samples} --threads 4 {params.vcf}| bcftools norm - -d snps -Ob -o {output.vcf} --threads 4 && bcftools index -f {output.vcf} && \
            bcftools convert --haplegendsample {params.prefix} {output.vcf} && \
            bcftools view -G {output.vcf} -Oz -o {output.sites} && tabix -f {output.sites} \
        )  2> {log}
        """
