OUTDIR_PANEL = "results/subrefs"


rule subset_refpanel:
    output:
        vcf=os.path.join(OUTDIR_PANEL, "{chrom}.size{size}.bcf"),
        hap=os.path.join(OUTDIR_PANEL, "{chrom}.size{size}.hap.gz"),
        leg=os.path.join(OUTDIR_PANEL, "{chrom}.size{size}.legend.gz"),
        sites=os.path.join(OUTDIR_PANEL, "{chrom}.size{size}.sites.vcf.gz"),
    params:
        outdir=OUTDIR_PANEL,
        vcf=lambda wildcards: REFPANEL[wildcards.chrom]["vcf"],
        samples=get_samples_list_comma,
    log:
        os.path.join(OUTDIR_PANEL, "{chrom}.subrefs{size}.llog"),
    conda:
        "../envs/pandas.yaml"
    threads: 1
    shell:
        """
        ./workflow/scripts/prep-refs.sh {wildcards.chrom} {params.vcf} {params.outdir} {params.samples} &> {log}
        """
