OUTDIR_PANEL = "results/subrefs"


def get_samples_list_comma(exclude=True):
    if exclude:
        return "^" + ",".join(SAMPLES.keys())
    else:
        return ",".join(SAMPLES.keys())


rule subset_refpanel:
    output:
        vcf=os.path.join(OUTDIR_PANEL, "{chrom}.bcf"),
        hap=os.path.join(OUTDIR_PANEL, "{chrom}.hap.gz"),
        leg=os.path.join(OUTDIR_PANEL, "{chrom}.legend.gz"),
        sites=os.path.join(OUTDIR_PANEL, "{chrom}.sites.vcf.gz"),
    params:
        outdir=OUTDIR_PANEL,
        vcf=lambda wildcards: REFPANEL[wildcards.chrom]["vcf"],
        samples=get_samples_list_comma(),
    log:
        os.path.join(OUTDIR_PANEL, "{chrom}.subrefs.llog"),
    conda:
        "../envs/pandas.yaml"
    threads: 1
    shell:
        """
        ./workflow/scripts/prep-refs.sh {wildcards.chrom} {params.vcf} {params.outdir} {params.samples} &> {log}
        """
