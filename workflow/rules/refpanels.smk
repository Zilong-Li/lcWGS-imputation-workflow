OUTDIR_PANEL = "results/subrefs"


def get_samples_list_comma(exclude=True):
    if exclude:
        return "^" + ",".join(SAMPLES.keys())
    else:
        return ",".join(SAMPLES.keys())


def get_ref_vcf(wildcards):
    return REFPANEL[wildcards.chrom]["vcf"]


rule subset_refpanel:
    output:
        vcf=os.path.join(OUTDIR_PANEL, "{chrom}.bcf"),
        hap=os.path.join(OUTDIR_PANEL, "{chrom}.hap.gz"),
        leg=os.path.join(OUTDIR_PANEL, "{chrom}.legend.gz"),
        sites=os.path.join(OUTDIR_PANEL, "{chrom}.sites.vcf.gz"),
    params:
        samples=get_samples_list_comma(),
        vcf=lambda wildcards: REFPANEL[wildcards.chrom]["vcf"],
    threads: 1
    shell:
        """
        ./workflow/scripts/prep-refs.sh {wildcards.chrom} {params.vcf} {OUTDIR_PANEL} {params.samples}
        """
