OUTDIR_DOWNSAMPLE = "results/downsample"


def get_bams_per_depth_chrom(wildcards):
    return expand(
        rules.downsample_bam.output, sample=SAMPLES.keys(), allow_missing=True
    )


rule downsample_bam:
    output:
        os.path.join(OUTDIR_DOWNSAMPLE, "{sample}_{depth}x_{chrom}.bam"),
    params:
        depth=lambda wildcards: SAMPLES[wildcards.sample]["depth"],
        bam=lambda wildcards: SAMPLES[wildcards.sample]["bam"],
    log:
        os.path.join(OUTDIR_DOWNSAMPLE, "{sample}_{depth}x_{chrom}.bam.llog"),
    threads: 1
    shell:
        """
        (
        FRAC=$(echo "scale=4 ; {wildcards.depth} / {params.depth}" | bc -l)
        samtools view -s $FRAC -o {output} {params.bam} {wildcards.chrom} && samtools index {output}
        ) &> {log}
        """


rule bamlist:
    input:
        get_bams_per_depth_chrom,
    output:
        os.path.join(OUTDIR_DOWNSAMPLE, "{depth}x_{chrom}.bamlist"),
    shell:
        """ echo {input} | tr ' ' '\\n' > {output} """
