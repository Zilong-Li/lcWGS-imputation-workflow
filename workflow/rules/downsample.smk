
## samtools depth -a $bam | awk '{s+=\$3;} END{print s/NR}'
rule downsample_bam:
    output:
        os.path.join(OUTDIR_DOWNSAMPLE, "{sample}_{depth}x_{chrom}.bam"),
    params:
        N="downsample_bam",
        depth=lambda wildcards: SAMPLES[wildcards.sample]["depth"],
        bam=lambda wildcards: SAMPLES[wildcards.sample]["bam"],
    log:
        os.path.join(OUTDIR_DOWNSAMPLE, "{sample}_{depth}x_{chrom}.bam.llog"),
    threads: 1
    conda:
        "../envs/pandas.yaml"
    shell:
        """
        (
        if [ {wildcards.depth} == 0 ];then \
            samtools view -o {output} {params.bam} {wildcards.chrom} && samtools index {output} \
        ;else\
            FRAC=$(echo "scale=4 ; {wildcards.depth} / {params.depth}" | bc -l) && \
            samtools view -s $FRAC -o {output} {params.bam} {wildcards.chrom} && samtools index {output} \
        ; fi
        ) &> {log}
        """


rule bamlist:
    input:
        expand(rules.downsample_bam.output, sample=SAMPLES.keys(), allow_missing=True),
    output:
        os.path.join(OUTDIR_DOWNSAMPLE, "{depth}x_{chrom}.bamlist"),
    log:
        os.path.join(OUTDIR_DOWNSAMPLE, "{depth}x_{chrom}.bamlist.llog"),
    params:
        N="bamlist",
    conda:
        "../envs/pandas.yaml"
    shell:
        """ echo {input} | tr ' ' '\\n' > {output} """
