
rule subset_sample_list:
    """put sample to keep in a file in case of arguments too long"""
    output:
        samples=temp(os.path.join(OUTDIR_PANEL, "{chrom}.size{size}.kept.samples")),
    params:
        N="subset_sample_list",
        samples=SAMPLES.keys(),
        exclude=if_exclude_samples_in_refpanel,
        vcf=lambda wildcards: REFPANEL[wildcards.chrom]["vcf"],
    log:
        os.path.join(OUTDIR_PANEL, "{chrom}.size{size}.kept.samples.llog"),
    conda:
        "../envs/quilt.yaml"
    script:
        "../scripts/subset_samples.R"


rule subset_refpanel_by_chrom:
    input:
        rules.subset_sample_list.output.samples,
    output:
        vcf=os.path.join(OUTDIR_PANEL, "{chrom}.size{size}.vcf.gz"),
        csi=os.path.join(OUTDIR_PANEL, "{chrom}.size{size}.vcf.gz.csi"),
        # hap=os.path.join(OUTDIR_PANEL, "{chrom}.size{size}.vcf.hap.gz"),
        # leg=os.path.join(OUTDIR_PANEL, "{chrom}.size{size}.vcf.legend.gz"),
        sites=os.path.join(OUTDIR_PANEL, "{chrom}.size{size}.sites.vcf.gz"),
        tsv=os.path.join(OUTDIR_PANEL, "{chrom}.size{size}.sites.tsv.gz"),
    params:
        N="subset_refpanel_by_chrom",
        prefix=lambda wildcards, output: os.path.splitext(output[0])[0],
        vcf=lambda wildcards: REFPANEL[wildcards.chrom]["vcf"],
    log:
        os.path.join(OUTDIR_PANEL, "{chrom}.subrefs{size}.llog"),
    conda:
        "../envs/pandas.yaml"
    shell:
        """
        bcftools view -v snps -m2 -M2 -g het --samples-file {input} --threads 4 {params.vcf}| bcftools norm - -d snps -Ob -o {output.vcf} --threads 4 && bcftools index -f {output.vcf} && \
        touch -m {output.vcf}.csi && \
        bcftools view -G {output.vcf} -Oz -o {output.sites} --threads 4 && tabix -f {output.sites} && \
        bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' {output.sites} | bgzip -c > {output.tsv} && \
        tabix -s1 -b2 -e2 {output.tsv} &> {log}
        """


rule subset_refpanel_by_chunkid:
    input:
        rules.subset_sample_list.output.samples,
    output:
        vcf=os.path.join(
            OUTDIR_PANEL, "refsize{size}", "{chrom}", "chunk_{chunkid}.vcf.gz"
        ),
        csi=os.path.join(
            OUTDIR_PANEL, "refsize{size}", "{chrom}", "chunk_{chunkid}.vcf.gz.csi"
        ),
        sites=os.path.join(
            OUTDIR_PANEL,
            "refsize{size}",
            "{chrom}",
            "chunk_{chunkid}.sites.vcf.gz",
        ),
        tsv=os.path.join(
            OUTDIR_PANEL, "refsize{size}", "{chrom}", "chunk_{chunkid}.tsv.vcf.gz"
        ),
    params:
        prefix=lambda wildcards, output: os.path.splitext(output[0])[0],
        vcf=lambda wildcards: REFPANEL[wildcards.chrom]["vcf"],
        rg=lambda wildcards: get_refpanel_chunk_region(
            wildcards.chrom, wildcards.chunkid
        ),
    log:
        os.path.join(
            OUTDIR_PANEL, "refsize{size}", "{chrom}", "chunk_{chunkid}.vcf.gz.llog"
        ),
    conda:
        "../envs/pandas.yaml"
    shell:
        """
        bcftools view -v snps -m2 -M2 --samples-file {input} --threads 4 {params.vcf} {params.rg} | bcftools norm - -d snps -Ob -o {output.vcf} --threads 4 && bcftools index -f {output.vcf} && \
        touch -m {output.vcf}.csi && \
        bcftools view -G {output.vcf} -Oz -o {output.sites} --threads 4 && tabix -f {output.sites} && \
        bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' {output.sites} | bgzip -c > {output.tsv} && \
        tabix -s1 -b2 -e2 {output.tsv} &> {log}
        """


def get_sites_refpanel_by_chunks(wildcards):
    d = get_refpanel_chunks(wildcards.chrom)
    ids = list(map(str, d.keys()))
    return expand(
        rules.subset_refpanel_by_chunkid.output.sites,
        chunkid=ids,
        allow_missing=True,
    )


rule concat_refpanel_sites_by_chunks:
    input:
        get_sites_refpanel_by_chunks,
    output:
        sites=os.path.join(
            OUTDIR_PANEL,
            "refsize{size}",
            "{chrom}.sites.vcf.gz",
        ),
        tsv=os.path.join(OUTDIR_PANEL, "refsize{size}", "{chrom}.tsv.vcf.gz"),
    log:
        os.path.join(OUTDIR_PANEL, "refsize{size}", "{chrom}.sites.vcf.gz.llog"),
    conda:
        "../envs/pandas.yaml"
    shell:
        """
        echo {input} | tr ' ' '\n' > {output.sites}.list && \
        bcftools concat -f {output.sites}.list -Da --threads 4 -Oz -o {output.sites} && \
        bcftools index -f {output.sites} && \
        bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' {output.sites} | bgzip -c > {output.tsv} && \
        tabix -s1 -b2 -e2 {output.tsv} &> {log}
        """
