
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
        ( \
            bcftools view -v snps -m2 -M2 --samples-file {input} --threads 4 {params.vcf}| bcftools norm - -d snps -Ob -o {output.vcf} --threads 4 && bcftools index -f {output.vcf} && \
            touch -m {output.vcf}.csi && \
            bcftools view -G {output.vcf} -Oz -o {output.sites} --threads 4 && tabix -f {output.sites} && \
            bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' {output.sites} | bgzip -c > {output.tsv} && tabix -s1 -b2 -e2 {output.tsv}
        )  &> {log}
        """


rule subset_refpanel_by_region2:
    input:
        rules.subset_sample_list.output.samples,
    output:
        vcf=os.path.join(
            OUTDIR_PANEL, "panelsize{size}", "vcfs" "{chrom}.{start}.{end}.vcf.gz"
        ),
        csi=os.path.join(
            OUTDIR_PANEL, "panelsize{size}", "vcfs" "{chrom}.{start}.{end}.vcf.gz.csi"
        ),
        sites=os.path.join(
            OUTDIR_PANEL,
            "panelsize{size}",
            "vcfs" "{chrom}.{start}.{end}.sites.vcf.gz",
        ),
        tsv=os.path.join(
            OUTDIR_PANEL, "panelsize{size}", "vcfs" "{chrom}.{start}.{end}.tsv.vcf.gz"
        ),
    params:
        N="subset_refpanel_by_region2",
        prefix=lambda wildcards, output: os.path.splitext(output[0])[0],
        vcf=lambda wildcards: REFPANEL[wildcards.chrom]["vcf"],
        start=lambda wildcards: max(
            1, int(wildcards.start) - int(config["glimpse"]["buffer"])
        ),
        end=lambda wildcards: int(wildcards.end) + int(config["glimpse"]["buffer"]),
    log:
        os.path.join(
            OUTDIR_PANEL, "panelsize{size}", "vcfs" "{chrom}.{start}.{end}.vcf.gz.llog"
        ),
    conda:
        "../envs/pandas.yaml"
    shell:
        """
        ( \
            bcftools view -v snps -m2 -M2 --samples-file {input} --threads 4 {params.vcf} {wildcards.chrom}:{params.start}-{params.end}| bcftools norm - -d snps -Ob -o {output.vcf} --threads 4 && bcftools index -f {output.vcf} && \
            touch -m {output.vcf}.csi && \
            bcftools view -G {output.vcf} -Oz -o {output.sites} --threads 4 && tabix -f {output.sites} && \
            bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' {output.sites} | bgzip -c > {output.tsv} && tabix -s1 -b2 -e2 {output.tsv}
        )  &> {log}
        """


def get_subset_refpanel_by_region2(wildcards):
    starts, ends = get_regions_list_from_glimpse_chunk(wildcards.chrom)
    return expand(
        rules.subset_refpanel_by_region2.output,
        zip,
        start=starts,
        end=ends,
        allow_missing=True,
    )


rule concat_refpanel_sites_by_region2:
    input:
        get_subset_refpanel_by_region2,
    output:
        sites=os.path.join(
            OUTDIR_PANEL,
            "panelsize{size}",
            "vcfs" "{chrom}.sites.vcf.gz",
        ),
        tsv=os.path.join(OUTDIR_PANEL, "panelsize{size}", "vcfs" "{chrom}.tsv.vcf.gz"),
    log:
        os.path.join(
            OUTDIR_PANEL, "panelsize{size}", "vcfs" "{chrom}.sites.vcf.gz.llog"
        ),
    conda:
        "../envs/pandas.yaml"
    shell:
        """
        ( \
            bcftools concat -D --threads 4 -Oz -o {output.sites} {input} && bcftools index -f {output.sites} \
            bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' {output.sites} | bgzip -c > {output.tsv} && tabix -s1 -b2 -e2 {output.tsv}
        ) & > {log}
        """
