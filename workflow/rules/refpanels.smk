
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


rule subset_refpanel:
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
        N="subset_refpanel",
        prefix=lambda wildcards, output: os.path.splitext(output[0])[0],
        vcf=lambda wildcards: REFPANEL[wildcards.chrom]["vcf"],
        start=lambda wildcards: max(1, wildcards.start - config["glimpse"]["buffer"]),
        end=lambda wildcards: wildcards.end + config["glimpse"]["buffer"],
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
