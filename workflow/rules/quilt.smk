def get_quilt_output_regular(wildcards):
    regions = get_regions_list_per_chrom(wildcards.chrom, config["quilt"]["chunksize"])
    return [
        f"results/quilt/{wildcards.chrom}/quilt.{wildcards.depth}x.regular.{wildcards.chrom}.{start}.{end}.vcf.gz"
        for start, end in regions
    ]


def get_quilt_output_mspbwt(wildcards):
    regions = get_regions_list_per_chrom(wildcards.chrom, config["quilt"]["chunksize"])
    return [
        f"results/quilt/{wildcards.chrom}/quilt.{wildcards.depth}x.mspbwt.{wildcards.chrom}.{start}.{end}.vcf.gz"
        for start, end in regions
    ]


def get_quilt_output_zilong(wildcards):
    regions = get_regions_list_per_chrom(wildcards.chrom, config["quilt"]["chunksize"])
    return [
        f"results/quilt/{wildcards.chrom}/quilt.{wildcards.depth}x.zilong.{wildcards.chrom}.{start}.{end}.vcf.gz"
        for start, end in regions
    ]


rule quilt_prepare:
    input:
        vcf=rules.subset_refpanel.output.vcf,
        hap=rules.subset_refpanel.output.hap,
        leg=rules.subset_refpanel.output.leg,
    output:
        "results/quilt/{chrom}/RData/QUILT_prepared_reference.{chrom}.{start}.{end}.RData",
    log:
        "results/quilt/{chrom}/RData/QUILT_prepared_reference.{chrom}.{start}.{end}.RData.llog",
    threads: 1
    shell:
        """
        /usr/bin/time -v /gpfs3/users/davies/xxd908/local/pkgs/QUILT/QUILT_prepare_reference.R \
            --reference_vcf_file=results/subrefs/{wildcards.chrom}.bcf \
            --reference_haplotype_file=results/subrefs/{wildcards.chrom}.hap.gz \
            --reference_legend_file=results/subrefs/{wildcards.chrom}.legend.gz \
            --chr={wildcards.chrom} \
            --regionStart={wildcards.start} \
            --regionEnd={wildcards.end} \
            --buffer={config[quilt][buffer]} \
            --nGen={config[quilt][nGen]} \
            --use_pbwt_index=TRUE \
            --use_mspbwt=TRUE \
            --outputdir=results/quilt/{wildcards.chrom} &> {log}
        """


rule quilt_run_regular:
    input:
        vcf=rules.subset_refpanel.output.vcf,
        hap=rules.subset_refpanel.output.hap,
        leg=rules.subset_refpanel.output.leg,
        bams=rules.bamlist.output,
        rdata=rules.quilt_prepare.output,
    output:
        temp(
            "results/quilt/{chrom}/quilt.{depth}x.regular.{chrom}.{start}.{end}.vcf.gz"
        ),
    log:
        os.path.join(
            "results/quilt/{chrom}/quilt.{depth}x.regular.{chrom}.{start}.{end}.vcf.gz.llog"
        ),
    threads: 1
    shell:
        """
        /usr/bin/time -v /gpfs3/users/davies/xxd908/local/pkgs/QUILT/QUILT.R \
            --reference_vcf_file=results/subrefs/{wildcards.chrom}.bcf \
            --reference_haplotype_file=results/subrefs/{wildcards.chrom}.hap.gz \
            --reference_legend_file=results/subrefs/{wildcards.chrom}.legend.gz \
            --prepared_reference_filename={input.rdata} \
            --bamlist={input.bams} \
            --chr={wildcards.chrom} \
            --regionStart={wildcards.start} \
            --regionEnd={wildcards.end} \
            --buffer={config[quilt][buffer]} \
            --nGen={config[quilt][nGen]} \
            --zilong=FALSE \
            --use_mspbwt=FALSE \
            --output_filename={output} &> {log}
        """


rule quilt_ligate_regular:
    input:
        get_quilt_output_regular,
    output:
        vcf="results/quilt/{chrom}/quilt.{depth}x.regular.{chrom}.bcf.gz",
        lst=temp("results/quilt/{chrom}/quilt.{depth}x.regular.{chrom}.vcf.list"),
    log:
        "results/quilt/{chrom}/quilt.{depth}x.regular.{chrom}.bcf.gz.llog",
    shell:
        """
        ( \
           echo {input} | tr ' ' '\n' > {output.lst} && \
           bcftools concat --file-list {output.lst} --output-type b --threads 4 -o {output.vcf} && \
           bcftools index -f {output.vcf} \
        ) &> {log}
        """


rule quilt_run_mspbwt:
    input:
        vcf=rules.subset_refpanel.output.vcf,
        hap=rules.subset_refpanel.output.hap,
        leg=rules.subset_refpanel.output.leg,
        bams=rules.bamlist.output,
        rdata=rules.quilt_prepare.output,
    output:
        temp("results/quilt/{chrom}/quilt.{depth}x.mspbwt.{chrom}.{start}.{end}.vcf.gz"),
    log:
        "results/quilt/{chrom}/quilt.{depth}x.mspbwt.{chrom}.{start}.{end}.vcf.gz.llog",
    threads: 1
    shell:
        """
        /usr/bin/time -v /gpfs3/users/davies/xxd908/local/pkgs/QUILT/QUILT.R \
            --reference_vcf_file=results/subrefs/{wildcards.chrom}.bcf \
            --reference_haplotype_file=results/subrefs/{wildcards.chrom}.hap.gz \
            --reference_legend_file=results/subrefs/{wildcards.chrom}.legend.gz \
            --prepared_reference_filename={input.rdata} \
            --bamlist={input.bams} \
            --chr={wildcards.chrom} \
            --regionStart={wildcards.start} \
            --regionEnd={wildcards.end} \
            --buffer={config[quilt][buffer]} \
            --nGen={config[quilt][nGen]} \
            --zilong=FALSE \
            --use_mspbwt=TRUE \
            --output_filename={output} &> {log}
        """


rule quilt_ligate_mspbwt:
    input:
        get_quilt_output_mspbwt,
    output:
        vcf="results/quilt/{chrom}/quilt.{depth}x.mspbwt.{chrom}.bcf.gz",
        lst=temp("results/quilt/{chrom}/quilt.{depth}x.mspbwt.{chrom}.vcf.list"),
    log:
        "results/quilt/{chrom}/quilt.{depth}x.mspbwt.{chrom}.bcf.gz.llog",
    shell:
        """
        ( \
           echo {input} | tr ' ' '\n' > {output.lst} && \
           bcftools concat --file-list {output.lst} --output-type b --threads 4 -o {output.vcf} && \
           bcftools index -f {output.vcf} \
        ) &> {log}
        """


rule quilt_run_zilong:
    input:
        vcf=rules.subset_refpanel.output.vcf,
        hap=rules.subset_refpanel.output.hap,
        leg=rules.subset_refpanel.output.leg,
        bams=rules.bamlist.output,
        rdata=rules.quilt_prepare.output,
    output:
        temp("results/quilt/{chrom}/quilt.{depth}x.zilong.{chrom}.{start}.{end}.vcf.gz"),
    log:
        "results/quilt/{chrom}/quilt.{depth}x.zilong.{chrom}.{start}.{end}.vcf.gz.llog",
    threads: 1
    shell:
        """
        /usr/bin/time -v /gpfs3/users/davies/xxd908/local/pkgs/QUILT/QUILT.R \
            --reference_vcf_file=results/subrefs/{wildcards.chrom}.bcf \
            --reference_haplotype_file=results/subrefs/{wildcards.chrom}.hap.gz \
            --reference_legend_file=results/subrefs/{wildcards.chrom}.legend.gz \
            --prepared_reference_filename={input.rdata} \
            --bamlist={input.bams} \
            --chr={wildcards.chrom} \
            --regionStart={wildcards.start} \
            --regionEnd={wildcards.end} \
            --buffer={config[quilt][buffer]} \
            --nGen={config[quilt][nGen]} \
            --pbwtL={config[quilt][pbwtL]} \
            --pbwtS={config[quilt][pbwtS]} \
            --zilong=TRUE \
            --use_mspbwt=FALSE \
            --output_filename={output} &> {log}
        """


rule quilt_ligate_zilong:
    input:
        get_quilt_output_zilong,
    output:
        vcf="results/quilt/{chrom}/quilt.{depth}x.zilong.{chrom}.bcf.gz",
        lst=temp("results/quilt/{chrom}/quilt.{depth}x.zilong.{chrom}.vcf.list"),
    log:
        "results/quilt/{chrom}/quilt.{depth}x.zilong.{chrom}.bcf.gz.llog",
    shell:
        """
        ( \
           echo {input} | tr ' ' '\n' > {output.lst} && \
           bcftools concat --file-list {output.lst} --output-type b --threads 4 -o {output.vcf} && \
           bcftools index -f {output.vcf} \
        ) &> {log}
        """
