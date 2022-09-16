def get_quilt_output(wildcards):
    regions = get_regions_list_per_chrom(wildcards.chrom, config["quilt"]["chunksize"])
    res = []
    for start, end in regions:
        res.append(
            f"results/quilt/{wildcards.chrom}/quilt.{wildcards.depth}x.regular.{wildcards.chrom}.{start}.{end}.vcf.gz"
        )
        res.append(
            f"results/quilt/{wildcards.chrom}/quilt.{wildcards.depth}x.mspbwt.{wildcards.chrom}.{start}.{end}.vcf.gz"
        )
        res.append(
            f"results/quilt/{wildcards.chrom}/quilt.{wildcards.depth}x.zilong.{wildcards.chrom}.{start}.{end}.vcf.gz"
        )
    return res


rule quilt_prepare:
    input:
        vcf=rules.subset_refpanel.output.vcf,
        hap=rules.subset_refpanel.output.hap,
        leg=rules.subset_refpanel.output.leg,
    output:
        os.path.join(
            "results/quilt/{chrom}/RData/QUILT_prepared_reference.{chrom}.{start}.{end}.RData"
        ),
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
            --outputdir=results/quilt/ &> {output}.llog
        """


rule quilt_run_regular:
    input:
        vcf=rules.subset_refpanel.output.vcf,
        hap=rules.subset_refpanel.output.hap,
        leg=rules.subset_refpanel.output.leg,
        bams=rules.bamlist.output,
        rdata=rules.quilt_prepare.output,
    output:
        os.path.join(
            "results/quilt/{chrom}/quilt.{depth}x.regular.{chrom}.{start}.{end}.vcf.gz"
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
            --output_filename={output} &> {output}.llog
        """


rule quilt_run_mspbwt:
    input:
        vcf=rules.subset_refpanel.output.vcf,
        hap=rules.subset_refpanel.output.hap,
        leg=rules.subset_refpanel.output.leg,
        bams=rules.bamlist.output,
        rdata=rules.quilt_prepare.output,
    output:
        os.path.join(
            "results/quilt/{chrom}/quilt.{depth}x.mspbwt.{chrom}.{start}.{end}.vcf.gz"
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
            --use_mspbwt=TRUE \
            --output_filename={output} &> {output}.llog
        """


rule quilt_run_zilong:
    input:
        vcf=rules.subset_refpanel.output.vcf,
        hap=rules.subset_refpanel.output.hap,
        leg=rules.subset_refpanel.output.leg,
        bams=rules.bamlist.output,
        rdata=rules.quilt_prepare.output,
    output:
        os.path.join(
            "results/quilt/{chrom}/quilt.{depth}x.zilong.{chrom}.{start}.{end}.vcf.gz"
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
            --pbwtL={config[quilt][pbwtL]} \
            --pbwtS={config[quilt][pbwtS]} \
            --zilong=TRUE \
            --use_mspbwt=FALSE \
            --output_filename={output} &> {output}.llog
        """


rule quilt:
    input:
        get_quilt_output,
    output:
        os.path.join("results/quilt/{chrom}/results.{depth}x_{chrom}.lst"),
    shell:
        """ echo {input} | tr ' ' '\n' > {output} """
