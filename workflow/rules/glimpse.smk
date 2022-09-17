def get_quilt_prepare_output(wildcards):
    regions = get_regions_list_per_chrom(wildcards.chrom, config["quilt"]["chunksize"])
    return [
        f"results/quilt/RData/QUILT_prepared_reference.{wildcards.chrom}.{start}.{end}.RData"
        for start, end in regions
    ]


rule glimpse_prepare:
    input:
        vcf=os.path.join("results/subrefs/{chrom}.bcf"),
    output:
        os.path.join("results/glimpse/{chrom}.chunks.txt"),
    shell:
        """
        GLIMPSE_chunk \
          --input ${input} \
          --reference {input} \
          --region {wildcards.chrom} \
          --window-size {config[glimpse][chunksize]} \
          --buffer-size {config[glimpse][buffer]} \
          --output {output}
        """


rule glimpse:
    input:
        get_quilt_prepare_output,
    output:
        os.path.join("results/quilt/RData/QUILT_prepared_reference.{chrom}.lst"),
    shell:
        """ echo {input} > {output} """
