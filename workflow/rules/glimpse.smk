OUTDIR_GLIMPSE = "results/glimpse/"

rule glimpse_prepare:
    input:
        vcf=rules.subset_refpanel.output.vcf,
    output:
        OUTDIR_GLIMPSE + "panelsize{size}/{chrom}/prep/.{chrom}.{start}.{end}.RData",
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
