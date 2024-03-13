
rule quilt_prepare_regular:
    input:
        vcf=rules.subset_refpanel_by_chunkid.output.vcf,
    output:
        os.path.join(
            OUTDIR_QUILT1,
            "refsize{size}",
            "prep_regular",
            "QUILT_prepared_reference.{chrom}.chunk_{chunkid}.RData",
        ),
    params:
        time=config["time"],
        outdir=lambda wildcards, output: os.path.dirname(output[0]),
        nGen=config["quilt1"]["nGen"],
        buffer=config["quilt1"]["buffer"],
        lowram=config["quilt1"]["lowram"],
        impute_rare_common=config["quilt1"]["impute_rare_common"],
        rare_af_threshold=config["quilt1"]["rare_af_threshold"],
        start=get_quilt_chunk_region_start,
        end=get_quilt_chunk_region_end,
        gmap=if_use_quilt_map_in_refpanel,
    log:
        os.path.join(
            OUTDIR_QUILT1,
            "refsize{size}",
            "prep_regular",
            "QUILT_prepared_reference.{chrom}.chunk_{chunkid}.RData.llog",
        ),
    conda:
        "../envs/quilt.yaml"
    threads: 1
    shell:
        """
        (
        if [ -s {params.gmap} ];then \
        {params.time} -v QUILT_prepare_reference.R \
            --genetic_map_file='{params.gmap}' \
            --reference_vcf_file={input.vcf} \
            --chr={wildcards.chrom} \
            --regionStart={params.start} \
            --regionEnd={params.end} \
            --buffer={params.buffer} \
            --nGen={params.nGen} \
            --use_hapMatcherR={params.lowram} \
            --use_mspbwt=FALSE \
            --impute_rare_common={params.impute_rare_common} \
            --rare_af_threshold={params.rare_af_threshold} \
            --outputdir={params.outdir} \
            --output_file={output} \
        ; else \
        {params.time} -v QUILT_prepare_reference.R \
            --reference_vcf_file={input.vcf} \
            --chr={wildcards.chrom} \
            --regionStart={params.start} \
            --regionEnd={params.end} \
            --buffer={params.buffer} \
            --use_hapMatcherR={params.lowram} \
            --nGen={params.nGen} \
            --use_mspbwt=FALSE \
            --impute_rare_common={params.impute_rare_common} \
            --rare_af_threshold={params.rare_af_threshold} \
            --outputdir={params.outdir} \
            --output_file={output} \
        ; fi
        ) &> {log}
        """


rule quilt_run_regular:
    input:
        vcf=rules.subset_refpanel_by_chunkid.output.vcf,
        bams=rules.bamlist.output,
        rdata=rules.quilt_prepare_regular.output,
    output:
        os.path.join(
            OUTDIR_QUILT1,
            "refsize{size}",
            "output",
            "quilt.down{depth}x.regular.{chrom}.chunk_{chunkid}.vcf.gz",
        ),
    params:
        time=config["time"],
        nGen=config["quilt1"]["nGen"],
        buffer=config["quilt1"]["buffer"],
        start=get_quilt_chunk_region_start,
        end=get_quilt_chunk_region_end,
        Ksubset=config["quilt1"]["Ksubset"],
        Knew=config["quilt1"]["Knew"],
        nGibbsSamples=config["quilt1"]["nGibbsSamples"],
        lowram=config["quilt1"]["lowram"],
        impute_rare_common=config["quilt1"]["impute_rare_common"],
        rare_af_threshold=config["quilt1"]["rare_af_threshold"],
        n_seek_its=config["quilt1"]["n_seek_its"],
        block_gibbs=config["quilt1"]["small_ref_panel_block_gibbs_iterations"],
        gibbs_iters=config["quilt1"]["small_ref_panel_gibbs_iterations"],
    log:
        os.path.join(
            OUTDIR_QUILT1,
            "refsize{size}",
            "output",
            "quilt.down{depth}x.regular.{chrom}.chunk_{chunkid}.vcf.gz.llog",
        ),
    conda:
        "../envs/quilt.yaml"
    threads: 1
    shell:
        """
        {params.time} -v QUILT.R \
            --reference_vcf_file={input.vcf} \
            --prepared_reference_filename={input.rdata} \
            --bamlist={input.bams} \
            --chr={wildcards.chrom} \
            --regionStart={params.start} \
            --regionEnd={params.end} \
            --buffer={params.buffer} \
            --nGen={params.nGen} \
            --zilong=FALSE \
            --use_mspbwt=FALSE \
            --Ksubset={params.Ksubset} \
            --Knew={params.Knew} \
            --nGibbsSamples={params.nGibbsSamples} \
            --use_hapMatcherR={params.lowram} \
            --impute_rare_common={params.impute_rare_common} \
            --rare_af_threshold={params.rare_af_threshold} \
            --n_seek_its={params.n_seek_its} \
            --small_ref_panel_block_gibbs_iterations='{params.block_gibbs}' \
            --small_ref_panel_gibbs_iterations={params.gibbs_iters} \
            --output_filename={output} &> {log}
        """


rule quilt_ligate_regular:
    input:
        get_quilt_regular_outputs,
    output:
        vcf=os.path.join(
            OUTDIR_QUILT1,
            "refsize{size}",
            "{chrom}",
            "quilt.down{depth}x.regular.{chrom}.vcf.gz",
        ),
        lst=temp(
            os.path.join(
                OUTDIR_QUILT1,
                "refsize{size}",
                "{chrom}",
                "quilt.down{depth}x.regular.{chrom}.vcf.list",
            )
        ),
    log:
        os.path.join(
            OUTDIR_QUILT1,
            "refsize{size}",
            "{chrom}",
            "quilt.down{depth}x.regular.{chrom}.vcf.gz.llog",
        ),
    params:
        N="quilt_ligate_regular",
        extra=config["extra_buffer_in_quilt"],
    conda:
        "../envs/quilt.yaml"
    shell:
        """
        ( \
        if [ {params.extra} -gt 0 ];then \
           echo {input} | tr ' ' '\n' > {output.lst} && \
           bcftools concat --ligate --file-list {output.lst} --output-type z --threads 4 -o {output.vcf} && \
           bcftools index -f {output.vcf} \
        ; else \
           echo {input} | tr ' ' '\n' > {output.lst} && \
           bcftools concat --file-list {output.lst} --output-type z --threads 4 -o {output.vcf} && \
           bcftools index -f {output.vcf} \
        ; fi \
        ) &> {log}
        """


rule quilt_prepare_mspbwt:
    input:
        vcf=rules.subset_refpanel_by_chunkid.output.vcf,
    output:
        os.path.join(
            OUTDIR_QUILT2,
            "refsize{size}",
            "prep_mspbwt",
            "QUILT_prepared_reference.{chrom}.chunk_{chunkid}.RData",
        ),
    params:
        time=config["time"],
        N="quilt_prepare_mspbwt",
        outdir=lambda wildcards, output: os.path.dirname(output[0]),
        nGen=config["quilt2"]["nGen"],
        buffer=config["quilt2"]["buffer"],
        start=get_quilt_chunk_region_start,
        end=get_quilt_chunk_region_end,
        gmap=if_use_quilt_map_in_refpanel,
        lowram=config["quilt2"]["lowram"],
        impute_rare_common=config["quilt2"]["impute_rare_common"],
        rare_af_threshold=config["quilt2"]["rare_af_threshold"],
        nindices=config["quilt2"]["mspbwt-nindices"],
    log:
        os.path.join(
            OUTDIR_QUILT2,
            "refsize{size}",
            "prep_mspbwt",
            "QUILT_prepared_reference.{chrom}.chunk_{chunkid}.RData.llog",
        ),
    conda:
        "../envs/quilt.yaml"
    threads: 1
    shell:
        """
        (
        if [ -s {params.gmap} ];then \
        {params.time} -v QUILT_prepare_reference.R \
            --genetic_map_file='{params.gmap}' \
            --reference_vcf_file={input.vcf} \
            --chr={wildcards.chrom} \
            --regionStart={params.start} \
            --regionEnd={params.end} \
            --use_hapMatcherR={params.lowram} \
            --buffer={params.buffer} \
            --nGen={params.nGen} \
            --use_mspbwt=TRUE \
            --impute_rare_common={params.impute_rare_common} \
            --rare_af_threshold={params.rare_af_threshold} \
            --mspbwt_nindices={params.nindices} \
            --outputdir={params.outdir} \
            --output_file={output} \
        ; else \
        {params.time} -v QUILT_prepare_reference.R \
            --reference_vcf_file={input.vcf} \
            --chr={wildcards.chrom} \
            --regionStart={params.start} \
            --regionEnd={params.end} \
            --buffer={params.buffer} \
            --use_hapMatcherR={params.lowram} \
            --nGen={params.nGen} \
            --use_mspbwt=TRUE \
            --rare_af_threshold={params.rare_af_threshold} \
            --impute_rare_common={params.impute_rare_common} \
            --mspbwt_nindices={params.nindices} \
            --outputdir={params.outdir} \
            --output_file={output} \
        ; fi \
        ) &> {log}
        """


rule quilt_run_mspbwt:
    input:
        vcf=rules.subset_refpanel_by_chunkid.output.vcf,
        bams=rules.bamlist.output,
        rdata=rules.quilt_prepare_mspbwt.output,
    output:
        os.path.join(
            OUTDIR_QUILT2,
            "refsize{size}",
            "output",
            "quilt.down{depth}x.mspbwt.{chrom}.chunk_{chunkid}.vcf.gz",
        ),
    params:
        time=config["time"],
        N="quilt_run_mspbwt",
        nGen=config["quilt2"]["nGen"],
        buffer=config["quilt2"]["buffer"],
        start=get_quilt_chunk_region_start,
        end=get_quilt_chunk_region_end,
        Ksubset=config["quilt2"]["Ksubset"],
        Knew=config["quilt2"]["Knew"],
        nGibbsSamples=config["quilt2"]["nGibbsSamples"],
        n_seek_its=config["quilt2"]["n_seek_its"],
        lowram=config["quilt2"]["lowram"],
        rare_af_threshold=config["quilt2"]["rare_af_threshold"],
        impute_rare_common=config["quilt2"]["impute_rare_common"],
        block_gibbs=config["quilt2"]["small_ref_panel_block_gibbs_iterations"],
        gibbs_iters=config["quilt2"]["small_ref_panel_gibbs_iterations"],
        mspbwtM=config["quilt2"]["mspbwtM"],
        mspbwtL=config["quilt2"]["mspbwtL"],
    log:
        os.path.join(
            OUTDIR_QUILT2,
            "refsize{size}",
            "output",
            "quilt.down{depth}x.mspbwt.{chrom}.chunk_{chunkid}.vcf.gz.llog",
        ),
    conda:
        "../envs/quilt.yaml"
    threads: 1
    shell:
        """
        {params.time} -v QUILT.R \
            --reference_vcf_file={input.vcf} \
            --prepared_reference_filename={input.rdata} \
            --bamlist={input.bams} \
            --use_hapMatcherR={params.lowram} \
            --impute_rare_common={params.impute_rare_common} \
            --chr={wildcards.chrom} \
            --regionStart={params.start} \
            --regionEnd={params.end} \
            --buffer={params.buffer} \
            --nGen={params.nGen} \
            --zilong=FALSE \
            --use_mspbwt=TRUE \
            --mspbwtM={params.mspbwtM} \
            --mspbwtL={params.mspbwtL} \
            --Ksubset={params.Ksubset} \
            --Knew={params.Knew} \
            --nGibbsSamples={params.nGibbsSamples} \
            --n_seek_its={params.n_seek_its} \
            --rare_af_threshold={params.rare_af_threshold} \
            --small_ref_panel_block_gibbs_iterations='{params.block_gibbs}' \
            --small_ref_panel_gibbs_iterations={params.gibbs_iters} \
            --output_filename={output} &> {log}
        """


rule quilt_ligate_mspbwt:
    input:
        get_quilt_mspbwt_outputs,
    output:
        vcf=os.path.join(
            OUTDIR_QUILT2,
            "refsize{size}",
            "{chrom}",
            "quilt.down{depth}x.mspbwt.{chrom}.vcf.gz",
        ),
        lst=temp(
            os.path.join(
                OUTDIR_QUILT2,
                "refsize{size}",
                "{chrom}",
                "quilt.down{depth}x.mspbwt.{chrom}.vcf.list",
            )
        ),
    log:
        os.path.join(
            OUTDIR_QUILT2,
            "refsize{size}",
            "{chrom}",
            "quilt.down{depth}x.mspbwt.{chrom}.vcf.gz.llog",
        ),
    params:
        N="quilt_ligate_mspbwt",
        extra=config["extra_buffer_in_quilt"],
    conda:
        "../envs/quilt.yaml"
    shell:
        """
        ( \
        if [ {params.extra} -gt 0 ];then \
           echo {input} | tr ' ' '\n' > {output.lst} && \
           bcftools concat --ligate --file-list {output.lst} --output-type z --threads 4 -o {output.vcf} && \
           bcftools index -f {output.vcf} \
        ; else \
           echo {input} | tr ' ' '\n' > {output.lst} && \
           bcftools concat --file-list {output.lst} --output-type z --threads 4 -o {output.vcf} && \
           bcftools index -f {output.vcf} \
        ; fi \
        ) &> {log}
        """
