
rule quilt_prepare_regular:
    input:
        vcf=rules.subset_refpanel.output.vcf,
        hap=rules.subset_refpanel.output.hap,
        leg=rules.subset_refpanel.output.leg,
    output:
        os.path.join(
            OUTDIR_PANEL,
            "panelsize{size}",
            "{chrom}",
            "prep_regular",
            "RData",
            "QUILT_prepared_reference.{chrom}.{start}.{end}.RData",
        ),
    params:
        time=config["time"],
        N="quilt_prepare_regular",
        nGen=config["quilt"]["nGen"],
        buffer=config["quilt"]["buffer"],
        gmap=if_use_quilt_map_in_refpanel,
        outdir=lambda wildcards, output: os.path.dirname(output[0])[:-5],
    log:
        os.path.join(
            OUTDIR_PANEL,
            "panelsize{size}",
            "{chrom}",
            "prep_regular",
            "RData",
            "QUILT_prepared_reference.{chrom}.{start}.{end}.RData.llog",
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
            --reference_haplotype_file={input.hap} \
            --reference_legend_file={input.leg} \
            --chr={wildcards.chrom} \
            --regionStart={wildcards.start} \
            --regionEnd={wildcards.end} \
            --buffer={params.buffer} \
            --nGen={params.nGen} \
            --use_pbwt_index=FALSE \
            --use_mspbwt=FALSE \
            --outputdir={params.outdir} \
        ; else \
        {params.time} -v QUILT_prepare_reference.R \
            --reference_vcf_file={input.vcf} \
            --reference_haplotype_file={input.hap} \
            --reference_legend_file={input.leg} \
            --chr={wildcards.chrom} \
            --regionStart={wildcards.start} \
            --regionEnd={wildcards.end} \
            --buffer={params.buffer} \
            --nGen={params.nGen} \
            --use_pbwt_index=FALSE \
            --use_mspbwt=FALSE \
            --outputdir={params.outdir} \
        ; fi
        ) &> {log}
        """


rule quilt_prepare_mspbwt:
    input:
        vcf=rules.subset_refpanel.output.vcf,
        hap=rules.subset_refpanel.output.hap,
        leg=rules.subset_refpanel.output.leg,
    output:
        os.path.join(
            OUTDIR_PANEL,
            "panelsize{size}",
            "{chrom}",
            "prep_mspbwt",
            "RData",
            "QUILT_prepared_reference.{chrom}.{start}.{end}.RData",
        ),
    params:
        time=config["time"],
        N="quilt_prepare_mspbwt",
        nGen=config["quilt"]["nGen"],
        buffer=config["quilt"]["buffer"],
        gmap=if_use_quilt_map_in_refpanel,
        nindices=config["quilt"]["mspbwt-nindices"],
        outdir=lambda wildcards, output: os.path.dirname(output[0])[:-5],
    log:
        os.path.join(
            OUTDIR_PANEL,
            "panelsize{size}",
            "{chrom}",
            "prep_mspbwt",
            "RData",
            "QUILT_prepared_reference.{chrom}.{start}.{end}.RData.llog",
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
            --reference_haplotype_file={input.hap} \
            --reference_legend_file={input.leg} \
            --chr={wildcards.chrom} \
            --regionStart={wildcards.start} \
            --regionEnd={wildcards.end} \
            --buffer={params.buffer} \
            --nGen={params.nGen} \
            --use_pbwt_index=FALSE \
            --use_mspbwt=TRUE \
            --mspbwt_nindices={params.nindices} \
            --outputdir={params.outdir} \
        ; else \
        {params.time} -v QUILT_prepare_reference.R \
            --reference_vcf_file={input.vcf} \
            --reference_haplotype_file={input.hap} \
            --reference_legend_file={input.leg} \
            --chr={wildcards.chrom} \
            --regionStart={wildcards.start} \
            --regionEnd={wildcards.end} \
            --buffer={params.buffer} \
            --nGen={params.nGen} \
            --use_pbwt_index=FALSE \
            --use_mspbwt=TRUE \
            --mspbwt_nindices={params.nindices} \
            --outputdir={params.outdir} \
        ; fi \
        ) &> {log}
        """


rule quilt_prepare_zilong:
    input:
        vcf=rules.subset_refpanel.output.vcf,
        hap=rules.subset_refpanel.output.hap,
        leg=rules.subset_refpanel.output.leg,
    output:
        os.path.join(
            OUTDIR_PANEL,
            "panelsize{size}",
            "{chrom}",
            "prep_zilong",
            "RData",
            "QUILT_prepared_reference.{chrom}.{start}.{end}.RData",
        ),
    params:
        time=config["time"],
        N="quilt_prepare_zilong",
        nGen=config["quilt"]["nGen"],
        buffer=config["quilt"]["buffer"],
        gmap=if_use_quilt_map_in_refpanel,
        outdir=lambda wildcards, output: os.path.dirname(output[0])[:-5],
    log:
        os.path.join(
            OUTDIR_PANEL,
            "panelsize{size}",
            "{chrom}",
            "prep_zilong",
            "RData",
            "QUILT_prepared_reference.{chrom}.{start}.{end}.RData.llog",
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
            --reference_haplotype_file={input.hap} \
            --reference_legend_file={input.leg} \
            --chr={wildcards.chrom} \
            --regionStart={wildcards.start} \
            --regionEnd={wildcards.end} \
            --buffer={params.buffer} \
            --nGen={params.nGen} \
            --use_pbwt_index=TRUE \
            --use_mspbwt=FALSE \
            --outputdir={params.outdir} \
        ; else \
        {params.time} -v QUILT_prepare_reference.R \
            --reference_vcf_file={input.vcf} \
            --reference_haplotype_file={input.hap} \
            --reference_legend_file={input.leg} \
            --chr={wildcards.chrom} \
            --regionStart={wildcards.start} \
            --regionEnd={wildcards.end} \
            --buffer={params.buffer} \
            --nGen={params.nGen} \
            --use_pbwt_index=TRUE \
            --use_mspbwt=FALSE \
            --outputdir={params.outdir} \
        ; fi \
        ) &> {log}
        """


rule quilt_run_regular:
    input:
        vcf=rules.subset_refpanel.output.vcf,
        hap=rules.subset_refpanel.output.hap,
        leg=rules.subset_refpanel.output.leg,
        bams=rules.bamlist.output,
        rdata=rules.quilt_prepare_regular.output,
    output:
        temp(
            os.path.join(
                OUTDIR_QUILT,
                "panelsize{size}",
                "{chrom}",
                "quilt.down{depth}x.regular.{chrom}.{start}.{end}.vcf.gz",
            )
        ),
    params:
        time=config["time"],
        N="quilt_run_regular",
        nGen=config["quilt"]["nGen"],
        buffer=config["quilt"]["buffer"],
        Ksubset=config["quilt"]["Ksubset"],
        nGibbsSamples=config["quilt"]["nGibbsSamples"],
        n_seek_its=config["quilt"]["n_seek_its"],
        n_burnin_its=config["quilt"]["n_burn_in_seek_its"],
        block_gibbs=config["quilt"]["small_ref_panel_block_gibbs_iterations"],
        gibbs_iters=config["quilt"]["small_ref_panel_gibbs_iterations"],
    log:
        os.path.join(
            OUTDIR_QUILT,
            "panelsize{size}",
            "{chrom}",
            "quilt.down{depth}x.regular.{chrom}.{start}.{end}.vcf.gz.llog",
        ),
    conda:
        "../envs/quilt.yaml"
    threads: 1
    shell:
        """
        {params.time} -v QUILT.R \
            --reference_vcf_file={input.vcf} \
            --reference_haplotype_file={input.hap} \
            --reference_legend_file={input.leg} \
            --prepared_reference_filename={input.rdata} \
            --bamlist={input.bams} \
            --chr={wildcards.chrom} \
            --regionStart={wildcards.start} \
            --regionEnd={wildcards.end} \
            --buffer={params.buffer} \
            --nGen={params.nGen} \
            --zilong=FALSE \
            --use_mspbwt=FALSE \
            --Ksubset={params.Ksubset} \
            --Knew={params.Ksubset} \
            --nGibbsSamples={params.nGibbsSamples} \
            --n_seek_its={params.n_seek_its} \
            --n_burn_in_seek_its={params.n_burnin_its} \
            --small_ref_panel_block_gibbs_iterations='{params.block_gibbs}' \
            --small_ref_panel_gibbs_iterations={params.gibbs_iters} \
            --output_filename={output} &> {log}
        """


rule quilt_ligate_regular:
    input:
        get_quilt_output_regular,
    output:
        vcf=os.path.join(
            OUTDIR_QUILT,
            "panelsize{size}",
            "{chrom}",
            "quilt.down{depth}x.regular.{chrom}.bcf.gz",
        ),
        lst=temp(
            os.path.join(
                OUTDIR_QUILT,
                "panelsize{size}",
                "{chrom}",
                "quilt.down{depth}x.regular.{chrom}.vcf.list",
            )
        ),
    log:
        os.path.join(
            OUTDIR_QUILT,
            "panelsize{size}",
            "{chrom}",
            "quilt.down{depth}x.regular.{chrom}.bcf.gz.llog",
        ),
    params:
        N="quilt_ligate_regular",
    conda:
        "../envs/quilt.yaml"
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
        rdata=rules.quilt_prepare_mspbwt.output,
    output:
        temp(
            os.path.join(
                OUTDIR_QUILT,
                "panelsize{size}",
                "{chrom}",
                "quilt.down{depth}x.mspbwt.{chrom}.{start}.{end}.vcf.gz",
            )
        ),
    params:
        time=config["time"],
        N="quilt_run_mspbwt",
        nGen=config["quilt"]["nGen"],
        buffer=config["quilt"]["buffer"],
        Ksubset=config["quilt"]["Ksubset"],
        nGibbsSamples=config["quilt"]["nGibbsSamples"],
        n_seek_its=config["quilt"]["n_seek_its"],
        n_burnin_its=config["quilt"]["n_burn_in_seek_its"],
        block_gibbs=config["quilt"]["small_ref_panel_block_gibbs_iterations"],
        gibbs_iters=config["quilt"]["small_ref_panel_gibbs_iterations"],
    log:
        os.path.join(
            OUTDIR_QUILT,
            "panelsize{size}",
            "{chrom}",
            "quilt.down{depth}x.mspbwt.{chrom}.{start}.{end}.vcf.gz.llog",
        ),
    conda:
        "../envs/quilt.yaml"
    threads: 1
    shell:
        """
        {params.time} -v QUILT.R \
            --reference_vcf_file={input.vcf} \
            --reference_haplotype_file={input.hap} \
            --reference_legend_file={input.leg} \
            --prepared_reference_filename={input.rdata} \
            --bamlist={input.bams} \
            --chr={wildcards.chrom} \
            --regionStart={wildcards.start} \
            --regionEnd={wildcards.end} \
            --buffer={params.buffer} \
            --nGen={params.nGen} \
            --zilong=FALSE \
            --use_mspbwt=TRUE \
            --Ksubset={params.Ksubset} \
            --Knew={params.Ksubset} \
            --nGibbsSamples={params.nGibbsSamples} \
            --n_seek_its={params.n_seek_its} \
            --n_burn_in_seek_its={params.n_burnin_its} \
            --small_ref_panel_block_gibbs_iterations='{params.block_gibbs}' \
            --small_ref_panel_gibbs_iterations={params.gibbs_iters} \
            --output_filename={output} &> {log}
        """


rule quilt_ligate_mspbwt:
    input:
        get_quilt_output_mspbwt,
    output:
        vcf=os.path.join(
            OUTDIR_QUILT,
            "panelsize{size}",
            "{chrom}",
            "quilt.down{depth}x.mspbwt.{chrom}.bcf.gz",
        ),
        lst=temp(
            os.path.join(
                OUTDIR_QUILT,
                "panelsize{size}",
                "{chrom}",
                "quilt.down{depth}x.mspbwt.{chrom}.vcf.list",
            )
        ),
    log:
        os.path.join(
            OUTDIR_QUILT,
            "panelsize{size}",
            "{chrom}",
            "quilt.down{depth}x.mspbwt.{chrom}.bcf.gz.llog",
        ),
    params:
        N="quilt_ligate_mspbwt",
    conda:
        "../envs/quilt.yaml"
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
        rdata=rules.quilt_prepare_zilong.output,
    output:
        temp(
            os.path.join(
                OUTDIR_QUILT,
                "panelsize{size}",
                "{chrom}",
                "quilt.down{depth}x.zilong.{chrom}.{start}.{end}.vcf.gz",
            )
        ),
    params:
        time=config["time"],
        N="quilt_run_zilong",
        nGen=config["quilt"]["nGen"],
        buffer=config["quilt"]["buffer"],
        pbwtL=config["quilt"]["pbwtL"],
        pbwtS=config["quilt"]["pbwtS"],
        Ksubset=config["quilt"]["Ksubset"],
        nGibbsSamples=config["quilt"]["nGibbsSamples"],
        n_seek_its=config["quilt"]["n_seek_its"],
        n_burnin_its=config["quilt"]["n_burn_in_seek_its"],
        block_gibbs=config["quilt"]["small_ref_panel_block_gibbs_iterations"],
        gibbs_iters=config["quilt"]["small_ref_panel_gibbs_iterations"],
    log:
        os.path.join(
            OUTDIR_QUILT,
            "panelsize{size}",
            "{chrom}",
            "quilt.down{depth}x.zilong.{chrom}.{start}.{end}.vcf.gz.llog",
        ),
    conda:
        "../envs/quilt.yaml"
    threads: 1
    shell:
        """
        {params.time} -v QUILT.R \
            --reference_vcf_file={input.vcf} \
            --reference_haplotype_file={input.hap} \
            --reference_legend_file={input.leg} \
            --prepared_reference_filename={input.rdata} \
            --bamlist={input.bams} \
            --chr={wildcards.chrom} \
            --regionStart={wildcards.start} \
            --regionEnd={wildcards.end} \
            --buffer={params.buffer} \
            --nGen={params.nGen} \
            --pbwtL={params.pbwtL} \
            --pbwtS={params.pbwtS} \
            --zilong=TRUE \
            --use_mspbwt=FALSE \
            --Ksubset={params.Ksubset} \
            --Knew={params.Ksubset} \
            --nGibbsSamples={params.nGibbsSamples} \
            --n_seek_its={params.n_seek_its} \
            --n_burn_in_seek_its={params.n_burnin_its} \
            --small_ref_panel_block_gibbs_iterations='{params.block_gibbs}' \
            --small_ref_panel_gibbs_iterations={params.gibbs_iters} \
            --output_filename={output} &> {log}
        """


rule quilt_ligate_zilong:
    input:
        get_quilt_output_zilong,
    output:
        vcf=os.path.join(
            OUTDIR_QUILT,
            "panelsize{size}",
            "{chrom}",
            "quilt.down{depth}x.zilong.{chrom}.bcf.gz",
        ),
        lst=temp(
            os.path.join(
                OUTDIR_QUILT,
                "panelsize{size}",
                "{chrom}",
                "quilt.down{depth}x.zilong.{chrom}.vcf.list",
            )
        ),
    log:
        os.path.join(
            OUTDIR_QUILT,
            "panelsize{size}",
            "{chrom}",
            "quilt.down{depth}x.zilong.{chrom}.bcf.gz.llog",
        ),
    params:
        N="quilt_ligate_zilong",
    conda:
        "../envs/quilt.yaml"
    shell:
        """
        ( \
           echo {input} | tr ' ' '\n' > {output.lst} && \
           bcftools concat --file-list {output.lst} --output-type b --threads 4 -o {output.vcf} && \
           bcftools index -f {output.vcf} \
        ) &> {log}
        """
