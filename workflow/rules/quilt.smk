rule quilt_prepare_mspbwt:
    input:
        vcf=rules.subset_refpanel.output.vcf,
        hap=rules.subset_refpanel.output.hap,
        leg=rules.subset_refpanel.output.leg,
    output:
        "results/quilt/panelsize{size}/{chrom}/prep_mspbwt/RData/QUILT_prepared_reference.{chrom}.{start}.{end}.RData",
    params:
        nGen=config["quilt"]["nGen"],
        buffer=config["quilt"]["buffer"],
        outdir=lambda wildcards, output: os.path.dirname(output[0])[:-5],
    log:
        "results/quilt/panelsize{size}/{chrom}/prep_mspbwt/RData/QUILT_prepared_reference.{chrom}.{start}.{end}.RData.llog",
    conda:
        "../envs/quilt.yaml"
    threads: 1
    shell:
        """
        /usr/bin/time -v /gpfs3/users/davies/xxd908/local/pkgs/QUILT/QUILT_prepare_reference.R \
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
            --outputdir={params.outdir} &> {log}
        """

rule quilt_prepare_zilong:
    input:
        vcf=rules.subset_refpanel.output.vcf,
        hap=rules.subset_refpanel.output.hap,
        leg=rules.subset_refpanel.output.leg,
    output:
        "results/quilt/panelsize{size}/{chrom}/prep_zilong/RData/QUILT_prepared_reference.{chrom}.{start}.{end}.RData",
    params:
        nGen=config["quilt"]["nGen"],
        buffer=config["quilt"]["buffer"],
        outdir=lambda wildcards, output: os.path.dirname(output[0])[:-5],
    log:
        "results/quilt/panelsize{size}/{chrom}/prep_zilong/RData/QUILT_prepared_reference.{chrom}.{start}.{end}.RData.llog",
    conda:
        "../envs/quilt.yaml"
    threads: 1
    shell:
        """
        /usr/bin/time -v /gpfs3/users/davies/xxd908/local/pkgs/QUILT/QUILT_prepare_reference.R \
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
            --outputdir={params.outdir} &> {log}
        """

rule quilt_run_regular:
    input:
        vcf=rules.subset_refpanel.output.vcf,
        hap=rules.subset_refpanel.output.hap,
        leg=rules.subset_refpanel.output.leg,
        bams=rules.bamlist.output,
        rdata=rules.quilt_prepare_mspbwt.output,
    output:
        temp(
            "results/quilt/panelsize{size}/{chrom}/quilt.{depth}x.regular.{chrom}.{start}.{end}.vcf.gz"
        ),
    params:
        nGen=config["quilt"]["nGen"],
        buffer=config["quilt"]["buffer"],
    log:
        os.path.join(
            "results/quilt/panelsize{size}/{chrom}/quilt.{depth}x.regular.{chrom}.{start}.{end}.vcf.gz.llog"
        ),
    conda:
        "../envs/quilt.yaml"
    threads: 1
    shell:
        """
        /usr/bin/time -v /gpfs3/users/davies/xxd908/local/pkgs/QUILT/QUILT.R \
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
            --output_filename={output} &> {log}
        """


rule quilt_ligate_regular:
    input:
        get_quilt_output_regular,
    output:
        vcf="results/quilt/panelsize{size}/{chrom}/quilt.{depth}x.regular.{chrom}.bcf.gz",
        lst=temp(
            "results/quilt/panelsize{size}/{chrom}/quilt.{depth}x.regular.{chrom}.vcf.list"
        ),
    log:
        "results/quilt/panelsize{size}/{chrom}/quilt.{depth}x.regular.{chrom}.bcf.gz.llog",
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
            "results/quilt/panelsize{size}/{chrom}/quilt.{depth}x.mspbwt.{chrom}.{start}.{end}.vcf.gz"
        ),
    params:
        nGen=config["quilt"]["nGen"],
        buffer=config["quilt"]["buffer"],
    log:
        "results/quilt/panelsize{size}/{chrom}/quilt.{depth}x.mspbwt.{chrom}.{start}.{end}.vcf.gz.llog",
    conda:
        "../envs/quilt.yaml"
    threads: 1
    shell:
        """
        /usr/bin/time -v /gpfs3/users/davies/xxd908/local/pkgs/QUILT/QUILT.R \
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
            --output_filename={output} &> {log}
        """


rule quilt_ligate_mspbwt:
    input:
        get_quilt_output_mspbwt,
    output:
        vcf="results/quilt/panelsize{size}/{chrom}/quilt.{depth}x.mspbwt.{chrom}.bcf.gz",
        lst=temp(
            "results/quilt/panelsize{size}/{chrom}/quilt.{depth}x.mspbwt.{chrom}.vcf.list"
        ),
    log:
        "results/quilt/panelsize{size}/{chrom}/quilt.{depth}x.mspbwt.{chrom}.bcf.gz.llog",
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
            "results/quilt/panelsize{size}/{chrom}/quilt.{depth}x.zilong.{chrom}.{start}.{end}.vcf.gz"
        ),
    params:
        nGen=config["quilt"]["nGen"],
        buffer=config["quilt"]["buffer"],
        pbwtL=config["quilt"]["pbwtL"],
        pbwtS=config["quilt"]["pbwtS"],
    log:
        "results/quilt/panelsize{size}/{chrom}/quilt.{depth}x.zilong.{chrom}.{start}.{end}.vcf.gz.llog",
    conda:
        "../envs/quilt.yaml"
    threads: 1
    shell:
        """
        /usr/bin/time -v /gpfs3/users/davies/xxd908/local/pkgs/QUILT/QUILT.R \
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
            --output_filename={output} &> {log}
        """


rule quilt_ligate_zilong:
    input:
        get_quilt_output_zilong,
    output:
        vcf="results/quilt/panelsize{size}/{chrom}/quilt.{depth}x.zilong.{chrom}.bcf.gz",
        lst=temp(
            "results/quilt/panelsize{size}/{chrom}/quilt.{depth}x.zilong.{chrom}.vcf.list"
        ),
    log:
        "results/quilt/panelsize{size}/{chrom}/quilt.{depth}x.zilong.{chrom}.bcf.gz.llog",
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
