
rule glimpse2_prepare_panel:
    input:
        refvcf=rules.subset_refpanel_by_region2.output.vcf,
    output:
        os.path.join(OUTDIR_GLIMPSE2, "refsize{size}", "{chrom}.{start}.{end}.spbwt"),
    log:
        os.path.join(
            OUTDIR_GLIMPSE2, "refsize{size}", "{chrom}.{start}.{end}.spbwt.llog"
        ),
    params:
        N="glimpse2_prepare_panel",
        time=config["time"],
        gmap=if_use_glimpse_map_in_refpanel,
        irg=get_glimpse_chunki_irg,
        org=get_glimpse_chunki_org,
    conda:
        "../envs/pandas.yaml"
    shell:
        """
        (
        if [ -s {params.gmap} ];then \
            {params.time} -v GLIMPSE2_split_reference \
            --keep-monomorphic-ref-sites \
            --reference {input.refvcf} \
            --map '{params.gmap}' \
            --input-region {params.irg} \
            --output-region {params.org} \
            --output {output} \
            --threads 4 && \
            mv {output}_*.bin {output} \
        ; else \
            {params.time} -v GLIMPSE2_split_reference \
            --keep-monomorphic-ref-sites \
            --reference {input.refvcf} \
            --input-region {params.irg} \
            --output-region {params.org} \
            --output {output} \
            --threads 4 && \
            mv {output}_*.bin {output} \
        ; fi \
        ) &> {log}
        """


rule glimpse2_phase:
    input:
        refbin=rules.glimpse2_prepare_panel.output,
        bams=rules.bamlist.output,
    output:
        temp(
            os.path.join(
                OUTDIR_GLIMPSE2,
                "panelsize{size}",
                "{chrom}",
                "down{depth}x.{chrom}.{start}.{end}.bcf",
            )
        ),
    log:
        os.path.join(
            OUTDIR_GLIMPSE2,
            "panelsize{size}",
            "{chrom}",
            "down{depth}x.{chrom}.{start}.{end}.bcf.llog",
        ),
    params:
        N="glimpse2_phase",
        time=config["time"],
        irg=get_glimpse_chunki_irg,
        org=get_glimpse_chunki_org,
        burnin=config["glimpse"]["burnin"],
        main=config["glimpse"]["main"],
        pbwtL=config["glimpse"]["pbwt-depth"],
        pbwtS=config["glimpse"]["pbwt-modulo"],
        ne=config["glimpse"]["ne"],
    conda:
        "../envs/pandas.yaml"
    shell:
        """
        (
            {params.time} -v GLIMPSE2_phase \
            --bam-list {input.bams} \
            --reference {input.refbin} \
            --burnin {params.burnin} \
            --main {params.main} \
            --pbwt-depth {params.pbwtL} \
            --pbwt-modulo {params.pbwtS} \
            --ne {params.ne} \
            --output {output} \
        ) &> {log}
        """


rule glimpse2_ligate:
    input:
        get_glimpse2_phase_outputs,
    output:
        vcf=os.path.join(
            OUTDIR_GLIMPSE2,
            "panelsize{size}",
            "{chrom}",
            "down{depth}x.{chrom}.vcf.gz",
        ),
        sample=os.path.join(
            OUTDIR_GLIMPSE2,
            "panelsize{size}",
            "{chrom}",
            "down{depth}x.{chrom}.vcf.gz.sample",
        ),
        tmp=temp(
            os.path.join(
                OUTDIR_GLIMPSE2,
                "panelsize{size}",
                "{chrom}",
                "stupid.down{depth}x.{chrom}.bcf",
            )
        ),
        lst=temp(
            os.path.join(
                OUTDIR_GLIMPSE2,
                "panelsize{size}",
                "{chrom}",
                "down{depth}x.{chrom}.vcf.list",
            ),
        ),
    log:
        os.path.join(
            OUTDIR_GLIMPSE2, "panelsize{size}", "{chrom}", "down{depth}x.{chrom}.llog"
        ),
    params:
        N="glimpse2_ligate",
        sample=config["samples"],
    conda:
        "../envs/pandas.yaml"
    shell:
        """
        echo {input} | tr ' ' '\\n' > {output.lst}
        GLIMPSE2_ligate --input {output.lst} --output {output.tmp} --threads 2 && \
        awk 'NR>1 {{ print $1 }}' {params.sample} > {output.sample} && \
        bcftools reheader -s {output.sample} -o {output.vcf} {output.tmp} && \
        bcftools index -f {output.vcf}
        """


rule glimpse_prepare_glvcf:
    input:
        bams=rules.bamlist.output,
        sites=lambda wildcards: expand(
            rules.concat_refpanel_sites_by_region2.output.sites,
            size=config["refsize"],
            allow_missing=True,
        ),
        tsv=lambda wildcards: expand(
            rules.concat_refpanel_sites_by_region2.output.tsv,
            size=config["refsize"],
            allow_missing=True,
        ),
    output:
        vcf=os.path.join(OUTDIR_GLIMPSE, "glvcf", "{chrom}", "down{depth}x.{chrom}.bcf"),
        csi=os.path.join(
            OUTDIR_GLIMPSE, "glvcf", "{chrom}", "down{depth}x.{chrom}.bcf.csi"
        ),
    log:
        os.path.join(OUTDIR_GLIMPSE, "glvcf/{chrom}/down{depth}.{chrom}.bcf.llog"),
    params:
        N="glimpse_prepare_glvcf",
        time=config["time"],
        fasta=config["genome"]["fasta"],
        bq=config["glimpse"]["bq"],
        mq=config["glimpse"]["mq"],
    conda:
        "../envs/pandas.yaml"
    shell:
        """
        (
        {params.time} -v bcftools mpileup -q {params.bq} -Q {params.mq} -f {params.fasta} -I -E -A -a 'FORMAT/DP' -r {wildcards.chrom} -T {input.sites[0]} -b {input.bams} -Ou \
        | bcftools call -Aim -C alleles -T {input.tsv[0]} -Ob -o {output.vcf} && bcftools index -f {output.vcf} \
        ) &> {log}
        """


rule glimpse_phase:
    input:
        refvcf=rules.subset_refpanel_by_region2.output.vcf,
        glvcf=rules.glimpse_prepare_glvcf.output.vcf,
    output:
        temp(
            os.path.join(
                OUTDIR_GLIMPSE,
                "panelsize{size}",
                "{chrom}",
                "down{depth}x.{chrom}.{start}.{end}.bcf",
            )
        ),
    log:
        os.path.join(
            OUTDIR_GLIMPSE,
            "panelsize{size}",
            "{chrom}",
            "down{depth}x.{chrom}.{start}.{end}.bcf.llog",
        ),
    params:
        N="glimpse_phase",
        time=config["time"],
        gmap=if_use_glimpse_map_in_refpanel,
        irg=get_glimpse_chunki_irg,
        org=get_glimpse_chunki_org,
        burnin=config["glimpse"]["burnin"],
        main=config["glimpse"]["main"],
        pbwtL=config["glimpse"]["pbwt-depth"],
        pbwtS=config["glimpse"]["pbwt-modulo"],
        ne=config["glimpse"]["ne"],
    conda:
        "../envs/pandas.yaml"
    shell:
        """
        (
        if [ -s {params.gmap} ];then \
            {params.time} -v GLIMPSE_phase \
            --input {input.glvcf} \
            --reference {input.refvcf} \
            --map '{params.gmap}' \
            --input-region {params.irg} \
            --output-region {params.org} \
            --burnin {params.burnin} \
            --main {params.main} \
            --pbwt-depth {params.pbwtL} \
            --pbwt-modulo {params.pbwtS} \
            --ne {params.ne} \
            --output {output} && \
            bcftools index -f {output} \
        ; else \
            {params.time} -v GLIMPSE_phase \
            --input {input.glvcf} \
            --reference {input.refvcf} \
            --input-region {params.irg} \
            --output-region {params.org} \
            --burnin {params.burnin} \
            --main {params.main} \
            --pbwt-depth {params.pbwtL} \
            --pbwt-modulo {params.pbwtS} \
            --ne {params.ne} \
            --output {output} && \
            bcftools index -f {output} \
        ; fi \
        ) &> {log}
        """


rule glimpse_ligate:
    input:
        get_glimpse_phase_outputs,
    output:
        vcf=os.path.join(
            OUTDIR_GLIMPSE, "panelsize{size}", "{chrom}", "down{depth}x.{chrom}.bcf"
        ),
        lst=temp(
            os.path.join(
                OUTDIR_GLIMPSE,
                "panelsize{size}",
                "{chrom}",
                "down{depth}x.{chrom}.vcf.list",
            ),
        ),
    log:
        os.path.join(
            OUTDIR_GLIMPSE, "panelsize{size}", "{chrom}", "down{depth}x.{chrom}.llog"
        ),
    params:
        N="glimpse_ligate",
    conda:
        "../envs/pandas.yaml"
    shell:
        """
        echo {input} | tr ' ' '\\n' > {output.lst}
        GLIMPSE_ligate --input {output.lst} --output {output.vcf} && bcftools index -f {output.vcf}
        """
