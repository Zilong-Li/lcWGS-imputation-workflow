OUTDIR_REPORT = "results/reports"


rule summary_speed:
    input:
        regular=collect_quilt_log_regular,
        mspbwt=collect_quilt_log_mspbwt,
        zilong=collect_quilt_log_zilong,
    output:
        os.path.join(OUTDIR_REPORT, "quilt.panelsize{size}.down{depth}x.{chrom}.txt"),
    params:
        awk1="{split($NF,a,\":\"); s+=a[1]*60+a[2];} END{print s/NR}",
        awk2="{s+=$NF}END{print s/NR}",
    log:
        os.path.join(OUTDIR_REPORT, "quilt.panelsize{size}.down{depth}x.{chrom}.txt.log"),
    conda:
        "../envs/quilt.yaml"
    shell:
        """
        time1=$(echo {input.regular} | tr ' ' '\\n' | xargs grep Elap | awk '{params.awk1}')
        ram1=$(echo {input.regular} | tr ' ' '\\n' | xargs grep Maximum | awk '{params.awk2}')
        time2=$(echo {input.mspbwt} | tr ' ' '\\n' | xargs grep Elap | awk '{params.awk1}')
        ram2=$(echo {input.mspbwt} | tr ' ' '\\n' | xargs grep Maximum | awk '{params.awk2}')
        time3=$(echo {input.zilong} | tr ' ' '\\n' | xargs grep Elap | awk '{params.awk1}')
        ram3=$(echo {input.zilong} | tr ' ' '\\n' | xargs grep Maximum | awk '{params.awk2}')
        echo $time1 $ram1 > {output}
        echo $time2 $ram2 >> {output}
        echo $time3 $ram3 >> {output}
        """


rule plot_speed:
    input:
        expand(rules.summary_speed.output, depth=config["downsample"], allow_missing=True),
    output:
        os.path.join(OUTDIR_REPORT, "quilt.panelsize{size}.{chrom}.pdf"),
    log:
        os.path.join(OUTDIR_REPORT, "quilt.panelsize{size}.{chrom}.pdf.llog"),
    conda:
        "../envs/quilt.yaml"
    script:
        "../scripts/plot_speed.R"
