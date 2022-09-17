OUTDIR_REPORT = "results/reports"


rule summary_speed:
    input:
        regular=collect_quilt_log_regular,
        mspbwt=collect_quilt_log_mspbwt,
        zilong=collect_quilt_log_zilong,
    output:
        os.path.join(OUTDIR_PANEL, "panelsize{size}.down{depth}x.{chrom}.txt"),
    shell:
        """
        echo {input} > {output}
        """
