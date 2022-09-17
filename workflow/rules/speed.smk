OUTDIR_REPORT = "results/reports"


rule summary_speed:
    input:
        regular=collect_quilt_log_regular,
        mspbwt=collect_quilt_log_mspbwt,
        zilong=collect_quilt_log_zilong,
    output:
        regular=os.path.join(OUTDIR_REPORT, "quilt.regular.panelsize{size}.down{depth}x.{chrom}.txt"),
        mspbwt=os.path.join(OUTDIR_REPORT, "quilt.mspbwt.panelsize{size}.down{depth}x.{chrom}.txt"),
        zilong=os.path.join(OUTDIR_REPORT, "quilt.zilong.panelsize{size}.down{depth}x.{chrom}.txt"),
    shell:
        """
        echo {input.regular} | tr ' ' '\\n' > {output.regular}
        echo {input.mspbwt} | tr ' ' '\\n' > {output.mspbwt}
        echo {input.zilong} | tr ' ' '\\n' > {output.zilong}
        """
