process MULTIQC {
    label 'process_medium'
    time '2h'

        publishDir "$params.Report_Dir/MultiQC/", mode:params.SaveMode, overwrite:params.Overwrite

        input:
        val x

        output:
        path "*multiqc_report.html", emit: report
        path "versions.yml"        , emit: versions optional true

        script:
        """
        multiqc --ignore "*PREQC*" -f $params.Report_Dir/Read_depth/Mosdepth/ $params.QC_Dir/
        """
}