process SAMTOOLS_DEPTH {
label 'high_mem' 
time '48h'
errorStrategy 'retry'
maxRetries 3
        
        publishDir "$params.Report_Dir/Read_depth/Samtools_depth/${SampleID}", mode:params.SaveMode, overwrite:params.Overwrite

        input:
        path bam
        path bai
        tuple val(SampleID), val(Technology)

        output:
        path "mapping_summary.txt"
        
        script:
        """
        samtools depth -a ${bam} | awk '{print \$3}' -  >> raw.stats 
        datamash mean 1 q1 1 median 1 q3 1 iqr 1 < raw.stats | awk '{print "${SampleID}",\$0}' >> mapping_summary.txt
        """
}
