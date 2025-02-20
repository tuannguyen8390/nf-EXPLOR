process MOSDEPTH {
label 'medium_job' 
time '3h'
errorStrategy 'retry'
maxRetries 3

        publishDir "$params.Report_Dir/Read_depth/Mosdepth/${SampleID}_${Technology}", mode:params.SaveMode, overwrite:params.Overwrite

        input:
        path bam
        path bai
        path genome
        path genome_index
        tuple val(SampleID), val(Technology), val (Kit), val ( Sex )

        output:
        tuple val(SampleID), val(Technology), val (Kit), val ( Sex )  , emit: map_info  
        file "*"
        
        script:
        """
        mosdepth  --threads $task.cpus \
                  -n \
                  --by 1000 \
                  --fasta ${genome} \
                  ${SampleID}_${Technology} \
                  ${bam}
        """
}