process MOSDEPTH {
label 'medium_job' 
time '3h'

        publishDir "$params.Report_Dir/Read_depth/Mosdepth/${SampleID}", mode:params.SaveMode, overwrite:params.Overwrite

        input:
        path bam
        path bai
        path genome
        path genome_index
        tuple val(SampleID), val(Technology), val (Kit)

        output:
        tuple val(SampleID), val(Technology), val (Kit)  , emit: map_info  
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