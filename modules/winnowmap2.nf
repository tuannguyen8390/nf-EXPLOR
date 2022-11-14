process WINNOWMAP2 {
label 'big_job' 
queue 'batch'
time '48h'

        scratch true
        stageInMode = 'symlink'
        stageOutMode = 'rsync'
        publishDir "$params.Map_Dir/${SampleID}_${Technology}", mode: params.SaveMode, overwrite: params.Overwrite

        input:
        tuple val(SampleID), val(Technology)
        path genome
        path genome_index
        path wmm_index

        output:
        path "*.sorted.bam", emit: bam
        path "*.sorted.bam.bai", emit: bai
        tuple val(SampleID), val(Technology), emit: info
        path mm_index
        path wmm_index
        
        script:
        if( params.MapMethod == 'Winnowmap2' && "${Technology}" == 'ONT')
        """
                winnowmap -t $task.cpus -W ${wmm_index} --MD -ax map-ont ${genome} $params.QC_Dir/${SampleID}_${Technology}/${SampleID}.fastq.gz | samtools sort -@ $task.cpus - -o ${SampleID}.sorted.bam 

                samtools index -@ $task.cpus ${SampleID}.sorted.bam
        """
        else if( params.MapMethod == 'Winnowmap2' && "${Technology}" == 'PB')
                """
                winnowmap -t $task.cpus -W ${wmm_index} --MD -ax map-pb ${genome} $params.QC_Dir/${SampleID}_${Technology}/${SampleID}.fastq.gz | samtools sort -@ $task.cpus - -o ${SampleID}.sorted.bam 
                samtools index -@ $task.cpus ${SampleID}.sorted.bam
                """
        else
                """
                error "Invalid alignment mode :  ${MapMethod}"
                """
}
