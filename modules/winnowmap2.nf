process WINNOWMAP2 {
label 'big_job' 
time '48h'

        scratch true
        stageInMode = 'symlink'
        stageOutMode = 'rsync'
        publishDir "$params.Map_Dir/${SampleID}_${Technology}", mode: params.SaveMode, overwrite: params.Overwrite

        input:
        tuple val(SampleID), val(Technology), val (Kit), path (fastq), path (fasta)
        path genome
        path genome_index
        path wmm_index

        output:
        path "*.sorted.bam", emit: bam
        path "*.sorted.bam.bai", emit: bai
        tuple val(SampleID), val(Technology), val (Kit), emit: info
        path mm_index
        path wmm_index
        
        script:
        if( params.MapMethod == 'Winnowmap2' && "${Technology}" == 'ONT')
        """
                winnowmap -t $task.cpus -W ${wmm_index} --MD -ax map-ont ${genome} ${fastq} | samtools sort -@ $task.cpus - -o ${SampleID}.sorted.bam 

                samtools index -@ $task.cpus ${SampleID}.sorted.bam
        """
        else if( params.MapMethod == 'Winnowmap2' && "${Technology}" == 'PB')
                """
                winnowmap -t $task.cpus -W ${wmm_index} --MD -ax map-pb ${genome} ${fastq} | samtools sort -@ $task.cpus - -o ${SampleID}.sorted.bam 
                samtools index -@ $task.cpus ${SampleID}.sorted.bam
                """
        else
                """
                error "Invalid alignment mode :  ${MapMethod}"
                """
}
