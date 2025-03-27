process MINIMAP2 {
label 'big_job' 
time '48h'
errorStrategy 'retry'
maxRetries 3

        scratch true
        stageInMode = 'symlink'
        stageOutMode = 'rsync'
        publishDir "$params.Map_Dir/${SampleID}_${Technology}", mode: params.SaveMode, overwrite: params.Overwrite

        input:
        tuple val( SampleID ), val( Technology ), val ( Kit ), val (Sex), path ( fastq ), path ( fasta )
        path genome
        path genome_index
        path mm_index

        output:
        path "*.sorted.bam", emit: bam
        path "*.sorted.bam.bai", emit: bai
        tuple val ( SampleID ), val( Technology ), val ( Kit ), val ( Sex ) , emit: info

        script:
        
        if( params.MapMethod == 'Minimap2' && "${Technology}" == 'ONT' ) 
            """
            minimap2 -ax map-ont ARS-bov-ont.mmi -t $task.cpus ${fastq} --MD | samtools sort -@ $task.cpus - -o ${SampleID}_${Technology}.sorted.bam 
            samtools index -@ $task.cpus ${SampleID}_${Technology}.sorted.bam                       
            """
        else if( params.MapMethod == 'Minimap2' && "${Technology}" == 'PB') 
            """
            minimap2 -ax map-hifi ARS-bov-hifi.mmi -t $task.cpus ${fastq} --MD | samtools sort -@ $task.cpus - -o ${SampleID}_${Technology}.sorted.bam  
            samtools index -@ $task.cpus ${SampleID}_${Technology}.sorted.bam          
            """
        else
            """
            error "Invalid alignment mode :  ${MapMethod}"
            """
}
