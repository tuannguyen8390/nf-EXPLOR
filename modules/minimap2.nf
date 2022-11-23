process MINIMAP2 {
label 'big_job' 
queue 'batch'
time '48h'

        scratch true
        stageInMode = 'symlink'
        stageOutMode = 'rsync'
        publishDir "$params.Map_Dir/${SampleID}_${Technology}", mode: params.SaveMode, overwrite: params.Overwrite

        input:
        tuple val( SampleID ), val( Technology ), path ( fastq ), path ( fasta )
        path genome
        path genome_index
        path mm_index

        output:
        path "*.sorted.bam", emit: bam
        path "*.sorted.bam.bai", emit: bai
        tuple val (SampleID ), val( Technology ), emit: info

        script:
        
        if( params.MapMethod == 'Minimap2' && "${Technology}" == 'ONT' ) 
            """
            minimap2 -ax map-ont ARS-bov-ont.mmi -t $task.cpus ${fastq} --MD | samtools sort -@ $task.cpus - -o ${SampleID}.sorted.bam 
            samtools index -@ $task.cpus ${SampleID}.sorted.bam                       
            """
        else if( params.MapMethod == 'Minimap2' && "${Technology}" == 'PB') 
            """
            minimap2 -ax map-pb ARS-bov-pb.mmi -t $task.cpus ${fastq} --MD | samtools sort -@ $task.cpus - -o ${SampleID}.sorted.bam  
            samtools index -@ $task.cpus ${SampleID}.sorted.bam          
            """
        else
            """
            error "Invalid alignment mode :  ${MapMethod}"
            """
}
