
process NGMLR {
label 'big_job' 
time '48h'
errorStrategy 'retry'
maxRetries 3

        scratch true
        stageInMode = 'symlink'
        stageOutMode = 'rsync'
        publishDir "$params.Map_Dir/${SampleID}_${Technology}", mode: params.SaveMode, overwrite: params.Overwrite

        input:
        tuple val(SampleID), val(Technology), val (Kit), val ( Sex ), path (fastq), path (fasta)
        path genome
        path genome_index

        output:
        path "*.sorted.bam", emit: bam
        path "*.sorted.bam.bai", emit: bai
        tuple val(SampleID), val(Technology), val (Kit), val ( Sex ), emit: info

        script:
        if( params.MapMethod == 'NGMLR' && "${Technology}" == 'ONT')
            """
            zcat ${fastq} | ngmlr --presets ont -t $task.cpus -r ${genome} | samtools view -bS - | samtools sort -@ $task.cpus - -o ${SampleID}.sorted.bam

            samtools index -@ $task.cpus ${SampleID}.sorted.bam
            """
        else if( params.MapMethod == 'NGMLR' && "${Technology}" == 'PB')
            """
            zcat ${fastq} | ngmlr --presets pacbio -t $task.cpus -r ${genome} | samtools view -bS - | samtools sort -@ $task.cpus - -o ${SampleID}.sorted.bam
            
            samtools index -@ $task.cpus ${SampleID}.sorted.bam
            """
        else
            """
            error "Invalid alignment mode :  ${MapMethod}"
            """

}

