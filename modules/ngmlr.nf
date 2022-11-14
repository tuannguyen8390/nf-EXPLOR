
process NGMLR {
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

        output:
        path "*.sorted.bam", emit: bam
        path "*.sorted.bam.bai", emit: bai
        tuple val(SampleID), val(Technology), emit: info

        script:
        if( params.MapMethod == 'NGMLR' && "${Technology}" == 'ONT')
            """
            zcat $params.QC_Dir/${SampleID}_${Technology}/${SampleID}.fastq.gz | ngmlr --presets ont -t $task.cpus -r ${genome} | samtools view -bS - | samtools sort -@ $task.cpus - -o ${SampleID}.sorted.bam

            samtools index -@ $task.cpus ${SampleID}.sorted.bam
            """
        else if( params.MapMethod == 'NGMLR' && "${Technology}" == 'PB')
            """
            zcat $params.QC_Dir/${SampleID}_${Technology}/${SampleID}.fastq.gz | ngmlr --presets pacbio -t $task.cpus -r ${genome} | samtools view -bS - | samtools sort -@ $task.cpus - -o ${SampleID}.sorted.bam
            
            samtools index -@ $task.cpus ${SampleID}.sorted.bam
            """
        else
            """
            error "Invalid alignment mode :  ${MapMethod}"
            """

}

