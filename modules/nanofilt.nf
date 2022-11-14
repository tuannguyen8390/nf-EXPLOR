process NANOFILT {
label 'medium_job' 
queue 'batch'
time '24h'

        scratch true
        stageInMode = 'symlink'
        stageOutMode = 'rsync'
        publishDir "$params.QC_Dir/${SampleID_LR}_${Technology}", mode: params.SaveMode, overwrite: params.Overwrite
        
        input:
        tuple val( SampleID_LR ), val( Technology ), val( Shortread_Avail ), path( FASTQ_LR_Dir )
        path SR_sample

        output:
        path "${SampleID_LR}.fastq.gz" optional true
        path "${SampleID_LR}.fasta.gz" optional true
        tuple val( SampleID_LR ), val( Technology ) , emit: LR_sample optional true 
        
        script:
        if ("${Shortread_Avail}" == 'FALSE')
        """

                cat ${FASTQ_LR_Dir}/*.fastq.gz | bgzip -d -c - -@ $task.cpus | NanoFilt -q 10 > ${SampleID_LR}.fastq
                bgzip -@ $task.cpus ${SampleID_LR}.fastq 
                seqtk seq -A ${SampleID_LR}.fastq > ${SampleID_LR}.fasta
                bgzip -@ $task.cpus ${SampleID_LR}.fasta 
                """
        
        else if ("${Shortread_Avail}" == 'TRUE')
                """
                dir=`grep -w ${SampleID_LR} ${SR_sample} | awk -F',' '{print \$9}'` 
                bgzip -d -c \$dir/*R1*.fastq.gz \$dir/*R2*.fastq.gz -@ $task.cpus | awk 'NR % 4 == 2' | sort | tr NT TN | ropebwt2 -LR | tr NT TN | fmlrc2-convert ${SampleID_LR}.npy

                cat ${FASTQ_LR_Dir}/*.fastq.gz | bgzip -d -c - -@ $task.cpus | NanoFilt -q 10 | awk 'NR%4==1||NR%4==2' | tr "@" ">" > ${SampleID_LR}_raw.fasta

                fmlrc2 -t $task.cpus --cache_size 12 ${SampleID_LR}.npy ${SampleID_LR}_raw.fasta ${SampleID_LR}.fasta

                bgzip -@ $task.cpus ${SampleID_LR}.fasta
                """
}