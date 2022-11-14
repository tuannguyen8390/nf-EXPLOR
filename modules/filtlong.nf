
process FILTLONG {
label 'big_job' 
queue 'batch'
time '24h'
        //scratch true
        stageInMode = 'symlink'
        stageOutMode = 'rsync'

        publishDir "$params.QC_Dir/${SampleID_LR}_${Technology}", mode: params.SaveMode, overwrite: params.Overwrite

        input:
        tuple val( SampleID_LR ), val( Technology ), val( Shortread_Avail ), path( FASTQ_LR_Dir )
        path SR_sample

        output:
        file  "${SampleID_LR}.fastq.gz" optional true 
        file  "${SampleID_LR}.fasta.gz" optional true 
        tuple val(SampleID_LR), val(Technology) , emit: LR_sample optional true 
                
       script:  
        if ("${Shortread_Avail}" == 'TRUE')
                """
                dir=`grep -w ${SampleID_LR} ${SR_sample} | awk -F',' '{print \$9}'` 
                
                cat ${FASTQ_LR_Dir}/*.fastq.gz > ${SampleID_LR}_PREQC.fastq.gz

                filtlong -1 \$dir/*R1*.fastq.gz  -2 \$dir/*R2*.fastq.gz  --keep_percent 90 --trim --split 500 --mean_q_weight 10 ${SampleID_LR}_PREQC.fastq.gz | gzip > ${SampleID_LR}.fastq.gz
                
                seqtk seq -A ${SampleID_LR}.fastq.gz > ${SampleID_LR}.fasta

                bgzip -@ $task.cpus ${SampleID_LR}.fasta 
                """

        else if ("${Shortread_Avail}" == 'FALSE')
                """
                cat ${FASTQ_LR_Dir}/*.fastq.gz > ${SampleID_LR}_PREQC.fastq.gz

                filtlong --keep_percent 90 --mean_q_weight 10 ${SampleID_LR}_PREQC.fastq.gz | gzip > ${SampleID_LR}.fastq.gz
                
                seqtk seq -A ${SampleID_LR}.fastq.gz > ${SampleID_LR}.fasta
                
                bgzip -@ $task.cpus ${SampleID_LR}.fasta 
                """
}

