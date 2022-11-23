
process FILTLONG {
cpus = 24
memory { 128.GB * task.attempt }
time { 36.hour * task.attempt } //Very long if you doing short-read polish
errorStrategy 'retry'
maxRetries 3
queue 'batch'

        scratch true
        stageInMode = 'symlink'
        stageOutMode = 'rsync'

        publishDir "$params.QC_Dir/${SampleID_LR}_${Technology}", mode: params.SaveMode, overwrite: params.Overwrite

        input:
        tuple val( SampleID_LR ), val( Technology ), val( Shortread_Avail ), path( FASTQ_LR_Dir )
        path SR_sample

        output:
        path "${SampleID_LR}_PREQC_NanoPlot/*"
        path "${SampleID_LR}_POSTQC_NanoPlot/*"
        path "${SampleID_LR}.fastq.gz"
        path "${SampleID_LR}.fasta.gz"
        tuple val( SampleID_LR ), val( Technology ), path( "${SampleID_LR}.fastq.gz" ), path( "${SampleID_LR}.fasta.gz" ) , emit: LR_sample  
                
       script:  
        if ("${Shortread_Avail}" == 'TRUE')
                """
                dir=`grep -w ${SampleID_LR} ${SR_sample} | awk -F',' '{print \$9}'` 
                
                cat ${FASTQ_LR_Dir}/*.fastq.gz > ${SampleID_LR}_PREQC.fastq.gz

                NanoPlot -t $task.cpus --fastq ${SampleID_LR}_PREQC.fastq.gz --outdir ${SampleID_LR}_PREQC_NanoPlot --loglength --plots dot

                cat \$dir/*R1* > R1.fastq.gz
                cat \$dir/*R2* > R2.fastq.gz

                filtlong -1 R1.fastq.gz -2 R2.fastq.gz  --keep_percent 90 --trim --split 500 --mean_q_weight 10 ${SampleID_LR}_PREQC.fastq.gz | gzip > ${SampleID_LR}.fastq.gz
                
                rm R1.fastq.gz R2.fastq.gz
                rm -rf ${SampleID_LR}_PREQC.fastq.gz

                NanoPlot -t $task.cpus --fastq ${SampleID_LR}.fastq.gz --outdir ${SampleID_LR}_POSTQC_NanoPlot --loglength --plots dot

                seqtk seq -A ${SampleID_LR}.fastq.gz | bgzip -@ $task.cpus > ${SampleID_LR}.fasta.gz

                """

        else if ("${Shortread_Avail}" == 'FALSE')
                """
                cat ${FASTQ_LR_Dir}/*.fastq.gz > ${SampleID_LR}_PREQC.fastq.gz

                NanoPlot -t $task.cpus --fastq ${SampleID_LR}_PREQC.fastq.gz --outdir ${SampleID_LR}_PREQC_NanoPlot --loglength --plots dot

                filtlong --keep_percent 90 --mean_q_weight 10 ${SampleID_LR}_PREQC.fastq.gz | gzip > ${SampleID_LR}.fastq.gz
                
                rm -rf ${SampleID_LR}_PREQC.fastq.gz
                
                NanoPlot -t $task.cpus --fastq ${SampleID_LR}.fastq.gz --outdir ${SampleID_LR}_POSTQC_NanoPlot --loglength --plots dot

                seqtk seq -A ${SampleID_LR}.fastq.gz | bgzip -@ $task.cpus > ${SampleID_LR}.fasta.gz
                
                """
}

