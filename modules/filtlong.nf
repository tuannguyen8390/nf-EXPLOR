
process FILTLONG {
cpus = 24
memory { 128.GB * task.attempt }
time { 36.hour * task.attempt } //Very long if you doing short-read polish
errorStrategy 'retry'
maxRetries 3

        scratch true
        stageInMode = 'symlink'
        stageOutMode = 'rsync'

        publishDir "$params.QC_Dir/${SampleID_LR}_${Technology}", mode: params.SaveMode, overwrite: params.Overwrite

        input:
        tuple val( SampleID_LR ), val( Technology ), val (Kit), val( Shortread_Avail ), path( FASTQ_LR_Dir )
        path SR_sample

        output:
        path "${SampleID_LR}_PREQC_NanoPlot/*"
        path "${SampleID_LR}_POSTQC_NanoPlot/*"
        path "${SampleID_LR}.fastq.gz"
        path "${SampleID_LR}.fasta.gz"
        tuple val( SampleID_LR ), val( Technology ), val (Kit), path( "${SampleID_LR}.fastq.gz" ), path( "${SampleID_LR}.fasta.gz" ) , emit: LR_sample  
                
       script:  
        if ("${Shortread_Avail}" == 'TRUE' && "${Technology}" == 'ONT')
                """
                dir=`grep -w ${SampleID_LR} ${SR_sample} | awk -F',' '{print \$9}'` 
                
                cat ${FASTQ_LR_Dir}/*.fastq.gz > tmp.fastq.gz
                
                porechop_abi -abi -i tmp.fastq.gz -o ${SampleID_LR}_PREQC.fastq.gz

                NanoPlot -t $task.cpus --fastq ${SampleID_LR}_PREQC.fastq.gz --outdir ${SampleID_LR}_PREQC_NanoPlot --loglength --plots dot

                cat \$dir/*R1* > R1.fastq.gz
                cat \$dir/*R2* > R2.fastq.gz

                filtlong -1 R1.fastq.gz -2 R2.fastq.gz --min_length 200 --trim --split 1000 ${SampleID_LR}_PREQC.fastq.gz | gzip > ${SampleID_LR}.fastq.gz
                
                rm -rf ${SampleID_LR}_PREQC.fastq.gz R1.fastq.gz R2.fastq.gz

                NanoPlot -t $task.cpus --fastq ${SampleID_LR}.fastq.gz --outdir ${SampleID_LR}_POSTQC_NanoPlot --loglength --plots dot

                seqtk seq -A ${SampleID_LR}.fastq.gz | bgzip -@ $task.cpus > ${SampleID_LR}.fasta.gz

                """
        else if ("${Shortread_Avail}" == 'TRUE' && "${Technology}" == 'PB')
                """
                dir=`grep -w ${SampleID_LR} ${SR_sample} | awk -F',' '{print \$9}'` 
                
                cat ${FASTQ_LR_Dir}/*.fastq.gz > ${SampleID_LR}_PREQC.fastq.gz

                NanoPlot -t $task.cpus --fastq ${SampleID_LR}_PREQC.fastq.gz --outdir ${SampleID_LR}_PREQC_NanoPlot --loglength --plots dot

                cat \$dir/*R1* > R1.fastq.gz
                cat \$dir/*R2* > R2.fastq.gz

                filtlong -1 R1.fastq.gz -2 R2.fastq.gz --min_length 200 --trim --split 1000 ${SampleID_LR}_PREQC.fastq.gz | gzip > ${SampleID_LR}.fastq.gz
                
                rm -rf ${SampleID_LR}_PREQC.fastq.gz R1.fastq.gz R2.fastq.gz

                NanoPlot -t $task.cpus --fastq ${SampleID_LR}.fastq.gz --outdir ${SampleID_LR}_POSTQC_NanoPlot --loglength --plots dot

                seqtk seq -A ${SampleID_LR}.fastq.gz | bgzip -@ $task.cpus > ${SampleID_LR}.fasta.gz

                """
        else if ("${Shortread_Avail}" == 'FALSE' && "${Technology}" == 'ONT')
                """
                cat ${FASTQ_LR_Dir}/*.fastq.gz > tmp.fastq.gz

                porechop_abi -abi -i tmp.fastq.gz -o ${SampleID_LR}_PREQC.fastq.gz # -tmp

                gunzip -c ${SampleID_LR}_PREQC.fastq.gz  | awk 'BEGIN{FS="\t";OFS=FS;n=1;old=""};\$1!~/^@/{print;next};\$1==old{\$1 = \$1 "_" n;n=n+1;print;next};{old=\$1;n=1;print;next}' > ${SampleID_LR}_PREQC.fastq
                gzip -f ${SampleID_LR}_PREQC.fastq

                NanoPlot -t $task.cpus --fastq ${SampleID_LR}_PREQC.fastq.gz --outdir ${SampleID_LR}_PREQC_NanoPlot --loglength --plots dot

                filtlong --min_length 200 ${SampleID_LR}_PREQC.fastq.gz | gzip > ${SampleID_LR}.fastq.gz
                
                rm -rf ${SampleID_LR}_PREQC.fastq.gz
                
                NanoPlot -t $task.cpus --fastq ${SampleID_LR}.fastq.gz --outdir ${SampleID_LR}_POSTQC_NanoPlot --loglength --plots dot

                seqtk seq -A ${SampleID_LR}.fastq.gz | bgzip -@ $task.cpus > ${SampleID_LR}.fasta.gz
                """
        else if ("${Shortread_Avail}" == 'FALSE' && "${Technology}" == 'PB')
                """
                cat ${FASTQ_LR_Dir}/*.fastq.gz > ${SampleID_LR}_PREQC.fastq.gz

                NanoPlot -t $task.cpus --fastq ${SampleID_LR}_PREQC.fastq.gz --outdir ${SampleID_LR}_PREQC_NanoPlot --loglength --plots dot

                filtlong --min_length 200 ${SampleID_LR}_PREQC.fastq.gz | gzip > ${SampleID_LR}.fastq.gz
                
                rm -rf ${SampleID_LR}_PREQC.fastq.gz
                
                NanoPlot -t $task.cpus --fastq ${SampleID_LR}.fastq.gz --outdir ${SampleID_LR}_POSTQC_NanoPlot --loglength --plots dot

                seqtk seq -A ${SampleID_LR}.fastq.gz | bgzip -@ $task.cpus > ${SampleID_LR}.fasta.gz
                """
}

