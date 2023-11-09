process NANOFILT {
cpus = 24
memory '128 GB'
time '72h'
errorStrategy 'retry'
maxRetries 3

        //scratch true
        stageInMode = 'symlink'
        stageOutMode = 'rsync'
        publishDir "$params.QC_Dir/${SampleID_LR}_${Technology}", mode: params.SaveMode, overwrite: params.Overwrite
        
        input:
        tuple val( SampleID_LR ), val( Technology ), val (Kit), val( Shortread_Avail ), path( FASTQ_LR_Dir )
        path SR_sample

        output:
        path "${SampleID_LR}_PREQC/*"
        path "${SampleID_LR}_POSTQC/*"      
        path "${SampleID_LR}.fasta.gz"
        path "${SampleID_LR}.fastq.gz" 
        tuple val( SampleID_LR ), val( Technology ), val (Kit), path ( "${SampleID_LR}.fasta.gz" ), path ( "${SampleID_LR}.fastq.gz" ), emit: LR_sample
        // With FMLRC2, the output is a FASTA file, not a FASTQ file, hence we have to swap the position & use the PREQC.fastq.gz as a "dummy" file

        script:
        if ("${Shortread_Avail}" == 'FALSE' && "${Technology}" == 'ONT')
                """
                cat ${FASTQ_LR_Dir}/*.fastq.gz > ${SampleID_LR}_PREQC.fastq.gz
                        
                NanoPlot -t $task.cpus --fastq ${SampleID_LR}_PREQC.fastq.gz --outdir ${SampleID_LR}_PREQC --loglength --plots dot

                bgzip -d -c ${SampleID_LR}_PREQC.fastq.gz -@ $task.cpus | NanoFilt -l 200 | bgzip -@ $task.cpus > ${SampleID_LR}.fastq.gz
        
                NanoPlot -t $task.cpus --fastq ${SampleID_LR}.fastq.gz --outdir ${SampleID_LR}_POSTQC --loglength --plots dot

                seqtk seq -A ${SampleID_LR}.fastq.gz | bgzip -@ $task.cpus > ${SampleID_LR}.fasta.gz
                
                rm -rf ${SampleID_LR}_PREQC.fastq.gz
                """
        else if ("${Shortread_Avail}" == 'FALSE' && "${Technology}" == 'PB')
                """
                cat ${FASTQ_LR_Dir}/*.fastq.gz > ${SampleID_LR}_PREQC.fastq.gz
                        
                NanoPlot -t $task.cpus --fastq ${SampleID_LR}_PREQC.fastq.gz --outdir ${SampleID_LR}_PREQC --loglength --plots dot

                bgzip -d -c ${SampleID_LR}_PREQC.fastq.gz -@ $task.cpus | NanoFilt -l 200 | bgzip -@ $task.cpus > ${SampleID_LR}.fastq.gz
        
                NanoPlot -t $task.cpus --fastq ${SampleID_LR}.fastq.gz --outdir ${SampleID_LR}_POSTQC --loglength --plots dot

                seqtk seq -A ${SampleID_LR}.fastq.gz | bgzip -@ $task.cpus > ${SampleID_LR}.fasta.gz
                
                rm -rf ${SampleID_LR}_PREQC.fastq.gz
                """
        else if ("${Shortread_Avail}" == 'TRUE' && "${Technology}" == 'ONT')
                """
                dir=`grep -w ${SampleID_LR} ${SR_sample} | awk -F',' '{print \$9}'` 

                bgzip -d -c \$dir/*R1*.fastq.gz \$dir/*R2*.fastq.gz -@ $task.cpus | awk 'NR % 4 == 2' | sort | tr NT TN | ropebwt2 -LR | tr NT TN | fmlrc2-convert ${SampleID_LR}.npy

                cat ${FASTQ_LR_Dir}/*.fastq.gz > ${SampleID_LR}_PREQC.fastq.gz

                NanoPlot -t $task.cpus --fastq ${SampleID_LR}_PREQC.fastq.gz --outdir ${SampleID_LR}_PREQC --loglength --plots hex dot

                bgzip -d -c  ${SampleID_LR}_PREQC.fastq.gz -@ $task.cpus | NanoFilt -l 200 | awk 'NR%4==1||NR%4==2' | tr "@" ">" > ${SampleID_LR}_raw.fasta

                fmlrc2 -t $task.cpus --cache_size d12 ${SampleID_LR}.npy ${SampleID_LR}_raw.fasta ${SampleID_LR}.fasta

                bgzip -@ $task.cpus ${SampleID_LR}.fasta
                ln -s ${SampleID_LR}.fasta.gz ${SampleID_LR}.fastq.gz

                NanoPlot -t $task.cpus --fasta ${SampleID_LR}.fasta.gz --outdir ${SampleID_LR}_POSTQC --loglength --plots dot
                
                rm -rf ${SampleID_LR}_PREQC.fastq.gz
                """
        else if ("${Shortread_Avail}" == 'TRUE' && "${Technology}" == 'PB')
                """
                dir=`grep -w ${SampleID_LR} ${SR_sample} | awk -F',' '{print \$9}'` 

                bgzip -d -c \$dir/*R1*.fastq.gz \$dir/*R2*.fastq.gz -@ $task.cpus | awk 'NR % 4 == 2' | sort | tr NT TN | ropebwt2 -LR | tr NT TN | fmlrc2-convert ${SampleID_LR}.npy

                cat ${FASTQ_LR_Dir}/*.fastq.gz > ${SampleID_LR}_PREQC.fastq.gz

                NanoPlot -t $task.cpus --fastq ${SampleID_LR}_PREQC.fastq.gz --outdir ${SampleID_LR}_PREQC --loglength --plots hex dot

                bgzip -d -c  ${SampleID_LR}_PREQC.fastq.gz -@ $task.cpus | NanoFilt -l 200 | awk 'NR%4==1||NR%4==2' | tr "@" ">" > ${SampleID_LR}_raw.fasta

                fmlrc2 -t $task.cpus --cache_size 12 ${SampleID_LR}.npy ${SampleID_LR}_raw.fasta ${SampleID_LR}.fasta

                bgzip -@ $task.cpus ${SampleID_LR}.fasta
                ln -s ${SampleID_LR}.fasta.gz ${SampleID_LR}.fastq.gz

                NanoPlot -t $task.cpus --fasta ${SampleID_LR}.fasta.gz --outdir ${SampleID_LR}_POSTQC --loglength --plots dot
                
                rm -rf ${SampleID_LR}_PREQC.fastq.gz
                """
}