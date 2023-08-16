process PEPPER {
label 'big_job'
time { 36.hour * task.attempt }
        
        scratch true
        //stageInMode = 'copy'
        stageOutMode = 'rsync'

        publishDir "$params.SNP_Dir/PEPPER/${SampleID}_${Technology}", mode:params.SaveMode, overwrite:params.Overwrite


        input:
        each chr
        path bed
        path bam
        path bai
        path genome
        path genome_index
        tuple val( SampleID ), val( Technology ), val ( Kit )

        output : 

        script :
        if( "${Technology}" == 'ONT' && "${Kit}" < '1041' )
                """
                if [[ ${Kit} -eq 1040 ]] 
                then
                        pepper_flag='--ont_r10_q20'
                elif [[ ${Kit} -eq 941 ]] 
                then
                        pepper_flag='--ont_r9_guppy5_sup'
                fi

                pos=`awk \'\$1==${chr}\' ARS-1.2.bed | awk '{print \$1":"\$2"-"\$3}'`
                nodeDir=`mktemp -d /tmp/PEPPERXXXXX`
                echo \$nodeDir
                echo \$pos

                run_pepper_margin_deepvariant call_variant \
                        -b ${bam} \
                        -f "${genome}" \
                        -o "\${nodeDir}" \
                        -t $task.cpus \
                        --sample_name ${SampleID} \
                        --output_prefix	${SampleID}_${Technology} \
                        --gvcf \
                        --region "\${pos}" \
                        \${pepper_flag}
                mkdir -p $params.SNP_Dir/PEPPER/${SampleID}_${Technology}/${chr}
                cp -rf \${nodeDir}/* $params.SNP_Dir/PEPPER/${SampleID}_${Technology}/${chr}
                """
        else if( "${Technology}" == 'ONT' && "${Kit}" >= '1041' )
                """
                echo "proceed with Deep Variant"
                """
        else if( "${Technology}" == 'PB')
                """
                pos=`awk \'\$1==${chr}\' ARS-1.2.bed | awk '{print \$1":"\$2"-"\$3}'`
                nodeDir=`mktemp -d /tmp/PEPPERXXXXX`
                echo \$nodeDir
                echo \$pos
                run_pepper_margin_deepvariant call_variant \
                        -b ${bam} \
                        -f "${genome}" \
                        -o "\${nodeDir}" \
                        -t $task.cpus \
                        --sample_name ${SampleID} \
                        --output_prefix	${SampleID}_${Technology} \
                        --gvcf \
                        --region "\${pos}" \
                        --hifi
                mkdir -p $params.SNP_Dir/PEPPER/${SampleID}_${Technology}/${chr}
                cp -rf \${nodeDir}/* $params.SNP_Dir/PEPPER/${SampleID}_${Technology}/${chr}
                """
        else
                """
                error "Invalid technology : ${Technology}"
                """
}
