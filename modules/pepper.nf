process PEPPER {
label 'big_job'
queue 'batch'
time '168h'
        
        scratch true
        stageInMode = 'copy'
        stageOutMode = 'rsync'

        publishDir "$params.SNP_Dir/PEPPER/${SampleID}_${Technology}", mode:params.SaveMode, overwrite:params.Overwrite


        input:
        path bam
        path bai
        path genome
        path genome_index
        tuple val(SampleID), val(Technology)

        
        output : 

        script :
        if( "${Technology}" == 'ONT')
                """
                nodeDir=`mktemp -d /tmp/PEPPERXXXXX`
                echo \$nodeDir
                run_pepper_margin_deepvariant call_variant \
                        -b "$params.Map_Dir/${SampleID}_${Technology}/${SampleID}.sorted.bam" \
                        -f "${genome}" \
                        -o "\${nodeDir}" \
                        -t $task.cpus \
                        --gvcf \
                        --ont_r9_guppy5_sup
                mkdir -p $params.SNP_Dir/PEPPER/${SampleID}_${Technology}
                cp -rf \${nodeDir}/* $params.SNP_Dir/PEPPER/${SampleID}_${Technology}
                """
        else if( "${Technology}" == 'PB')
                """
                nodeDir=`mktemp -d /tmp/PEPPERXXXXX`
                echo \$nodeDir
                run_pepper_margin_deepvariant call_variant \
                        -b "$params.Map_Dir/${SampleID}_${Technology}/${SampleID}.sorted.bam" \
                        -f "${genome}" \
                        -o "\${nodeDir}" \
                        -t $task.cpus \
                        --gvcf \
                        --hifi
                mkdir -p $params.SNP_Dir/PEPPER/${SampleID}_${Technology}
                cp -rf \${nodeDir}/* $params.SNP_Dir/PEPPER/${SampleID}_${Technology}
                """
        else
                """
                error "Invalid technology : ${Technology}"
                """
}
