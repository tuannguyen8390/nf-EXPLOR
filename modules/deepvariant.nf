process DEEPVARIANT {
label 'medium_job'
time { 24.hour * task.attempt }
        
        scratch true
        stageInMode = 'copy'
        stageOutMode = 'rsync'

        publishDir "$params.SNP_Dir/DEEPVARIANT/${SampleID}_${Technology}", mode:params.SaveMode, overwrite:params.Overwrite


        input:
        each chr
        path bam
        path bai
        path genome
        path genome_index
        tuple val( SampleID ), val( Technology ), val ( Kit ), val ( Sex ) 

        output : 

        script :
        if( "${Technology}" == 'ONT' && "${Kit}" >= '1041')
                """
                deepvariant_flag='ONT_R104'
                nodeDir=`mktemp -d /tmp/DEEPVARXXXX`
                echo \$nodeDir


                /opt/deepvariant/bin/run_deepvariant \
                        --reads=${bam} \
                        --ref=${genome} \
                        --output_vcf="\${nodeDir}"/${SampleID}_${Technology}.vcf.gz \
                        --output_gvcf="\${nodeDir}"/${SampleID}_${Technology}.gvcf.gz \
                        --sample_name=${SampleID} \
                        --num_shards=$task.cpus \
                        --regions="${chr}" \
                        --model_type=\${deepvariant_flag}

                
                mkdir -p $params.SNP_Dir/DEEPVARIANT/${SampleID}_${Technology}/${chr}
                cp -rf \${nodeDir}/* $params.SNP_Dir/DEEPVARIANT/${SampleID}_${Technology}/${chr}        
                """
        else if( "${Technology}" == 'ONT' && "${Kit}" <= '1041')
                """
                echo "proceed with PEPPER"
                """
        else if( "${Technology}" == 'PB')
                """  
                deepvariant_flag='PACBIO'
                nodeDir=`mktemp -d /tmp/DEEPVARXXXX`
                echo \$nodeDir
                /opt/deepvariant/bin/run_deepvariant \
                        --reads=${bam} \
                        --ref=${genome} \
                        --output_vcf=\${nodeDir}/${SampleID}_${Technology}.vcf.gz \
                        --output_gvcf=\${nodeDir}/${SampleID}_${Technology}.gvcf.gz \
                        --sample_name=${SampleID} \
                        --num_shards=$task.cpus \
                        --regions="${chr}" \
                        --model_type=\${deepvariant_flag}


                mkdir -p $params.SNP_Dir/DEEPVARIANT/${SampleID}_${Technology}/${chr}
                cp -rf \${nodeDir}/* $params.SNP_Dir/DEEPVARIANT/${SampleID}_${Technology}/${chr}


                """
        else
                """
                error "Invalid technology : ${Technology}"
                """
}

process DEEPVARIANT_Y {
label 'medium_job'
time { 24.hour * task.attempt }
        
        scratch true
        stageInMode = 'copy'
        stageOutMode = 'rsync'

        publishDir "$params.SNP_Dir/DEEPVARIANT/${SampleID}_${Technology}", mode:params.SaveMode, overwrite:params.Overwrite


        input:
        path bam
        path bai
        path genome
        path genome_index
        tuple val( SampleID ), val( Technology ), val ( Kit ), val ( Sex ) 

        output : 

        script :
        if( "${Technology}" == 'ONT' && "${Kit}" >= '1041' && "${Sex}" == 'M')
                """
                deepvariant_flag='ONT_R104'
                nodeDir=`mktemp -d /tmp/DEEPVARXXXX`
                echo \$nodeDir

                /opt/deepvariant/bin/run_deepvariant \
                        --reads=${bam} \
                        --ref=${genome} \
                        --output_vcf="\${nodeDir}"/${SampleID}_${Technology}.vcf.gz \
                        --output_gvcf="\${nodeDir}"/${SampleID}_${Technology}.gvcf.gz \
                        --sample_name=${SampleID} \
                        --num_shards=$task.cpus \
                        --regions=$params.chrY \
                        --model_type=\${deepvariant_flag}

                
                mkdir -p $params.SNP_Dir/DEEPVARIANT/${SampleID}_${Technology}/$params.chrY
                cp -rf \${nodeDir}/* $params.SNP_Dir/DEEPVARIANT/${SampleID}_${Technology}/$params.chrY   
                """
        else if( "${Technology}" == 'ONT' && "${Kit}" <= '1041' && "${Sex}" == 'M')
                """
                echo "proceed with PEPPER"
                """
        else if( "${Technology}" == 'PB')
                """  
                deepvariant_flag='PACBIO'
                nodeDir=`mktemp -d /tmp/DEEPVARXXXX`
                echo \$nodeDir
                /opt/deepvariant/bin/run_deepvariant \
                        --reads=${bam} \
                        --ref=${genome} \
                        --output_vcf=\${nodeDir}/${SampleID}_${Technology}.vcf.gz \
                        --output_gvcf=\${nodeDir}/${SampleID}_${Technology}.gvcf.gz \
                        --sample_name=${SampleID} \
                        --num_shards=$task.cpus \
                        --regions=$params.chrY \
                        --model_type=\${deepvariant_flag}


                mkdir -p $params.SNP_Dir/DEEPVARIANT/${SampleID}_${Technology}/$params.chrY
                cp -rf \${nodeDir}/* $params.SNP_Dir/DEEPVARIANT/${SampleID}_${Technology}/$params.chrY


                """
        else
                """
                echo "Anim is ${Sex}, skip variant calling on Y"
                """
}
