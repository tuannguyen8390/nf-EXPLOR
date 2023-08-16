process DEEPVARIANT {
label 'big_job'
time { 36.hour * task.attempt }
        
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
        tuple val( SampleID ), val( Technology ), val ( Kit )

        output : 

        script :
        if( "${Technology}" == 'ONT' && "${Kit}" >= '1041')
                """
                deepvariant_flag='ONT_R104'
                nodeDir=`mktemp -d /tmp/DEEPVARXXXX`
                echo \$nodeDir
                zcat ${genome} > \$nodeDir/ARS-UCD1.2_Btau5.0.1Y.fa
                cat ${genome_index} > \$nodeDir/ARS-UCD1.2_Btau5.0.1Y.fa.fai

                /opt/deepvariant/bin/run_deepvariant \
                        --reads=${bam} \
                        --ref=\$nodeDir/ARS-UCD1.2_Btau5.0.1Y.fa \
                        --output_vcf="\${nodeDir}"/${SampleID}_${Technology}.vcf.gz \
                        --output_gvcf="\${nodeDir}"/${SampleID}_${Technology}.gvcf.gz \
                        --sample_name=${SampleID} \
                        --num_shards=$task.cpus \
                        --regions="${chr}" \
                        --model_type=\${deepvariant_flag}
                rm -rf \$nodeDir/ARS-UCD1.2_Btau5.0.1Y.fa  
                rm -rf \$nodeDir/ARS-UCD1.2_Btau5.0.1Y.fa.fai
                
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
                zcat ${genome} > \$nodeDir/ARS-UCD1.2_Btau5.0.1Y.fa
                cat ${genome_index} > \$nodeDir/ARS-UCD1.2_Btau5.0.1Y.fa.fai
                /opt/deepvariant/bin/run_deepvariant \
                        --reads=${bam} \
                        --ref=\$nodeDir/ARS-UCD1.2_Btau5.0.1Y.fa \
                        --output_vcf=\${nodeDir}/${SampleID}_${Technology}.vcf.gz \
                        --output_gvcf=\${nodeDir}/${SampleID}_${Technology}.gvcf.gz \
                        --sample_name=${SampleID} \
                        --num_shards=$task.cpus \
                        --regions="${chr}" \
                        --model_type=\${deepvariant_flag}
                rm -rf \$nodeDir/ARS-UCD1.2_Btau5.0.1Y.fa  
                rm -rf \$nodeDir/ARS-UCD1.2_Btau5.0.1Y.fa.fai

                mkdir -p $params.SNP_Dir/DEEPVARIANT/${SampleID}_${Technology}/${chr}
                cp -rf \${nodeDir}/* $params.SNP_Dir/DEEPVARIANT/${SampleID}_${Technology}/${chr}


                """
        else
                """
                error "Invalid technology : ${Technology}"
                """
}
