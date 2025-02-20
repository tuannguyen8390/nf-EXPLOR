process CLAIR3 {
label 'medium_job'
time { 24.hour * task.attempt }


        scratch true
        stageInMode = 'copy'
        stageOutMode = 'rsync'
        
        publishDir "$params.SNP_Dir/Clair3/${SampleID}_${Technology}", mode:params.SaveMode, overwrite:params.Overwrite

        input:
        each chr
        path bam
        path bai
        path genome
        path genome_index
        tuple val(SampleID), val(Technology), val (Kit) , val ( Sex )
        path clair3_model_path

        output : 

        script :
        if( "${Technology}" == 'ONT')
                """
                if [[ ${Kit} -eq 1041 ]] 
                then
                        ONT_model_path='models/r1041_e82_400bps_sup_g615'
                elif [[ ${Kit} -eq 1040 ]] 
                then
                        ONT_model_path='models/r104_e81_sup_g5015'
                elif [[ ${Kit} -eq 941 ]] 
                then
                        ONT_model_path='models/r941_prom_sup_g5014'
                fi

                nodeDir=`mktemp -d /tmp/CLAIRXXXXXX`
                echo \$ONT_model_path
                echo \$nodeDir

                run_clair3.sh --bam_fn=${bam} \
                        --ref_fn=${genome} \
                        --threads=$task.cpus \
                        --platform=ont\
                        --sample_name=${SampleID} \
                        --model_path=\$ONT_model_path \
                        --ctg_name=${chr} \
                        --remove_intermediate_dir \
                        --gvcf \
                        --output="\${nodeDir}"

                mkdir -p $params.SNP_Dir/Clair3/${SampleID}_${Technology}/${chr}
                cp -rf \$nodeDir/* $params.SNP_Dir/Clair3/${SampleID}_${Technology}/${chr}
                """
        else if( "${Technology}" == 'PB')
                """
                PB_model_path='models/hifi'

                nodeDir=`mktemp -d /tmp/CLAIRXXXXXX`
                echo \$PB_model_path
                echo \$nodeDir

                run_clair3.sh --bam_fn=${bam} \
                        --ref_fn=${genome} \
                        --threads=$task.cpus \
                        --platform=hifi \
                        --sample_name=${SampleID} \
                        --model_path=\$PB_model_path \
                        --ctg_name=${chr} \
                        --remove_intermediate_dir \
                        --gvcf \
                        --output=\${nodeDir}
                rm -rf ARS-2.0.fa  

                mkdir -p $params.SNP_Dir/Clair3/${SampleID}_${Technology}/${chr}
                cp -rf \${nodeDir}/* $params.SNP_Dir/Clair3/${SampleID}_${Technology}/${chr}
                """
        else
                """
                error "Invalid technology : ${Technology}"
                """
}


process CLAIR3_Y {
label 'medium_job'
time { 24.hour * task.attempt }


        scratch true
        stageInMode = 'copy'
        stageOutMode = 'rsync'
        
        publishDir "$params.SNP_Dir/Clair3/${SampleID}_${Technology}", mode:params.SaveMode, overwrite:params.Overwrite

        input:
        path bam
        path bai
        path genome
        path genome_index
        tuple val(SampleID), val(Technology), val (Kit) , val ( Sex )
        path clair3_model_path

        output : 

        script :
        if( "${Technology}" == 'ONT' && "${Sex}" == 'M') 
                """
                if [[ ${Kit} -eq 1041 ]] 
                then
                        ONT_model_path='models/r1041_e82_400bps_sup_g615'
                elif [[ ${Kit} -eq 1040 ]] 
                then
                        ONT_model_path='models/r104_e81_sup_g5015'
                elif [[ ${Kit} -eq 941 ]] 
                then
                        ONT_model_path='models/r941_prom_sup_g5014'
                fi

                nodeDir=`mktemp -d /tmp/CLAIRXXXXXX`
                echo \$ONT_model_path
                echo \$nodeDir

                run_clair3.sh --bam_fn=${bam} \
                        --ref_fn=${genome} \
                        --threads=$task.cpus \
                        --platform=ont\
                        --sample_name=${SampleID} \
                        --model_path=\$ONT_model_path \
                        --ctg_name=$params.chrY \
                        --remove_intermediate_dir \
                        --gvcf \
                        --output="\${nodeDir}"

                mkdir -p $params.SNP_Dir/Clair3/${SampleID}_${Technology}/$params.chrY
                cp -rf \$nodeDir/* $params.SNP_Dir/Clair3/${SampleID}_${Technology}/$params.chrY
                """
        else if( "${Technology}" == 'PB' && "${Sex}" == 'M')
                """
                PB_model_path='models/hifi'

                nodeDir=`mktemp -d /tmp/CLAIRXXXXXX`
                echo \$PB_model_path
                echo \$nodeDir

                run_clair3.sh --bam_fn=${bam} \
                        --ref_fn=${genome} \
                        --threads=$task.cpus \
                        --platform=hifi \
                        --sample_name=${SampleID} \
                        --model_path=\$PB_model_path \
                        --ctg_name=$params.chrY \
                        --remove_intermediate_dir \
                        --gvcf \
                        --output=\${nodeDir}
                rm -rf ARS-2.0.fa  

                mkdir -p $params.SNP_Dir/Clair3/${SampleID}_${Technology}/$params.chrY
                cp -rf \${nodeDir}/* $params.SNP_Dir/Clair3/${SampleID}_${Technology}/$params.chrY
                """
        else
                """
                echo "Anim is ${Sex}, skip variant calling on Y"
                """
}
