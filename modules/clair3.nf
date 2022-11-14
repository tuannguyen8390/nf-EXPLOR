process CLAIR3 {
label 'big_job'
queue 'batch'
time '400h'

        scratch true
        stageInMode = 'copy'
        stageOutMode = 'rsync'
        
        publishDir "$params.SNP_Dir/Clair3/${SampleID}_${Technology}", mode:params.SaveMode, overwrite:params.Overwrite

        input:
        path bam
        path bai
        path genome
        path genome_index
        tuple val(SampleID), val(Technology)
        path ONT_model_path
        path PB_model_path

        output : 

        script :
        if( "${Technology}" == 'ONT')
                """
                nodeDir=`mktemp -d /tmp/CLAIRXXXXXX`
                echo \$nodeDir
                run_clair3.sh --bam_fn=$params.Map_Dir/${SampleID}_${Technology}/${SampleID}.sorted.bam \
                        --ref_fn=${genome}  \
                        --threads=$task.cpus \
                        --platform="ont"\
                        --sample_name=${SampleID} \
                        --model_path=${ONT_model_path} \
                        --include_all_ctgs \
                        --remove_intermediate_dir \
                        --gvcf \
                        --output=\${nodeDir}
                mkdir -p $params.SNP_Dir/Clair3/${SampleID}_${Technology}
                cp -rf \$nodeDir/* $params.SNP_Dir/Clair3/${SampleID}_${Technology}
                """
        else if( "${Technology}" == 'PB')
                """
                nodeDir=`mktemp -d /tmp/CLAIRXXXXXX`
                echo \$nodeDir
                run_clair3.sh --bam_fn=$params.Map_Dir/${SampleID}_${Technology}/${SampleID}.sorted.bam \
                        --ref_fn=${genome}  \
                        --threads=$task.cpus \
                        --platform="hifi" \
                        --sample_name=${SampleID} \
                        --model_path=${PB_model_path} \
                        --include_all_ctgs \
                        --remove_intermediate_dir \
                        --gvcf \
                        --output=\${nodeDir}
                mkdir -p $params.SNP_Dir/Clair3/${SampleID}_${Technology}
                cp -rf \${nodeDir}/* $params.SNP_Dir/Clair3/${SampleID}_${Technology}
                """
        else
                """
                error "Invalid alignment mode :  ${SNPCallMethod}"
                """
}
