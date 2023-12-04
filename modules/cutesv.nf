process CUTESV {
label 'medium_job' 
time { 12.hour * task.attempt } 

        scratch true
        stageInMode = 'symlink'
        stageOutMode = 'rsync'

        publishDir "$params.SV_Dir/CuteSV/${SampleID}_${Technology}", mode: params.SaveMode, overwrite: params.Overwrite

        input:
        path bam
        path bai
        path genome
        path genome_index
        tuple val( SampleID ), val ( Technology ), val ( Kit )
        
        output:
        file "*.vcf.gz"

        script:
        """
        
        cuteSV --threads $task.cpus \
                --genotype \
                --max_cluster_bias_INS 100 \
                --diff_ratio_merging_INS 0.3 \
	        --max_cluster_bias_DEL 100 \
	        --diff_ratio_merging_DEL 0.3 \
                --sample ${SampleID} \
                --report_readid ${bam} \
                ${genome} \
                ${SampleID}_${Technology}.vcf \
                .

        bgzip -@ $task.cpu ${SampleID}_${Technology}.vcf
        """
}