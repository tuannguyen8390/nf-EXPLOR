process CUTESV {
label 'big_job' 
queue 'batch'
time '24h'
clusterOptions = "--account='dbioanim6'"

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
        zcat ${genome} > genome.fa
        cuteSV --threads $task.cpus \
                --genotype \
                --max_cluster_bias_INS 100 \
                --diff_ratio_merging_INS 0.3 \
	        --max_cluster_bias_DEL 100 \
	        --diff_ratio_merging_DEL 0.3 \
                --sample ${SampleID} \
                --report_readid ${bam} \
                genome.fa \
                ${SampleID}_${Technology}.vcf \
                .

        gzip ${SampleID}_${Technology}.vcf
        """
}