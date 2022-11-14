process CUTESV {
label 'medium_job' 
queue 'batch'
time '48h'

        scratch true
        stageInMode = 'symlink'
        stageOutMode = 'rsync'

        publishDir "$params.SV_Dir/CuteSV/${SampleID}_${Technology}", mode: params.SaveMode, overwrite: params.Overwrite

        input:
        path bam
        path bai
        path genome
        path genome_index
        tuple val(SampleID), val(Technology)
        

        output:
        file "*"

        script:
        """
        zcat ${genome} > genome.fa
        cuteSV  --threads $task.cpus --genotype --retain_work_dir --sample ${SampleID} --report_readid ${bam} genome.fa ${SampleID}_SV.vcf .
        gzip ${SampleID}_SV.vcf
        rm -rf genome.fa
        """
}