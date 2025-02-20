process LONGSHOT {
time { 24.hour * task.attempt }
label 'medium_job'
errorStrategy 'retry'
maxRetries 3

        //scratch true
        scratch true
        stageInMode = 'copy'
        stageOutMode = 'rsync'
        
        publishDir "$params.SNP_Dir/Longshot/${SampleID}_${Technology}", mode:params.SaveMode, overwrite:params.Overwrite

        input:
        each chr
        path bam
        path bai
        path genome
        path genome_index
        tuple val( SampleID ), val( Technology ), val ( Kit ), val ( Sex ) 

        output : 
        path "${SampleID}_${Technology}_${chr}.vcf.gz" 

        script:
        """
        longshot --bam $bam --ref ${genome} --region ${chr} --out ${SampleID}_${Technology}_${chr}.vcf
        bgzip ${SampleID}_${Technology}_${chr}.vcf
        """
}

