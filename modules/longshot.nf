// SNVs calling settings 

process LONGSHOT {
// queue 'batch' // problems in slurm, as transferred to partition argument --> check
time { 8.hour * task.attempt }

        //scratch true
        stageInMode = 'copy'
        stageOutMode = 'rsync'
        
        publishDir "$params.SNP_Dir/Longshot/${SampleID}_${Technology}", mode:params.SaveMode, overwrite:params.Overwrite

        input:
        each chr
        path bam
        path bai
        path genome
        path genome_index
        tuple val( SampleID ), val( Technology ), val ( Kit )

        output : 
        path "${SampleID}_${Technology}_${chr}.vcf.gz" 

        script:
        """
        source activate ONT
        zcat $genome > ARS-UCD1.2_Btau5.0.1Y.fa
        longshot --bam $bam --ref ARS-UCD1.2_Btau5.0.1Y.fa --region ${chr} --out ${SampleID}_${Technology}_${chr}.vcf
        bgzip ${SampleID}_${Technology}_${chr}.vcf
        """
}

