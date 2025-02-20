process SNIFFLES2 {
label 'medium_job'
time { 24.hour * task.attempt } 
errorStrategy 'retry'
maxRetries 3

        scratch true
        stageInMode = 'symlink'
        stageOutMode = 'rsync'


        publishDir "$params.SV_Dir/Sniffles2/${SampleID}_${Technology}", mode: params.SaveMode, overwrite: params.Overwrite

        input:
        path bam
        path bai
        path genome
        path genome_index
        tuple val(SampleID), val(Technology), val (Kit), val ( Sex ) 

        output:
        file "*"

        script:
        """
        sniffles --input ${bam} --vcf ${SampleID}_${Technology}_SV.vcf.gz --snf ${SampleID}_${Technology}.snf --reference ${genome} --output-rnames --threads $task.cpus
        
        bgzip ${SampleID}_${Technology}.snf
        """
}

// sniffles --version | head -n 1 | sed 's/ Version //' >> versions.txt