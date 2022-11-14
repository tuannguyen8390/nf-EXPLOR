process SNIFFLES2 {
label 'medium_job' 
queue 'batch'
time '6h'

        scratch true
        stageInMode = 'symlink'
        stageOutMode = 'rsync'


        publishDir "$params.SV_Dir/Sniffles2/${SampleID}_${Technology}", mode: params.SaveMode, overwrite: params.Overwrite

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
        sniffles --input ${bam} --vcf ${SampleID}_SV.vcf --snf ${SampleID}.snf --output-rnames --threads $task.cpus
        
        gzip ${SampleID}_SV.vcf
        gzip ${SampleID}.snf
        """
}

// sniffles --version | head -n 1 | sed 's/ Version //' >> versions.txt