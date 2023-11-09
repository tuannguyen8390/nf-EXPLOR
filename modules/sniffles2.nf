process SNIFFLES2 {
cpus = 24
memory { 64.GB * task.attempt }
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
        tuple val(SampleID), val(Technology), val (Kit)

        output:
        file "*"

        script:
        """
        sniffles --input ${bam} --vcf ${SampleID}_SV.vcf --snf ${SampleID}.snf --output-rnames --threads $task.cpus
        
        bgzip -@ $task.cpus ${SampleID}_SV.vcf
        bgzip -@ $task.cpus ${SampleID}.snf
        """
}

// sniffles --version | head -n 1 | sed 's/ Version //' >> versions.txt