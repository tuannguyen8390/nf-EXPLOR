process DYSGU {
label 'small_job' 
queue 'batch'
time '24h'
container 'kcleal/dysgu'

        scratch true
        stageInMode = 'symlink'
        stageOutMode = 'rsync'


        publishDir "$params.SV_Dir/Dysgu/${SampleID}_${Technology}", mode: params.SaveMode, overwrite: params.Overwrite

        input:
        path bam
        path bai
        path genome
        path genome_index
        tuple val(SampleID), val(Technology)

        output:
        file "*.vcf.gz"

        script:
        if( "${Technology}" == 'ONT' ) 
        """
        temp_file="/tmp/\$RANDOM"

        dysgu call --mode nanopore ${genome} \$temp_file ${bam} | gzip - > ${SampleID}_${Technology}.vcf.gz
        """
        else if( "${Technology}" == 'PB') 
        """
        temp_file="/tmp/\$RANDOM"

        dysgu call --mode pacbio ${genome} \$temp_file ${bam} | gzip - > ${SampleID}.vcf.gz
        """
}
