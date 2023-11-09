process CREATE_LR_QC{
        errorStrategy 'terminate'
        scratch true
        stageInMode = 'symlink'
        stageOutMode = 'rsync'

        input:
        tuple val( SampleID_LR ), val( Technology ), val ( Kit ), val( Shortread_Avail ), path( FASTQ_LR_Dir )
        path SR_sample
        path params.QC_Dir

        output:
        tuple val( SampleID_LR ), val( Technology ) , val ( Kit ), path ( "${SampleID_LR}.fastq.gz" ),  path ( "${SampleID_LR}.fasta.gz" ), emit: LR_sample optional true

        script:
        """
        if [ -f "$params.QC_Dir/${SampleID_LR}_${Technology}/${SampleID_LR}.fastq.gz" ]; then
            ln -s $params.QC_Dir/${SampleID_LR}_${Technology}/${SampleID_LR}.fastq.gz . 
            ln -s $params.QC_Dir/${SampleID_LR}_${Technology}/${SampleID_LR}.fasta.gz . 
        elif [ -f "$params.QC_Dir/${SampleID_LR}_${Technology}/${SampleID_LR}.fasta" ]; then
            echo "File not found!"
            exit 127
        fi
        """
}
