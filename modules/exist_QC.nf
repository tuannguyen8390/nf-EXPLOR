process CREATE_LR_QC{
    errorStrategy 'terminate'

        input:
        tuple val( SampleID_LR ), val( Technology ), val( Shortread_Avail ), path( FASTQ_LR_Dir )
        path SR_sample

        output:
        file  "${SampleID_LR}.fasta.gz" optional true 
        tuple val(SampleID_LR), val(Technology) , emit: LR_sample optional true

        script:
        """
        if [ -f "$params.QC_Dir/${SampleID_LR}_${Technology}/${SampleID_LR}.fasta.gz" ]; then
            echo "Fasta files available, assuming already QC-ed"
        else 
            echo "File not found!"
            exit 127
        fi
        """
}
