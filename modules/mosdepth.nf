process MOSDEPTH {
label 'medium_job' 
time '3h'

        publishDir "$params.Report_Dir/Read_depth/Mosdepth/${SampleID}", mode:params.SaveMode, overwrite:params.Overwrite

        input:
        path bam
        path bai
        path genome
        path genome_index
        tuple val(SampleID), val(Technology), val (Kit)

        output:
        file "*.txt"

        script:
        """
        mosdepth  --threads $task.cpus -n --fasta params.GenomeDir/ARS-UCD1.2_Btau5.0.1Y.fa ${SampleID}_${Technology} ${bam}
        """
}