nextflow.enable.dsl=2


process make_sample_meta {
label 'single_cpu_job' 
queue 'batch'
time '1h'  
    publishDir "meta/test", mode: 'copy'

    input:

    output:
    path "test_metadata_SR.csv" 
    path "test_metadata_LR.csv"

    script:
    """
    sed "s|BASEDIR|$baseDir|g" $baseDir/meta/test/test_metadata_SR.csv > test_metadata_SR.csv 
    sed "s|BASEDIR|$baseDir|g" $baseDir/meta/test/test_metadata_LR.csv > test_metadata_LR.csv
    """
} 

process download_genome {
    publishDir "asset/genome_compact", mode: 'copy'

    input:

    output:
    path "*", emit:fagz

    script:
    """
    wget https://bovlrctest.s3.amazonaws.com/ARS-UCD1.2_Btau5.0.1Y.fa.gz
    wget https://bovlrctest.s3.amazonaws.com/ARS-UCD1.2_Btau5.0.1Y.fa.fai
    """
}

process download_asset {
    publishDir "$baseDir", mode: 'copy'

    input:

    output:

    script:
    """
    wget https://bovlrctest.s3.amazonaws.com/asset.zip
    unzip asset.zip    
    """
}


process download_clair3_model {
    publishDir "asset/models", mode: 'copy'

    input:

    output:
    path "*"

    script:
    """
    wget http://www.bio8.cs.hku.hk/clair3/clair3_models/clair3_models.tar.gz 
    tar -zxvf clair3_models.tar.gz 
    """
}

process make_map_index {
container 'tuannguyen90/clrc_mapping:1.0'
label 'single_cpu_job' 
queue 'batch'
time '3h'  
    publishDir "asset/genome_compact", mode: 'copy'

    input:
    path x

    output:
    path "*"

    script:
    """
    minimap2 -x map-ont -d ARS-bov-ont.mmi ARS-UCD1.2_Btau5.0.1Y.fa.gz
    minimap2 -x map-pb -d ARS-bov-pb.mmi ARS-UCD1.2_Btau5.0.1Y.fa.gz
    meryl count k=15 output merylDB ARS-UCD1.2_Btau5.0.1Y.fa.gz
    meryl print greater-than distinct=0.9998 merylDB > repetitive_k15.txt
    """
}

process make_gatk_index {
container 'broadinstitute/gatk:4.3.0.0'
label 'single_cpu_job' 
queue 'batch'
time '1h'  
    publishDir "asset/genome_compact", mode: 'copy'

    input:
    path x

    output:
    path "*"

    script:
    """
    gatk CreateSequenceDictionary -R ARS-UCD1.2_Btau5.0.1Y.fa.gz
    """
} 


workflow {
    make_sample_meta()
    download_genome()
    download_asset()
    download_clair3_model()
    make_map_index(download_genome.out.fagz)
    make_gatk_index(download_genome.out.fagz)
}
