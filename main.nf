#!/usr/bin/env nextflow
/*
====================================================================================
Nextflow main workflow for BovLRC process
Tuan Nguyen | tuan.nguyen@agriculture.vic.gov.au
04.02.2022  | Version 0.1
====================================================================================

Usage

nextflow run BovLRC.nf
/
=====================================================

*/

nextflow.enable.dsl = 2

include { LR_QC         }       from './workflows/LR_QC' 
include { CREATE_LR_QC  }       from './modules/exist_QC' addParams(options: 
[
    QC_Dir:params.QC_Dir,
])
include { SNV_SV_LR     }       from './workflows/SNV_SV_LR' addParams(options: [
    LR_MetaDir:params.LR_MetaDir,
    GenomeDir: params.GenomeDir,
    ONT_SNP_model:params.ONT_SNP_model,
    PB_SNP_model:params.PB_SNP_model,
    MapMethod:params.MapMethod,
])

// Primary channels
        Channel.fromPath(params.LR_MetaDir)
                 .splitCsv(header:true)
                 .map{ row-> tuple("$row.SampleID"),("$row.Technology"),("$row.Shortread_Avail"), file("$row.FASTQ_LR_Dir")}
                 .set {LR_sample_ch}

        Channel.fromPath(params.SR_MetaDir)
                 .collect()
                 .set {SR_sample_ch}

// Main workflow
workflow {
    if (!params.enable_QC){
        CREATE_LR_QC(LR_sample_ch,SR_sample_ch)
        SNV_SV_LR(CREATE_LR_QC.out.LR_sample)
    }
    else if (!params.enable_SNV_SV){
        LR_QC(LR_sample_ch,SR_sample_ch)
    }
    else
    {
        LR_QC(LR_sample_ch,SR_sample_ch)
        SNV_SV_LR(LR_QC.out.LR_sample)
    }
}
