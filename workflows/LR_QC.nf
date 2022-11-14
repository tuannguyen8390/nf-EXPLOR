#!/usr/bin/env nextflow
/*
====================================================================================
Nextflow wrapper for Long read quality control (Using Nextflow DSL2)
Tuan Nguyen | tuan.nguyen@agriculture.vic.gov.au
08.02.2022
====================================================================================

///////////////////////////////////////////////////////////////
IMPORTANT* : Currently, this script is only for testing purpose
ONLY WORK WITH FILE IN FASTQ.GZ FORMAT
Make sure your file is formatted with bgzip and named as follow
"SAMPLE_ID.fastq.gz"


Weakness : Possible scaling up for large number of samples as FASTQ files can be split into multiple parts
///////////////////////////////////////////////////////////////
*/



nextflow.enable.dsl=2

include {NANOFILT       } from '../modules/nanofilt.nf'     
include {FILTLONG       } from '../modules/filtlong.nf'

/*
#==============================================
Define main workflow
#==============================================
*/

workflow LR_QC {
        take : 
        LR_sample_ch
        SR_sample_ch
        main :
        // Execution

        if( params.enable_nanofilt ){
                NANOFILT(LR_sample_ch, SR_sample_ch)
                LR_sample = NANOFILT.out.LR_sample
                }
        if( params.enable_filtlong ){
                FILTLONG(LR_sample_ch, SR_sample_ch)
                LR_sample = FILTLONG.out.LR_sample
                }

        emit : 
        LR_sample
}