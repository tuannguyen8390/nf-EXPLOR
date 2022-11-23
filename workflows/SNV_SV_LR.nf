#!/usr/bin/env nextflow
/*
====================================================================================
Nextflow wrapper for Long read analysis of structural variants & small variants discovery (Using Nextflow DSL2)
Tuan Nguyen | tuan.nguyen@agriculture.vic.gov.au
04.02.2022  | Version 0.1
====================================================================================


===This was written in DSL2 for future scalability===
=====================================================

*/
nextflow.enable.dsl=2
/*

#==============================================
Import modules
#==============================================
*/
include { MINIMAP2       } from '../modules/minimap2.nf'
include { NGMLR          } from '../modules/ngmlr.nf'
include { WINNOWMAP2     } from '../modules/winnowmap2.nf'
include { CLAIR3         } from '../modules/clair3'
include { PEPPER         } from '../modules/pepper'
include { SNIFFLES2      } from '../modules/sniffles2'
include { CUTESV         } from '../modules/cutesv'
include { MOSDEPTH       } from '../modules/mosdepth' 
include { SAMTOOLS_DEPTH } from '../modules/samtools_depth'
include { DYSGU          } from '../modules/dysgu'
/*
#==============================================
Define main workflow
#==============================================
*/

workflow SNV_SV_LR {
        take : 
        LR_sample

        main :
        // Channels
        mm_index = Channel.fromPath("$params.GenomeDir/ARS-bov*.mmi")
                        .collect()

        wmm_index = Channel.fromPath("$params.GenomeDir/repetitive_k15.txt")
                        .collect()

        genome = Channel.fromPath("$params.GenomeDir/ARS-UCD1.2_Btau5.0.1Y.fa.gz")
                        .collect()

        genome_index = Channel.fromPath("$params.GenomeDir/ARS-UCD1.2_Btau5.0.1Y.fa.fai")
                        .collect()

        ONT_model_path = Channel.fromPath("$params.ONT_SNP_model")
                        .collect()

        PB_model_path = Channel.fromPath("$params.PB_SNP_model")
                        .collect()

        // Execution
        // Invoke mapping
        if (params.enable_minimap2){
                MINIMAP2(LR_sample,genome,genome_index,mm_index)
                map_bam = MINIMAP2.out.bam
                map_bai = MINIMAP2.out.bai
                map_info =  MINIMAP2.out.info 
        }
        else if (params.enable_winnowmap2){
                WINNOWMAP2(LR_sample,genome,genome_index,wmm_index)
                map_bam = WINNOWMAP2.out.bam
                map_bai = WINNOWMAP2.out.bai
                map_info =  WINNOWMAP2.out.info 
        }
        else if (params.enable_ngmlr){
                NGMLR(LR_sample,genome,genome_index)
                map_bam = NGMLR.out.bam
                map_bai = NGMLR.out.bai
                map_info =  NGMLR.out.info 
        }
        
        // Invoke Clair3
        if (params.enable_clair3){
                CLAIR3(map_bam, map_bai, genome, genome_index, map_info, ONT_model_path, PB_model_path)
                }
        if (params.enable_pepper && !params.enable_nanofilt){
                PEPPER(map_bam, map_bai, genome, genome_index, map_info)
                } 
        // P.E.P.P.E.R won't work with Fasta inputs - hence the exclusion of Nanofilt

        // Invoke SVs 
        if (params.enable_sniffles2) {
                SNIFFLES2(map_bam, map_bai, genome, genome_index, map_info)   
                }
        if (params.enable_cutesv) {
                CUTESV(map_bam, map_bai, genome, genome_index, map_info)   
                }
        if (params.enable_dysgu) {
                DYSGU(map_bam, map_bai, genome, genome_index, map_info)   
                }
        // Invoke calculate coverage  
        MOSDEPTH(map_bam,map_bai,genome,genome_index,map_info)
        //SAMTOOLS_DEPTH(map_bam,map_bai,map_info) // Currently disabled due to long running time & large memory requirements for large alignments                
}
