# Project nf-EXPLOR | Nextflow pipeline for EXPloring variation in LOng REad Sequencing 
v.0.0.2 - [Tuan Nguyen](tuan.nguyen@agriculture.vic.gov.au) 

#### Currently freely available for usage in the Bovine Long-Read Consortium (BovLRC) :cow:

 
## 1. Clone this Github

```
git clone https://github.com/tuannguyen8390/nf-EXPLOR.git
```

The pipeline deployed multiple bioinformatics software for the detection of Single Nucldeotide Polymorphism (SNPs) & Structural Variants (SVs). The pipeline (version 0.0.3) currently freely available & it was designed to deal with data from both Oxford Nanopore as well as PacBio (However we only test at the moment with ONT). Written with Nextflow DSL2.


## 2. Obtain & install Docker/Shifter/Singularity 

Installation guide for Docker can be found [here](https://docs.docker.com/get-docker/)

Installation guide for Shifter can be found [here](https://www.nersc.gov/users/software/nersc-software/shifter/)

Installation guide for Singularity can be found [here](https://sylabs.io/guides/3.5/user-guide/quick_start.html)

## 2. Edit the config file

Nextflow should operate on any system you installing it on (whether it is PBS, SLURM, AWS, Google Cloud...), all you need to do is open the `nextflow.config` file & edit a few things based on your own configuration (marked as "BASED PARAMETERS", these including things like analysis directory, executor used, adjusting computational resource...  

:triangular_flag_on_post: I suggest backing up the original `nextflow.config` so you have a reference later on. 

## 3. Pull assets (genome - ARS2.0, we suggest using this genome for reproducibility across partner of the consortium), then perform some initial setup

Run the following command to pull assets (genome) and perform some initial setup (choose 1 among Shifter/Docker/Singularity only)

```
nextflow run setup.nf -profile shifter/docker/singularity
```

## 4. Test run the pipeline (choose 1 among Shifter/Docker/Singularity only)

Edit the nextflow.config files to suit your local environment 

```
nextflow run setup.nf -profile shifter/docker/singularity,test
```

## 5. :rocket: Run the pipeline. The pipeline works using 2 metadata spreadsheet in the `meta` folder, in which:

:triangular_flag_on_post: `metadata_SR.csv` : metadata for short-read data

:triangular_flag_on_post: `metadata_LR.csv` : metadata for long-read data

```
nextflow run main.nf -profile shifter/docker/singularity
```

## 5*. If you run AWS, you can use the following command to run the pipeline

```
nextflow run main.nf -profile shifter/docker/singularity,awsbatch
```



# Pipeline overview

## 1. QC :

- [FiltLong](https://github.com/rrwick/Filtlong) : QC for both LongReads and ShortReads ( **DEFAULT + RECOMMENDED**)

- [NanoFilt](https://github.com/wdecoster/nanofilt) + [FMLRC2](https://github.com/HudsonAlpha/fmlrc2) : NanoFilt for QC of Long-Read samples, and FMLRC2 + NanoFilt for QC of Short-Read samples ( Currently **NOT _COMPATIBLE_** with PEPPER & DEEPVARIANT - use with caution !!!)

## 2. Mapping:

- [Minimap2](https://github.com/lh3/minimap2) : ( **DEFAULT for BovLRC participants**)

- [Winnowmap2](https://github.com/marbl/Winnowmap)

- [NGMLR](https://github.com/philres/ngmlr)

## 3. SNP Caller: All callers can be run in parallel & deploy per chromosome ( Chr 1 - 29 & X & Y as the pipe currently deployed in cattle )
- [Clair3](https://github.com/HKU-BAL/Clair3) : ( **DEFAULT for BovLRC participants**) - Please note that extra ONT models can be found on [Clair3_rerio_models](https://github.com/nanoporetech/rerio/tree/master/clair3_models)
 
- [PEPPER](https://github.com/kishwarshafin/pepper) - By default, Flowcell < 10.4 will be analyzed with PEPPER

- [DEEPVARIANT](https://github.com/google/deepvariant) - By default, Flowcell >= 10.4 will be analyzed with DEEPVARIANT & HIFI 

- [Longshot](https://github.com/pjedge/longshot) 

## 4. SV Caller: All callers can be run in parallel

- [Sniffles2](https://github.com/fritzsedlazeck/Sniffles) ( **DEFAULT for BovLRC participants**)

- [DYSGU](https://github.com/kcleal/dysgu)

- [CuteSV2](https://github.com/tjiangHIT/cuteSV) 

## 5. Reporting

- PRE/POST QC : NanoPlot
- Alignment Depth : Mosdepth
- MultiQC
  
---

I've absolutely no doubt that there should be some problems :). It runs on my end, but perhaps not yours. If that is the case, please email to [Tuan Nguyen](mailto:tuan.nguyen@agriculture.vic.gov.au) :star: 
