# Project nf-EXPLOR | Nextflow pipeline for EXPloring variation in LOng REad Sequencing

*Currently freely available for usage in the Bovine Long-Read Consortium (BovLRC)
*

[Tuan Nguyen](tuan.nguyen@agriculture.vic.gov.au)

###

Initial setup: 
1. Clone this Github
```
git clone https://github.com/tuannguyen8390/AgVic_CLRC.git
```

The pipeline deployed multiple bioinformatics software for the detection of Single Nucldeotide Polymorphism (SNPs) & Structural Variants (SVs). The pipeline (version 0.0.3) currently freely available & it was designed to deal with data from both Oxford Nanopore as well as PacBio (However we only test at the moment with ONT). Written with Nextflow DSL2.


2. Obtain & install Docker/Shifter/Singularity 
- Installation guide for Docker can be found [here](https://docs.docker.com/get-docker/)
- Installation guide for Shifter can be found [here](https://www.nersc.gov/users/software/nersc-software/shifter/)
- Installation guide for Singularity can be found [here](https://sylabs.io/guides/3.5/user-guide/quick_start.html)

3. Pull assets (genome) and perform some initial setup
- Run the following command to pull assets (genome) and perform some initial setup (choose 1 among Shifter/Docker/Singularity only)
```
nextflow run setup.nf -profile shifter/docker/singularity
```

4. Test run the pipeline ((choose 1 among Shifter/Docker/Singularity only)

Edit the nextflow.config files to suit your local environment 

```
nextflow run setup.nf -profile shifter/docker/singularity,test
```

5. Run the pipeline
```
nextflow run main.nf -profile shifter/docker/singularity
```

5*. If you run AWS, you can use the following command to run the pipeline
```
nextflow run main.nf -profile shifter/docker/singularity,awsbatch
```

The pipeline works using 2 metadata spreadsheet in the `meta` folder, in which:

`metadata_SR.csv` : metadata for short-read data

`metadata_LR.csv` : metadata for long-read data

Please refer to these files in the `nextflow.config`. Simply forking/cloning the github & editing the file per your own configuration - as we cannot customize the pipeline for everyone :).

---

## Pipeline overview

1. QC :

- [FiltLong](https://github.com/rrwick/Filtlong) : QC for both LongReads and ShortReads (**DEFAULT + RECOMMENDED**)

- [NanoFilt](https://github.com/wdecoster/nanofilt) + [FMLRC2](https://github.com/HudsonAlpha/fmlrc2) : NanoFilt for QC of Long-Read samples, and FMLRC2 + NanoFilt for QC of Short-Read samples ( Currently **NOT _COMPATIBLE_** with PEPPER & DEEPVARIANT - use with caution !!!)

2. Mapping:

- [Minimap2](https://github.com/lh3/minimap2) : (**DEFAULT for BovLRC participants**)

- [Winnowmap2](https://github.com/marbl/Winnowmap)

- [NGMLR](https://github.com/philres/ngmlr)

3. SNP Caller: All callers can be run in parallel & deploy per chromosome ( Chr 1 - 29 & X as the pipe currently deployed in cattle )
- [Clair3](https://github.com/HKU-BAL/Clair3) : (**DEFAULT for BovLRC participants**)

- [PEPPER](https://github.com/kishwarshafin/pepper) - By default, Flowcell < 10.4 will be analyzed with PEPPER

- [DEEPVARIANT](https://github.com/google/deepvariant) - By default, Flowcell >= 10.4 will be analyzed with DEEPVARIANT & HIFI 

- [Longshot](https://github.com/pjedge/longshot) 

4. SV Caller: All callers can be run in parallel

- [Sniffles2](https://github.com/fritzsedlazeck/Sniffles) (**DEFAULT for BovLRC participants**)

- [DYSGU](https://github.com/kcleal/dysgu)

- [CuteSV2](https://github.com/tjiangHIT/cuteSV) 

5. Reporting
- PRE/POST QC : NanoPlot
- Alignment Depth : Mosdepth
- MultiQC
  
---

I've absolutely no doubt that there should be some problems :). It runs on my end, but perhaps not yours. If that is the case, please email to [Tuan Nguyen](mailto:tuan.nguyen@agriculture.vic.gov.au)
