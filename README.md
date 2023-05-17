# Project VALOR (**V**ariant **A**nalysis using **LO**ng **R**ead sequencing) 

*Currently deployed for Cattle SVs & SNPs Discovery in the Bovine Long-Read Consortium (BovLRC)
*
[Tuan Nguyen](tuan.nguyen@agriculture.vic.gov.au)
##

Initial setup: 
1. Clone this Github
```
git clone https://github.com/tuannguyen8390/AgVic_CLRC.git
```

Pipeline developed for usage in the Bovine Long-Read Consortium (BovLRC). The pipeline deployed multiple bioinformatics software for the detection of Single Nucldeotide Polymorphism (SNPs) & Structural Variants (SV). The pipeline (version 0.0.2) currently deployed. Written in Nextflow DSL2.


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

Please refer to these files for editing your own. You can run with your own files deploying `--LR_MetaDir` AND/OR `--SR_MetaDir`

---

## Pipeline overview

1. QC :
[FiltLong](https://github.com/rrwick/Filtlong) : QC for both LongReads and ShortReads (**DEFAULT**)

[NanoFilt](https://github.com/wdecoster/nanofilt) + [FMLRC2](https://github.com/HudsonAlpha/fmlrc2) : NanoFilt for QC of Long-Read samples, and FMLRC2 + NanoFilt for QC of Short-Read samples .

2. Mapping:
[Minimap2](https://github.com/lh3/minimap2) : (**DEFAULT**)

[Winnowmap2](https://github.com/marbl/Winnowmap)

[NGMLR](https://github.com/philres/ngmlr)

3. SNP Caller: All callers are run in parallel
[Clair3](https://github.com/HKU-BAL/Clair3) : (**RECOMMEND FOR DOWNSTREAM ANALYSIS**)

[PEPPER](https://github.com/kishwarshafin/pepper) - By default, Flowcell < 10.4 will be analyzed with PEPPER

[DEEPVARIANT](https://github.com/google/deepvariant) - By default, Flowcell >= 10.4 will be analyzed with DEEPVARIANT & HIFI (**RECOMMEND FOR DOWNSTREAM ANALYSIS**)

[Longshot](https://github.com/pjedge/longshot) 

4. SV Caller: All callers are run in parallel
[Sniffles2](https://github.com/fritzsedlazeck/Sniffles) (**RECOMMEND FOR DOWNSTREAM ANALYSIS**)

[DYSGU](https://github.com/kcleal/dysgu) (**RECOMMEND FOR DOWNSTREAM ANALYSIS**)

[CuteSV2](https://github.com/tjiangHIT/cuteSV) (**RECOMMEND FOR DOWNSTREAM ANALYSIS**)

---

If you have any queries, please email to [Tuan Nguyen](mailto:tuan.nguyen@agriculture.vic.gov.au)
