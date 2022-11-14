# AgVic_CLRC 

Cattle SVs & SNPs Discovery for Bovine Long-Read Consortium

Written in Nextflow DSL2

Pipeline developed for usage in the Bovine Long-Read Consortium (BovLRC). The pipeline deployed multiple bioinformatics software for the detection of Single Nucldeotide Polymorphism (SNPs) & Structural Variants (SV). The pipeline (version 0.0.1) currently deployed.

1 - QC :
[LongFilt]() : QC for both LongReads and ShortReads (**DEFAULT**)
[NanoFilt]() + [FMLRC2]() : NanoFilt for QC of Long-Read samples, and FMLRC2 + NanoFilt for QC of Short-Read samples .
2 - Mapping:
[Minimap2]() : (**DEFAULT**)
[Winnowmap2]()
[NGMLR](https://github.com/philres/ngmlr)
3 - SNP Caller:
[Clair3]() 
[PEPPER]()
4 - SV Caller:
[Sniffles]()
[DYSGU]()
[SVCute]()

If you have any queries, please email to [Tuan Nguyen](mailto:tuan.nguyen@agriculture.vic.gov.au)