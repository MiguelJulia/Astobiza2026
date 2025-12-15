# Astobiza2026
Scripts for analysing RNAseq data, scRNAseq data, and plotting figures for the paper:

Fibroblastic aspartoacylase suppresses TGFβ-mediated responses and cancer progression, Astobiza et al 2026. 

The repository is organised in the following manner:
```
.
├── RNAseq
│   ├── 01_RNAseq_pipeline.sh
│   ├── postanalysis
│   │   └── 02_Ranalysis_exploratory_analysis_aspa.Rmd
│   └── sample_sheet_2.csv
└── scRNAseq
    ├── BC and BLC
    │   ├── 01_Breast_and_Bladder_Signatures_and_Genes.Rmd
    │   ├── BC preprocessing
    │   │   └── BC_Analysis_FINAL.Rmd
    │   └── BLC preprocessing
    │       ├── 01_scRNAseq_pipeline.sh
    │       └── postanalysis
    │           ├── 00_Ranalisys.R
    │           └── 01_Ranalisys_allsamples.Rmd
    └── PC
        ├── 00_Ranalysis.Rmd
        └── markers
            ├── Bladder.fibroblast.table.rethinkingclustering.csv
            ├── Bladder.table.rethinkingclustering.csv
            ├── markers_annotation.tsv
            ├── markers_fer_classification.tsv
            ├── markers_fibroblast_clusters_onlyCluster10.tsv
            ├── markers_fibroblast.tsv
            └── markers_mycelltype.tsv
```
# RNAseq
- Here you have the code and sample sheed to replicate the sample analysis and postanalysis of our RNAseq samples.
- Inside the subfolder *postanalysis* you have the script to replicate figures.

# scRNAseq
- Here you have the code and sample sheed to replicate the sample analysis and postanalysis of the scRNAseq used in the paper. Our own dataset is inside *PC* and public datasets under *BC and BLC* and inside their corresponding subfolders.
- Inside the subfolder *postanalysis* you have the script to replicate figures.

