# Analysis of publicly available bulk RNA sequencing datasets

### Folder Structure:
* `src`: contains all scripts necessary to generate figures and tables for bulk rna-seq analysis
* `annotations`: contains dataset annotation files required to execute scripts
* `raw_datasets`: contains raw datasets downloaded from GREIN or from GEO entries (some sample names were modified slightly, like inclusion of extra underscores, for compatibility with pipeline)
* `deseq_files`: folder where DESeq2 files from `runDEseq.py` will be deposited
* `paralog_data`: contains paralog data for each KO target in each dataset (downloaded Aug 11 2023). Can be regenerated using `findParalogs.py`
* `de_analysis`: folder where all plots and tables of this analysis will be deposited after running `main.py`


### Setup Instructions:
* download [author-provided DESeq2 dataset for GSE151825](https://drive.google.com/drive/folders/1Oz9S5Ydvu5K6QPuBMKXlM0BPKcfvaT0N) and deposit in `deseq_files` folder
* setup environment using using `nitc-env.yaml` (example `conda env create -f nitc-env.yml`)
* navigate to `rnaseq` folder before running scripts to prevent issues with local paths


### Reproducing Data:
1. Run `src/deDeseq.py` to generate all DESeq files and populate the `deseq_files` folder
2. Run `src/main.py` to generate plots and tables for each dataset 


`DESeq_summary_plots.R` in the parent directory contains code to generate all other plots related to this data
