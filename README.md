## Environment

Analysis was performed in a Conda environment. See **NSCLCL.yml** for details.

## Experimental conditions

### RNA-Seq data

RNA-Seq with 90 samples:
  4 cell lines: CUTO8, CUTO9, CUTO29 and YU1077
3 inhibitors plus control: alectinib, brigatinib, lorlatinib (2 replicates each)
2 time points: 6h, 24h
(YU1077 wasn't treated with alectinib)

## Data
Raw RNA-Seq data have been deposited in Arrayexpress: E-MTAB-11342

* Processed data are available in data/NSCLCL_data.RData. This file containes the following objects:
* normalized_counts: DESeq2-normalized counts from cell line RNA-Seq data
* res_diff_expr: DESeq2 output from cell line RNA-Seq data
* sample_info: sample information cell line RNA-Seq data
* geneset_ls: genesets downloaded from MSigDB v7.2

## Manuscript
The main analysis, as reported in the manuscript

### RNA-Seq 
```{r}
source("scripts/manuscript.R")
```

