# Week 2

``` transcript_count = args[1]
in_cols = unlist(strsplit(args[2], ","))
res_filename = args[3]
normalized_filename = args[4]
tmp_dir = args[5] 
```
Here, the argument "transcript_count" accepts the transcript count matrix and "in_cols" accepts the columns for subsetting. "res_filename" and "normalized_filename" accepts the location for the storing the differentially expressed transcript.

## Design Matrix .
``` 
coldata <- data.frame(
  sample = c( "SEEM.WT.CTL1A", "SEEM.WT.CTL2A","SEEM.WT.CTL3A", 
              "SEDM.WT.CTL1A", "SEDM.WT.CTL2A", "SEDM.WT.CTL3A"),
  condition = c( "control", "control",  "control",
                 "treated","treated","treated"), 
  row.names = "sample" )
coldata$condition <- as.factor(coldata$condition)
``` 
This creates the meta data for the DGE analysis. This data frame assigns SEEM columns as control and SEDM columns as treated samples.

## Differential Gene Expression Analysis
``` 
dds <- DESeqDataSetFromMatrix(countData = transcript_count_matrix, colData = coldata, design = ~ condition)
dds <- dds[rowSums(counts(dds)) >= 10,] #prefilters to remove genes with less than 10 reads
# Control is set as reference
dds$condition <- relevel(dds$condition, ref = "control")
dds <- DESeq(dds)
```
Here, dds is the DESeqDataSet which is constructed for DGE analysis.

## Example Rscript
``` 
Rscript ./week-2/src/DESeq.R ./week-2/data/transcript_count_matrix.csv "SEEM.WT.CTL1A,SEEM.WT.CTL2A,SEEM.WT.CTL3A,SEDM.WT.CTL1A,SEDM.WT.CTL2A,SEDM.WT.CTL3A" ./week-2/results/condition_infected_vs_control_dge.csv ./week-2/results/normalized_counts_dge.csv ./week-2/data
``` 
