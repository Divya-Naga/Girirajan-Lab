# Week-5
## Intersection of between the significantly DE genes from SEDM vs SEEM analysis and A687 (24hr) vs Nalm6 (0hr) analysis
<img src="https://github.com/Divya-Naga/Girirajan-Lab/blob/main/week-5/results/intersection_significant_DE_genes.png">
The venn diagram shows that 6 genes are shared by both the analysis of genes with significant Differential Expression.

Here, the argument "transcript_count" accepts the transcript count matrix and "in_cols" accepts the columns for subsetting. "res_filename" and "normalized_filename" accepts the location for the storing the differentially expressed transcript.

```
gtf_file_input = as.character(args[1])
SEDM_vs_SEEM_input = as.character(args[2])
A687_vs_Nalm6_input = as.character(args[3])
Venn_Out = args[4]
tmp_dir = args[5]
```

Here, the argument "gtf_file_input" accepts the zip file containing the gtf file and "SEDM_vs_SEEM_input" accepts the SEDM vs SEEM analysis and "A687_vs_Nalm6_input" accepts the A687 (24hr) vs Nalm6 (0hr) Analysis. "Venn_Out" file location for storing the venn diagram.
```
Rscript ./week-5/src/gtf_DESeq.R ./week-5/data/merged.annotated.gtf.gz ./week-2/results/condition_infected_vs_control_dge.csv ./week-4/results/condition_infected_vs_control_dge.csv ./week-5/results/intersection_significant_DE_genes.svg ./week-5/data
```
