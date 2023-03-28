library("biomaRt")
library("dplyr")
library("DESeq2")

args <- commandArgs(trailingOnly=TRUE)
# check to see that one argument is given
if (length(args)!=5) {
  stop("Please provide a valid transcript_count, in_cols, res_filename, normalized_filename and tmp_dir", call.=FALSE)
}

transcript_count = args[1]
in_cols = unlist(strsplit(args[2], ","))
res_filename = args[3]
normalized_filename = args[4]
tmp_dir = args[5]

transcript_count_data <- read.csv(transcript_count, row.names = "transcript_id")

Sys.setenv(BIOMART_CACHE = tmp_dir)
print(biomartCacheInfo())

#Count Matrix with the desired columns
transcript_count_matrix <- as.matrix(transcript_count_data[, c("SEEM.WT.CTL1A", "SEEM.WT.CTL2A",
                                                               "SEEM.WT.CTL3A", "SEDM.WT.CTL1A",
                                                              "SEDM.WT.CTL2A", "SEDM.WT.CTL3A")])

# Sample information 
coldata <- data.frame(
  sample = c( "SEEM.WT.CTL1A", "SEEM.WT.CTL2A","SEEM.WT.CTL3A", 
              "SEDM.WT.CTL1A", "SEDM.WT.CTL2A", "SEDM.WT.CTL3A"),
  condition = c( "control", "control",  "control",
                 "treated","treated","treated"), 
  row.names = "sample" )
coldata$condition <- as.factor(coldata$condition)

# DESeqDataSet for Differential gene expression analysis 
dds <- DESeqDataSetFromMatrix(countData = transcript_count_matrix, colData = coldata, 
                              design = ~ condition)
dds <- dds[rowSums(counts(dds)) >= 10,] #prefilters to remove genes with less than 10 reads

# Control is set as reference
dds$condition <- relevel(dds$condition, ref = "control")

dds <- DESeq(dds) # performs Differential Gene Expression Analysis

# Result 
res <- results(dds,alpha=0.05)  
res <- res[order(res$padj),] # ordered w.r.t p values
#Output result file
write.csv(as.data.frame(res), file=res_filename)

#Normalized counts
normalized_counts <- counts(dds, normalized=TRUE)
write.csv(as.data.frame(normalized_counts), file= normalized_filename)

