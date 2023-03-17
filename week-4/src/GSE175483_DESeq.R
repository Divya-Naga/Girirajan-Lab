library("biomaRt")
library("dplyr")
library("DESeq2")

args <- commandArgs(trailingOnly=TRUE)
# check to see that one argument is given
if (length(args)!=10) {
  stop("Please provide a valid individual_count_files_control, individual_count_files_replicate, meta_counts_filename, res_filename, normalized_filename and tmp_dir", call.=FALSE)
}

individual_count_files_control = as.character(c(args[1],args[2],args[3]))
individual_count_files_replicate = as.character(c(args[4],args[5],args[6]))
meta_counts_filename = args[7]
res_filename = args[8]
normalized_filename = args[9]
tmp_dir = args[10]

#Extracting the treated sample replicates 
for (i in 1:length(individual_count_files_control)){
  metadata <- read.table(gzfile(individual_count_files_control[i]), header= TRUE) 
  n<-dim(metadata)[1]
  metadata <- metadata[1:(n-5),]
  assign(paste0("Nalm6_0hr_Pred_Rep",i),metadata)
}

#Extracting the treated sample replicates 
for (i in 1:length(individual_count_files_replicate)){
  metadata <- read.table(gzfile(individual_count_files_replicate[i]), header= TRUE) 
  n<-dim(metadata)[1]
  metadata <- metadata[1:(n-5),]
  assign(paste0("A687_24hr_Pred_Rep",i),metadata)
}

# Transcript count matrix for meta counts
transcript_count_matrix_NEW <- cbind(Nalm6_0hr_Pred_Rep1 = Nalm6_0hr_Pred_Rep1[,2],
                                Nalm6_0hr_Pred_Rep2 = Nalm6_0hr_Pred_Rep2[,2],
                                Nalm6_0hr_Pred_Rep3 = Nalm6_0hr_Pred_Rep3[,2],
                                A687_24hr_Pred_Rep1 = A687_24hr_Pred_Rep1[,2],
                                A687_24hr_Pred_Rep2 = A687_24hr_Pred_Rep2[,2],
                                A687_24hr_Pred_Rep3 = A687_24hr_Pred_Rep3[,2]) 
row.names(transcript_count_matrix_NEW) <- Nalm6_0hr_Pred_Rep1[,1]
transcript_count_matrix_M <- as.matrix(transcript_count_matrix_NEW)

#Output file for Meta counts matrix
write.csv(as.data.frame(transcript_count_matrix_M), file=meta_counts_filename)

# Sample information 
coldata <- data.frame(
  sample = c( "Nalm6_0hr_Pred_Rep1", "Nalm6_0hr_Pred_Rep2","Nalm6_0hr_Pred_Rep3", 
              "A687_24hr_Pred_Rep1", "A687_24hr_Pred_Rep2", "A687_24hr_Pred_Rep3"),
  condition = c( "control", "control",  "control",
                 "treated","treated","treated"), 
  row.names = "sample" )
coldata$condition <- as.factor(coldata$condition)

# DESeqDataSet for Differential gene expression analysis 
dds <- DESeqDataSetFromMatrix(countData = transcript_count_matrix_M, colData = coldata, 
                              design = ~ condition)
dds <- dds[rowSums(counts(dds)) >= 10,] #prefilters to remove genes with less than 10 reads

# Control is set as reference
dds$condition <- relevel(dds$condition, ref = "control")

dds <- DESeq(dds) # performs Differential Gene Expression Analysis

# Result 
res <- results(dds, alpha = 0.05) # cut-off p-value < 0.05  
res <- res[order(res$padj),] # ordered w.r.t p values

#Output result file
write.csv(as.data.frame(res), file=res_filename)

#Normalized counts
normalized_counts <- counts(dds, normalized=TRUE)
write.csv(as.data.frame(normalized_counts), file= normalized_filename)
