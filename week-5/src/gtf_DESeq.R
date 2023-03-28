library(ggvenn)
library("biomaRt")
library("dplyr")
library("DESeq2")

args <- commandArgs(trailingOnly=TRUE)
# check to see that one argument is given
if (length(args)!=5) {
  stop("Please provide a valid gtf_file_input as a gz file, SEDM_vs_SEEM_input, A687_vs_Nalm6_input, Venn_Out and tmp_dir", call.=FALSE)
}

gtf_file_input = as.character(args[1])
SEDM_vs_SEEM_input = as.character(args[2])
A687_vs_Nalm6_input = as.character(args[3])
Venn_Out = args[4]
tmp_dir = args[5]

# gtf file gene mapping
gtf <- rtracklayer::import(gzfile(gtf_file_input))
gtf_type <- data.frame(type = as.character(gtf@elementMetadata@listData[["type"]]),
                       Transcript_id = gtf@elementMetadata@listData[["transcript_id"]],
                       gene_name = gtf@elementMetadata@listData[["gene_name"]])
gtf_type <- dplyr::filter(gtf_type, type == "transcript")

# gene mapping FOR SEDM vs SEEM analysis
SEDM_vs_SEEM_res <- read.csv(SEDM_vs_SEEM_input)
gene_mapping_week2 <- data.frame(Transcript_id = SEDM_vs_SEEM_res[["X"]])
gene_mapping_SEDM_vs_SEEM <- merge(gene_mapping_week2, gtf_type, by = "Transcript_id")
gene_mapping_SEDM_vs_SEEM <- gene_mapping_SEDM_vs_SEEM[c("Transcript_id","gene_name")]

#gene names A687_24hr vs Nalm6_0hr
A687_24hr_vs_Nalm6_0hr_res <- read.csv(A687_vs_Nalm6_input)
gene_names_A687_24hr_vs_Nalm6_0hr <- data.frame(gene_name = A687_24hr_vs_Nalm6_0hr_res[["X"]])

#Venn Diagram
png(filename=Venn_Out) #output graph

  X = list(
    "SEDM vs SEEM" = gene_mapping_SEDM_vs_SEEM[["gene_name"]], 
    "A687_24hr vs Nalm6_0hr" = gene_names_A687_24hr_vs_Nalm6_0hr[["gene_name"]])
ggvenn(
  X, 
  fill_color = c("#0073C2FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
)
dev.off()



