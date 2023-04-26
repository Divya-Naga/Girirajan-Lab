library("biomaRt")
library("clusterProfiler")
library("enrichplot")
library("msigdbr")
library("dplyr")
library("DESeq2")
args <- commandArgs(trailingOnly=TRUE)
# check to see that one argument is given
if (length(args)!=4) {
  stop("Please provide a valid gene_file1, gene_file2, gseaout_file, category,dot_Out and tmp_dir", call.=FALSE)
}

gene_file1 = args[1]
gene_file2 = args[2]
gseaout_file = args[3]
tmp_dir = args[4]

gene_table1 = read.csv(gene_file1)
gene_table2 = read.csv(gene_file2)

Sys.setenv(BIOMART_CACHE = tmp_dir)
print(biomartCacheInfo())

#Overlapping SDE genes between the two files (padj <0.05)
gene_table1 <- na.omit(gene_table1[gene_table1$padj < 0.05,])
gene_table2 <- na.omit(gene_table2[gene_table2$padj < 0.05,])

Overlapping_SDE_gene <- as.data.frame(intersect(row.names(gene_table1), row.names(gene_table2)))
colnames(Overlapping_SDE_gene) <- c("ensembl_names")
Overlapping_SDE_gene$`ensembl_names`<-gsub("\\..*","",as.character(Overlapping_SDE_gene$`ensembl_names`)) 

# get the genes and store as a vector
genes <- as.data.frame(Overlapping_SDE_gene$`ensembl_names`)

# convert ensemble id to entrez id :: link: https://support.bioconductor.org/p/114325/
mart <- useMart("ensembl","hsapiens_gene_ensembl")
entrez_genes <- getBM(c("ensembl_gene_id", "entrezgene_id"), "ensembl_gene_id", genes, mart)

msig_list <- list("H", "C1", "C2", "C3", "C4", "C5", "C6","C7")

for (x in msig_list) {
#Term to gene for specified category of gene set 
m_t2g <- msigdbr(species = "Homo sapiens", category = x) %>%
  dplyr::select(gs_name, entrez_gene)

#Over-presentation analysis
em <- enricher(gene = entrez_genes[, 2],
               TERM2GENE = m_t2g)


  
# save to file .. 
  write.table(em, file = paste0(gseaout_file,"GSEAout_",x,".csv"), sep=",", row.names=TRUE, col.names=TRUE)
}

# use enrichGO for go term analysis :: link: http://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-go.html
goenrich = enrichGO(
  gene=entrez_genes[, 2],
  OrgDb='org.Hs.eg.db',
  pAdjustMethod="BH",
  pvalueCutoff=0.05,
  ont="BP"
)

pdf(paste0(gseaout_file,"Dot_Plot.pdf"))
myplot <- dotplot(goenrich, showCategory=20, font.size=8) + ggtitle("dotplot for GSEA")
print(myplot)
dev.off()




