library("biomaRt")
library("clusterProfiler")
library("enrichplot")


args = commandArgs(trailingOnly=TRUE)
# check to see that one argument is given
if (length(args)!=3) {
  stop("gene_list, gsea table, tmp_dir must be given", call.=FALSE)
}

gene_file = args[1]
gseaout_file = args[2]
tmp_dir = args[6]

# get the genes and store as a vector
genes = read.table(gene_file, header=FALSE)

Sys.setenv(BIOMART_CACHE = tmp_dir)
print(biomartCacheInfo())

# convert ensemble id to entrez id :: link: https://support.bioconductor.org/p/114325/
mart <- useMart("ensembl","hsapiens_gene_ensembl")
entrez_genes <- getBM(c("ensembl_gene_id", "entrezgene_id"), "ensembl_gene_id", genes, mart)

# use enrichGO for go term analysis :: link: http://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-go.html
goenrich = enrichGO(
    gene=entrez_genes[, 2],
    OrgDb='org.Hs.eg.db',
    pAdjustMethod="BH",
    pvalueCutoff=0.05,
    ont="BP"
)

# save to file .. 
write.table(goenrich, file=gseaout_file, sep=",", row.names=TRUE, col.names=TRUE)
