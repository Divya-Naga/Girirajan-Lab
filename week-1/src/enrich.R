library("biomaRt")
library("clusterProfiler")
library("enrichplot")
library("msigdbr")
library("dplyr")

args = commandArgs(trailingOnly=TRUE)
# check to see that one argument is given
if (length(args)!=4) {
  stop("gene_list, gsea table, tmp_dir must be given", call.=FALSE)
}

gene_file = args[1]
gseaout_file = args[2]
category = as.character(args[3])
tmp_dir = args[4]

# get the genes and store as a vector
genes = read.table(gene_file, header=FALSE)

Sys.setenv(BIOMART_CACHE = tmp_dir)
print(biomartCacheInfo())

# convert ensemble id to entrez id :: link: https://support.bioconductor.org/p/114325/
mart <- useMart("ensembl","hsapiens_gene_ensembl")
entrez_genes <- getBM(c("ensembl_gene_id", "entrezgene_id"), "ensembl_gene_id", genes, mart)

#Term to gene for specified category of gene set 
m_t2g <- msigdbr(species = "Homo sapiens", category = category) %>%
  dplyr::select(gs_name, entrez_gene)

#Over-presentation analysis
em <- enricher(gene = entrez_genes[, 2],
               TERM2GENE = m_t2g)

# save to file .. 
write.table(em, file=gseaout_file, sep=",", row.names=TRUE, col.names=TRUE)
