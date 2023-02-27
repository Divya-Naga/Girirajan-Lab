# Week-1
## Changed the code from enrich.R to involve addiotional database support from msigdatabase

Molecular Signatures Database contains 8 major collections of gene sets:

H: hallmark gene sets <br> 
C1: positional gene sets <br> 
C2: curated gene sets <br> 
C3: motif gene sets <br> 
C4: computational gene sets <br> 
C5: GO gene sets <br> 
C6: oncogenic signatures <br> 
C7: immunologic signatures <br> 

Here, "Category" accepts an argument which specifies the type of collection.
``` gene_file = args[1]
gseaout_file = args[2]
category = as.character(args[3])
tmp_dir = args[4] 
```
To retrieve all human gene sets for the specified category
```
#Term to gene for specified category of gene set 
m_t2g <- msigdbr(species = "Homo sapiens", category = category) %>%
  dplyr::select(gs_name, entrez_gene)
```
Example of an R script input:
```
Divyas-MacBook-Pro:Girirajan-Lab DivyaNaga$ Rscript ./week-1/src/enrich.R ./week-1/data/trtmntvsctrl_sig_genes.txt ./week-1/data/gseaout.csv C5 week-1/data/
```
