#-----------------------------------------------------
#Umschreiben der Hallmarkpathways von Gennamen in Ensemble IDs
#-----------------------------------------------------

#BiocManager::install("biomaRt")
library(biomaRt)
genesets = readRDS("data/hallmarks_genesets.rds")

mart = useEnsembl(dataset = "hsapiens_gene_ensembl", biomart='ensembl')
genesets_ids = lapply(genesets[[1]], FUN = function(x){getBM(attributes = "ensembl_gene_id",
                                                             filters = "external_gene_name",
                                                             values = x,
                                                             mart = mart )})
genesets_ids = sapply(genesets_ids, FUN = function(genesets_ids){return(as.vector(genesets_ids))})
names(genesets_ids) = names(genesets[[1]])

save(genesets_ids, file = 'data/geneset_ids.RData')