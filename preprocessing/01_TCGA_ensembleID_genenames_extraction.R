#--------------------------------
#Extraktion aller Ensemble ids und Gennamen aus den Expressionsdaten des TCGA
#--------------------------------
tcga_exp = readRDS("data/tcga_tumor_log2TPM.RDS")

#extrahieren aller Gene, die in den Expressionsdaten vorkommen
tcga_genes = rownames(tcga_exp)

#dieser Vektor enth?lt sowohl Ensemble ids, als auch Gennamen und muss daher gespalten werden
tcga_genes = strsplit(tcga_genes, split = '|', fixed = TRUE)

#speichern der Ensemble ids ohne Versionsnummer als eigenen Vektor
tcga_geneids = sapply(tcga_genes, function(tcga_genes){return(tcga_genes[1])})
tcga_geneids = strsplit(tcga_geneids, split = '.', fixed = TRUE)
tcga_geneids = sapply(tcga_geneids, function(tcga_geneids){return(tcga_geneids[1])})

#speichern der Gennamen ohne Versionsnummer als eigenen Vektor
tcga_genenames = sapply(tcga_genes, function(tcga_genes){return(tcga_genes[2])})
tcga_genenames = strsplit(tcga_genenames, split = '.', fixed = TRUE)
tcga_genenames = sapply(tcga_genenames, function(tcga_genenames){return(tcga_genenames[1])})

tcga_genes = cbind.data.frame(tcga_geneids,tcga_genenames)

#speichern aller Ensemble ids und Gennamen in einem Dataframe
save(tcga_genes, file = 'data/tcga_genes.RData')

    


