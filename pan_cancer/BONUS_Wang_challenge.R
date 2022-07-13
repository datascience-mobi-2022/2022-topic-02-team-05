#--------------------------
#challange von Wang: Erstellen einer Funktion, die
#einen Krebstypen dataframe und genesets als input nimmt und einen df mit p-Werten ausgibt
#--------------------------

#Daten laden:
load('data/tcga_exp_cleaned.RData')
tcga_anno = readRDS("data/tcga_tumor_annotation.RDS")
load('data/our_genesets.RData')
load('data/geneset_ids.RData')

#umbennen der rownames zu reinen EnsembleIDs
tcga_genes_cleaned = rownames(tcga_exp_cleaned)

#dieser Vektor enthält sowohl EnsembleIDs, als auch Gennamen und muss daher gespalten werden
tcga_genes_cleaned = strsplit(tcga_genes_cleaned, split = '|', fixed = TRUE)

#speichern der EnsembleIDs ohne Versionsnummer als eigenen Vektor
tcga_geneids = sapply(tcga_genes_cleaned, function(tcga_genes_cleaned){return(tcga_genes_cleaned[1])})
tcga_geneids = strsplit(tcga_geneids, split = '.', fixed = TRUE)
tcga_geneids = sapply(tcga_geneids, function(tcga_geneids){return(tcga_geneids[1])})
rownames(tcga_exp_cleaned) = tcga_geneids

#createn einer Liste mit allen Patienten in dataframes sortiert nach krebs
cancers = list();cancers = vector('list',length(table(tcga_anno$cancer_type_abbreviation)))
names(cancers) = names(table(tcga_anno$cancer_type_abbreviation))
i=1
for (i in 1:length(cancers)){
  cancers[[i]] = tcga_exp_cleaned[,tcga_anno$cancer_type_abbreviation == names(cancers)[i]]
}

#erstellen der finalen Funktion
enrichment = function(expressiondata, genesets = genesets_ids){
  ESmatrix = sapply(genesets, FUN = function(x){
    ins = na.omit(match(x,rownames(expressiondata)))#indices der gene im aktuellen set
    outs = -ins#indices der gene nicht im aktuellen set
    
    #gibt einen vektor aus, der für jeden Patienten den p-Wert für das aktuelle Gen enthält
    res = NULL
    for (i in 1:ncol(expressiondata)){#testet für jeden patienten
      res[i] = wilcox.test(expressiondata[ins,i],expressiondata[outs,i],'two.sided')$p.value
    }
    return(res)
  })
  row.names(ESmatrix) = colnames(expressiondata); return(ESmatrix)
}

pvalueslist = lapply(cancers, enrichment)#für die tests für jeden krebstypen durch

