#challange von Wang
#Daten laden:
load('~/GitHub/2022-topic-02-team-05/data/tcga_exp_cleaned.RData')
tcga_anno = readRDS("~/GitHub/2022-topic-02-team-05/data/tcga_tumor_annotation.RDS")
load('data/our_genesets.RData')
load('data/geneset_ids.RData')

#umbennen der rowanems zu reinen ensembl ids
tcga_genes_cleaned = rownames(tcga_exp_cleaned)
#dieser vetor enth?lt sowohl enseblm id als auch genenamen und muss daher gespalten werden
tcga_genes_cleaned = strsplit(tcga_genes_cleaned, split = '|', fixed = TRUE)
#speicher der ensembl ids ohne versionsnummer als eigenen vektor
tcga_geneids = sapply(tcga_genes_cleaned, function(tcga_genes_cleaned){return(tcga_genes_cleaned[1])})
tcga_geneids = strsplit(tcga_geneids, split = '.', fixed = TRUE)
tcga_geneids = sapply(tcga_geneids, function(tcga_geneids){return(tcga_geneids[1])})
rownames(tcga_exp_cleaned) = tcga_geneids


#createn einer liste mit allen patienten in dfs sortiert nach krebs
cancers = list();cancers = vector('list',length(table(tcga_anno$cancer_type_abbreviation)))
names(cancers) = names(table(tcga_anno$cancer_type_abbreviation))
i=1
for (i in 1:length(cancers)){
  cancers[[i]] = tcga_exp_cleaned[,tcga_anno$cancer_type_abbreviation == names(cancers)[i]]
}

#function die einen krebstypen df und genesets als input nimmt und ein df mit pvalues ausgibt
enrichment = function(expressiondata, genesets = genesets_ids){
  ESmatrix = sapply(genesets, FUN = function(x){
    ins = na.omit(match(x,rownames(expressiondata)))#indices der gene im aktuellen set
    outs = -ins#indices der gene nicht im aktuellen set
    
    #gibt einen vektor der für jeden patienten den pval für das aktuelle gene enthält
    res = NULL
    for (i in 1:ncol(expressiondata)){#testet für jeden patienten
      res[i] = wilcox.test(expressiondata[ins,i],expressiondata[outs,i],'two.sided')$p.value
    }
    return(res)
  })
  row.names(ESmatrix) = colnames(expressiondata); return(ESmatrix)
}

pvalueslist = lapply(cancers, enrichment)#für die tests für jeden krebstypen durch


      
#test mit dummydaten
# set.seed(100)
# 
# # create example data
# 
# df1 <- data.frame(matrix(rnorm(50*50, 0, 2), 50, 50), row.names = 1:50)
# df2 <- data.frame(matrix(rnorm(50*50, 0, 1.5), 50, 50), row.names = 1:50)
# 
# df.list = list("LUAD" = df1, "BRCA" = df2)
# 
# pathway1 <- c(1, 4, 6, 34, 27, 14)
# pathway2 <- c(5, 8, 9, 19, 5, 19, 20, 50)
# 
# pw.list = list("pathway1" = pathway1, "pathway2"=pathway2)
# pvalueslist = lapply(df.list, FUN = function(x){return(enrichment(x,pw.list))})
