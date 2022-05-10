if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("gage")
BiocManager::install("fgsea")
library(dplyr)
library(gage)
library(fgsea)

load('~/GitHub/2022-topic-02-team-05/data/tcga_exp_cleaned.RData')
load('~/GitHub/2022-topic-02-team-05/data/our_genesets.RData')

#----------------------------------------------------------
#bennen der gene in tcga exp nur mit namen ohne ensembl id
#--------------------------------------------------------
#gennamen extrahieren
tcga_genes = rownames(tcga_exp_cleaned)

#dieser vetor enth?lt sowohl enseblm id als auch genenamen und muss daher gespalten werden
tcga_genes = strsplit(tcga_genes, split = '|', fixed = TRUE)

#speicher der gennamen als eigenen vektor
tcga_genenames = sapply(tcga_genes, function(tcga_genes){return(tcga_genes[2])})

#entfernen der Versionsnummern f?r ensembl id und gennamen
tcga_genenames = strsplit(tcga_genenames, split = '.', fixed = TRUE)
tcga_genenames = sapply(tcga_genenames, function(tcga_genenames){return(tcga_genenames[1])})

#----------------------
#versuch eine gsea für einen patienten zu implementieren
#------------------------
probe_exp = tcga_exp_cleaned[,1]
names(probe_exp) = tcga_genenames
probe_exp = sort(probe_exp, decreasing = TRUE)

pval = 0.05
pathways = gmtPathways("C:/Users/jakob/Desktop/c2.cp.v7.5.1.symbols.gmt")

probe_gsea_3 <- fgseaMultilevel(pathways = our_genesets, 
                                stats = probe_exp,
                                minSize=3,
                                maxSize=600) #%>% 
  
  #dplyr::filter(padj < !!pval)











tcga_gsea = gseGO(geneList=our_genesets_sorted[[1]], 
      ont ="ALL", 
      keyType = "SYMBOL", 
      nPerm = 100, #the higher the more accurrate but longer
      minGSSize = 3, 
      maxGSSize = 800, 
      pvalueCutoff = 0.05, 
      verbose = TRUE, 
      OrgDb = "org.Hs.eg.db", 
      pAdjustMethod = "none")
