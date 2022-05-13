#--------------------------------------------------
# In diesem Dukoment versuchen wir die Gereinigten expressionsdaten mittels gsea zu
#pathway activitäten zusammenzufassen un din einer matrix mit patienten als spalten und 
# enrichment scores als zeilen zusammenzufassen
#--------------------------------------------------

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("gage")
BiocManager::install("fgsea")
library(dplyr)
library(gage)
library(fgsea)
library(pheatmap)

load('~/GitHub/2022-topic-02-team-05/data/tcga_exp_cleaned.RData')
load('~/GitHub/2022-topic-02-team-05/data/our_genesets.RData')
genesets = readRDS("~/GitHub/2022-topic-02-team-05/data/hallmarks_genesets.rds")

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
#versuch eine gsea zu implementieren
#------------------------

#function die eine gsea für einen beliebigen patienten durchführt und nur den NES benannt ausgibt
#alle pathways die nocht siginifikant von der base level abweichen werden 0 gesetzt

GSEA_NES = function(patient, pathways = c(genesets[[1]], our_genesets), pvalue = 0.05){
  
  #sorting and naming the exp data
  names(patient) = tcga_genenames
  #hiervor z-transforamtion einfügen
  patient = sort(patient, decreasing = TRUE) 
  
  #duarchführen der GSEA
  res = fgseaMultilevel(pathways = pathways, 
                        stats = patient,
                        minSize=3)
  #extrahieren der NES werte & nullsetzten aller nichtsiginifikanten
  ret = res$NES; names(ret) = res$pathway
  #ret[res$padj > pvalue] = 0
  ret[is.na(ret)] = 0
  
  message('I´m still standing')
  return(ret)
}

#-------------------------------------------
#Durchführung der GSEA
#-----------------------------------------------

#pathways = gmtPathways("C:/Users/jakob/Desktop/c2.cp.v7.5.1.symbols.gmt") alle canonischen pathways msigr vllt als alternative
GSEA_matrix = apply(tcga_exp_cleaned[,1:20], 2, GSEA_NES) %>% as.data.frame()

#darstellen der GSEA matrix als heatmap
#rote sind überexpremierte pathways, blue unterexpremierte, weiß keine abweichung
pheatmap(as.matrix(GSEA_matrix),
         breaks = seq(-max(GSEA_matrix), max(GSEA_matrix), length.out = 201),
         color = colorRampPalette(c('blue','light blue', 'white','yellow', 'red'),
                                  bias = 1,
                                  space = 'rgb',
                                  interpolate = 'linear'
         )(200),
         clustering_method = 'average',
         treeheight_row = 25, treeheight_col = 20, cellwidth = 20,cellheight = 20,
         show_colnames = FALSE,
         legend_breaks = c(-max(GSEA_matrix),0, max(GSEA_matrix)),
         legend_labels = c('underexpressed', 'normal expression', 'overexpressed')
        )


save(GSEA_matrix, file = '~/GitHub/2022-topic-02-team-05/data/GSEA_matrix.RData')



  

