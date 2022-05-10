if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("gage")
BiocManager::install("fgsea")
library(dplyr)
library(gage)
library(fgsea)

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

#pathways = gmtPathways("C:/Users/jakob/Desktop/c2.cp.v7.5.1.symbols.gmt") alle canonischen pathways msigr vllt als alternative


#function die eine gsea für einen beliebigen patienten durchführt und den NES ausgibt
GSEA = function(patient){
  pathways = c(genesets[[1]], our_genesets)
  names(patient) = tcga_genenames
  patient = sort(patient, decreasing = TRUE)
  
  res = fgseaMultilevel(pathways = pathways, 
                        stats = patient,
                        minSize=3
                        )
  message('I´m still standing')
  return(res)
}

#liste mit allen infos (ES NES pWERT) zu den GSEAS der einzelnen patienten
GSEA_list = apply(tcga_exp_cleaned[,1:10], 2, GSEA) #erstmal nur die ersten zehn patienten weil das sonst exig dauert
#hier muss jetzt eine fuktion hin die die NEM spalte aller dataframes aus der GSEA liste pulled und als df speichert


save(GSEA_list, file = '~/GitHub/2022-topic-02-team-05/data/GSEA_matrix.RData')


