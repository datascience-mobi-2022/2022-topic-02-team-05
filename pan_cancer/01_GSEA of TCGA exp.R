#--------------------------------------------------
# In diesem Dukoment versuchen wir die Gereinigten expressionsdaten mittels gsea zu
#pathway activit?ten zusammenzufassen un din einer matrix mit patienten als spalten und 
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
load("~/GitHub/2022-topic-02-team-05/data/geneset_ids.RData")
tcga_anno = readRDS('~/GitHub/2022-topic-02-team-05/data/tcga_tumor_annotation.RDS')


#--------------------------------
#Extraktion aller ensembl ids und gennamen aus den cleanen exp daten
#--------------------------------
tcga_genes_cleaned = rownames(tcga_exp_cleaned)
#dieser vetor enth?lt sowohl enseblm id als auch genenamen und muss daher gespalten werden
tcga_genes_cleaned = strsplit(tcga_genes_cleaned, split = '|', fixed = TRUE)
#speicher der ensembl ids ohne versionsnummer als eigenen vektor
tcga_geneids = sapply(tcga_genes_cleaned, function(tcga_genes_cleaned){return(tcga_genes_cleaned[1])})
tcga_geneids = strsplit(tcga_geneids, split = '.', fixed = TRUE)
tcga_geneids = sapply(tcga_geneids, function(tcga_geneids){return(tcga_geneids[1])})

#speicher der genenames ohne versionsnummer als eigenen vektor
tcga_genenames = sapply(tcga_genes_cleaned, function(tcga_genes_cleaned){return(tcga_genes_cleaned[2])})
tcga_genenames = strsplit(tcga_genenames, split = '.', fixed = TRUE)
tcga_genenames = sapply(tcga_genenames, function(tcga_genenames){return(tcga_genenames[1])})

tcga_genes_cleaned = cbind.data.frame(tcga_geneids,tcga_genenames)
#speichern eines datframes der der die Ensembl ids und genenamen aller genen der exp daten enth?lt
save(tcga_genes_cleaned, file = '~/GitHub/2022-topic-02-team-05/data/tcga_genes_cleaned.RData')

#----------------------
#versuch eine gsea zu implementieren
#------------------------

#function die eine gsea f?r einen beliebigen patienten durchf?hrt und nur den NES benannt ausgibt
GSEA_NES = function(patient){
  
  #sorting and naming the exp data
  names(patient) = tcga_genes_cleaned$tcga_geneids
  patient = sort(patient, decreasing = TRUE) 
  
  #duarchf?hren der GSEA
  res = fgseaMultilevel(pathways = pathways, 
                        stats = patient,
                        minSize=3)
  #extrahieren der NES werte & nullsetzten aller nichtsiginifikanten
  ret = res$NES; names(ret) = res$pathway
  #ret[res$padj > pvalue] = 0
  ret[is.na(ret)] = 0
  
  message('I?m still standing')
  return(ret)
}

#-------------------------------------------
#Durchf?hrung der GSEA
#-----------------------------------------------
#pathways = gmtPathways("C:/Users/jakob/Desktop/c2.cp.v7.5.1.symbols.gmt") alle canonischen pathways msigr vllt als alternative
pathways = c(genesets_ids, our_genesets)
#pvalue = 0.05

#ztransforamtion der daten zur gsea analyse nach Peng et al.
tcga_exp_ztrans = apply(as.matrix(tcga_exp_cleaned), 2, scale)


#verbesserte Version
# extraktion aller patienten mit einem tumortype
tcga_anno = tcga_anno[order(tcga_anno$cancer_type_abbreviation),]

cancers = list();cancers = vector('list',length(table(tcga_anno$cancer_type_abbreviation)))
names(cancers) = names(table(tcga_anno$cancer_type_abbreviation))

y = 1 ;z = 1;type = 'ACC'
for (j in 1:length(tcga_anno$cancer_type_abbreviation)){
  
  if(tcga_anno$cancer_type_abbreviation[j] == type){
      cancers[[y]][[z]] = tcga_anno$sample[j]
      z = z+1
  } else {
      type = tcga_anno$cancer_type_abbreviation[j]
      y = y+1
      z = 1
      cancers[[y]][[z]] = tcga_anno$sample[j]
      
  }
  message(j)
}
#list mit allen patient mit einem cancertype
canc = list()
for(i in 1:length(cancers)){
    canc[[i]] = (sapply(cancers[[i]], FUN = function(x){return(x)}))
}
names(canc) = names(table(tcga_anno$cancer_type_abbreviation))


tcga_cancers = lapply(canc, FUN = function(x){
  return(tcga_exp_cleaned[,tcga_exp_cleaned == x])
})

test = tcga_exp_cleaned[3,colnames(tcga_exp_cleaned) == canc[['ACC']]]

GSEA_matrix = apply(tcga_exp_ztrans[,1:100], 2, GSEA_NES) %>% as.data.frame()
save(GSEA_matrix, file = '~/GitHub/2022-topic-02-team-05/data/GSEA_matrix.RData')

#darstellen der GSEA matrix als heatmap
#rote sind ?berexpremierte pathways, blue unterexpremierte, wei? keine abweichung
pheatmap(as.matrix(GSEA_matrix),
         breaks = seq(-max(GSEA_matrix), max(GSEA_matrix), length.out = 201),
         color = colorRampPalette(c('blue','lightskyblue','lightblue1', 'white','lightyellow','yellow', 'red'),
                                  bias = 1,
                                  space = 'rgb',
                                  interpolate = 'linear'
         )(200),
         clustering_method = 'average',
         treeheight_row = 25, treeheight_col = 20, cellwidth = 7,cellheight = 8,
         show_colnames = FALSE, fontsize_row = 8, border_color = NA,
         legend_breaks = c(-max(GSEA_matrix),0, max(GSEA_matrix)),
         legend_labels = c('underexpressed', 'normal expression', 'overexpressed')
        )






  

