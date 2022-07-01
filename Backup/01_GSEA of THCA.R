#-------------------------------
#In diesem Dokument werden wir eine GSEA für THCA mit allen unseren Pathways durchführen
#-------------------------------

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("phenoTest")
# BiocManager::install("gage")
# BiocManager::install("fgsea")
library(dplyr)
library(fgsea)
library(phenoTest) #alternatives package
# library(gage)

load('data/thca_tumor_exp_cleaned.RData')
load('data/thca_normal_exp_cleaned.RData')
load('data/pathways.RData')

#GSEA Funktion definieren
GSEA = function(patientsorted, pathways){
  res = fgseaMultilevel(pathways = pathways, 
                        stats = patientsorted)
  #extrahieren der NES werte & nullsetzten aller NAs
  ret = res$ES; names(ret) = res$pathway
  message('I?m still standing')
  return(ret)
}
#Liste mit nach Sortierten patienten erstellen
#hierzu wird zuerst der FC zw tumor und normalgewebe berrechnet und dann jeder patient
#danach sortiert
log2FC = thca_tumor_exp_cleaned - thca_normal_exp_cleaned
patient_sorted = apply(log2FC, 2, FUN = function(x){sort(x, decreasing = TRUE)}, simplify = F)

#Durchführen der GSEA
thca_gsea = sapply(patient_sorted, FUN = function(x){GSEA(x, pathways)})
save(thca_gsea, file = 'data/thca_gsea.RData')


#Darstellung als Heatmap
library(pheatmap)
#darstellen der GSEA matrix als heatmap
pheatmap(thca_gsea,
         breaks = seq(-max(thca_gsea), max(thca_gsea), length.out = 201),
         color = colorRampPalette(c('blue','lightskyblue','lightblue1', 'white','lightyellow','yellow', 'red'),
                                  bias = 1,
                                  space = 'rgb',
                                  interpolate = 'linear'
         )(200),
         clustering_method = 'average',
         treeheight_row = 40, treeheight_col = 20, cellwidth = 9,cellheight = 0.75,
         show_colnames = TRUE, show_rownames = FALSE, fontsize_col = 8, border_color = NA,
         legend_breaks = c(-max(thca_gsea),0, max(thca_gsea)),
         legend_labels = c('underexpressed', 'normal expression', 'overexpressed')
)
