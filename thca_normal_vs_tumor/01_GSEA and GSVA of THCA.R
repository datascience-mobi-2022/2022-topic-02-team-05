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
library(GSVA)
library(ComplexHeatmap)

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
# library(pheatmap)
# #darstellen der GSEA matrix als heatmap
# pheatmap(thca_gsea,
#          breaks = seq(-max(thca_gsea), max(thca_gsea), length.out = 201),
#          color = colorRampPalette(c('blue','lightskyblue','lightblue1', 'white','lightyellow','yellow', 'red'),
#                                   bias = 1,
#                                   space = 'rgb',
#                                   interpolate = 'linear'
#          )(200),
#          clustering_method = 'average',
#          treeheight_row = 40, treeheight_col = 20, cellwidth = 9,cellheight = 0.75,
#          show_colnames = TRUE, show_rownames = FALSE, fontsize_col = 8, border_color = NA,
#          legend_breaks = c(-max(thca_gsea),0, max(thca_gsea)),
#          legend_labels = c('underexpressed', 'normal expression', 'overexpressed')
# )

#---------------------------------------------
#Analyse der Daten mittels GSVA
#---------------------------------------------
#Durchführen eine GSVA zuerst auf Normal und Tumorgewebe
thca_norm_gsva = gsva(as.matrix(thca_normal_exp_cleaned), pathways,
                 method = 'gsva',
                 kcdf = 'Gaussian'  #Da wir kontinuierliche Daten haben
)
thca_tumor_gsva = gsva(as.matrix(thca_tumor_exp_cleaned), pathways,
                      method = 'gsva',
                      kcdf = 'Gaussian'  #Da wir kontinuierliche Daten haben
)

#--------------------------------------
#Volcano plot
#--------------------------------------

#Berrechnen des Foldchanges zwischne beiden Daten
thca_logFC_gsva = apply(thca_tumor_gsva, 1, mean) - apply(thca_norm_gsva, 1, mean)
#pvalue berechen
thca_pval_gsva = c()
for (i in (1:nrow(thca_norm_gsva))){
  res = wilcox.test(thca_tumor_gsva[i,], thca_norm_gsva[i,], alternative = 'two.sided')$p.value
  thca_pval_gsva = append(thca_pval_gsva, res)
};rm(i);rm(res)

#Volcanoplot
alpha.kor = 0.1 #Siginfikanzniveau
load('data/thca_genes_cleaned.RData')
thca_volcano = data.frame(thca_logFC_gsva, thca_pval_gsva) #die 2 Vektoren f?r unseren Volcano PLot werden in einen df gepackt, damit daraus ein plot erstellt werden kann
with(thca_volcano, plot(thca_volcano$thca_logFC_gsva, -log10(thca_volcano$thca_pval_gsva), main = "Volcano plot of THCA pathway activity", xlab='log2(foldchange)',
                     ylab = '-log10(Pvalues)', pch = 20))

with(subset(thca_volcano, thca_logFC_gsva>alpha.kor & thca_pval_gsva < alpha.kor), points(thca_logFC_gsva, -log10(thca_pval_gsva), col="orange", pch = 20))
with(subset(thca_volcano, thca_logFC_gsva< -alpha.kor & thca_pval_gsva < alpha.kor), points(thca_logFC_gsva, -log10(thca_pval_gsva), col="green", pch = 20))
with(subset(thca_volcano, abs(thca_logFC_gsva) < alpha.kor | thca_pval_gsva > alpha.kor), points(thca_logFC_gsva, -log10(thca_pval_gsva), pch=19, col="gray"))
