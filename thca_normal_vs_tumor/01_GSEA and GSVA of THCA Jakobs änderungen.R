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
load('data/thca_anno.RData')

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

#-----------------------------------------------
#Darstellung der GSEA Ergebnisse als Heatmap
#-----------------------------------------------
#Annotationen für die Pathways                              
path.anno = rowAnnotation(Pathway = c(rep('met',611),rep('hall', 46)),
                          col = list(Pathway = c("met" = "deepskyblue", "hall" = "blue4")),
                          annotation_legend_param = list(title = 'Pathway type',
                                                         at = c('met', 'hall'),
                                                         labels = c('Metabolic', 'Hallmark')))
#Annotation der kmeans cluster und den Histological type
top.anno = HeatmapAnnotation(Type = thca_anno$histological_type,
                             col = list(Type = c('Classical/usual' = 'olivedrab4',
                                                 'Follicular (>= 99% follicular patterned)' = 'darkolivegreen',
                                                 'Tall Cell (>= 50% tall cell features)' = 'green3',
                                                 'other' = 'honeydew3')),
                             annotation_legend_param = list(
                                title = 'Histological type',
                                at = names(table(thca_anno$histological_type)),
                                lables = names(table(thca_anno$histological_type)))
                             )
Heatmap(thca_gsea,
        show_row_names = F, show_column_names = F, width = unit(21, 'cm'), height = unit(16, 'cm'),
        heatmap_legend_param = list(
          title = "Thyroid cancer pathway activity", at = c(-2, 2), 
          labels = c("underexpressed", "overexpressed")),
        row_dend_reorder = T, column_dend_reorder = T,
        left_annotation = path.anno,
        top_annotation = top.anno,
        column_km = 3 #Splitted unserer Heatmap in die 3 Subtypen die wir in der PanCancer gefunden haben
)



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

# #Volcanoplot
# alpha.kor = 0.1 #Siginfikanzniveau
# load('data/thca_genes_cleaned.RData')
# thca_volcano = data.frame(thca_logFC_gsva, thca_pval_gsva) #die 2 Vektoren f?r unseren Volcano PLot werden in einen df gepackt, damit daraus ein plot erstellt werden kann
# with(thca_volcano, plot(thca_volcano$thca_logFC_gsva, -log10(thca_volcano$thca_pval_gsva), main = "Volcano plot of THCA pathway activity", xlab='log2(foldchange)',
#                      ylab = '-log10(Pvalues)', pch = 20))
# 
# with(subset(thca_volcano, thca_logFC_gsva>alpha.kor & thca_pval_gsva < alpha.kor), points(thca_logFC_gsva, -log10(thca_pval_gsva), col="orange", pch = 20))
# with(subset(thca_volcano, thca_logFC_gsva< -alpha.kor & thca_pval_gsva < alpha.kor), points(thca_logFC_gsva, -log10(thca_pval_gsva), col="green", pch = 20))
# with(subset(thca_volcano, abs(thca_logFC_gsva) < alpha.kor | thca_pval_gsva > alpha.kor), points(thca_logFC_gsva, -log10(thca_pval_gsva), pch=19, col="gray"))
