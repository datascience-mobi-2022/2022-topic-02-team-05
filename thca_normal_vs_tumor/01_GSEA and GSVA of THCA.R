#-------------------------------
#In diesem Dokument werden wir eine GSEA für THCA mit allen unseren Pathways durchführen
#-------------------------------

#if (!require("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")
#BiocManager::install("phenoTest")
#BiocManager::install("gage")
#BiocManager::install("fgsea")
#BiocManager::install("GSVA")

library(dplyr)
library(fgsea)
library(GSVA)
library(ComplexHeatmap)
library(ggplot2)

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

#Berrechnen des Foldchanges zwischen beiden Daten
thca_logFC_gsva = apply(thca_tumor_gsva, 1, mean) - apply(thca_norm_gsva, 1, mean)

#pvalue berechen
thca_pval_gsva = c()
for (i in (1:nrow(thca_norm_gsva))){
  res = wilcox.test(thca_tumor_gsva[i,], thca_norm_gsva[i,], alternative = 'two.sided')$p.value
  thca_pval_gsva = append(thca_pval_gsva, res)
};rm(i);rm(res)

#signifikanzlevel ohne bonferroni 
alpha = 0.025

#Erstellen dataframe
data.thca = data.frame(thca_logFC_gsva, thca_pval_gsva)
cbind(data.thca, rownames(data.thca)) -> data.thca
colnames(data.thca) <- c("logFC", "Pvalues", "pathway_names")

#hinzufügen einer Spalte, die sagt, ob das Gen up- oder downregulated wird
#hinzufügen einer Spalte diffexpressed mit NOs 
data.thca$diffexpressed <- "NO"
#wenn log2Foldchange > 0.1 and pvalue < alpha.kor, set as "UP" 
data.thca$diffexpressed[data.thca$logFC > 0.00001 & data.thca$Pvalue > alpha] <- "UP"
# if log2Foldchange < -0.1 and pvalue < 0.05, set as "DOWN"
data.thca$diffexpressed[data.thca$logFC < -0.00001 & data.thca$Pvalue > alpha] <- "DOWN"

volcano2 = ggplot(data = data.thca, aes(x = logFC, y = -log10(Pvalues), color = diffexpressed, label = pathway_names)) +
  geom_point() +
  theme_minimal() +
  ggtitle("Volcanoplot - GSVA of tumor tissue vs GSVA of normal tissue") +
  geom_text(data = subset(data.thca, -log10(thca_pval_gsva) > 0.5) , size = 1.5, nudge_y = 0.03, check_overlap = TRUE)


volcano2


#---------------------------------------------
#Plotten des P-Werts aller überexprimierten und aller unterexprimierten Gene 
#---------------------------------------------

upregulated <- as.data.frame(data.thca$Pvalues[data.thca$diffexpressed == "UP"])
rownames(upregulated) <- data.thca$pathway_names[data.thca$diffexpressed == "UP"]
colnames(upregulated) <- c("Pvalues_up")
normalized_up <- (upregulated$Pvalues_up - mean(upregulated$Pvalues_up)) / var(upregulated$Pvalues_up)
cbind(upregulated, normalized_up) -> upregulated

downregulated <- as.data.frame(data.thca$Pvalues[data.thca$diffexpressed == "DOWN"])
rownames(downregulated) <- data.thca$pathway_names[data.thca$diffexpressed == "DOWN"]
colnames(downregulated) <- c("Pvalues_down") 
normalized_down <- (downregulated$Pvalues_down - mean(downregulated$Pvalues_down)) / var(downregulated$Pvalues_down)
cbind(downregulated, normalized_down) -> downregulated

upregulated_geranked <- rank(upregulated$normalized_up)
downregulated_geranked <- rank(downregulated$normalized_down)

up_plot <- ggplot(data = upregulated, aes(y = -log10(Pvalues_up), x = rank(log10(Pvalues_up)), label = rownames(upregulated))) + 
  geom_point() + 
  geom_text(size = 1, hjust = -0.1, check_overlap = TRUE)  

up_plot

down_plot <- ggplot(data = downregulated, aes(y = -log10(Pvalues_down), x = rank(log10(Pvalues_down)), label = rownames(downregulated))) + 
  geom_point() + 
  geom_text(size = 1, hjust = -0.1, check_overlap = TRUE)  

down_plot


#volcano <- ggplot(data = data.thca, aes(x = log2fc.thca, y = -log10(p.values), color = diffexpressed, label = thca_genenames)) +
  #geom_point() + 
  #theme_minimal() + 
  #ggtitle("Volcanoplot") + 
  #geom_hline(yintercept = -log10(alpha.kor), col = "black", show.legend = TRUE) + 
  #geom_text(data = subset(data.thca, -log10(p.values) > 22), size = 2, check_overlap = TRUE, nudge_x = 0.1, nudge_y = 0.8, color = "black")

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
