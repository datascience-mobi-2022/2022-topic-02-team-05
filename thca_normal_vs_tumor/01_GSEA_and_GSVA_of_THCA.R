#-------------------------------
#In diesem Dokument werden wir eine GSEA und GSVA für THCA mit allen unseren pathways durchführen
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
load('data/thca_anno.RData')

#GSEA Funktion definieren
GSEA = function(patientsorted, pathways){
  res = fgseaMultilevel(pathways = pathways, 
                        stats = patientsorted)
  #extrahieren der NES werte & Nullsetzten aller NAs
  ret = res$ES; names(ret) = res$pathway
  message('I?m still standing')
  return(ret)
}

#Liste mit sortierten Patienten erstellen
#hierzu wird zuerst der foldchange zwischen Tumor- und Normalgewebe berechnet und dann jeder Patient
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

#Zusammenfügen der Tumor und Normal Exp Daten in einen Df für die GSVA
thca_exp_data = cbind.data.frame(thca_tumor_exp_cleaned, thca_normal_exp_cleaned)

#Durchführen der GSVA
thca_gsva_data = gsva(as.matrix(thca_exp_data), pathways,
                      method = 'gsva',
                      kcdf = 'Gaussian'  #Da wir kontinuierliche Daten haben
)

#Splitten des GSVA Dfs in Tumor und Normalgewebe
thca_norm_gsva  = thca_gsva_data[, colnames(thca_normal_exp_cleaned) == colnames(thca_gsva_data)]
thca_tumor_gsva = thca_gsva_data[, colnames(thca_tumor_exp_cleaned) == colnames(thca_gsva_data)]


#--------------------------------------
#Volcano plot erstellen
#--------------------------------------

#Berechnen des foldchanges zwischen beiden Daten
thca_logFC_gsva = apply(thca_tumor_gsva, 1, mean) - apply(thca_norm_gsva, 1, mean)

#p-Value berechen
thca_pval_gsva = c()
for (i in (1:nrow(thca_norm_gsva))){
  res = wilcox.test(thca_tumor_gsva[i,], thca_norm_gsva[i,], alternative = 'two.sided')$p.value
  thca_pval_gsva = append(thca_pval_gsva, res)
};rm(i);rm(res)

#Speichern der p-Values für die Regression zur Selektion eines pathways
names(thca_pval_gsva) = rownames(thca_gsva_data)
save(thca_pval_gsva, file = 'data/regression/thca_pval_gsva.RData')

#Signifikanzlevel ohne bonferroni 
alpha = 0.025 #,da beidseitiger Test

#Erstellen dataframe
data.thca = data.frame(thca_logFC_gsva, thca_pval_gsva)
cbind(data.thca, rownames(data.thca)) -> data.thca
colnames(data.thca) <- c("logFC", "Pvalues", "pathway_names")

#hinzufügen einer Spalte, die sagt, ob das Gen up- oder downregulated wird
#hinzufügen einer Spalte diffexpressed mit NOs 
data.thca$diffexpressed <- "NO"

#wenn log2Foldchange > 0 und pvalue < alpha.kor, als "UP" definieren
data.thca$diffexpressed[data.thca$logFC > 0 & data.thca$Pvalue < alpha] <- "UP"
# if log2Foldchange < 0 and pvalue < 0.05, als "DOWN" definieren
data.thca$diffexpressed[data.thca$logFC < 0 & data.thca$Pvalue < alpha] <- "DOWN"

volcano2 = ggplot(data = data.thca, aes(x = logFC, y = -log10(Pvalues), color = diffexpressed, label = pathway_names)) +
  geom_point() +
  theme_minimal() +
  xlab("log2 foldchange") +
  ylab("-log(P-value)") +
  ggtitle("Volcanoplot - GSVA of tumor tissue vs normal tissue") +
  geom_text(data = subset(data.thca, -log10(thca_pval_gsva) > 11) , size = 1.5, nudge_y = 0.03, check_overlap = TRUE)

volcano2


#---------------------------------------------
#Plotten der P-Werte aller überexprimierten und aller unterexprimierten Gene 
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

#plotten der pathways, die in THCA upregulated sind
up_plot <- ggplot(data = upregulated, aes(y = -log10(Pvalues_up), x = rank(log10(Pvalues_up)), label = rownames(upregulated))) + 
  geom_point() + 
  geom_text(size = 1, hjust = -0.1, check_overlap = TRUE) +
  xlab("P-value ranks") +
  ylab("-log(P-value)") +
  ggtitle("Pathways upregulated in THCA") +
  theme_light()

up_plot

#plotten der pathways, die in THCA downregulated sind
down_plot <- ggplot(data = downregulated, aes(y = -log10(Pvalues_down), x = rank(log10(Pvalues_down)), label = rownames(downregulated))) + 
  geom_point() + 
  geom_text(data = downregulated, size = 1, hjust = -0.1, check_overlap = TRUE) +
  xlab("P-value ranks") +
  ylab("-log(P-value)") +
  ggtitle("Pathways downregulated in THCA") +
  theme_light()

down_plot

#Ausgabe der Namen der 10 am signifikantesten Hoch/runterregulierten pathways
names_dn = rownames(downregulated)[match(sort(downregulated$Pvalues_down), downregulated$Pvalues_down)]
names_dn[1:10]

names_up = rownames(upregulated)[match(sort(upregulated$Pvalues_up), upregulated$Pvalues_up)]
names_up[1:10]

#-------------------------------------------------
# Vergleich der GSEA und GSVA ergebnisse
#Dazu wird die Corealtion zwischen beiden Methoden für jeden Pathway bestimmt
#-------------------------------------------------
GSVA_act = thca_tumor_gsva - thca_norm_gsva
cor = vector(length = nrow(GSVA_act))
i = 1; for (i in 1:nrow(GSVA_act)){
  res = cor(thca_gsea[i,], GSVA_act[i,])
  cor[i] = res
 
}; rm(i, res); names(cor) = rownames(GSVA_act)
boxplot(cor)
