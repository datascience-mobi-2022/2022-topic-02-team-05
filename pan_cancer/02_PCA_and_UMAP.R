#--------------------------------------
#In diesem Dokument nehmen wir die GSEA Daten und führen eine PCA, sowie 
#ein clustering der transformierten Daten durch. So können wir verschiedene
#unterschiedliche Pathwayaktivitäten zuordnen
#------------------------------------------

library(ggplot2)
library(metaplot)
library(gridExtra)
library(umap)

#install.packages("umap")
load('data/tcga_gsva.RData')
load('data/tcga_anno.RData')

#--------------------------------------------
#Durchführen der PCA und UMAP
#---------------------------------------------

set.seed(123)
PCA = prcomp(t(tcga_gsva))
PCA_data = as.data.frame(PCA$x)
PCA_data$type = tcga_anno$cancer_type_abbreviation
PCA_data$form = tcga_anno$cancer_form

UMAP = umap(PCA$x)
UMAP_data = as.data.frame(UMAP$layout)  #UMAP der PCA-Daten
UMAP_data$type = tcga_anno$cancer_type_abbreviation
UMAP_data$form = tcga_anno$cancer_form

colours = c('blue4','dodgerblue4','deepskyblue','cyan',
            'lightblue1','aquamarine','chartreuse4',
            'aquamarine3','darkolivegreen','darkolivegreen4',
            'chartreuse2','darkolivegreen2','lemonchiffon',
            'yellow','peachpuff1','gold','orange','red',
            'indianred3','orangered4','sienna4','tan3','salmon2',
            'plum','rosybrown1','violetred','magenta','magenta4',
            'maroon','wheat1','snow3','gray29','black')

#Plottet die PCA mit Krebsarten nach Farbe
ggplot(PCA_data, aes(x = PC1, y = PC2, color = type)) + geom_point(size = 2) +
  scale_color_manual(values = colours) +
  labs(title = 'PCA of TCGA pathway activity') +
  theme(legend.key = element_rect(fill = 'white')) +
  theme_minimal() +
  guides(color = guide_legend(override.aes = list(size = 5)))

#Plottet die UMAP mit Krebsarten nach Farbe
ggplot(UMAP_data, aes(x = V1, y = V2, color = type)) + geom_point(size = 1) +
  scale_color_manual(values = colours) +
  labs(title = 'UMAP of TCGA pathway activity') +
  theme(legend.key = element_rect(fill = 'white')) +
  theme_minimal() +
  guides(color = guide_legend(override.aes = list(size = 5)))

form.col = c('Adenocarcinoma' = 'deepskyblue3',
             'Squamous cell carcinoma' = 'dodgerblue4',
             'Transitional cell carcinoma' = 'deepskyblue',
             'Melanoma'='cyan',
             'Carcinoma'='azure3',
             'Sarcoma'='darkolivegreen',
             'Glioblastoma'='gold',
             'Leukemia'='red')

#Plottet die PCA mit Krebsüberformen nach Farbe
ggplot(PCA_data, aes(x = PC1, y = PC2, color = form)) + geom_point(size = 2) +
  scale_color_manual(values = form.col) +
  labs(title = 'PCA of TCGA pathway activity') +
  theme(legend.key = element_rect(fill = 'white')) +
  theme_minimal() +
  guides(color = guide_legend(override.aes = list(size = 5)))
#Plottet unserer UMAP mit Krebsüberformen nach Farbe
ggplot(UMAP_data, aes(x = V1, y = V2, color = form)) + geom_point(size = 1) +
  scale_color_manual(values = form.col) +
  labs(title = 'UMAP of TCGA pathway activity') +
  theme(legend.key = element_rect(fill = 'white')) +
  theme_minimal() +
  guides(color = guide_legend(override.aes = list(size = 5)))


#------------------------------------------------------------
#zum Vergelich führen wir nun PCA und UMAP auf den Genen statt auf den Pathway
#Aktivitäten durch
#------------------------------------------------------------
library(Seurat) #package für schnelle PCA
set.seed(123)
load('data/tcga_exp_cleaned.RData') 

PCA_genes = RunPCA(t(tcga_exp_cleaned), npcs = 10)

#UMAP
UMAP_genes = umap(PCA_genes@feature.loadings)
UMAP_genes_data = as.data.frame(UMAP_genes$layout)
UMAP_genes_data$type = tcga_anno$cancer_type_abbreviation
UMAP_genes_data$form = tcga_anno$cancer_form

#Plottet unserer UMAP mit Krebsarten nach Farbe
ggplot(UMAP_genes_data, aes(x = V1, y = V2, color = type)) + geom_point(size = 1) +
  scale_color_manual(values = colours) +
  labs(title = 'UMAP of TCGA gene expression') +
  theme(legend.key = element_rect(fill = 'white')) +
  theme_minimal() +
  guides(color = guide_legend(override.aes = list(size = 5)))

#Plottet unserer UMAP mit Krebsüberformen nach Farbe
ggplot(UMAP_genes_data, aes(x = V1, y = V2, color = form)) + geom_point(size = 1) +
  scale_color_manual(values = form.col) +
  labs(title = 'UMAP of TCGA gene expression') +
  theme(legend.key = element_rect(fill = 'white')) +
  theme_minimal() +
  guides(color = guide_legend(override.aes = list(size = 5)))


#--------------------------------------------
#Visualisierung von der durch die PCs erklärten Varianz und den Pathwaykomponenten
#----------------------------------------------

#Varianz die erklärt wird
PCs = 1:658; var_prop = PCA$sdev^2/sum(PCA$sdev^2); var_cum = cumsum(var_prop)
Var = cbind.data.frame(PCs, var_prop, var_cum)
v1 = ggplot(Var[1:10,], aes(PCs, var_prop)) + geom_point()+
  labs(y = "Proportion of Variance Explained")+
  scale_x_continuous(name = 'Principle Components', breaks = 1:10)

v2 = ggplot(Var[1:10,], aes(PCs, var_cum)) + geom_point()+
  labs(y = "Cummulative Proportion of Variance Explained")+
  scale_x_continuous(name = 'Principle Components', breaks = 1:10)
grid.arrange(grobs = list(v1,v2), ncol= 2)

#visualisierung der pathways, die in einer PC stecken
pathways = as.data.frame(PCA$rotation)
pathways$names = row.names(pathways)

pathways = pathways[order(abs(pathways$PC1), decreasing = TRUE),]
t1 = ggplot(pathways[1:10,], aes(names, PC1)) + geom_col()+
  labs(x = 'Pathway', y= 'Contribution to PC1')+
  theme(axis.text = element_text(angle = 90,))

pathways = pathways[order(abs(pathways$PC2), decreasing = TRUE),]
t2 = ggplot(pathways[1:10,], aes(names, PC2)) + geom_col()+
  labs(x = 'Pathway', y= 'Contribution to PC2')+
  theme(axis.text = element_text(angle = 90,))

pathways = pathways[order(abs(pathways$PC3), decreasing = TRUE),]
t3 = ggplot(pathways[1:10,], aes(names, PC3)) + geom_col()+
  labs(x = 'Pathway', y= 'Contribution to PC3')+
  theme(axis.text = element_text(angle = 90,))
grid.arrange(grobs = list(t1, t2, t3), ncol = 3)
