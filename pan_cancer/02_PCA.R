#--------------------------------------
# In diesem Dkument nehemn wir die GSEA Daten und f?hren eine PCA sowie wie
# ein clustering der transformierten Daten durch so k?nnen wir hoffentlich verschiedenen
# unterschiedliche Pathway activt?ten zuordnen
#------------------------------------------
library(ggplot2)
library(metaplot)
library(gridExtra)
library(umap)

#install.packages("umap")
load('~/GitHub/2022-topic-02-team-05/data/tcga_gsva.RData')
tcga_anno = readRDS('~/GitHub/2022-topic-02-team-05/data/tcga_tumor_annotation.RDS')

#--------------------------------------------
#Durchf?hren der PCA und UMAP
#---------------------------------------------
set.seed(123)
PCA = prcomp(t(tcga_gsva))
PCA_data = as.data.frame(PCA$x)
PCA_data$type = sapply(rownames(PCA_data), FUN = function(x){
  return(tcga_anno$cancer_type_abbreviation[x == tcga_anno$sample])})

UMAP = umap(PCA$x)
UMAP_data = as.data.frame(UMAP$layout)  #UMAP der PCA-Daten
UMAP_data$type = sapply(rownames(UMAP_data), FUN = function(x){
  return(tcga_anno$cancer_type_abbreviation[x == tcga_anno$sample])})

colours = c('blue4','dodgerblue4','deepskyblue','cyan',
            'lightblue1','aquamarine','chartreuse4',
            'aquamarine3','darkolivegreen','darkolivegreen4',
            'chartreuse2','darkolivegreen2','lemonchiffon',
            'yellow','peachpuff1','gold','orange','red',
            'indianred3','orangered4','sienna4','tan3','salmon2',
            'plum','rosybrown1','violetred','magenta','magenta4',
            'maroon','wheat1','snow3','gray29','black')

#Plottet die PCA
ggplot(PCA_data, aes(x = PC1, y = PC2, color = type)) + geom_point(size = 2) +
  scale_color_manual(values = colours) +
  labs(title = 'UMAP of TCGA expression data') +
  theme(legend.key = element_rect(fill = 'white')) +
  theme_minimal() +
  guides(color = guide_legend(override.aes = list(size = 5)))
#Plottet unserer UMAP
ggplot(UMAP_data, aes(x = V1, y = V2, color = type)) + geom_point(size = 1) +
  scale_color_manual(values = colours) +
  labs(title = 'UMAP of TCGA expression data') +
  theme(legend.key = element_rect(fill = 'white')) +
  theme_minimal() +
  guides(color = guide_legend(override.aes = list(size = 5)))





#--------------------------------------------
#Visualisierung von erklärt varianz und den Pathwaykomponenten
#----------------------------------------------

#--------------------------------
#Varianz die erkl?rt wird
#--------------------------------
PCs = 1:658; var_prop = PCA$sdev^2/sum(PCA$sdev^2); var_cum = cumsum(var_prop)
Var = cbind.data.frame(PCs, var_prop, var_cum)
v1 = ggplot(Var[1:10,], aes(PCs, var_prop)) + geom_point()+
  labs(y = "Proportion of Variance Explained")+
  scale_x_continuous(name = 'Principle Components', breaks = 1:10)

v2 = ggplot(Var[1:10,], aes(PCs, var_cum)) + geom_point()+
  labs(y = "Cummulative Proportion of Variance Explained")+
  scale_x_continuous(name = 'Principle Components', breaks = 1:10)
grid.arrange(grobs = list(v1,v2), ncol= 2)

#visualisierung der pathways die in einer PC stecken
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

grid.arrange(p1, arrangeGrob(p2,p3, ncol=2), heights=c(2.5/4, 2.5/4), ncol=1,top = 'PCA of GSEA data')

grid.arrange(arrangeGrob(v1, v2, ncol = 2),
             arrangeGrob(t1,t2,t3, ncol = 3),
             arrangeGrob(p1,p2,p3, ncol = 3),
             ncol=1)
