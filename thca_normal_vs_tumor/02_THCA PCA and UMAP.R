library(ggplot2)
library(metaplot)
library(gridExtra)
library(umap)

#install.packages("umap")
load('data/thca_gsea.RData')
load('data/thca_anno.RData')

#--------------------------------------------
#Durchf?hren der PCA und UMAP
#---------------------------------------------
set.seed(123)
PCA = prcomp(t(thca_gsea))
PCA_data = as.data.frame(PCA$x)
PCA_data$stage = thca_anno$ajcc_pathologic_tumor_stage
PCA_data$gender = thca_anno$gender
PCA_data$type = sapply(thca_anno$histological_type, FUN = function(x){return(strsplit(x, split = '-', fixed = TRUE)[[1]][2])})


UMAP = umap(PCA$x)
UMAP_data = as.data.frame(UMAP$layout)  #UMAP der PCA-Daten
UMAP_data$stage = thca_anno$ajcc_pathologic_tumor_stage
UMAP_data$gender = thca_anno$gender
UMAP_data$type = sapply(thca_anno$histological_type, FUN = function(x)
                    {return(strsplit(x, split = '-', fixed = TRUE)[[1]][2])})

colours = c('darkgreen','yellow3','yellow2','dodgerblue','blue4')

#Plottet die PCA
ggplot(PCA_data, aes(x = PC1, y = PC2, color = stage, shape = type)) + geom_point(size = 2) +
  scale_color_manual(values = colours) +
  labs(title = 'PCA of THCA expression data') +
  theme(legend.key = element_rect(fill = 'white')) +
  theme_minimal() +
  guides(color = guide_legend(override.aes = list(size = 5)))
#Plottet unserer UMAP
ggplot(UMAP_data, aes(x = V1, y = V2, color = stage, shape = type)) + geom_point(size = 2) +
  scale_color_manual(values = colours) +
  labs(title = 'UMAP of THCA expression data') +
  theme(legend.key = element_rect(fill = 'white')) +
  theme_minimal() +
  guides(color = guide_legend(override.aes = list(size = 5)))

#--------------------------------------
#Zum Vergleichen führen wir nun auch PCA und UMAP auf die reinen Exp Daten durch
#Dau verwenden wir den log2FC
#--------------------------------------
set.seed(123)
load('data/thca_normal_exp_cleaned.RData')
load('data/thca_tumor_exp_cleaned.RData')

#log2FC für jeden Expressionwert berrechnen
geneFC = thca_tumor_exp_cleaned - thca_normal_exp_cleaned

#PCA
PCA_genes = prcomp(t(geneFC))
PCA_genes_data = as.data.frame(PCA_genes$x)
PCA_genes_data$type = thca_anno$histological_type

#UMAP
UMAP_genes = umap(PCA_genes$x)
UMAP_genes_data = as.data.frame(UMAP_genes$layout)
UMAP_genes_data$type = thca_anno$histological_type

type.col = c('Classical/usual' = 'blue4',
             'Follicular (>= 99% follicular patterned)' = 'darkolivegreen3',
             'Tall Cell (>= 50% tall cell features)' = 'orange',
             'other' = 'honeydew3')

#Plottet die PCA
ggplot(PCA_genes_data, aes(x = PC1, y = PC2, color = type)) + geom_point(size = 2) +
  scale_color_manual(values = type.col) +
  labs(title = 'PCA of THCA expression data') +
  theme(legend.key = element_rect(fill = 'white')) +
  theme_minimal() +
  guides(color = guide_legend(override.aes = list(size = 5)))
#Plottet unserer UMAP
ggplot(UMAP_genes_data, aes(x = V1, y = V2, color = type)) + geom_point(size = 2) +
  scale_color_manual(values = type.col) +
  labs(title = 'UMAP of THCA expression data') +
  theme(legend.key = element_rect(fill = 'white')) +
  theme_minimal() +
  guides(color = guide_legend(override.aes = list(size = 5)))



