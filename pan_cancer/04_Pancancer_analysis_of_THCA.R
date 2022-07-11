#-----------------------------------------
#im Folgenden extrahieren wir aus der heatmap nur den THCA Datensatz um eventuell Subtypen in 
#unserem Tumor zu finden
#Außerdem PCA und UMAP und dann beide als heatmap
#-----------------------------------------

library(ggplot2)
library(metaplot)
library(gridExtra)
library(umap)
library(ComplexHeatmap)
library(cluster)

#install.packages("umap")
load('data/tcga_gsva.RData')
load('data/tcga_anno.RData')
set.seed(123)


#----------------------------
#slicen der Daten
#----------------------------

thca_gsva = tcga_gsva[,tcga_anno$cancer_type_abbreviation == 'THCA']
thca_pan_anno = tcga_anno[tcga_anno$cancer_type_abbreviation == 'THCA',]

#Umbennen der Histologischen stufe (entfernen der unnötigen vorsatzes THCA)
thca_pan_anno$histological_type = sapply(thca_pan_anno$histological_type, FUN = function(x){
  return(strsplit(x, split = '- ', fixed = TRUE)[[1]][2])})

#Alle NAs werden other genannt
thca_pan_anno$histological_type[is.na(thca_pan_anno$histological_type)] = 'other'


#---------------------------
#PCA
#---------------------------

PCA_thca = prcomp(t(thca_gsva))
PCA_data_thca = as.data.frame(PCA_thca$x)
PCA_data_thca$hist = thca_pan_anno$histological_type


#---------------------------
#UMAP
#---------------------------

UMAP_thca = umap(PCA_thca$x)
UMAP_data_thca = as.data.frame(UMAP_thca$layout)  #UMAP der PCA-Daten
UMAP_data_thca$hist = thca_pan_anno$histological_type

#In der großen pancancer Daten sehen wir drei Cluster
km_data = kmeans(t(thca_gsva), centers = 3)
cluster = sapply(km_data$cluster, FUN = function(x){
  if (x == 1) return('A')
  else if (x == 2) return('B')
       else return('C')
})
PCA_data_thca$Cluster = cluster
UMAP_data_thca$Cluster = cluster

#Wir sehen auch drei histologisch unterschiedliche Subtypen in den Annotationen

#Plottet unsere PCA
ggplot(PCA_data_thca, aes(x = PC1, y = PC2, col = Cluster, shape = hist)) + 
  geom_point(size = 2) +
  labs(title = 'PCA of THCA cancer pathway activity') +
  theme_minimal()

#Plottet unsere UMAP
ggplot(UMAP_data_thca, aes(x = V1, y = V2, col = Cluster, shape = hist)) + 
  geom_point(size = 1) +
  labs(title = 'UMAP of THCA cancer pathway activity') +
  theme_minimal()

#Plotten der Heatmap mit Patienten in Drei Clustern
#Annotationen für die Pathways                              
path.anno = rowAnnotation(Pathway = c(rep('met',612),rep('hall', 46)),
                          col = list(Pathway = c("met" = "deepskyblue", "hall" = "blue4")),
                          annotation_legend_param = list(title = 'Pathway type',
                                                         at = c('met', 'hall'),
                                                         labels = c('Metabolic', 'Hallmark')))
#Annotation der k-means cluster und des histological types
top.anno = HeatmapAnnotation(Kmeans = cluster,
                             Type = thca_pan_anno$histological_type,
                             col = list(Kmeans = c('A' = 'darkgoldenrod1', 'B' = 'khaki1', 'C' = 'khaki3'),
                                        Type = c('Classical/usual' = 'olivedrab4',
                                                 'Follicular (>= 99% follicular patterned)' = 'darkolivegreen',
                                                 'Tall Cell (>= 50% tall cell features)' = 'green3',
                                                 'other' = 'honeydew3')),
                             annotation_legend_param = list(
                               list(title = 'Kmeans Cluster',
                                    at = c('A', 'B', 'C'),
                                    labels = c('Cluster A','Cluster B', 'Cluster C')),
                               list(title = 'Histological type',
                                    at = names(table(thca_pan_anno$histological_type)),
                                    lables = names(table(thca_pan_anno$histological_type)))
                               ))
Heatmap(thca_gsva,
        show_row_names = F, show_column_names = F, width = unit(21, 'cm'), height = unit(16, 'cm'),
        heatmap_legend_param = list(
          title = "Thyroid cancer pathway activity", at = c(-2, 2), 
          labels = c("underexpressed", "overexpressed")),
        row_dend_reorder = T, column_dend_reorder = T,
        left_annotation = path.anno,
        top_annotation = top.anno,
        column_split = km_data$cluster #Splitted unserer Heatmap in die 2 Subtypen
)
