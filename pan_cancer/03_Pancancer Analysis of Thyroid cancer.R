#-----------------------------------------
#im folgenden extrahieren wir aus der heatmap nur den THCA Datensatz um eventuell subtypen in 
#unserem tuor zu finden dazu wie sonst auch PCA UMAP und dann Heatmap
#-----------------------------------------
library(ggplot2)
library(metaplot)
library(gridExtra)
library(umap)
library(ComplexHeatmap)

#install.packages("umap")
load('data/tcga_gsva.RData')
tcga_anno = readRDS('data/tcga_tumor_annotation.RDS')
set.seed(123)

#slicen der Daten
thca_gsva = tcga_gsva[,tcga_anno$cancer_type_abbreviation == 'THCA']

#PCA
PCA_thca = prcomp(t(thca_gsva))
PCA_data_thca = as.data.frame(PCA_thca$x)

#UMAP
UMAP_thca = umap(PCA_thca$x)
UMAP_data_thca = as.data.frame(UMAP_thca$layout)  #UMAP der PCA-Daten

#----------------------------------
#kmeans clustering der Daten
#----------------------------------
#Beste clustering form:
set.seed(123)
elbow_data = c()#leerer vektor für die Daten
for (i in 1:10){
  res = kmeans(UMAP_data_thca, centers = i)$tot.withinss
  elbow_data = append(elbow_data, res)
}; rm(i);rm(res)
plot(elbow_data)
#=> Der elbowplot zeigt, 3 Cluster sind die beste Anzahl für die Daten
#Kmeans clustern mit 3 Zentren und zuordnug der Punkte
km_data = kmeans(UMAP_data_thca, centers = 3)
PCA_data_thca$Cluster = km_data$cluster
UMAP_data_thca$Cluster = km_data$cluster

#Plottet unserer PCA
ggplot(PCA_data_thca, aes(x = PC1, y = PC2, col = Cluster)) + geom_point(size = 2) +
  labs(title = 'PCA of THCA cacncer pathway activity') +
  theme_minimal()
#Plottet unserer UMAP
ggplot(UMAP_data_thca, aes(x = V1, y = V2, col = Cluster)) + geom_point(size = 1) +
  labs(title = 'UMAP of THCA cacncer pathway activity') +
  theme_minimal()

#Plotten der Heatmap mit patienten in Drei Clustern
#Annotationen für die Pathways                              
path.anno = rowAnnotation(Pathway = c(rep('met',612),rep('hall', 46)),
                          col = list(Pathway = c("met" = "deepskyblue", "hall" = "blue4")),
                          annotation_legend_param = list(title = 'Pathway type',
                                                         at = c('met', 'hall'),
                                                         labels = c('Metabolic', 'Hallmark')))
#Annotation der kmeans cluster
patient.anno = HeatmapAnnotation(Kmeans = km_data$cluster,
                                 annotation_legend_param = list(title = 'Kmeans Cluster',
                                                                at = c(1, 2, 3),
                                                                labels = c('Cluster 1','Cluster 2','Cluster 3')))
Heatmap(thca_gsva,
        show_row_names = F, show_column_names = F, width = unit(25, 'cm'), height = unit(18, 'cm'),
        heatmap_legend_param = list(
          title = "Thyroid cancer pathway activity", at = c(-2, 2), 
          labels = c("underexpressed", "overexpressed")),
        row_dend_reorder = T, column_dend_reorder = T,
        left_annotation = path.anno,
        top_annotation = patient.anno,
        column_split = km_data$cluster #Splitted unserer Heatmap in die drei Subtypen
)
