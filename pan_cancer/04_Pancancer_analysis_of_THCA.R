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
load("data/pathways.RData")
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

#---------------------------------------------------------------
#Plotten der Heatmap mit Patienten in Drei Clustern
#----------------------------------------------------------------

#markieren der Hallmark pathways             
path.type = c(rep('other',612),rep('hall', 46))
#Markieren der pathways für beta 4 intergrin abh carcinogenese
path.type[match(c("REACTOME_INTERLEUKIN_36_PATHWAY",
                  "REACTOME_TYPE_I_HEMIDESMOSOME_ASSEMBLY",
                  "PID_INTEGRIN4_PATHWAY"),
          names(pathways))] = 'beta4'
#Markieren der sig. geänderteten DNA Reperatur und instabilität
path.type[match(c("KAUFFMANN_MELANOMA_RELAPSE_DN",
                  "SESTO_RESPONSE_TO_UV_C4"),
                names(pathways))] = 'CIN'
#Markieren von THCA assozierten TSGs und Oncogenes plus prolif signaling
path.type[match(c("RAMPON_ENRICHED_LEARNING_ENVIRONMENT_EARLY_UP",
          "ROVERSI_GLIOMA_LOH_REGIONS",
          "REACTOME_GLI_PROTEINS_BIND_PROMOTERS_OF_HH_RESPONSIVE_GENES_TO_PROMOTE_TRANSCRIPTION",
          "REACTOME_PROTEIN_METHYLATION",    
          "JONES_TCOF1_TARGETS",  
          "TORCHIA_TARGETS_OF_EWSR1_FLI1_FUSION_TOP20_DN",
          "FARMER_BREAST_CANCER_CLUSTER_8",      
          "REACTOME_SIGNALING_BY_MST1",   
          "OKAWA_NEUROBLASTOMA_1P36_31_DELETION"),
          names(pathways))] = 'prolif' 
#Markieren der sig. geänderteten metabolischen pathways
path.type[match(c("REACTOME_ETHANOL_OXIDATION",            
                  "REACTOME_ABACAVIR_METABOLISM",                     
                  "REACTOME_SYNTHESIS_OF_BILE_ACIDS_AND_BILE_SALTS_VIA_27_HYDROXYCHOLESTEROL",  
                  "REACTOME_PYRIMIDINE_CATABOLISM"),
                names(pathways))] = 'met'

#Annotationen für die Pathways
path.anno = rowAnnotation(Pathway = path.type,
              col = list(Pathway = c("other" = "deepskyblue",
                                     "hall" = "blue4",
                                     "beta4" = 'gold',
                                     'CIN' = 'green',
                                     'prolif' = 'red',
                                     'met' = 'black'
                                     )),
              annotation_legend_param = list(title = 'Pathway type',
                 at = c('other', 'hall', 'beta4', 'CIN', 'prolif', 'met'),
                 labels = c('other', 'Hallmark', 'Integrin carcinogenesis', 
                            'genomic instability', 'altered proliferative signaling',
                            'altered metabolism')))
#Annotation der k-means cluster und des histological types
top.anno = HeatmapAnnotation(Kmeans = cluster,
             Type = thca_pan_anno$histological_type,
             col = list(Kmeans = c('A' = 'darkgoldenrod1', 'B' = 'khaki1', 'C' = 'khaki3'),
                        Type = c('Classical/usual' = 'olivedrab4',
                                 'Follicular (>= 99% follicular patterned)' = 'darkolivegreen',
                                 'Tall Cell (>= 50% tall cell features)' = 'green3',
                                 'other' = 'honeydew3')),
             annotation_legend_param = list(
               list(title = 'Kmeans Cluster',nrow = 1,
                    at = c('A', 'B', 'C'),
                    labels = c('Cluster A','Cluster B', 'Cluster C')),
               list(title = 'Histological type',nrow = 1,
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

#---------------------------------------------------------------
#Plotten der Heatmap mit Patienten in Drei Clustern
#mit nur den top 50 var genen plus hallmarks und signifikant unterschieliche 
#aus focussed GSVA tum vs normal
#----------------------------------------------------------------
path.var = apply(thca_gsva[which(path.type == 'other'),], 1, var)
thca_gsva_topvar = thca_gsva[order(path.var, decreasing = T)[1:50],]
thca_gsva_topvar = rbind(thca_gsva[which(!path.type == 'other'),], thca_gsva_topvar)

#markieren der Hallmark pathways    
load('data/geneset_ids.RData')
path.top = rep('other',114)
path.top[match(names(genesets_ids),rownames(thca_gsva_topvar))] = 'hall'
#Markieren der pathways für beta 4 intergrin abh carcinogenese
path.top[match(c("REACTOME_INTERLEUKIN_36_PATHWAY",
                  "REACTOME_TYPE_I_HEMIDESMOSOME_ASSEMBLY",
                  "PID_INTEGRIN4_PATHWAY"),
               rownames(thca_gsva_topvar))] = 'beta4'
#Markieren der sig. geänderteten DNA Reperatur und instabilität
path.top[match(c("KAUFFMANN_MELANOMA_RELAPSE_DN",
                  "SESTO_RESPONSE_TO_UV_C4"),
               rownames(thca_gsva_topvar))] = 'CIN'
#Markieren von THCA assozierten TSGs und Oncogenes plus prolif signaling
path.top[match(c("RAMPON_ENRICHED_LEARNING_ENVIRONMENT_EARLY_UP",
                  "ROVERSI_GLIOMA_LOH_REGIONS",
                  "REACTOME_GLI_PROTEINS_BIND_PROMOTERS_OF_HH_RESPONSIVE_GENES_TO_PROMOTE_TRANSCRIPTION",
                  "REACTOME_PROTEIN_METHYLATION",    
                  "JONES_TCOF1_TARGETS",  
                  "TORCHIA_TARGETS_OF_EWSR1_FLI1_FUSION_TOP20_DN",
                  "FARMER_BREAST_CANCER_CLUSTER_8",      
                  "REACTOME_SIGNALING_BY_MST1",   
                  "OKAWA_NEUROBLASTOMA_1P36_31_DELETION"),
               rownames(thca_gsva_topvar))] = 'prolif' 
#Markieren der sig. geänderteten metabolischen pathways
path.top[match(c("REACTOME_ETHANOL_OXIDATION",            
                  "REACTOME_ABACAVIR_METABOLISM",                     
                  "REACTOME_SYNTHESIS_OF_BILE_ACIDS_AND_BILE_SALTS_VIA_27_HYDROXYCHOLESTEROL",  
                  "REACTOME_PYRIMIDINE_CATABOLISM"),
               rownames(thca_gsva_topvar))] = 'met'

#Annotationen für die Pathways
path.anno = rowAnnotation(Pathway = path.top,
                col = list(Pathway = c("other" = "deepskyblue",
                                       "hall" = "blue3",
                                       "beta4" = 'gold',
                                       'CIN' = 'green',
                                       'prolif' = 'red',
                                       'met' = 'black'
                )),
                annotation_legend_param = list(title = 'Pathway type', nrow = 1,
                                               at = c('other', 'hall', 'beta4', 'CIN', 'prolif', 'met'),
                                               labels = c('other', 'Hallmark', 'Integrin carcinogenesis', 
                                                          'genomic instability', 'altered proliferative signaling',
                                                          'altered metabolism')))

hm = Heatmap(thca_gsva_topvar,
        show_row_names = T, show_column_names = F, width = unit(22, 'cm'), height = unit(14, 'cm'),
        row_names_gp = gpar(fontsize = 4),
        heatmap_legend_param = list(direction = 'horizontal',
          title = "Thyroid cancer pathway activity", at = c(-2, 2), 
          labels = c("underexpressed", "overexpressed")),
        row_dend_reorder = T, column_dend_reorder = T,
        left_annotation = path.anno,
        top_annotation = top.anno,
        column_split = km_data$cluster #Splitted unserer Heatmap in die 2 Subtypen
)
draw(hm,merge_legend = T, heatmap_legend_side = "bottom",
     annotation_legend_side = "bottom")
