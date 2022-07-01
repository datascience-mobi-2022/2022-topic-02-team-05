#In diesem Dokument werden wir einen GSVA für unserer Pancancer Analyse
#durchführen um diese Daten dann mit der Z-Transformierten GSEA zu vergleichen
#--------------------------------------------------------
#--------------------------------------------------------

#BiocManager::install("GSVA")
#BiocManager::install("ComplexHeatmap")
library(GSVA)
library(ComplexHeatmap)
library(pheatmap)

load('data/tcga_exp_cleaned.RData')     #Exp Daten
load('data/pathways.RData')             #pathways
tcga_anno = readRDS('~/GitHub/2022-topic-02-team-05/data/tcga_tumor_annotation.RDS')

tcga_gsva = gsva(as.matrix(tcga_exp_cleaned), pathways,
                 method = 'gsva',
                 kcdf = 'Gaussian'  #Da wir kontinuierliche Daten haben
                 )
#darstellen der GSVA matrix als heatmap
# pheatmap(tcga_gsva,
#          breaks = seq(-max(tcga_gsva), max(tcga_gsva), length.out = 201),
#          color = colorRampPalette(c('blue','lightskyblue','lightblue1', 'white','lightyellow','yellow', 'red'),
#                                   bias = 1,
#                                   space = 'rgb',
#                                   interpolate = 'linear'
#          )(200),
#          clustering_method = 'complete',
#          treeheight_row = 0, treeheight_col = 0, cellwidth = 0.07, cellheight = 0.5,
#          show_colnames = F, show_rownames = F, border_color = NA,
#          legend_breaks = c(-max(tcga_gsva),0, max(tcga_gsva)),
#          legend_labels = c('underexpressed', 'normal expression', 'overexpressed')
# )

save(tcga_gsva, file = 'data/tcga_gsva.RData')

#-------------------------------------
#Plotten der GSVA pathway activity mit complex Heatmap
#-------------------------------------
#Farben für Krebstypen definieren
colours = c('blue4','dodgerblue4','deepskyblue','cyan',
            'lightblue1','aquamarine','chartreuse4',
            'aquamarine3','darkolivegreen','darkolivegreen4',
            'chartreuse2','darkolivegreen2','lemonchiffon',
            'yellow','peachpuff1','gold','orange','red',
            'indianred3','orangered4','sienna4','tan3','salmon2',
            'plum','rosybrown1','violetred','magenta','magenta4',
            'maroon','wheat1','snow3','gray29','black')
names(colours) = names(table(tcga_anno$cancer_type_abbreviation))

#Annotationen für die Krebstypen
patient.anno = HeatmapAnnotation(Cancer = tcga_anno$cancer_type_abbreviation,
                                 col = list(Cancer = colours),
                                 annotation_legend_param = list(title = 'Cancer type',
                                   at = names(table(tcga_anno$cancer_type_abbreviation)),
                                   labels = names(table(tcga_anno$cancer_type_abbreviation))))
#Annotationen für die Pathways                              
path.anno = rowAnnotation(Pathway = c(rep('met',612),rep('hall', 46)),
                            col = list(Pathway = c("met" = "deepskyblue", "hall" = "blue4")),
                            annotation_legend_param = list(title = 'Pathway type',
                                                           at = c('met', 'hall'),
                                                           labels = c('Metabolic', 'Hallmark')))
#Plotten der Heatmap
Heatmap(tcga_gsva,
        show_row_names = F, show_column_names = F, width = unit(25, 'cm'), height = unit(18, 'cm'),
        heatmap_legend_param = list(
          title = "TCGA pathway activity", at = c(-2, 2), 
          labels = c("underexpressed", "overexpressed")),
        row_dend_reorder = T, column_dend_reorder = T,
        top_annotation = patient.anno,
        left_annotation = path.anno
)

#--------------------------
#Diese letzte Heatmap vergleicht nun die Durchschnittlichen expressionwerte für
#einen Krebstyp mit dem Rest
#--------------------------

#Liste mit Df mit GSVA daten für jeden Krebs
cancers_gsva = list(); cancers_gsva = vector('list',length(table(tcga_anno$cancer_type_abbreviation)))
names(cancers_gsva) = names(table(tcga_anno$cancer_type_abbreviation))
i=1; for (i in 1:length(cancers_gsva)){
  cancers_gsva[[i]] = tcga_gsva[,tcga_anno$cancer_type_abbreviation == names(cancers_gsva)[i]]
};rm(i)
cancers_gsva_means = lapply(cancers_gsva, FUN = function(x){
  apply(x, 1, median)
})
#Matrix mit den Mittlewerten der Pathway activity für jden Krebs
means_data = matrix(unlist(cancers_gsva_means), ncol = 33, dimnames = list(
  names(pathways), names(table(tcga_anno$cancer_type_abbreviation))
))

#Plotten der Heatmap
Heatmap(means_data,
        show_row_names = F, show_column_names = T, width = unit(20, 'cm'), height = unit(18, 'cm'),
        heatmap_legend_param = list(
          title = "Mean cancer pathway activity", at = c(-2, 2), 
          labels = c("underexpressed", "overexpressed")),
        row_dend_reorder = T, column_dend_reorder = T,
        left_annotation = path.anno
)
