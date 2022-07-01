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


#Plotten mit complex Heatmap
patient.anno = HeatmapAnnotation(Cancer = tcga_anno$cancer_type_abbreviation[1:400],
                                 annotation_legend_param = list(title = 'Cancer type',
                                   at = names(table(tcga_anno$cancer_type_abbreviation[1:400])),
                                   labels = names(table(tcga_anno$cancer_type_abbreviation[1:400]))))
                              
path.anno = rowAnnotation(Pathway = c(rep('met',612),rep('hall', 46)),
                            col = list(Pathway = c("met" = "deepskyblue", "hall" = "blue4")),
                            annotation_legend_param = list(title = 'Pathway type',
                                                           at = c('met', 'hall'),
                                                           labels = c('Metabolic', 'Hallmark')))

Heatmap(tcga_gsva[, 1:400],
        show_row_names = F, show_column_names = F, width = unit(25, 'cm'), height = unit(18, 'cm'),
        heatmap_legend_param = list(
          title = "TCGA pathway activity", at = c(-2, 2), 
          labels = c("underexpressed", "overexpressed")),
        row_dend_reorder = T, column_dend_reorder = T,
        top_annotation = patient.anno,
        left_annotation = path.anno
)


