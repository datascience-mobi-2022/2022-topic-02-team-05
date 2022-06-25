#In diesem Dokument werden wir einen GSVA für unserer Pancancer Analyse
#durchführen um diese Daten dann mit der Z-Transformierten GSEA zu vergleichen
#--------------------------------------------------------
#--------------------------------------------------------

#BiocManager::install("GSVA")
library(GSVA)
library(pheatmap)

load('data/tcga_exp_cleaned.RData')     #Exp Daten
load('data/our_genesets_final.RData')   #MSigDB pathways
load('data/geneset_ids.RData')          #Hallmarkpathways

pathways = c(our_genesets_final, genesets_ids)
rm(our_genesets_final); rm(genesets_ids)

tcga_gsva = gsva(as.matrix(tcga_exp_cleaned), pathways,
                 method = 'gsva',
                 kcdf = 'Gaussian'  #Da wir kontinuierliche Daten haben
                 )
#darstellen der GSVA matrix als heatmap
pheatmap(tcga_gsva,
         breaks = seq(-max(tcga_gsva), max(tcga_gsva), length.out = 201),
         color = colorRampPalette(c('blue','lightskyblue','lightblue1', 'white','lightyellow','yellow', 'red'),
                                  bias = 1,
                                  space = 'rgb',
                                  interpolate = 'linear'
         )(200),
         clustering_method = 'complete',
         treeheight_row = 0, treeheight_col = 0, cellwidth = 0.25, cellheight = 0.5,
         show_colnames = F, show_rownames = F, border_color = NA,
         legend_breaks = c(-max(tcga_gsva),0, max(tcga_gsva)),
         legend_labels = c('underexpressed', 'normal expression', 'overexpressed')
)

save(tcga_gsva, file = 'data/tcga_gsva.RData')
