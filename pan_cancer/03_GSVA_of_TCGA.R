#--------------------------------------------------------
#In diesem Dokument werden wir eine GSVA für unserer Pancancer Analyse
#durchführen um diese Daten dann mit der z-transformierten GSEA zu vergleichen
#--------------------------------------------------------

#BiocManager::install("GSVA")
#BiocManager::install("ComplexHeatmap")
library(GSVA)
library(ComplexHeatmap)

load('data/tcga_exp_cleaned.RData')     #Exp Daten
load('data/pathways.RData')             #pathways
load('data/tcga_anno.RData')            #Annotationen für alle Daten
load('data/Cancer_form_anno.RData')     #Annotation für die Krebsarten
#load('data/tcga_gsva.RData') für den Fall das man den Code oben nicht nochmal ausführen möchte

tcga_gsva = gsva(as.matrix(tcga_exp_cleaned), pathways,
                 method = 'gsva',
                 kcdf = 'Gaussian'  #Da wir kontinuierliche Daten haben
                 )

save(tcga_gsva, file = 'data/tcga_gsva.RData')


#-------------------------------------
#Plotten der GSVA pathway activity mit complex Heatmap
#-------------------------------------

#Farben für Krebstypen definieren
patient.col = c('blue4','dodgerblue4','deepskyblue','cyan',
            'lightblue1','aquamarine','chartreuse4',
            'aquamarine3','darkolivegreen','darkolivegreen4',
            'chartreuse2','darkolivegreen2','lemonchiffon',
            'yellow','peachpuff1','gold','orange','red',
            'indianred3','orangered4','sienna4','tan3','salmon2',
            'plum','rosybrown1','violetred','magenta','magenta4',
            'maroon','wheat1','snow3','gray29','black')
names(patient.col) = names(table(tcga_anno$cancer_type_abbreviation))

#Farben für die Art des Krebses definieren
form.col = c('Adenocarcinoma' = 'deepskyblue3',
             'Squamous cell carcinoma' = 'dodgerblue4',
             'Transitional cell carcinoma' = 'deepskyblue',
             'Melanoma'='cyan',
             'Carcinoma'='azure3',
             'Sarcoma'='darkolivegreen',
             'Glioblastoma'='gold',
             'Leukemia'='red')

#Annotationen für die Krebstypen
top.anno = HeatmapAnnotation(Type = tcga_anno$cancer_type_abbreviation, #Krebsabkürzung als farbe
                             Form = tcga_anno$cancer_form, #Krebsform als farbe
                                 col = list(Type = patient.col, #Def der Farben
                                            Form = form.col),
                                 annotation_legend_param = list( #Def der Legende
                                   Type = list(nrow = 2,
                                        title = 'Cancer type',
                                        at = names(table(tcga_anno$cancer_type_abbreviation)),
                                        labels = names(table(tcga_anno$cancer_type_abbreviation))),
                                   Form = list(nrow = 1,
                                        title = 'Cancer form',
                                        at = names(table(tcga_anno$cancer_form)),
                                        labels = names(table(tcga_anno$cancer_form)))
                                   ))

#Annotationen für die Pathways                              
path.anno = rowAnnotation(Pathway = c(rep('met',612),rep('hall', 46)),
                            col = list(Pathway = c("met" = "deepskyblue", "hall" = "blue4")),
                            annotation_legend_param = list(title = 'Pathway type',
                              at = c('met', 'hall'),
                              labels = c('Metabolic', 'Hallmark')))
#Plotten der Heatmap
Heatmap(tcga_gsva,
        show_row_names = F, show_column_names = F, width = unit(22, 'cm'), height = unit(18, 'cm'),
        heatmap_legend_param = list(
          title = "TCGA pathway activity", at = c(-2, 2), 
          labels = c("underexpressed", "overexpressed")),
        row_dend_reorder = T, column_dend_reorder = T,
        top_annotation = top.anno,
        left_annotation = path.anno
)

#-------------------------------------
#Plotten der GSVA pathway activity mit complex Heatmap 
#allerdings nur mit den 100 höchstvarainten pathways
#-------------------------------------
path.var = apply(tcga_gsva, 1, var)
tcga_gsva_topvar = tcga_gsva[order(path.var, decreasing = T)[1:100],]
hm2 = Heatmap(tcga_gsva_topvar,
        show_row_names = T, show_column_names = F, width = unit(20, 'cm'), height = unit(15, 'cm'),
        row_names_gp = gpar(fontsize = 4),
        heatmap_legend_param = list(
          direction = "horizontal",
          title = "TCGA pathway activity", at = c(-2, 2), 
          labels = c("underexpressed", "overexpressed")),
        row_dend_reorder = T, column_dend_reorder = T,
        top_annotation = top.anno
)
draw(hm2,merge_legend = T, heatmap_legend_side = "bottom",
     annotation_legend_side = "bottom")
#--------------------------
#Diese letzte heatmap vergleicht nun die durchschnittlichen Expressionwerte für
#einen Krebstyp mit dem Rest
#--------------------------

#Liste mit dataframe mit GSVA daten für jeden Krebs
cancers_gsva = list(); cancers_gsva = vector('list',length(table(tcga_anno$cancer_type_abbreviation)))
names(cancers_gsva) = names(table(tcga_anno$cancer_type_abbreviation))
i=1; for (i in 1:length(cancers_gsva)){
  cancers_gsva[[i]] = tcga_gsva[,tcga_anno$cancer_type_abbreviation == names(cancers_gsva)[i]]
};rm(i)
cancers_gsva_means = lapply(cancers_gsva, FUN = function(x){
  apply(x, 1, median)
})
#Matrix mit den Mittelwerten der Pathwayaktivität für jeden Krebs
means_data = matrix(unlist(cancers_gsva_means), ncol = 33, dimnames = list(
  names(pathways), names(table(tcga_anno$cancer_type_abbreviation))
))
rm(cancers_gsva);rm(cancers_gsva_means) #Environment aufräumen

#Neue Annotation für die Krebsform definieren
form.anno = HeatmapAnnotation(Type = Cancer_form_anno$Cancer_form,
                              col = list(Type = form.col),
                              annotation_legend_param = list(title = 'Cancer type',
                                nrow = 1, 
                                at = names(table(Cancer_form_anno$Cancer_form)),
                                labels = names(table(Cancer_form_anno$Cancer_form))))

#Clustern der Daten in vier Cluster
cluster_data = hclust(dist(t(means_data)), method = 'complete')

#Plotten der Heatmap
Heatmap(means_data,
        show_row_names = F, show_column_names = T, width = unit(20, 'cm'), height = unit(18, 'cm'),
        heatmap_legend_param = list(
          title = "Mean cancer pathway activity", at = c(-2, 2), 
          labels = c("underexpressed", "overexpressed")),
        row_dend_reorder = T, column_dend_reorder = T,
        left_annotation = path.anno,
        top_annotation = form.anno,
        column_split = cutree(cluster_data, k = 4 )
)

#--------------------------
#Diese heatmap vergleicht nun die durchschnittlichen Expressionwerte für
#einen Krebstyp mit dem Rest aber ernuet nehemn wir dann hierfür nur die 
#100 Variantesen pathways
#--------------------------

means_data_topvar = means_data[order(path.var, decreasing = T)[1:100],]
hm4 = Heatmap(means_data_topvar,
        show_row_names = T, show_column_names = T, width = unit(22, 'cm'), height = unit(15, 'cm'),
        row_names_gp = gpar(fontsize = 4),
        heatmap_legend_param = list(
          direction = 'horizontal',
          title = "Mean cancer pathway activity", at = c(-2, 2), 
          labels = c("underexpressed", "overexpressed")),
        row_dend_reorder = T, column_dend_reorder = T,
        #left_annotation = path.anno,
        top_annotation = form.anno,
        column_split = cutree(cluster_data, k = 4 )
)
draw(hm4,merge_legend = T, heatmap_legend_side = "bottom",
     annotation_legend_side = "bottom")
