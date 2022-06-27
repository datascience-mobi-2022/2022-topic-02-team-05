#-----------------------------
#DIESES DOKUMENT ERSETZT DIE 01 METABOLIC PATHWAY SELECTION
# Hier nehemen wir alle canonischen pathways von msigdbr und überprüfen diese auf kompatibilität mit unseren daten
# Sollten pathways zu 90% in unseren gecleanten genen enthalen sein werden sie weiter verwendet
#-----------------------------
#install.packages("msigdbr")
#if (!require("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install(version = "3.14")
library(msigdbr)

#---------------------------
#Laden der MSigDBr Pathways und Formatierung
#---------------------------

#oeffnen der Datenbank
msigdbr_pathways = msigdbr(species = "Homo sapiens", category = "C2")
pathway_list = list(); pathway_list = vector('list',length(table(msigdbr_pathways$gs_name)))
names(pathway_list) = names(table(msigdbr_pathways$gs_name))
i=1; for (i in 1:length(pathway_list)){
  pathway_list[[i]] = msigdbr_pathways[msigdbr_pathways$gs_name == names(pathway_list)[i],6]
};rm(i)

#neuformatieren in die From der hallmarkpathways 
pathway_list = lapply(pathway_list, FUN = function(x){res = unlist(unname(x))
  return(res)})
rm(msigdbr_pathways)

#----------------------------
#Vergleichen zu wie viele der Genen in diesen Canonical Pathways letzendlich in den gecleanten
#tcga genen enthalten sind und gibt einen Prozentsatz aus
#----------------------------

load('data/tcga_genes_cleaned.RData')
coverage = sapply(pathway_list, FUN = function(x){
              enthalten = x %in% tcga_genes_cleaned$tcga_geneids
              res = sum(enthalten)/length(enthalten)
              return(res)
            })
hist(coverage, breaks = 100)
#Aus dem histogramm schließen wir dass wir nur pathways mit 99% abdeckung haben wollen
our_genesets = pathway_list[which(unname(coverage >= 0.99))]

#-----------------------------
#Vergleichen der Pathways mit Jaccard untereinander um dopplungen zu vermeiden
#-----------------------------

jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}
# leeren Data-frame erstellen:
data.jaccard = data.frame(matrix(ncol=length(our_genesets), nrow = length(our_genesets)),
                          row.names = names(our_genesets))
colnames(data.jaccard) <- names(our_genesets) #n
#jaccard_index:
i = 1; j = 1; for (j in 1:length(our_genesets)){
  for (i in 1:length(our_genesets)){
    jacindex <- jaccard(our_genesets[[i]], our_genesets[[j]])
    data.jaccard[i,j] <- jacindex
    i = i+1
  }
  i = 1
  j = j+1
};rm(i);rm(j)
#Darstellen als Heatmap
library(pheatmap)
pheatmap(as.matrix(data.jaccard),
         breaks = seq(0, max(data.jaccard), length.out = 51),
         color = colorRampPalette(c('lightskyblue','lightcyan','white','yellow','orange', 'red'),
                                  bias = 1,
                                  space = 'rgb',
                                  interpolate = 'linear'
         )(50),
         clustering_method = 'average', treeheight_row = 0, treeheight_col = 0,
         cellwidth = 1, cellheight = 1,
         show_colnames = FALSE,show_rownames = FALSE, border_color = 'lightcyan2',
         legend_breaks = c(0, 1),
         legend_labels = c('unlike','same')
)
#Die Heatmap zeigt das einige Pathways sich überschneiden
#um diese doppelten Pathways zu entfernen vergleiche wir wie hoch die Summe aller
#Jaccardindices eines pathays ist (je höher desto mehr doppelt er sich)
#Anschließend wähle wir nur die aus die im 1s bereich liegen und entfernen die anderen
Sumjac = apply(data.jaccard, 1, sum)
selected = which(median(Sumjac)+0.5*sd(Sumjac) >= Sumjac)#alle pathways im 1s bereich

#Neue Heatmap zum Test obs was gebracht hat
data.jaccard.cleaned = data.jaccard[selected,selected]
pheatmap(as.matrix(data.jaccard.cleaned),
         breaks = seq(0, max(data.jaccard), length.out = 51),
         color = colorRampPalette(c('lightskyblue','lightcyan','white','yellow','orange', 'red'),
                                  bias = 1,
                                  space = 'rgb',
                                  interpolate = 'linear'
         )(50),
         clustering_method = 'average', treeheight_row = 0, treeheight_col = 0,
         cellwidth = 1, cellheight = 1,
         show_colnames = FALSE,show_rownames = FALSE, border_color = 'lightcyan2',
         legend_breaks = c(0, 1),
         legend_labels = c('unlike','same')
)
#Auswahl der Pathways die unterschieldich genug sind (ins. 874 stueck)
our_genesets_new = our_genesets[selected]

#-------------------------------------------
#Nun werden die Pathways noch mit den Hallmark pathways abgegelichen um dort dopplung zu vermeiden
#Dazu benutzen wir wieder einen Jaccard index
#-------------------------------------------

load('data/geneset_ids.RData')
# leeren Data-frame erstellen:
hallm.jaccard = data.frame(matrix(ncol=length(genesets_ids), nrow = length(our_genesets_new)),
                          row.names = names(our_genesets_new))
colnames(hallm.jaccard) <- names(genesets_ids)
#jaccard_index:
i = 1; j = 1; for (j in 1:length(genesets_ids)){
  for (i in 1:length(our_genesets_new)){
    jacindex <- jaccard(our_genesets_new[[i]], genesets_ids[[j]])
    hallm.jaccard[i,j] <- jacindex
    i = i+1
  }
  i = 1
  j = j+1
};rm(i);rm(j);rm(jacindex)
#Heatmap
pheatmap(as.matrix(hallm.jaccard),
         breaks = seq(0, max(hallm.jaccard), length.out = 51),
         color = colorRampPalette(c('lightskyblue','lightcyan','white','yellow','orange', 'red'),
                                  bias = 1,
                                  space = 'rgb',
                                  interpolate = 'linear'
         )(50),
         clustering_method = 'average', treeheight_row = 0, treeheight_col = 0,
         cellwidth = 8, cellheight = 0.5,
         show_colnames = TRUE,show_rownames = FALSE, border_color = 'lightcyan2',
         legend_breaks = c(0, 1),
         legend_labels = c('unlike','same')
)
#Auch hier gibt es wieder ähnlichkeiten die entfernt werden müssen
similarity = apply(hallm.jaccard, 1 , sum)
hist(similarity, breaks = 100)
#Auch hier werden wir wieder nur die pathways behalten deren Aehnlichkeit nicht groeßer
#als der 1s Berrecih ist
hallm.selected = which(median(similarity)+0.5*sd(similarity) >= similarity)#alle pathways im 1s bereich

#Überprüfeung mit Heatmap
pheatmap(as.matrix(hallm.jaccard[hallm.selected,]),
         breaks = seq(0, max(hallm.jaccard), length.out = 51),
         color = colorRampPalette(c('lightskyblue','lightcyan','white','yellow','orange', 'red'),
                                  bias = 1,
                                  space = 'rgb',
                                  interpolate = 'linear'
         )(50),
         clustering_method = 'average', treeheight_row = 0, treeheight_col = 0,
         cellwidth = 8, cellheight = 0.5,
         show_colnames = TRUE,show_rownames = FALSE, border_color = 'lightcyan2',
         legend_breaks = c(0, 1),
         legend_labels = c('unlike','same')
)
#Keine Dopplungen mehr vorhanden
#Pathways sind hiermit fertig selektiert
our_genesets_final = our_genesets_new[hallm.selected]
save(our_genesets_final, file = 'data/our_genesets_final.RData')

#Zusammenführen in ein Pathwayliste
pathways = c(our_genesets_final, genesets_ids)
save(pathways, file = 'data/pathways.RData')
