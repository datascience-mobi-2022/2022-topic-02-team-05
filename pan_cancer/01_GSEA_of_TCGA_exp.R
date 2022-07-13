#--------------------------------------------------
#In diesem Dokument versuchen wir die gereinigten Expressionsdaten mittels GSEA zu
#Pathway-Aktivitätswerten zusammenzufassen und in einer Matrix mit patienten als Spalten und 
#enrichment scores als zeilen aufzutragen
#--------------------------------------------------

#if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#BiocManager::install("gage")
#BiocManager::install("fgsea")
library(dplyr)
library(gage)
library(fgsea)
library(pheatmap)

load('data/tcga_exp_cleaned.RData')
tcga_anno = readRDS('data/tcga_tumor_annotation.RDS')
load('data/geneset_ids.RData')
load('data/our_genesets.RData')


#-----------------------------------------------
#Extrahieren aller Patienten mit dem gleichem Krebstyp in einen Dataframe
#-----------------------------------------------

cancers_exp = list(); cancers_exp = vector('list',length(table(tcga_anno$cancer_type_abbreviation)))
names(cancers_exp) = names(table(tcga_anno$cancer_type_abbreviation))
i=1; for (i in 1:length(cancers_exp)){
  cancers_exp[[i]] = tcga_exp_cleaned[,tcga_anno$cancer_type_abbreviation == names(cancers_exp)[i]]
}

#------------------------------------------------
#Durchführung einer einfachen Enrichment Analysis auf Basis von Wangs Vorschlag
#------------------------------------------------

#Funktion, die den Krebstypen dataframe und genesets als input nimmt und einen dataframe mit p-Werten ausgibt
enrichment = function(expressiondata, genesets = genesets_ids){
  ESmatrix = sapply(genesets, FUN = function(x){
    ins = na.omit(match(x,rownames(expressiondata)))#indices der Gene im aktuellen set
    outs = -ins#indices der ene nicht im aktuellen set
    
    #gibt einen Vektor, der für jeden Patienten den p-Wert für das aktuelle Gen enthält
    res = NULL
    for (i in 1:ncol(expressiondata)){#testet für jeden patienten
      res[i] = wilcox.test(expressiondata[ins,i],expressiondata[outs,i],'two.sided')$p.value
    }
    return(res)
  })
  row.names(ESmatrix) = colnames(expressiondata); return(ESmatrix)
}

#Anwendung der Analyse auf alle dataframes der einzelen Krebstypen 
enrichment_matrix = lapply(cancers_exp, FUN = function(x){
                      return(enrichment(x, c(geneset_ids, our_genesets)))
                    })
save(enrichment_matrix, file = 'data/GSEA_matrix.RData')


#-------------------------------------------
#Z-Transformation der Expressionsdaten für die GSEA für jedes Gen in einem Krebs über alle Patienten
#-------------------------------------------

cancers_exp_scaled = lapply(cancers_exp, FUN = function(x){
                        res = apply(x, 1, scale)
                        res[is.nan(res)] = 0 # Setzt alle Gene null die eine sd von nul haben, und deswegen bei scalen unendlich werden
                        return(t(res))
                      })
i=1; for (i in 1:length(cancers_exp_scaled)){
  colnames(cancers_exp_scaled[[i]]) = colnames(cancers_exp[[i]])}


#--------------------------------------------
#Sortieren der Gene der einzelnen Patienten nach den z-transformierten Daten
#jeder Krebstyp wird als Liste aller Patienten ausgegeben mit den sortierten und benannten Expressionsdaten
#--------------------------------------------

sorting = function(expressiondata, scaleddata){
  liste = list()
  liste = vector('list', ncol(scaleddata))
  names(liste) = colnames(scaleddata)
  i=1; for (i in 1:length(liste)){
    ord = scaleddata[,i]
    ord = order(abs(ord), decreasing = TRUE)
    res = expressiondata[,i]
    names(res) = rownames(expressiondata)
    patient_sorted = res[ord]
    liste[[i]] = patient_sorted
  }
  return(liste)
}  
#ausführen des sortings für alle Krebstypen
cancers_sorted = list(); cancers_sorted = vector('list',length(cancers_exp))
names(cancers_sorted) = names(cancers_exp)
i=1; for (i in 1:length(cancers_sorted)){
  cancers_sorted[[i]] = sorting(cancers_exp[[i]], cancers_exp_scaled[[i]])
}
save(cancers_sorted, file = 'data/cancers_sorted.RData')


#----------------------
#Funktion die unsere GSEA durchführt
#Input sind die Expressionsdaten der sortierten Patienten,
#sowie alle zu analysierenden pathways als Liste
#------------------------

GSEA = function(patientsorted, pathways = pathways){
  res = fgseaMultilevel(pathways = pathways, 
                        stats = patientsorted,
                        minSize=3)
  #extrahieren der NES werte & nullsetzten aller NAs
  ret = res$NES; names(ret) = res$pathway
  ret[is.na(ret)] = 0
  
  message('I?m still standing')
  return(ret)
}


#-------------------------------------------
#Durchführung der GSEA
#-------------------------------------------

load('data/cancers_sorted.RData')
load('data/our_genesets.RData')
load('data/geneset_ids.RData')
pathways = c(our_genesets, genesets_ids);rm(genesets_ids);rm(our_genesets)

#Durchführung der GSEA für den ersten krebstyp
GSEA_ACC = sapply(cancers_sorted[['ACC']], FUN = function(x){
  return(GSEA(x, pathways))})

save(GSEA_ACC, file = 'data/GSEA/GSEA_ACC.RData')

#darstellen der GSEA Matrix als heatmap
pheatmap(as.matrix(GSEA_matrix),
         breaks = seq(-max(GSEA_matrix), max(GSEA_matrix), length.out = 201),
         color = colorRampPalette(c('blue','lightskyblue','lightblue1', 'white','lightyellow','yellow', 'red'),
                                  bias = 1,
                                  space = 'rgb',
                                  interpolate = 'linear'
         )(200),
         clustering_method = 'average',
         treeheight_row = 25, treeheight_col = 20, cellwidth = 7,cellheight = 8,
         show_colnames = FALSE, fontsize_row = 8, border_color = NA,
         legend_breaks = c(-max(GSEA_matrix),0, max(GSEA_matrix)),
         legend_labels = c('underexpressed', 'normal expression', 'overexpressed')
        )






  

