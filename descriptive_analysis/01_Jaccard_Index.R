#----------------------------------------------
#In diesem Dokument suchen wir Gemeinsamkeiten zwischen den Hallmarkgenen und unseren Pathways
#Verwendung des Jaccard-Indexes (Schnittmenge im Verh?ltnis zur Vereinigungsmeng)
#----------------------------------------------

#loading data:

load("data/our_genesets_final.RData")
load("data/geneset_ids.RData")


#----------------------------------------------
#Jaccard Index zwischen den Genen
#----------------------------------------------

#für jedes Listenelement aus our_genesets(Vektor mit Gennamen) schauen, ob es in Hallmark_geneset 1 drinnen ist
#wenn es drinnen ist: count+1 in Matrix; dann n?chstes Gen

# leeren Data-frame erstellen:

data.jaccard <- data.frame(matrix(ncol=length(our_genesets_final), nrow = length(genesets_ids)))

colnames(data.jaccard) <- names(our_genesets_final) #n
rownames(data.jaccard) <- names(genesets_ids) #m

#jaccard_index berechnen
n = 1
m = 1

# Jaccar-Fkt. definieren: zwei Vektoren a und b als input und gibt die Schnittmenge der Vektoren 
# (intersection) geteilt durch die Vereinigungsmenge (union) als output. 

jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}


#---------------------------------------------------------------------
#Jaccardindex für all unsere Pathways vs die Hallmarkpathways einzeln darstellen
#um die Ähnlichkeit der Pathways zu bestimmen
#---------------------------------------------------------------------

for (set_our in our_genesets_final){
  for (set_hm in genesets_ids){
    jacindex <- jaccard(set_our, set_hm)
    print(jacindex)
    data.jaccard[m,n] <- jacindex
    m <- m+1
  }
  m <- 1
  n <- n+1
}

# Darstellung des Jaccardindices (Zahl zwischen 0 und 1) mithilfe einer Heatmap
pheatmap(as.matrix(data.jaccard),
         color = colorRampPalette(c("lightskyblue", "lightcyan", "white", "yellow", "orange", "red"),
                                  bias = 1,
                                  space = "rgb",
                                  interpolate = "linear"
         ) (50), 
         clustering_method = "average", treeheight_row = 20, treeheight_col = 20, 
         labels_col = FALSE,  
         fontsize = 5, border_color = "black", 
         legend_breaks = c(0, 1), 
         legend_labels = c("unlike", "same"),
         main = "Jaccardindex - our pathways vs hallmark pathways",
)


#performing the same (jaccard index) for our-genesets against themselves and the hallmark genesets
data.jaccard.our <- data.frame(matrix(ncol=length(our_genesets_final), nrow = length(our_genesets_final)))
data.jaccard.hm <- data.frame(matrix(ncol=length(genesets_ids), nrow = length(genesets_ids)))

colnames(data.jaccard.our) <- names(our_genesets_final) #n
rownames(data.jaccard.our) <- names(our_genesets_final) #m

colnames(data.jaccard.hm) <- names(genesets_ids) #n
rownames(data.jaccard.hm) <- names(genesets_ids) #m

n=1
m=1

for (set_our in our_genesets_final){
  for (set_our2 in our_genesets_final){
    jacindex <- jaccard(set_our, set_our2)
    data.jaccard.our[m,n] <- jacindex
    m <- m+1
  }
  m <- 1
  n <- n+1
}

n=1
m=1

for (set_hm1 in genesets_ids){
  for (set_hm in genesets_ids){
    jacindex <- jaccard(set_hm1, set_hm)
    data.jaccard.hm[m,n] <- jacindex
    m <- m+1
  }
  m <- 1
  n <- n+1
}


#--------------
#heatmaps
#--------------

#our genesets vs themselves
pheatmap(as.matrix(data.jaccard.our),
         color = colorRampPalette(c('lightskyblue','lightcyan','white','yellow','orange', 'red'),
                                  bias = 1,
                                  space = 'rgb',
                                  interpolate = 'linear'
         )(50),
         clustering_method = 'average', treeheight_row = 20, treeheight_col = 20,
         legend_breaks = c(0, 1),
         legend_labels = c('unlike','same'),
         fontsize = 1, border_color = 'lightcyan2'
)

#hallmark genesets vs themselves
pheatmap(as.matrix(data.jaccard.hm),
         color = colorRampPalette(c('lightskyblue','lightcyan','white','yellow','orange', 'red'),
                                  bias = 1,
                                  space = 'rgb',
                                  interpolate = 'linear'
         )(50),
         clustering_method = 'average', treeheight_row = 20, treeheight_col = 20,
         cellwidth = 9, cellheight = 9,
         fontsize = 8, border_color = 'lightcyan2',
         legend_breaks = c(0, 1),
         legend_labels = c('unlike','same')
)       
