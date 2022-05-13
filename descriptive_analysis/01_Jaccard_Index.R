#----------------------------------------------
# In diesem Dokument suchen wir Gemeinsamkeiten zwischen den Hallmarkgenen und unseren Pathways
# Verwendung des Jaccard-Indexes (Schnittmenge im Verh?ltnis zur Vereinigungsmeng)
#----------------------------------------------

#loading data:

load("~/GitHub/2022-topic-02-team-05/data/our_genesets.RData")
load("~/GitHub/2022-topic-02-team-05/data/geneset_ids.RData")


# Jaccard-Index zwischen Genen:
# f?r jedes Listenelement aus our_genesets(Vektor mit Gennamen) schauen ob es in Hallmark_geneset 1 drinnen ist
# wenn es drinnen ist: count+1 in Matrix; dann n?chstes Gen

# leeren Data-frame erstellen:

data.jaccard <- data.frame(matrix(ncol=length(our_genesets), nrow = length(genesets_ids)))

colnames(data.jaccard) <- names(our_genesets) #n
rownames(data.jaccard) <- names(genesets_ids) #m

#jaccard_index:

n = 1
m = 1

jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}

for (set_our in our_genesets){
  for (set_hm in genesets_ids){
    jacindex <- jaccard(set_our, set_hm)
    print(jacindex)
    data.jaccard[m,n] <- jacindex
    m <- m+1
  }
  m <- 1
  n <- n+1
}

heatmap(as.matrix(data.jaccard), mar = c(12.5,10), main= 'Similarity of hallmark and our genesets', xlab = 'our genesets', ylab = 'hallmark genesets')
library(pheatmap)
pheatmap(as.matrix(data.jaccard),
         breaks = seq(0, max(data.jaccard), length.out = 51),
         color = colorRampPalette(c('lightskyblue','lightcyan','white','yellow','orange', 'red'),
                                  bias = 1,
                                  space = 'rgb',
                                  interpolate = 'linear'
         )(50),
         clustering_method = 'average', treeheight_row = 20, treeheight_col = 20,
         cellwidth = 10, cellheight = 10,
         fontsize = 8, border_color = 'lightcyan2',
         legend_breaks = c(0, max(data.jaccard)),
         legend_labels = c('unlike','same')
)

       