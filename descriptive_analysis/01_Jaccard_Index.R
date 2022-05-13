#----------------------------------------------
# In diesem Dokument suchen wir Gemeinsamkeiten zwischen den Hallmarkgenen und unseren Pathways
# Verwendung des Jaccard-Indexes (Schnittmenge im Verhältnis zur Vereinigungsmeng)
#----------------------------------------------

#loading data:

load("~/GitHub/2022-topic-02-team-05/data/our_genesets.RData")
hallmarks_genesets = readRDS("~/GitHub/2022-topic-02-team-05/data/hallmarks_genesets.rds")


# Jaccard-Index zwischen Genen:
# für jedes Listenelement aus our_genesets(Vektor mit Gennamen) schauen ob es in Hallmark_geneset 1 drinnen ist
# wenn es drinnen ist: count+1 in Matrix; dann nächstes Gen

# leeren Data-frame erstellen:

data.jaccard <- data.frame(matrix(ncol=length(our_genesets), nrow = length(hallmarks_genesets$genesets)))

colnames(data.jaccard) <- names(our_genesets) #n
rownames(data.jaccard) <- names(hallmarks_genesets$genesets) #m

#jaccard_index:

n = 1
m = 1

jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}

for (set_our in our_genesets){
  for (set_hm in hallmarks_genesets$genesets){
    jacindex <- jaccard(set_our, set_hm)
    print(jacindex)
    data.jaccard[m,n] <- jacindex
    m <- m+1
  }
  m <- 1
  n <- n+1
}

heatmap(as.matrix(data.jaccard), mar = c(12.5,10), main= 'Similarity of hallmark and our genesets', xlab = 'our genesets', ylab = 'hallmark genesets')


       