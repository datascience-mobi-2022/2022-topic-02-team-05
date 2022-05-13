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
for (geneset.our in our_genesets){
  for (geneset.hm in hallmarks_genesets$genesets) { 
    sapply(geneset.our, setequal(geneset.our,geneset.hm))
  }
}

n = 1
m = 1
for (geneset.our in our_genesets){
  for (gene in geneset.our){
    for (geneset.hm in hallmarks_genesets$genesets){ 
      count[m,n] <- gene %in% geneset.hm
      n <- n+1
    }
  }
}




       