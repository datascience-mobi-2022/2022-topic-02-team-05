#----------------------------------------------
# In diesem Dokument suchen wir Gemeinsamkeiten zwischen den Hallmarkgenen und unseren Pathways
# Verwendung des Jaccard-Indexes (Schnittmenge im Verhältnis zur Vereinigungsmeng)
#----------------------------------------------

#loading data:

load("~/GitHub/2022-topic-02-team-05/data/our_genesets.RData")
hallmarks_genesets = readRDS("~/GitHub/2022-topic-02-team-05/data/hallmarks_genesets.rds")


# zwischen Genen: