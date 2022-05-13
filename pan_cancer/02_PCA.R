#--------------------------------------
# In diesem Dkument nehemn wir die GSEA Daten und führen eine PCA sowie wie
# ein clustering der transformierten Daten durch so können wir hoffentlich verschiedenen
# unterschiedliche Pathway activtäten zuordnen
#------------------------------------------

load('~/GitHub/2022-topic-02-team-05/data/GSEA_matrix.RData')

#--------------------------------------------
#Durchführen der PCA
#---------------------------------------------

PCA = prcomp(GSEA_matrix)
