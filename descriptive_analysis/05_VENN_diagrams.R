#----------------------------------------------------------------------
# In this document, VENN diagramms scripted. They show how large the portion of 
# all genes (19200x) is that is in our pathways (600x)
#-----------------------------------------------------------------------

install.packages("VennDiagram")                       # Install VennDiagram package
library("VennDiagram")                                # Load VennDiagram package

load("~/GitHub/2022-topic-02-team-05/data/tcga_genes_cleaned.RData")
load("~/GitHub/2022-topic-02-team-05/data/tcga_exp_cleaned.RData")

#------------------------------------------------------------------------
# unlisting der Gennamen von Expressionsdaten
#------------------------------------------------------------------------

list <- rownames(tcga_exp_cleaned)
ensemble_vector_tcga <- unlist(list)

#-----------------------------------------------------------------------
# Gennamen aus unseren Pathways extrahieren und in einen Vektor packen
#-----------------------------------------------------------------------


ensemble_vector_pathways <- 
  
#-----------------------------------------------------------------------
# Gennamen finden, die in Pathways und in Expressionsdaten sind
#-----------------------------------------------------------------------  
ensemble_vector_both <- c()

for i in length(ensemble_vector_tcga){
  a <- is.element(i, ensemble_vector_pathways)
  if a == TRUE {
    append(ensemble_vector_both, ensemble_vector[i])
  }
}

#-----------------------------------------------------------------------
# VENN-Diagramm machen
#-----------------------------------------------------------------------

grid.newpage()                                        # Move to new plotting page
draw.pairwise.venn(area1 = ensemble_vector_tcga,      # Create pairwise venn diagram
                   area2 = ensemble_vector_pathways,
                   cross.area = ensemble_vector_both)

