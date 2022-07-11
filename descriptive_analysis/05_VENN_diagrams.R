#----------------------------------------------------------------------
#VENN Diagramm wird erstellt, um zu bestimmen, wie viele der Gene der 
#TCGA-Expressionsmatrix (19200+) in unseren pathways sind (600+)
#-----------------------------------------------------------------------

#install.packages("VennDiagram")                       # Install VennDiagram package
library("VennDiagram")                                # Load VennDiagram package

load("data/tcga_genes_cleaned.RData")
load("data/tcga_exp_cleaned.RData")
load("data/pathways.RData")

#unlisting der Gennamen von Expressionsdaten
list <- rownames(tcga_exp_cleaned)
ensemble_vector_tcga <- unlist(list)

#identische Ensemble IDs rauswerfen
ensemble_vector_tcga <- unique(ensemble_vector_tcga) #kein Gen wird rausgeworfen! (Alle unterschiedlich)

#Gennamen aus unseren Pathways extrahieren und in einen Vektor packen
pathw <- unlist(pathways, use.names = FALSE)
ensemble_vector_pathways <- pathw
  
#Gennamen finden, die in Pathways und in Expressionsdaten sind
ensemble_vector_both <- intersect(ensemble_vector_pathways, ensemble_vector_tcga)

#VENN-Diagramm erstellen
grid.newpage()                                        # Move to new plotting page
draw.pairwise.venn(area1 = length(ensemble_vector_tcga),      # Create pairwise venn diagram
                   area2 = length(ensemble_vector_pathways),
                   cross.area = length(ensemble_vector_both),
                   fill = c("pink", "orange"),
                   lty = "blank",
                   category = c("TCGA", "Pathways"))

