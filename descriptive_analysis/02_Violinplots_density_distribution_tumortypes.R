#-------------------------------------------------------------------------------
# Darstellung der Verteilung der Genexpressions-Verteilung (density plot) für 
# jeden unserer 5 Krebstypen: Darstellung über violin plots
#-------------------------------------------------------------------------------

#Daten laden:
load('~/GitHub/2022-topic-02-team-05/data/tcga_exp_cleaned.RData')
tcga_anno = readRDS("~/GitHub/2022-topic-02-team-05/data/tcga_tumor_annotation.RDS")

#Packages:
install.packages("vioplot")
library("vioplot")

#Patienten herauspicken nach den einzelnen Tumortypes:
BRCA_genes <- tcga_exp_cleaned[,tcga_anno$cancer_type_abbreviation == 'BRCA']
  
KIRC_genes <- tcga_exp_cleaned[,tcga_anno$cancer_type_abbreviation == 'KIRC']

LUAD_genes <- tcga_exp_cleaned[,tcga_anno$cancer_type_abbreviation == 'LUAD']

PRAD_genes <- tcga_exp_cleaned[,tcga_anno$cancer_type_abbreviation == 'PRAD']

THCA_genes <- tcga_exp_cleaned[,tcga_anno$cancer_type_abbreviation == 'THCA']

#rbind(BRCA_genes, KIRC_genes,LUAD_genes,PRAD_genes,THCA_genes)

#Mittelwerte der einzelnen Gene:
BRCA_mean <- apply(BRCA_genes,1, mean)
KIRC_mean <- apply(KIRC_genes,1, mean)
LUAD_mean <- apply(LUAD_genes,1, mean)
PRAD_mean <- apply(PRAD_genes,1, mean)
THCA_mean <- apply(THCA_genes,1, mean)

#Violinplots mit einzelnen Tumortypes:

vioplot(BRCA_mean)

#alle Violinplots in einem Plot! 
#Achsen noch beschriften


  