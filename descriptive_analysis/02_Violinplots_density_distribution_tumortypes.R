#-------------------------------------------------------------------------------
# Darstellung der Verteilung der Genexpressions-Verteilung (density plot) f?r 
# jeden unserer 5 Krebstypen: Darstellung ?ber violin plots
#-------------------------------------------------------------------------------

#Daten laden:
load('~/GitHub/2022-topic-02-team-05/data/tcga_exp_cleaned.RData')
tcga_anno = readRDS("~/GitHub/2022-topic-02-team-05/data/tcga_tumor_annotation.RDS")

#Packages:
#install.packages("vioplot")
library("vioplot")

#Patienten herauspicken nach den einzelnen Tumortypes:
BRCA_genes <- tcga_exp_cleaned[,tcga_anno$cancer_type_abbreviation == 'BRCA']
  
KIRC_genes <- tcga_exp_cleaned[,tcga_anno$cancer_type_abbreviation == 'KIRC']

LUAD_genes <- tcga_exp_cleaned[,tcga_anno$cancer_type_abbreviation == 'LUAD']

PRAD_genes <- tcga_exp_cleaned[,tcga_anno$cancer_type_abbreviation == 'PRAD']

THCA_genes <- tcga_exp_cleaned[,tcga_anno$cancer_type_abbreviation == 'THCA']


#Median der einzelnen Gene:
BRCA_mean <- apply(BRCA_genes,1, median)
KIRC_mean <- apply(KIRC_genes,1, median)
LUAD_mean <- apply(LUAD_genes,1, median)
PRAD_mean <- apply(PRAD_genes,1, median)
THCA_mean <- apply(THCA_genes,1, median)

#BRCA_mean <- apply(BRCA_genes,1, var)
#KIRC_mean <- apply(KIRC_genes,1, var)
#LUAD_mean <- apply(LUAD_genes,1, var)
#PRAD_mean <- apply(PRAD_genes,1, var)
#THCA_mean <- apply(THCA_genes,1, var)

#Violinplots mit einzelnen Tumortypes:
vioplot(BRCA_mean, xlab = "Thyorid cancer: Geneexpression Density Plot", ylab = "gene expression")
vioplot(KIRC_mean)
vioplot(LUAD_mean)
vioplot(PRAD_mean)
vioplot(THCA_mean)

#alle Violinplots in einem Plot; daf?r einen Dataframe mit allen mean-Expressionsdaten
#erstellen mit zugeh?rigem Cancer-Type

BRCA_mframe <- data.frame('gene_mean' = BRCA_mean, 'cancer_type' = 'BRCA')
KIRC_mframe <- data.frame('gene_mean' = KIRC_mean, 'cancer_type' = 'KIRC')
LUAD_mframe <- data.frame('gene_mean' = LUAD_mean, 'cancer_type' = 'LUAD')
PRAD_mframe <- data.frame('gene_mean' = PRAD_mean, 'cancer_type' = 'PRAD')
THCA_mframe <- data.frame('gene_mean' = THCA_mean, 'cancer_type' = 'THCA')

five_cancers_frame <- rbind(BRCA_mframe, KIRC_mframe, LUAD_mframe, PRAD_mframe, THCA_mframe)

#Nun Violinplot erstellen, sodass alle plots nebeneinander sind

vioplot(five_cancers_frame$gene_mean ~ five_cancers_frame$cancer_type, col = 2:length(levels(five_cancers_frame$cancer_type)),
        xlab = "cancer type", ylab = "gene expression")

  