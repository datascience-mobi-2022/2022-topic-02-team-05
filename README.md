---
editor_options: 
  markdown: 
    wrap: 72
---

# 2022-topic-02-team-05

# Thyroid cancer (THCA)

# Maren Schneider, Anna Lange, Jakob Then, David Matuschek

#folder structur

## preprocessing

### 01_TCGA_ensembleID_genenames_extraction: extraction of ensembleIDs
and genenames from the expression data from the TCGA-dataframe

### 02_Renaming_hallmarkpathways: renaming the genenames from the
hallmarkpathways-dataframe into ensembleIDs

### 03_metabolic_pathway_selection: selection of the metabolic pathways
from the msigdbr database

### 04_TCGA_exp_cleaning_and_biotype_analysis: Cleaning of the TCGA
expression data based on variance and biotype

### 05_THCA_exp_cleaning: Cleaning of the THCA expressiondata based on
variance and biotypes

### 06_Additional_annotations: adding the full name of the abbreviated
cancer type in the TCGA annotations and the histological grade to the
annotations of the THCA matrix

## descriptive analysis

### 01_Jaccard_Index: shows the similarity of the genes from our pathways
vs the hallmarkpathways

### 02_Violinplots_density_distribution_tumortypes: shows the
distribution of the geneexpression in 5 different cancer types from the
TCGA data

### 03_Plotting_Mean_over_Variance: plotting the mean of the
geneexpression from the TCGA data over the the variance

### 04_VolcanoPlots: comparing mean geneexpression of tumor tissue vs
normal tissue

### 05_VENN_diagrams: showing how many genes from the TCGA expressiondata
are in the pathways we picked to analyze

## pan_cancer

### 01_GSEA_of_TCGA_exp: computing of enrichment scores for the
expressiondata from the TCGA based on GSEA

### 02_PCA_and_UMAP: PCA and UMAP of the expressiondata from the TCGA and
the enrichment scores of the pathways

### 03_GSVA_of_TCGA: GSVA of the THCA expressiondata and comparison with
GSEA

### 04_Pancancer_analysis_of_THCA: clustering of the GSVA values for the
TCGA expressiondata and visualising of THCA subtypes also UMAP and PCA
of those GSVA values

### BONUS_Wang_challenge: function that takes a cancertype dataframe and
genesets as input and computes the p-value

## thca_normal_vs_tumor

### 01_GSEA_and_GSVA_of_THCA: GSEA and GSVA of the THCA expressiondata
and volcanoplots based on GSVA of THCA tumor tissue vs normal tissue

### 02_THCA_PCA_and_UMAP: PCA and UMAP of the expression data of THCA
expressiondata (via foldchange) and of the GSEA values for the THCA data
