## 2022-topic-02-team-05

# Thyroid cancer (THCA)

# Maren Schneider, Anna Lange, Jakob Then, David Matuschek


## abstract
In recent years bioinformatic methods became a tool of utmost importance in medical research. To detect specific gene and pathway patterns in different cancer types or histological types, here, we perform pan-cancer analysis. Furthermore, we carry out a focused analysis to specify different subcategories within a certain cancer type and to identify targets for therapy. As main methods in identifying up- or downregulated pathways, Gene Set Enrichment Analysis (GSEA) and Gene Set Variation Analysis (GSVA) are used. GSVA of The Cancer Genome Atlas (TCGA) expression data, containing 33 different tumor types, reveals four clusters of cancer types, which are defined by different histological types like glioblastoma and adenocarcinoma. Thus, the histological types seem to correlate with a specific set of pathways that are especially enriched in certain types. Furthermore, GSVA results of Thyroid cancer (THCA) expression data show that thyroid carcinogenesis is associated with the upregulation of proliferation signaling pathways like the hedgehog pathway, alpha6beta4 integrin signaling pathway, and associated pathways such as interleukin 36 signaling. Based on those proliferation signaling pathways, three subclusters form inside the THCA patients from the pan-cancer data. One THCA subcluster could be linked to the follicular histological subtype and is defined by increased mTOR and MAPK activity while having low alpha6beta4 activity. In contrast, another THCA subcluster is defined by a low mTOR and MAPK activity, but a high alpha6beta4 activity. The third THCA subtype we identified shows both ways of proliferation signaling. These results promise more success in THCA treatment, as a more precise diagnosis of the distinct THCA subtypes is possible allowing for better treatment. To improve the understanding of THCA, and thereby bettering patients’ prognosis, this project focuses on finding genes that have a significantly altered expression in THCA compared to other cancers and especially to normal tissue.

## folder structure
### preprocessing

01_TCGA_ensembleID_genenames_extraction: extraction of ensembleIDs and
genenames from the expression data from the TCGA-dataframe

02_Renaming_hallmarkpathways: renaming the genenames from the
hallmarkpathways-dataframe into ensembleIDs

03_metabolic_pathway_selection: selection of the metabolic pathways from
the msigdbr database

04_TCGA_exp_cleaning_and_biotype_analysis: Cleaning of the TCGA
expression data based on variance and biotype

05_THCA_exp_cleaning: Cleaning of the THCA expressiondata based on
variance and biotypes

06_Additional_annotations: adding the full name of the abbreviated
cancer type in the TCGA annotations and the histological grade to the
annotations of the THCA matrix

### descriptive analysis

01_Jaccard_Index: shows the similarity of the genes from our pathways vs the hallmarkpathways

02_Violinplots_density_distribution_tumortypes: shows the distribution of the geneexpression in 5 different cancer types from the TCGA data

03_Plotting_Mean_over_Variance: plotting the mean of the geneexpression from the TCGA data over the the variance

04_VolcanoPlots: comparing mean geneexpression of tumor tissue vs normal tissue

05_VENN_diagrams: showing how many genes from the TCGA expressiondata are in the pathways we picked to analyze

### pan_cancer

01_GSEA_of_TCGA_exp: computing of enrichment scores for the expressiondata from the TCGA based on GSEA

02_PCA_and_UMAP: PCA and UMAP of the expressiondata from the TCGA and the enrichment scores of the pathways

03_GSVA_of_TCGA: GSVA of the THCA expressiondata and comparison with GSEA

04_Pancancer_analysis_of_THCA: clustering of the GSVA values for the TCGA expressiondata and visualising of THCA subtypes also UMAP and PCA of those GSVA values

BONUS_Wang_challenge: function that takes a cancertype dataframe and genesets as input and computes the p-value

### thca_normal_vs_tumor

01_GSEA_and_GSVA_of_THCA: GSEA and GSVA of the THCA expressiondata and volcanoplots based on GSVA of THCA tumor tissue vs normal tissue

02_THCA_PCA_and_UMAP: PCA and UMAP of the expression data of THCA expressiondata (via foldchange) and of the GSEA values for the THCA data

### neuronal_network_and_linear_regression: 
contains code for the linear regression to predict pathway activity and the neuronal network as an alternative

### output: 
contains image files of the plots from or R code

### markdown: 
contains the final report as a Rmarkdown ("Topic02_Team05.Rmd.pdf") and bib folders for references, citations and pictures in the report 
