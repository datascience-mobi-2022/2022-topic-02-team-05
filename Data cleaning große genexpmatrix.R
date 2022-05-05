#----------------------------------------
#preprocessing of tcga gene exp data
#----------------------------------------


#checking for NAs
#erstmal nur in den ersten tausend genen weil ich sonst fehler bekomm ):
tcga_exp_narm = na.omit(tcga_exp[1:1000,])

#sorting out low variance genes
tcga_exp_var = apply(tcga_exp_narm, 1, var)
hist(log(tcga_exp_var), breaks = 50, probability = TRUE)

#cutting der gene mit sehr niedriger exression d.h. log(var) < -1
tcga_exp_hvar = tcga_exp_narm[log(tcga_exp_var) > -1, ]
                                                
#------------------------------
#auswahl nicer pathways
#------------------------------

#install.packages("msigdbr")
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("biomaRt")

#loading our selected pathways
gsea_pathways = msigdbr(species = "Homo sapiens")
TERT_pathway = gsea_pathways[gsea_pathways$gs_name == "BIOCARTA_TEL_PATHWAY", ]
MAPK_pathway = gsea_pathways[gsea_pathways$gs_name == "BIOCARTA_MAPK_PATHWAY", ]
thyroidhormone_pathway = gsea_pathways[gsea_pathways$gs_name == "GOBP_THYROID_HORMONE_METABOLIC_PROCESS", ]

our_genesets = list(TERT_pathway$gene_symbol,
                    MAPK_pathway$gene_symbol,
                    thyroidhormone_pathway$gene_symbol
                    )
names(our_genesets) = c('TERT_pathway', 'MAPK_pathway', 'thyroidhormone_pathway')


#---------------------------------------------------
#genetypen der pathways checken
#alles zur biomart function hab ich hier her:
#https://cran.r-project.org/web/packages/biomartr/vignettes/Functional_Annotation.html
#---------------------------------------------------

#function die Gennamen nimmt und dafür den gentypen gibt
checkgenetypes = function(pathway){biomart(pathway,
                                           mart = "ENSEMBL_MART_ENSEMBL",
                                           dataset = "hsapiens_gene_ensembl",
                                           attributes = c("description","gene_biotype","transcript_biotype"),
                                           filters = "external_gene_name")
                                   message('Let me check that for you :D')   
                                  }
TERT_genetypes = checkgenetypes(TERT_pathway$gene_symbol)




