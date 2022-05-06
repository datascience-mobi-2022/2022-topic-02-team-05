#wichtige packages
#install.packages("msigdbr")
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("biomaRt")
library(biomaRt)
library(msigdbr)
library(dplyr)
library(ggplot2)
#laden der daten von carl
#tcga_exp = readRDS("C:/Users/jakob/Documents/R/Data Analysis Seminar/tcga_tumor_log2TPM.RDS")
#tcga_annot = readRDS("C:/Users/jakob/Documents/R/Data Analysis Seminar/tcga_tumor_annotation.RDS")
#tcga_tumor_norm = readRDS("C:/Users/jakob/Documents/R/Data Analysis Seminar/tcga_tumor_normal_datascience_proj_2020.RDS")
#genesets = readRDS("C:/Users/jakob/Documents/R/Data Analysis Seminar/hallmarks_genesets.rds")


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

#loading our selected pathways
gsea_pathways = msigdbr(species = "Homo sapiens")
TERT_pathway = gsea_pathways[gsea_pathways$gs_name == "BIOCARTA_TEL_PATHWAY", ]
MAPK_pathway = gsea_pathways[gsea_pathways$gs_name == "BIOCARTA_MAPK_PATHWAY", ]
thyroidhormone_pathway = gsea_pathways[gsea_pathways$gs_name == "GOBP_THYROID_HORMONE_METABOLIC_PROCESS", ]

#liste aller pathways mit allen inofs zum jeweiligen pathway
our_genesets = list(TERT_pathway$gene_symbol,
                    MAPK_pathway$gene_symbol,
                    thyroidhormone_pathway$gene_symbol
                    )
names(our_genesets) = c('TERT_pathway', 'MAPK_pathway', 'thyroidhormone_pathway')



#---------------------------------------------------
#genetypen der pathways checken
#---------------------------------------------------

#function die Gennamen nimmt und dafür den gentypen gibt
checkbiotypes = function(pathway){mart = useEnsembl(dataset = "hsapiens_gene_ensembl", biomart='ensembl')
  
                                  res = getBM(attributes = c("ensembl_gene_id", "gene_biotype"),
                                        filters = "external_gene_name",
                                        values = pathway,
                                        mart = mart )
                                  message('Let me check that for you :D')
                                  
                                  return(res)
                                  }
#-----------------------------------------
#gene biotypes unserer pathways bestimmen
#-----------------------------------------
our_genesets_biotypes = sapply(our_genesets, checkbiotypes)
#entfernen der ids
our_genesets_biotypes = our_genesets_biotypes[2,]




#function die alle biotypes eines pathways als barplot plottet
#genelist ist dabei die liste unserer pathways, number die stelle des pathways in der liste
plotbiotypes = function(genelist, number){
  x = as.data.frame(table(as.vector(genelist[[number]])))
  colnames(x) = c('Biotype','Ammount')
  
  ggplot(x, aes(Biotype, Ammount)) + geom_bar(stat = 'identity') +  
    labs(title = names(genelist[number])) +
    theme(axis.text = element_text(angle = 0,))
}

#plotten aller unserer pathways
#funktioniert auch noch nicht!!!
par(mfrow=c(1,3))
for (i in 1:length(our_genesets_biotypes)) {
  plot = plotbiotypes(our_genesets_biotypes, i)
  
}


#-------------------------------
#gene biotypes von carls pathways bestimmen
#----------------------------------

#DAUERT EWIG!!!!!!!
genesets_biotypes = sapply(genesets[[1]], checkbiotypes)
#entfernen der ids
genesets_biotypes = genesets_biotypes[2,]


#plotten des ersten von carls pathways
plotbiotypes(genesets_biotypes, 1)




#---------------------------------------------------------
#analyse der biotypen der gene aus den tcga expressionsdaten
#---------------------------------------------------------

#gennamen extrahieren
tcga_genes = rownames(tcga_exp_hvar)

#dieser vetor enthält sowohl enseblm id als auch genenamen und muss daher gespalten werden
tcga_genes = strsplit(tcga_genes, split = '|', fixed = TRUE)

#speicher der gennamen als eigenen vektor
tcga_genenames = sapply(tcga_genes, function(tcga_genes){return(tcga_genes[2])})

#entfernen der Versionsnummern für ensembl id und gennamen
tcga_genenames = strsplit(tcga_genenames, split = '.', fixed = TRUE)
tcga_genenames = sapply(tcga_genenames, function(tcga_genenames){return(tcga_genenames[1])})

#biotypes der tcga gene bestimmen
tcga_biotypes = checkbiotypes(tcga_genenames)

#darstellen aller biotypes im set
x = as.data.frame(table(as.vector(tcga_biotypes$gene_biotype)))
colnames(x) = c('Biotype','Ammount')

ggplot(x, aes(Biotype, Ammount)) + geom_bar(stat = 'identity') +  
  labs(title = 'TCGA expression data biotypes') +
  theme(axis.text = element_text(angle = 90,))


#--------------------------------------------------------------------------
#entfernen aller nicht proteincodierenden gene aus der tcga expressions matrix
#--------------------------------------------------------------------------
tcga_exp_protcod = tcga_exp_hvar[tcga_biotypes$gene_biotype == 'protein_coding', ]

