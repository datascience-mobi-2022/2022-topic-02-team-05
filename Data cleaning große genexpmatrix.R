#wichtige packages
#install.packages("msigdbr")
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
#BiocManager::install("biomaRt")
library(biomaRt)
library(msigdbr)
library(dplyr)
library(ggplot2)
#laden der daten von carl
tcga_exp = readRDS("data/tcga_tumor_log2TPM.RDS")
tcga_annot = readRDS("data/tcga_tumor_annotation.RDS")
tcga_tumor_norm = readRDS("data/tcga_tumor_normal_datascience_proj_2022.RDS")
genesets = readRDS("data/hallmarks_genesets.rds")


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
JAKSTAT_pathway = gsea_pathways[gsea_pathways$gs_name == "KEGG_JAK_STAT_SIGNALING_PATHWAY", ]
ferroptosis_pathway = gsea_pathways[gsea_pathways$gs_name == "WP_FERROPTOSIS", ]
translationInitiation_pathway = gsea_pathways[gsea_pathways$gs_name == "BIOCARTA_EIF_PATHWAY", ]
mehrEnergie_pathway = gsea_pathways[gsea_pathways$gs_name == "BIOCARTA_ETC_PATHWAY", ]
glycolysis_pathway = gsea_pathways[gsea_pathways$gs_name == "BIOCARTA_GLYCOLYSIS_PATHWAY", ]
citric_cycle_pathway = gsea_pathways[gsea_pathways$gs_name == "BIOCARTA_KREB_PATHWAY", ]
synapsis_pathway = gsea_pathways[gsea_pathways$gs_name == "BIOCARTA_PDZS_PATHWAY", ]
aminoacid_pathway = gsea_pathways[gsea_pathways$gs_name == "REACTOME_AMINO_ACID_SYNTHESIS_AND_INTERCONVERSION_TRANSAMINATION", ]
estrogen_pathway = gsea_pathways[gsea_pathways$gs_name == "WP_ESTROGEN_METABOLISM", ]
ACE_pathway = gsea_pathways[gsea_pathways$gs_name == "BIOCARTA_ACE2_PATHWAY", ] #Angiotensin converting enzyme: Niere 
inflam_response_pathway = gsea_pathways[gsea_pathways$gs_name == "BIOCARTA_INFLAM_PATHWAY", ] #Cytokines and Inflammatory Response
CaMK_pathway = gsea_pathways[gsea_pathways$gs_name == "KEGG_HEMATOPOIETIC_CELL_LINEAGE", ] #Hematopoietic cell lineage
stemmcell_reg_pathway = gsea_pathways[gsea_pathways$gs_name == "REACTOME_TRANSCRIPTIONAL_REGULATION_OF_PLURIPOTENT_STEM_CELLS", ] #Transcriptional regulation of pluripotent stem cells
mammarygland1_embryo_pathway = gsea_pathways[gsea_pathways$gs_name == "WP_MAMMARY_GLAND_DEVELOPMENT_PATHWAY_EMBRYONIC_DEVELOPMENT_STAGE_1_OF_4", ]
mammarygland2_puberty_pathway = gsea_pathways[gsea_pathways$gs_name == "WP_MAMMARY_GLAND_DEVELOPMENT_PATHWAY_PUBERTY_STAGE_2_OF_4", ]
mammarygland3_pregnancy_pathway = gsea_pathways[gsea_pathways$gs_name == "WP_MAMMARY_GLAND_DEVELOPMENT_PATHWAY_PREGNANCY_AND_LACTATION_STAGE_3_OF_4", ]
mammarygland4_abbau_pathway = gsea_pathways[gsea_pathways$gs_name == "WP_MAMMARY_GLAND_DEVELOPMENT_PATHWAY_INVOLUTION_STAGE_4_OF_4", ]



#liste aller pathways mit allen inofs zum jeweiligen pathway
our_genesets = list(TERT_pathway$gene_symbol,
                    MAPK_pathway$gene_symbol,
                    thyroidhormone_pathway$gene_symbol,
                    JAKSTAT_pathway$gene_symbol,
                    ferroptosis_pathway$gene_symbol,
                    translationInitiation_pathway$gene_symbol,
                    mehrEnergie_pathway$gene_symbol
                    )
names(our_genesets) = c('TERT_pathway',
                        'MAPK_pathway',
                        'thyroidhormone_pathway',
                        'JAKSTAT_pathway',
                        'ferroptosis_pathway',
                        'translationInitiation_pathway',
                        'mehrEnergie_pathway'
                        )



#---------------------------------------------------
#genetypen der pathways checken
#---------------------------------------------------

#function die Gennamen nimmt und daf?r den gentypen gibt
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

#plotten aller unserer pathways einzeln
#funktioniert auch noch nicht!!!
par(mfrow=c(1,3))
for (i in 1:length(our_genesets_biotypes)) {
  plot = plotbiotypes(our_genesets_biotypes, i)
  
}

#plotten einer gesamt?bersicht ?ber die biotypes aller unserer pathways
res = NULL
for (i in 1:length(our_genesets_biotypes)){
  res = c(res, our_genesets_biotypes[[i]])
}
x = as.data.frame(table(res))
colnames(x) = c('Biotype','Ammount')

ggplot(x, aes(Biotype, Ammount)) + geom_bar(stat = 'identity') +  
  labs(title = 'Our pathways biotypes') +
  theme(axis.text = element_text(angle = 90,))


#-------------------------------
#gene biotypes von carls pathways bestimmen
#----------------------------------

#DAUERT EWIG!!!!!!!
genesets_biotypes = sapply(genesets[[1]], checkbiotypes)
#entfernen der ids
genesets_biotypes = genesets_biotypes[2,]


#plotten des ersten von carls pathways
plotbiotypes(genesets_biotypes, 1)
#besser hier ein funktion die alle pathways nacheinander in eine Grafik plottet

#plotten einer ?bersicht ?ber die biotypes aller von carls pathways
res = NULL
for (i in 1:length(genesets_biotypes)){
  res = c(res, genesets_biotypes[[i]])
}
x = as.data.frame(table(res))
colnames(x) = c('Biotype','Ammount')

ggplot(x, aes(Biotype, Ammount)) + geom_bar(stat = 'identity') +  
  labs(title = 'Carls pathways biotypes') +
  theme(axis.text = element_text(angle = 90,))


#---------------------------------------------------------
#analyse der biotypen der gene aus den tcga expressionsdaten
#---------------------------------------------------------

#gennamen extrahieren
tcga_genes = rownames(tcga_exp_hvar)

#dieser vetor enth?lt sowohl enseblm id als auch genenamen und muss daher gespalten werden
tcga_genes = strsplit(tcga_genes, split = '|', fixed = TRUE)

#speicher der gennamen als eigenen vektor
tcga_genenames = sapply(tcga_genes, function(tcga_genes){return(tcga_genes[2])})

#entfernen der Versionsnummern f?r ensembl id und gennamen
tcga_genenames = strsplit(tcga_genenames, split = '.', fixed = TRUE)
tcga_genenames = sapply(tcga_genenames, function(tcga_genenames){return(tcga_genenames[1])})

#hier am besten noch doppelte namen mit eventuell anderen versionsnummern ?berpr?fen und die neuere nehemen

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
tcga_exp_protcod = na.omit(tcga_exp_protcod)
##hier geht sich das iwie nicht mit den zeilen aus wieso kommen wir von 1000 auf 800 auf 300 auf 140

#----------------------------------------------------
#pca der gereinigten tcga expressionsdaten
#----------------------------------------------------

tcga_exp_pca = prcomp(tcga_exp_protcod)

#plotten der prozentualen var die durch die ersten 10 PCs erkl?rt wird
barplot(sqrt(tcga_exp_pca$sdev[1:10]/sum(tcga_exp_pca$sdev)),
        names.arg = colnames(tcga_exp_pca$x)[1:10])  


