#-----------------------------------------------------
#Dokument in dem wir alle Metabolischen pathways laden und die Hallmarkpathways
#in ENsembl IDs umschreiben
#-----------------------------------------------------

#install.packages("msigdbr")
#if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")

library(msigdbr)

#?ffnen der Datenbank
gsea_pathways = msigdbr(species = "Homo sapiens")

#Auswahl der einzelnen pathways
TERT_pathway = gsea_pathways[gsea_pathways$gs_name == "BIOCARTA_TEL_PATHWAY", ]
MAPK_pathway = gsea_pathways[gsea_pathways$gs_name == "BIOCARTA_MAPK_PATHWAY", ]
thyroidhormone_pathway = gsea_pathways[gsea_pathways$gs_name == "GOBP_THYROID_HORMONE_METABOLIC_PROCESS", ]
JAKSTAT_pathway = gsea_pathways[gsea_pathways$gs_name == "KEGG_JAK_STAT_SIGNALING_PATHWAY", ]
ferroptosis_pathway = gsea_pathways[gsea_pathways$gs_name == "WP_FERROPTOSIS", ]
translationInitiation_pathway = gsea_pathways[gsea_pathways$gs_name == "BIOCARTA_EIF_PATHWAY", ]
mehrEnergie_pathway = gsea_pathways[gsea_pathways$gs_name == "BIOCARTA_ETC_PATHWAY", ]
#glycolysis_pathway = gsea_pathways[gsea_pathways$gs_name == "BIOCARTA_GLYCOLYSIS_PATHWAY", ]
citric_cycle_pathway = gsea_pathways[gsea_pathways$gs_name == "BIOCARTA_KREB_PATHWAY", ]
synapsis_pathway = gsea_pathways[gsea_pathways$gs_name == "BIOCARTA_PDZS_PATHWAY", ]
#aminoacid_pathway = gsea_pathways[gsea_pathways$gs_name == "REACTOME_AMINO_ACID_SYNTHESIS_AND_INTERCONVERSION_TRANSAMINATION", ]
estrogen_pathway = gsea_pathways[gsea_pathways$gs_name == "WP_ESTROGEN_METABOLISM", ]
ACE_pathway = gsea_pathways[gsea_pathways$gs_name == "BIOCARTA_ACE2_PATHWAY", ] #Angiotensin converting enzyme: Niere 
inflam_response_pathway = gsea_pathways[gsea_pathways$gs_name == "BIOCARTA_INFLAM_PATHWAY", ] #Cytokines and Inflammatory Response
CaMK_pathway = gsea_pathways[gsea_pathways$gs_name == "KEGG_HEMATOPOIETIC_CELL_LINEAGE", ] #Hematopoietic cell lineage
stemmcell_reg_pathway = gsea_pathways[gsea_pathways$gs_name == "REACTOME_TRANSCRIPTIONAL_REGULATION_OF_PLURIPOTENT_STEM_CELLS", ] #Transcriptional regulation of pluripotent stem cells
mammarygland1_embryo_pathway = gsea_pathways[gsea_pathways$gs_name == "WP_MAMMARY_GLAND_DEVELOPMENT_PATHWAY_EMBRYONIC_DEVELOPMENT_STAGE_1_OF_4", ]
mammarygland2_puberty_pathway = gsea_pathways[gsea_pathways$gs_name == "WP_MAMMARY_GLAND_DEVELOPMENT_PATHWAY_PUBERTY_STAGE_2_OF_4", ]
mammarygland3_pregnancy_pathway = gsea_pathways[gsea_pathways$gs_name == "WP_MAMMARY_GLAND_DEVELOPMENT_PATHWAY_PREGNANCY_AND_LACTATION_STAGE_3_OF_4", ]
mammarygland4_abbau_pathway = gsea_pathways[gsea_pathways$gs_name == "WP_MAMMARY_GLAND_DEVELOPMENT_PATHWAY_INVOLUTION_STAGE_4_OF_4", ]
TGFBEMT_pathway = gsea_pathways[gsea_pathways$gs_name == "FOROUTAN_TGFB_EMT_UP", ]
CaspaseCascade_pathway = gsea_pathways[gsea_pathways$gs_name == "BIOCARTA_CASPASE_PATHWAY", ]
NKcytotoxicity = gsea_pathways[gsea_pathways$gs_name == "KEGG_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY", ]
Tcellactivation = gsea_pathways[gsea_pathways$gs_name == "BIOCARTA_TCR_PATHWAY", ]
Bcellactivation = gsea_pathways[gsea_pathways$gs_name == "BIOCARTA_BCR_PATHWAY", ]
purin_denovo_pathway = gsea_pathways[gsea_pathways$gs_name == "KEGG_PURINE_METABOLISM", ]
Interleukin_pathway = gsea_pathways[gsea_pathways$gs_name == "BIOCARTA_IL1R_PATHWAY", ]
wnt_pathway = gsea_pathways[gsea_pathways$gs_name == "BIOCARTA_WNT_PATHWAY", ]
oxStress_pathway = gsea_pathways[gsea_pathways$gs_name == "BIOCARTA_ARENRF2_PATHWAY", ]
IGF1_pathway = gsea_pathways[gsea_pathways$gs_name == "BIOCARTA_IGF1_PATHWAY", ]
ACHR_pathway = gsea_pathways[gsea_pathways$gs_name == "BIOCARTA_ACH_PATHWAY", ]
HSP27_pathway = gsea_pathways[gsea_pathways$gs_name == "BIOCARTA_HSP27_PATHWAY", ]
RAN_pathway = gsea_pathways[gsea_pathways$gs_name == "BIOCARTA_RAN_PATHWAY", ]
glycosylation_pathway = gsea_pathways[gsea_pathways$gs_name == "BIOCARTA_AMAN_PATHWAY", ]
TATAregulated_genes_pathway = gsea_pathways[gsea_pathways$gs_name == "TATAAA_TATA_01", ] #TATA-regulated genes
TGGAAAregulated_genes_pathway = gsea_pathways[gsea_pathways$gs_name == "TGGAAA_NFAT_Q4_01", ] #TGGAAA-regulated genes
Retinoate_pathway = gsea_pathways[gsea_pathways$gs_name == "GOMF_RETINAL_DEHYDROGENASE_ACTIVITY", ] #Retinoate from Retinol via Dehydrogenase
Retinol_pathway = gsea_pathways[gsea_pathways$gs_name == "KEGG_RETINOL_METABOLISM", ] #Retinol metabolism
Farnesly_diphosphate_pathway = gsea_pathways[gsea_pathways$gs_name == "GOBP_FARNESYL_DIPHOSPHATE_METABOLIC_PROCESS", ] #for example squalene synthesis 
cholesterol_pathway = gsea_pathways[gsea_pathways$gs_name == "BIOCARTA_FXR_PATHWAY", ] #cholesterol regulation
cortisol_pathway = gsea_pathways[gsea_pathways$gs_name == "GOBP_CORTISOL_BIOSYNTHETIC_PROCESS", ] #cortisol synthesis
steroid_hormone_pathway = gsea_pathways[gsea_pathways$gs_name == "KEGG_STEROID_HORMONE_BIOSYNTHESIS", ] #steroid hormone synthesis in general
hypogonadism_pathway = gsea_pathways[gsea_pathways$gs_name == "HP_MALE_HYPOGONADISM", ] #Decreased functionality of the male gonad, i.e., of the testis, with reduced spermatogenesis or testosterone synthesis.
estradiol_responsive_genes_pathway = gsea_pathways[gsea_pathways$gs_name == "SYED_ESTRADIOL_RESPONSE", ] #genes responsive to estradiol regulation 
Oocyte_maturation_pathway = gsea_pathways[gsea_pathways$gs_name == "BIOCARTA_MPR_PATHWAY", ] #initiation of oocyte maturation via progesterone
Vitamin_D_pathway = gsea_pathways[gsea_pathways$gs_name == "GOBP_VITAMIN_D_BIOSYNTHETIC_PROCESS", ] #VitaminD synthesis
urea_pathway = gsea_pathways[gsea_pathways$gs_name == "GOBP_UREA_CYCLE", ] #urea synthesis

#liste aller pathways mit allen inofs zum jeweiligen pathway
our_genesets = list(TERT_pathway$human_ensembl_gene,
                    MAPK_pathway$human_ensembl_gene,
                    thyroidhormone_pathway$human_ensembl_gene,
                    JAKSTAT_pathway$human_ensembl_gene,
                    ferroptosis_pathway$human_ensembl_gene,
                    translationInitiation_pathway$human_ensembl_gene,
                    mehrEnergie_pathway$human_ensembl_gene,
                    #glycolysis_pathway$human_ensembl_gene,
                    citric_cycle_pathway$human_ensembl_gene,
                    synapsis_pathway$human_ensembl_gene,
                    #aminoacid_pathway$human_ensembl_gene,
                    estrogen_pathway$human_ensembl_gene,
                    ACE_pathway$human_ensembl_gene,
                    inflam_response_pathway$human_ensembl_gene,
                    CaMK_pathway$human_ensembl_gene,
                    stemmcell_reg_pathway$human_ensembl_gene,
                    mammarygland1_embryo_pathway$human_ensembl_gene,
                    mammarygland2_puberty_pathway$human_ensembl_gene,
                    mammarygland3_pregnancy_pathway$human_ensembl_gene,
                    mammarygland4_abbau_pathway$human_ensembl_gene,
                    TGFBEMT_pathway$human_ensembl_gene,
                    CaspaseCascade_pathway$human_ensembl_gene,
                    NKcytotoxicity$human_ensembl_gene,
                    Tcellactivation$human_ensembl_gene,
                    Bcellactivation$human_ensembl_gene,
                    purin_denovo_pathway$human_ensembl_gene,
                    Interleukin_pathway$human_ensembl_gene,
                    wnt_pathway$human_ensembl_gene,
                    oxStress_pathway$human_ensembl_gene,
                    IGF1_pathway$human_ensembl_gene,
                    ACHR_pathway$human_ensembl_gene,
                    HSP27_pathway$human_ensembl_gene,
                    RAN_pathway$human_ensembl_gene,
                    glycosylation_pathway$human_ensembl_gene,
                    TATAregulated_genes_pathway$human_ensembl_gene,
                    TGGAAAregulated_genes_pathway$human_ensembl_gene,
                    Retinoate_pathway$human_ensembl_gene,
                    Retinol_pathway$human_ensembl_gene,
                    Farnesly_diphosphate_pathway$human_ensembl_gene,
                    cholesterol_pathway$human_ensembl_gene,
                    cortisol_pathway$human_ensembl_gene,
                    steroid_hormone_pathway$human_ensembl_gene,
                    hypogonadism_pathway$human_ensembl_gene,
                    estradiol_responsive_genes_pathway$human_ensembl_gene,
                    Oocyte_maturation_pathway$human_ensembl_gene,
                    Vitamin_D_pathway$human_ensembl_gene,
                    urea_pathway$human_ensembl_gene
                    )

names(our_genesets) = c('TERT_pathway',
                        'MAPK_pathway',
                        'thyroidhormone_pathway',
                        'JAKSTAT_pathway',
                        'ferroptosis_pathway',
                        'translationInitiation_pathway',
                        'mehrEnergie_pathway',
                        #'glycolysis_pathway',
                        'citric_cycle_pathway',
                        'synapsis_pathway',
                        #'aminoacid_pathway',
                        'estrogen_pathway',
                        'ACE_pathway',
                        'inflam_response_pathway',
                        'CaMK_pathway',
                        'stemmcell_reg_pathway',
                        'mammarygland1_embryo_pathway',
                        'mammarygland2_puberty_pathway',
                        'mammarygland3_pregnancy_pathway',
                        'mammarygland4_abbau_pathway',
                        'TGFb EMT Genes UP',
                        'Caspase Cascade Pathway',
                        'NK cytotoxicity',
                        'Tcell activation',
                        'Bcell activation',
                        'purin_denovo_pathway',
                        'Interleukin_pathway',
                        'wnt_pathway',
                        'oxStress_pathway',
                        'IGF1_pathway',
                        'ACHR_pathway',
                        'HSP27_pathway',
                        'RAN_pathway',
                        'glycosylation_pathway',
                        'TATAregulated_genes_pathway',
                        'TGGAAAregulated_genes_pathway',
                        'Retinoate_pathway',
                        'Retinol_pathway',
                        'Farnesly_diphosphate_pathway',
                        'cholesterol_pathway',
                        'cortisol_pathway',
                        'steroid_hormone_pathway',
                        'hypogonadism_pathway',
                        'estradiol_responsive_genes_pathway', 
                        'Oocyte_maturation_pathway',
                        'Vitamin_D_pathway',
                        'urea_pathway')

#saving the list with all our genesets in our github folder data
save(our_genesets, file = '~/GitHub/2022-topic-02-team-05/data/our_genesets.RData')

#--------------------------------
#Extraktion aller ensembl ids und gennamen aus den exp daten
#--------------------------------
tcga_exp = readRDS("~/GitHub/2022-topic-02-team-05/data/tcga_tumor_log2TPM.RDS")

#extrhiren aller gene die in den Expressionsdaten vorkommen
tcga_genes = rownames(tcga_exp)

#dieser vetor enth?lt sowohl enseblm id als auch genenamen und muss daher gespalten werden
tcga_genes = strsplit(tcga_genes, split = '|', fixed = TRUE)

#speicher der ensembl ids ohne versionsnummer als eigenen vektor
tcga_geneids = sapply(tcga_genes, function(tcga_genes){return(tcga_genes[1])})
tcga_geneids = strsplit(tcga_geneids, split = '.', fixed = TRUE)
tcga_geneids = sapply(tcga_geneids, function(tcga_geneids){return(tcga_geneids[1])})

#speicher der genenames ohne versionsnummer als eigenen vektor
tcga_genenames = sapply(tcga_genes, function(tcga_genes){return(tcga_genes[2])})
tcga_genenames = strsplit(tcga_genenames, split = '.', fixed = TRUE)
tcga_genenames = sapply(tcga_genenames, function(tcga_genenames){return(tcga_genenames[1])})

tcga_genes = cbind.data.frame(tcga_geneids,tcga_genenames)
#speichern eines datframes der der die Ensembl ids und genenamen aller genen der exp daten enth?lt
save(tcga_genes, file = '~/GitHub/2022-topic-02-team-05/data/tcga_genes.RData')

#---------------------------------------------------
#Umschreiben der Hallmarkpathways von Gennamen in Ensembl ids 
#-----------------------------------------------------
#BiocManager::install("biomaRt")
library(biomaRt)
genesets = readRDS("~/GitHub/2022-topic-02-team-05/data/hallmarks_genesets.rds")

mart = useEnsembl(dataset = "hsapiens_gene_ensembl", biomart='ensembl')
genesets_ids = lapply(genesets[[1]], FUN = function(x){getBM(attributes = "ensembl_gene_id",
                                           filters = "external_gene_name",
                                           values = x,
                                           mart = mart )})
genesets_ids = sapply(genesets_ids, FUN = function(genesets_ids){return(as.vector(genesets_ids))})
save(genesets_ids, file = '~/GitHub/2022-topic-02-team-05/data/geneset_ids.RData')
    


