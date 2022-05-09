#-----------------------------------------------------
#Dokument in dem wir alle Metabolischen pathways laden
#-----------------------------------------------------

#install.packages("msigdbr")
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")

library(msigdbr)

#öffnen der Datenbank
gsea_pathways = msigdbr(species = "Homo sapiens")

#Auswahl der einzelnen pathways
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
RAN_pathway = gsea_pathways[gsea_pathways$gs_name == "BIOCARTA_RAN_PATHWAY", ]
glycosylation_pathway = gsea_pathways[gsea_pathways$gs_name == "BIOCARTA_AMAN_PATHWAY", ]




#liste aller pathways mit allen inofs zum jeweiligen pathway
our_genesets = list(TERT_pathway$gene_symbol,
                    MAPK_pathway$gene_symbol,
                    thyroidhormone_pathway$gene_symbol,
                    JAKSTAT_pathway$gene_symbol,
                    ferroptosis_pathway$gene_symbol,
                    translationInitiation_pathway$gene_symbol,
                    mehrEnergie_pathway$gene_symbol,
                    glycolysis_pathway$gene_symbol,
                    citric_cycle_pathway$gene_symbol,
                    synapsis_pathway$gene_symbol,
                    aminoacid_pathway$gene_symbol,
                    estrogen_pathway$gene_symbol,
                    ACE_pathway$gene_symbol,
                    inflam_response_pathway$gene_symbol,
                    CaMK_pathway$gene_symbol,
                    stemmcell_reg_pathway$gene_symbol,
                    mammarygland1_embryo_pathway$gene_symbol,
                    mammarygland2_puberty_pathway$gene_symbol,
                    mammarygland3_pregnancy_pathway$gene_symbol,
                    mammarygland4_abbau_pathway$gene_symbol,
                    TGFBEMT_pathway$gene_symbol,
                    CaspaseCascade_pathway$gene_symbol,
                    NKcytotoxicity$gene_symbol,
                    Tcellactivation$gene_symbol,
                    Bcellactivation$gene_symbol,
                    purin_denovo_pathway$gene_symbol,
                    Interleukin_pathway$gene_symbol,
                    wnt_pathway$gene_symbol,
                    oxStress_pathway$gene_symbol,
                    IGF1_pathway$gene_symbol,
                    ACHR_pathway$gene_symbol,
                    HSP27_pathway$gene_symbol,
                    RAN_pathway$gene_symbol,
                    glycosylation_pathway$gene_symbol
)
names(our_genesets) = c('TERT_pathway',
                        'MAPK_pathway',
                        'thyroidhormone_pathway',
                        'JAKSTAT_pathway',
                        'ferroptosis_pathway',
                        'translationInitiation_pathway',
                        'mehrEnergie_pathway',
                        'glycolysis_pathway',
                        'citric_cycle_pathway',
                        'synapsis_pathway',
                        'aminoacid_pathway',
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
                        'glycosylation_pathway'
)

#saving the list with all our genesets in our github folder data
save(our_genesets, file = 'data/our_genesets.RData')
