#----------------------------------------------
#In diesem Dokument werden die THCA Expressionsdaten nach hoher Varianz gereinigt und nach relevanten Biotypes filtriert
#----------------------------------------------

library(biomaRt)
library(dplyr)
library(ggplot2)


#------------------------------------------------
#Extrahieren, der für uns relevanten Daten in jeweils einen Vektor
#------------------------------------------------

tcga_tumor_normal_datascience_proj_2022 = readRDS("data/tcga_tumor_normal_datascience_proj_2022.RDS")

thca_tumor_exp <- tcga_tumor_normal_datascience_proj_2022$THCA$tumor
thca_normal_exp <- tcga_tumor_normal_datascience_proj_2022$THCA$normal
thca_anno <- tcga_tumor_normal_datascience_proj_2022$THCA$clinical
rm(tcga_tumor_normal_datascience_proj_2022)

save(thca_anno, file = 'data/thca_anno.RData')


#------------------------------------------------
#Generelles Preprocessing + nur die hochvarianten Gene werden behalten
#------------------------------------------------

thca_tumor_exp = na.omit(thca_tumor_exp) #keine NAs
thca_normal_exp = na.omit(thca_normal_exp) #keine NAs

#Berechnen der Varianz aller Gene
thca_tumor_exp_var = apply(thca_tumor_exp, 1, var)

#Plotten als Histogramm
hist(log(thca_tumor_exp_var), breaks = 50, probability = TRUE)

thca_tumor_exp_hvar = thca_tumor_exp[log(thca_tumor_exp_var) > -1.25, ] #Expressionsdaten der hochvarianten Gene
save(thca_tumor_exp_hvar, file = 'data/thca_tumor_exp_hvar.RData') #noch 15402 Gene

#gleiche Gene beim normal tissue löschen
thca_normal_exp_hvar = thca_normal_exp[log(thca_tumor_exp_var) > -1.25, ]
save(thca_normal_exp_hvar, file = 'data/thca_normal_exp_hvar.RData') #noch 15402 Gene


#-----------------------------------------------
#Rownames in Ensemble IDs umschreiben und Analyse der Biotypes der Gene
#-----------------------------------------------

#Funktion, die Gennamen nimmt und dafÜr den Gentypen/Biotyp gibt
mart = useEnsembl(dataset = "hsapiens_gene_ensembl", biomart='ensembl')
checkbiotypes = function(pathway){
  
  res = getBM(attributes = c("ensembl_gene_id", "gene_biotype"),
              filters = "ensembl_gene_id",
              values = pathway,
              mart = mart )
  message('Let me check that for you :D')
  
  return(res)
}

#Rownames von hoch varianten Genen in Ensemble IDs umschreiben
thca_tumor_genes_hvar = rownames(thca_tumor_exp_hvar)
thca_tumor_genes_hvar = strsplit(thca_tumor_genes_hvar, split = '|', fixed = TRUE)

#speichern der Ensemble IDs ohne Versionsnummer als eigenen Vektor
thca_geneids_hvar = sapply(thca_tumor_genes_hvar, function(thca_tumor_genes_hvar){return(thca_tumor_genes_hvar[1])})
thca_geneids_hvar = strsplit(thca_geneids_hvar, split = '.', fixed = TRUE)
thca_geneids_hvar = sapply(thca_geneids_hvar, function(thca_geneids_hvar){return(thca_geneids_hvar[1])})

#speichern der Gennamen ohne Versionsnummer als eigenen Vektor
thca_genenames_hvar = sapply(thca_tumor_genes_hvar, function(thca_tumor_genes_hvar){return(thca_tumor_genes_hvar[2])})
thca_genenames_hvar = strsplit(thca_genenames_hvar, split = '.', fixed = TRUE)
thca_genenames_hvar = sapply(thca_genenames_hvar, function(thca_genenames_hvar){return(thca_genenames_hvar[1])})

#Zusammenfügen in einen Dataframe von Gennamen und IDs
thca_genes_hvar = cbind.data.frame(thca_geneids_hvar,thca_genenames_hvar)

#IDs als rownames der Expressionsdaten von tumor und normal tissue
rownames(thca_tumor_exp_hvar) <- thca_geneids_hvar
rownames(thca_normal_exp_hvar) <- thca_geneids_hvar


#---------------------------------------------------------------------
#Löschen aller nicht-proteincodierenden pathways
#---------------------------------------------------------------------

#biotype-checking von den hochvarianten Genen über Ensemble IDs:
thca_biotypes_hvar = checkbiotypes(rownames(thca_normal_exp_hvar))
thca_protein_hvar = thca_biotypes_hvar$ensembl_gene_id[thca_biotypes_hvar$gene_biotype == 'protein_coding']

#Tumor tissue: biotypes, die nicht proteincodierend sind, rauswerfen:
thca_tumor_exp_cleaned = thca_tumor_exp_hvar[rownames(thca_tumor_exp_hvar) %in% thca_protein_hvar,]

#Normal tissue: biotypes genauso rauswerfen:
thca_normal_exp_cleaned <- thca_normal_exp_hvar[rownames(thca_tumor_exp_hvar) %in% thca_protein_hvar,]

thca_tumor_exp_hvar[sapply(rownames(thca_tumor_exp_hvar), function(thca_tumor_exp_hvar){isin(thca_protein_hvar, thca_tumor_exp_hvar)})] #Was macht diese Zeile(Anna)

#prüfen auf NAs
sum(is.na(thca_tumor_exp_cleaned))
sum(is.na(thca_normal_exp_cleaned))

#Speichern der gecleanten Expressionsdaten: mit Ensemble IDs als Zeilennamen
save(thca_tumor_exp_cleaned, file = 'data/thca_tumor_exp_cleaned.RData')
save(thca_normal_exp_cleaned, file = 'data/thca_normal_exp_cleaned.RData')


#--------------------------------
#Extraktion aller Ensemble IDs und Gennamen aus den gecleanten Expressionsdaten
#--------------------------------

genes_cleaned_indexes =  thca_genes_hvar$thca_geneids_hvar %in% rownames(thca_tumor_exp_cleaned)
thca_genes_cleaned = thca_genes_hvar[genes_cleaned_indexes == TRUE,]

#speichern eines dataframes, der die Ensemble IDS und Gennamen aller Gene der Expressionsdaten enthält
save(thca_genes_cleaned, file = 'data/tcga_genes_cleaned.RData')













