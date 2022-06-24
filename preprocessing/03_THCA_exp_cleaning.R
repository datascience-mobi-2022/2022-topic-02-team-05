#----------------------------------------------
#In diesem Dokument werden die THCA Expressionsdaten nach hoher Varianz gereinigt und nach relevanten Biotypes filtriert
#----------------------------------------------

library(biomaRt)
library(dplyr)
library(ggplot2)

#------------------------------------------------
#Extrahieren nur der für uns relevanten Daten in jeweils einem Vektor
#------------------------------------------------

load("data/tcga_tumor_normal_datascience_proj_2022")

thca_tumor_exp <- tcga_tumor_normal_datascience_proj_2022$THCA$tumor
thca_normal_exp <- tcga_tumor_normal_datascience_proj_2022$THCA$normal
thca_anno <- tcga_tumor_normal_datascience_proj_2022$THCA$clinical

save(thca_anno, file = 'data/thca_anno')

#------------------------------------------------
# Generelles Preprocessing und nur die hochvarianten Gene behalten!
#------------------------------------------------

thca_tumor_exp = na.omit(thca_tumor_exp) #keine NA's
thca_normal_exp = na.omit(thca_normal_exp) #keine NA's

#Berrechnen der Varianz aller Gene
thca_tumor_exp_var = apply(thca_tumor_exp, 1, var)

#PLotten als Histogramm
#hist(thca_tumor_exp_var, breaks = 50, probability = TRUE) #man sieht keine schöne Verteilung!
hist(log(thca_tumor_exp_var), breaks = 50, probability = TRUE)
#---------------------------------------
thca_tumor_exp_hvar = thca_tumor_exp[log(thca_tumor_exp_var) > -1.25, ] #Expressionsdaten der hochvarianten Gene
save(thca_tumor_exp_hvar, file = 'data/thca_tumor_exp_hvar.RData') #noch 15402 Gene

#gleiche Gene Beim normal tissue rauswerfen:
thca_normal_exp_hvar = thca_normal_exp[log(thca_tumor_exp_var) > -1.25, ]
save(thca_normal_exp_hvar, file = 'data/thca_normal_exp_hvar.RData') #noch 15402 Gene

#-----------------------------------------------
# Rownames in Ensemble/Ids umschreiben und Analyse der Biotypes der Gene:
#-----------------------------------------------
#function die Gennamen nimmt und dafÜr den gentypen/Biotype gibt
mart = useEnsembl(dataset = "hsapiens_gene_ensembl", biomart='ensembl')
checkbiotypes = function(pathway){
  
  res = getBM(attributes = c("ensembl_gene_id", "gene_biotype"),
              filters = "ensembl_gene_id",
              values = pathway,
              mart = mart )
  message('Let me check that for you :D')
  
  return(res)
}


#Rownames von high variant genes in ensemble IDs umschreiben:
thca_tumor_genes_hvar = rownames(thca_tumor_exp_hvar)
thca_tumor_genes_hvar = strsplit(thca_tumor_genes_hvar, split = '|', fixed = TRUE)

#speicher der ensembl ids ohne versionsnummer als eigenen vektor
thca_geneids_hvar = sapply(thca_tumor_genes_hvar, function(thca_tumor_genes_hvar){return(thca_tumor_genes_hvar[1])})
thca_geneids_hvar = strsplit(thca_geneids_hvar, split = '.', fixed = TRUE)
thca_geneids_hvar = sapply(thca_geneids_hvar, function(thca_geneids_hvar){return(thca_geneids_hvar[1])})

#speicher der genenames ohne versionsnummer als eigenen vektor
thca_genenames_hvar = sapply(thca_tumor_genes_hvar, function(thca_tumor_genes_hvar){return(thca_tumor_genes_hvar[2])})
thca_genenames_hvar = strsplit(thca_genenames_hvar, split = '.', fixed = TRUE)
thca_genenames_hvar = sapply(thca_genenames_hvar, function(thca_genenames_hvar){return(thca_genenames_hvar[1])})

#Zusammenfügen in einen Dataframe von genenames und IDs
thca_genes_hvar = cbind.data.frame(thca_geneids_hvar,thca_genenames_hvar)

#IDs als rownames der Expressionsdaten von tumor und normal tissue
rownames(thca_tumor_exp_hvar) <- thca_geneids_hvar
rownames(thca_normal_exp_hvar) <- thca_geneids_hvar

#---------------------------------------------------------------------
# da wir nur proteincodierende Gene in unseren Pathways haben, werden alle anderen rausgeschmissen

#biotype-checking von den hochvarianten Genen über ID:
thca_biotypes_hvar = checkbiotypes(rownames(thca_normal_exp_hvar))
thca_protein_hvar = thca_biotypes_hvar$ensembl_gene_id[thca_biotypes_hvar$gene_biotype == 'protein_coding']

#Tumor tissue: biotypes, die nicht proteincodierend sind, rauswerfen:
thca_tumor_exp_cleaned = thca_tumor_exp_hvar[rownames(thca_tumor_exp_hvar) %in% thca_protein_hvar,]

#Normal tissue Biotypes genauso rauswerfen:
thca_normal_exp_cleaned <- thca_normal_exp_hvar[rownames(thca_tumor_exp_hvar) %in% thca_protein_hvar,]

thca_tumor_exp_hvar[sapply(rownames(thca_tumor_exp_hvar), function(thca_tumor_exp_hvar){isin(thca_protein_hvar, thca_tumor_exp_hvar)})] #Was macht diese Zeile(Anna)

#checking for NA's
sum(is.na(thca_tumor_exp_cleaned)) #keine NAs waren vorhanden!??????????????????????????
sum(is.na(thca_normal_exp_cleaned))

#Speichern der gecleanten Expressionsdaten: mit Ensemble IDs als Zeilennamen
save(thca_tumor_exp_cleaned, file = 'data/thca_tumor_exp_cleaned.RData')
save(thca_normal_exp_cleaned, file = 'data/thca_normal_exp_cleaned.RData')

#--------------------------------
# Extraktion aller ensembl ids und gennamen aus den gecleanen Expressionsdaten
#--------------------------------

genes_cleaned_indexes =  thca_genes_hvar$thca_geneids_hvar %in% rownames(thca_tumor_exp_cleaned)
thca_genes_cleaned = thca_genes_hvar[genes_cleaned_indexes == TRUE,]

#speichern eines dataframes, der die Ensembl ids und genenamen aller gene der exp daten enth?lt
save(thca_genes_cleaned, file = 'data/tcga_genes_cleaned.RData')













