#----------------------------------------------
#In diesem Dokument werden die TCGA Expressionsdaten nach hoher Varianz gereinigt und die Biotypes analysiert
#----------------------------------------------

#BiocManager::install("biomaRt")
library(biomaRt)
library(dplyr)
library(ggplot2)

#laden der daten
tcga_exp = readRDS("data/tcga_tumor_log2TPM.RDS")
load("data/geneset_ids.RData")
load("data/our_genesets.RData")

#----------------------------------------------
#Cleaning der TCGA expressionsdaten
#----------------------------------------------
#checking for NAs
#erstmal nur in den ersten dreitausend genen weil ich sonst fehler bei der varianz bekomm ):
gc() #gibt arbeitsspeciher frei der f?r die gro?en datenmengen gebraucht wird
tcga_exp_narm = na.omit(tcga_exp) #keine na's vorhanden!

#Berrechnen der Varianz aller Gene
tcga_exp_var = apply(tcga_exp_narm, 1, var)

#PLotten als Histogramm
hist(log(tcga_exp_var), breaks = 50, probability = TRUE)

#speichern als plot
#funktioniert noch nicht
# savePlot(filename = "output/tcga_exp_genevariance",
#          device = dev.cur(),
#          type = "jpg")

#cutting der gene mit sehr niedriger exression d.h. log(var) < -1
#erstmal willk?rlich festgestzt
tcga_exp_hvar = tcga_exp_narm[log(tcga_exp_var) > -1, ] #Expressionsdaten der hochvarianten Gene!
save(tcga_exp_hvar, file = 'data/tcga_exp_hvar.RData')

#-----------------------------------------------
# Analyse der Biotypes von Unseren pathways, Dr. Herrmanns pathways und den tcga expressionsdaten
#-----------------------------------------------
#function die Gennamen nimmt und daf?r den gentypen gibt
mart = useEnsembl(dataset = "hsapiens_gene_ensembl", biomart='ensembl')
checkbiotypes = function(pathway){

res = getBM(attributes = c("ensembl_gene_id", "gene_biotype"),
            filters = "ensembl_gene_id",
            values = pathway,
            mart = mart )
message('Let me check that for you :D')

return(res)
}

#-----------------------------------------
#gene biotypes unserer Metabolic pathways bestimmen
#-----------------------------------------
our_genesets_biotypes = sapply(our_genesets, checkbiotypes)
our_genesets_biotypes = our_genesets_biotypes[2,]

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
genesets_biotypes = sapply(genesets_ids, checkbiotypes)
genesets_biotypes = genesets_biotypes[2,]


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
load('data/tcga_genes.RData')
tcga_biotypes = checkbiotypes(tcga_genes$tcga_geneids)
tcga_biotypes_only = tcga_biotypes$gene_biotype

#darstellen aller biotypes im set
x = as.data.frame(table(as.vector(tcga_biotypes_only)))
colnames(x) = c('Biotype','Ammount')

ggplot(x, aes(Biotype, Ammount)) + geom_bar(stat = 'identity') +  
  labs(title = 'TCGA expression data biotypes') +
  theme(axis.text = element_text(angle = 90,))


#--------------------------------------------------------------------------
#entfernen aller nicht proteincodierenden gene aus der tcga expressions matrix
#--------------------------------------------------------------------------
load('data/tcga_exp_hvar.RData')

#Rownames von high variant genes in ensemble IDs umschreiben:
tcga_genes_hvar = rownames(tcga_exp_hvar)
tcga_genes_hvar = strsplit(tcga_genes_hvar, split = '|', fixed = TRUE)

#speicher der ensembl ids ohne versionsnummer als eigenen vektor
tcga_geneids_hvar = sapply(tcga_genes_hvar, function(tcga_genes_hvar){return(tcga_genes_hvar[1])})
tcga_geneids_hvar = strsplit(tcga_geneids_hvar, split = '.', fixed = TRUE)
tcga_geneids_hvar = sapply(tcga_geneids_hvar, function(tcga_geneids_hvar){return(tcga_geneids_hvar[1])})

#IDs als rownames:
rownames(tcga_exp_hvar) <- tcga_geneids_hvar

#biotype-checking von den hochvarianten Genen über ID:
tcga_biotypes_hvar = checkbiotypes(rownames(tcga_exp_hvar))
tcga_protein_hvar = tcga_biotypes_hvar$ensembl_gene_id[tcga_biotypes_hvar$gene_biotype == 'protein_coding']

#biotypes, die nicht proteincodierend sind, rauswerfen:
tcga_exp_cleaned = tcga_exp_hvar[rownames(tcga_exp_hvar) %in% tcga_protein_hvar,]
tcga_exp_hvar[sapply(rownames(tcga_exp_hvar), function(tcga_exp_hvar){isin(tcga_protein_hvar, tcga_exp_hvar)})]

#checking for NA's
sum(is.na(tcga_exp_cleaned)) #keine NAs waren vorhanden!

#--------------------------------
#Extraktion aller ensembl ids und gennamen aus den cleanen exp daten und
#Umbennen der expressionsdaten zeilennamen in nur ids
#--------------------------------
genes_cleaned_indexes =  tcga_genes$tcga_geneids %in% rownames(tcga_exp_cleaned)
tcga_genes_cleaned = tcga_genes[genes_cleaned_indexes == TRUE,]

#speichern eines datframes der der die Ensembl ids und genenamen aller genen der exp daten enth?lt
save(tcga_genes_cleaned, file = 'data/tcga_genes_cleaned.RData')

#Speichern der gecleanten Expressionsdaten: mit Ensemble IDs als Zeilennamen
save(tcga_exp_cleaned, file = 'data/tcga_exp_cleaned.RData')




