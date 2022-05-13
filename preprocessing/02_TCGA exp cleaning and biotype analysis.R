#----------------------------------------------
#In diesem Dokument werden die TCGA Expressionsdaten nach hoher Varianz gereinigt und die Biotypes analysiert
#----------------------------------------------

#BiocManager::install("biomaRt")
library(biomaRt)
library(dplyr)
library(ggplot2)

#laden der daten
tcga_exp = readRDS("~/GitHub/2022-topic-02-team-05/data/tcga_tumor_log2TPM.RDS")
genesets_ids = readRDS("~/GitHub/2022-topic-02-team-05/data/genesets_ids.rds")
load("~/GitHub/2022-topic-02-team-05/data/our_genesets.RData")

#----------------------------------------------
#Cleaning der TCGA expressionsdaten
#----------------------------------------------
#checking for NAs
#erstmal nur in den ersten dreitausend genen weil ich sonst fehler bei der varianz bekomm ):
gc() #gibt arbeitsspeciher frei der für die großen datenmengen gebraucht wird
tcga_exp_narm = na.omit(tcga_exp)

#Berrechnen der Varianz aller Gene
tcga_exp_var = apply(tcga_exp_narm, 1, var)

#PLotten als Histogramm
hist(log(tcga_exp_var), breaks = 50, probability = TRUE)

#speichern als plot
#funktioniert noch nicht
savePlot(filename = "~/GitHub/2022-topic-02-team-05/output/tcga_exp_genevariance",
         device = dev.cur(),
         type = "jpg")

#cutting der gene mit sehr niedriger exression d.h. log(var) < -1
#erstmal willkürlich festgestzt
tcga_exp_hvar = tcga_exp_narm[log(tcga_exp_var) > -1, ]

#wdh für die nächsten gene
tcga_exp_narm_2 = na.omit(tcga_exp[3001:6000,])
tcga_exp_var_2 = apply(tcga_exp_narm_2, 1, var)
tcga_exp_hvar_2 = tcga_exp_narm_2[log(tcga_exp_var_2) > -1, ]

tcga_exp_hvar = rbind(tcga_exp_hvar, tcga_exp_hvar_2)

#-----------------------------------------------
# Analyse der Biotypes von Unseren pathways, Dr. Herrmanns pathways und den tcga expressionsdaten
#-----------------------------------------------
#function die Gennamen nimmt und dafür den gentypen gibt
checkbiotypes = function(pathway){mart = useEnsembl(dataset = "hsapiens_gene_ensembl", biomart='ensembl')

res = getBM(attributes = c("ensembl_gene_id", "gene_biotype"),
            filters = "external_gene_name",
            values = pathway,
            mart = mart )
message('Let me check that for you :D')

return(res)
}

#function die alle biotypes eines pathways als barplot plottet
#genelist ist dabei die liste unserer pathways, number die stelle des pathways in der liste
plotbiotypes = function(genelist, number){
  x = as.data.frame(table(as.vector(genelist[[number]])))
  colnames(x) = c('Biotype','Ammount')
  
  ggplot(x, aes(Biotype, Ammount)) + geom_bar(stat = 'identity') +  
    labs(title = names(genelist[number])) +
    theme(axis.text = element_text(angle = 0,))
}


#-----------------------------------------
#gene biotypes unserer Metabolic pathways bestimmen
#-----------------------------------------
our_genesets_biotypes = sapply(our_genesets, checkbiotypes)
our_genesets_biotypes = our_genesets_biotypes[2,]#entfernen der ensembl ids

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
genesets_biotypes = sapply(genesets[[1]], checkbiotypes)
genesets_biotypes = genesets_biotypes[2,]#entfernen der ids


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
tcga_exp_cleaned = tcga_exp_hvar[tcga_biotypes$gene_biotype == 'protein_coding', ]
tcga_exp_cleaned = na.omit(tcga_exp_cleaned)
##hier geht sich das iwie nicht mit den zeilen aus wieso kommen wir von 1000 auf 800 auf 300 auf 140

save(our_genesets_biotypes, file = '~/GitHub/2022-topic-02-team-05/data/our_genesets_biotypes.RData')
save(genesets_biotypes, file = '~/GitHub/2022-topic-02-team-05/data/genesets_biotypes.RData')
save(tcga_biotypes, file = '~/GitHub/2022-topic-02-team-05/data/tcga_biotypes.RData')
save(tcga_exp_cleaned, file = '~/GitHub/2022-topic-02-team-05/data/tcga_exp_cleaned.RData')

