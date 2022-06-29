#----------------------------------------------
#In diesem Dokument wird die Varianz und der Mean für jedes Gen für jeden Tumortype geplotted
#----------------------------------------------

#Laden der bereits gecleanten Daten aus der großen tcga matrix 
load('data/tcga_exp_cleaned.RData')

#load ggplot2 
library(ggplot2)

#Variance of every gene
Var_clean <- apply(tcga_exp_cleaned, 1, var)
y1 = as.vector(Var_clean)
#mean of every gene 
Mea_clean <- apply(tcga_exp_cleaned, 1, mean)
x1 = as.vector(Mea_clean)

#Standardisieren der Varianz!
stand_var_clean <- scale(Var_clean)
#y = as.vector(stand_var_clean)

#data frame out of Mea_clean and stand_var_clean
as.vector(row.names(tcga_exp_cleaned)) -> names1
data.frame(x1, y1, names1) -> df_VoverM

#plotten der Werte mit Varianz über Mittelwert 
Var_over_Mean_plot = 
  ggplot(df_VoverM, aes(x1, y1, label = names1)) + 
  geom_point(shape = 21, size = 1.5, fill = "red") +
  geom_text(data = subset(df_VoverM, y1 > 33), size = 2, check_overlap = TRUE, nudge_x = 2, nudge_y = 1) + 
  ggtitle("Variance over Mean cleaned matrix") + 
  xlab("Mean") + 
  ylab("Variance") 
 
Var_over_Mean_plot

#Gennamen der hochvarianten (>33) Gene in einer Liste 
highVariance_cleaned <- subset(df_VoverM, y1 > 33)
highVariance_cleaned <- highVariance_cleaned$names1

#BiocManager::install("AnnotationDbi")
#BiocManager::install("org.Hs.eg.db")
library("AnnotationDbi")
library("org.Hs.eg.db")

#columns(org.Hs.eg.db)

MAPIDS1 <- mapIds(org.Hs.eg.db, keys = highVariance_cleaned, keytype = "ENSEMBL", column = "ALIAS")
as.data.frame(MAPIDS1) -> highVariancegenes_cleaned


#-------------------------------------------------------------
#Laden der  ungecleanten Daten aus der großen tcga matrix 
#-------------------------------------------------------------

load('data/tcga_tumor_log2TPM.RDS')
load("data/tcga_genes")

tcga_genes$tcga_geneids -> col1
rownames(tcga_tumor_log2TPM) <- col1
tcga_tumor_log2TPM

#Variance of every gene
Var_clean_big <- apply(tcga_tumor_log2TPM, 1, var)
y2 = as.vector(Var_clean_big)
#mean of every gene 
Mea_clean_big <- apply(tcga_tumor_log2TPM, 1, mean)
x2 = as.vector(Mea_clean_big)

#Standardisieren der Varianz
stand_var_clean <- scale(Var_clean_big)
#y = as.vector(stand_var_clean)

#data frame out of Mea_clean and stand_var_clean
as.vector(row.names(tcga_tumor_log2TPM)) -> names2
data.frame(x2, y2, names2) -> df_VoverM_big

#plotten der Werte mit Varianz über Mittelwert 
Var_over_Mean_plot_big = 
  ggplot(df_VoverM_big, aes(x2, y2, label = names2)) + 
  geom_point(shape = 21, size = 1.5, fill = "green") +
  geom_text(data = subset(df_VoverM_big, y2 > 33), size = 2, check_overlap = TRUE, nudge_x = 2, nudge_y = 1) + 
  ggtitle("Variance over Mean big matrix") + 
  xlab("Mean") + 
  ylab("Variance") 

Var_over_Mean_plot_big

#Gennamen der hochvarianten (>33) Gene in einer Liste 
highVariance_uncleaned <- subset(df_VoverM_big, y1 > 33)
highVariance_uncleaned <- highVariance_uncleaned$names2

#BiocManager::install("AnnotationDbi")
#BiocManager::install("org.Hs.eg.db")
library("AnnotationDbi")
library("org.Hs.eg.db")

#columns(org.Hs.eg.db)

MAPIDS2 <- mapIds(org.Hs.eg.db, keys = highVariance_uncleaned, keytype = "ENSEMBL", column = "ALIAS")
as.data.frame(MAPIDS2) -> highVariancegenes_uncleaned
