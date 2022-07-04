#-----------------------------------
#f?r jedes gen p-ewrt berechnen vorher foldchange, log von foldchange xachse, y achse pwert, twosided ttest
#---------------------
library(ggplot2)

tcga_tumor_norm = readRDS("data/tcga_tumor_normal_datascience_proj_2022.RDS")
thca.tumor <- tcga_tumor_norm$THCA$tumor
thca.norm <- tcga_tumor_norm$THCA$normal


#--------------------------------------------------------------------------------
#Preproscessing
#--------------------------------------------------------------------------------

#Na#s rauswerfen
gc() #gibt arbeitsspeciher frei der f?r die gro?en datenmengen gebraucht wird
thca_tumor_normal <- tcga_tumor_norm$THCA
thca.norm <- thca_tumor_normal$normal
thca.tumor <- thca_tumor_normal$tumor
thca.norm <- na.omit(thca.norm)
thca.tumor <- na.omit(thca.tumor)

#Berrechnen der Varianz
thca.tumor.var = apply(thca.tumor, 1, var)
thca.norm.var = apply(thca.norm, 1, var)

#PLotten als Histogramm: Logarythmus, damit man was sieht: sonst alle sehr kleine Varianzen!
hist(log(thca.tumor.var), breaks = 50, probability = TRUE)
hist(log(thca.norm.var), breaks = 50, probability = TRUE)

#cutting der gene mit sehr niedriger exression d.h. log(var) < -1
#erstmal willk?rlich festgestzt 
#die Werte d?rfen ja nur gel?scht werden wenn sie in beiden eine niedrige Varianz haben

thca.tumor.v = thca.tumor[log(thca.tumor.var) > -1 | log(thca.norm.var) > -1, ]
thca.norm.v = thca.norm[log(thca.tumor.var) > -1 | log(thca.norm.var) > -1, ]
thca.tumor.va = na.omit(thca.tumor.v)
thca.norm.va = na.omit(thca.norm.v)


#--------------------------------
#Mittelwert berechnen
#--------------------------------

mean.thca.tumor = apply(thca.tumor.va, 1, mean)
mean.thca.norm = apply(thca.norm.va, 1, mean)


#-----------------------------------------------------------------------------------
#Foldchange(FC) besrechnen!
#den log2FC_gene bekommen wir durch mean(condition1) - mean(condition2) 
#-----------------------------------------------------------------------------------

log2fc.thca = mean.thca.norm - mean.thca.tumor

#Leeren Vektor definieren:
p.values <- c()

#for-loop: ttest-paired durchführen über alle Gene zwischen normal und tumor tissue
for (i in (1:nrow(thca.norm.va))){
  x <- t.test(thca.norm.va[i,], thca.tumor.va[i,], alternative = 'two.sided')$p.value
  p.values <- append(p.values, x)
}

#-----------------------------
#Bonferroni-Adjustment, um alphafehlerkummulation zu vermeiden
#-----------------------------

n = nrow(thca.norm.va)
bf = 1/n

#adjusted significanzlevel
alpha = 0.025
alpha.kor = alpha*bf

#gene, die sich signifikant von alpha unterscheiden in tumor von normal tissue
fc.sig <- log2fc.thca[p.values < alpha.kor]
pv.sig <- p.values[p.values < alpha.kor]


#-----------------------------------------
#Erstellen eines dataframes
#-----------------------------------------

#Erstellen dataframe
data.thca = data.frame(log2fc.thca, p.values)

#rownames in Ensemble IDs 
thca_genes = rownames(data.thca)
thca_genes = strsplit(thca_genes, split = '|', fixed = TRUE)

thca_geneids = sapply(thca_genes, function(thca_genes){return(thca_genes[1])})
thca_geneids = strsplit(thca_geneids, split = '.', fixed = TRUE)
thca_geneids = sapply(thca_geneids, function(thca_geneids){return(thca_geneids[1])})

thca_genenames = sapply(thca_genes, function(thca_genes){return(thca_genes[2])})
thca_genenames = strsplit(thca_genenames, split = '.', fixed = TRUE)
thca_genenames = sapply(thca_genenames, function(thca_genenames){return(thca_genenames[1])})

thca_genes = cbind.data.frame(thca_geneids,thca_genenames)

rownames(data.thca) <- thca_geneids 

#neue spalte mit Genenames
cbind.data.frame(data.thca, thca_genenames) -> data.thca

#hinzufügen einer Spalte, die sagt, ob das Gen up- oder downregulated wird
#hinzufügen einer Spalte diffexpressed mit NOs 
data.thca$diffexpressed <- "NO"
#wenn log2Foldchange > 0.1 and pvalue < alpha.kor, set as "UP" 
data.thca$diffexpressed[data.thca$log2fc.thca > 0.1 & data.thca$p.value < alpha.kor] <- "UP"
# if log2Foldchange < -0.1 and pvalue < 0.05, set as "DOWN"
data.thca$diffexpressed[data.thca$log2fc.thca < -0.1 & data.thca$p.value < alpha.kor] <- "DOWN"


#-----------------------------------------------
#Erstellen des Volcanoplots
#-----------------------------------------------

volcano <- ggplot(data = data.thca, aes(x = log2fc.thca, y = -log10(p.values), color = diffexpressed, label = thca_genenames)) +
  geom_point() + 
  theme_minimal() + 
  ggtitle("Volcanoplot") + 
  geom_hline(yintercept = -log10(alpha.kor), col = "black", show.legend = TRUE) + 
  geom_text(data = subset(data.thca, -log10(p.values) > 22), size = 2, check_overlap = TRUE, nudge_x = 0.1, nudge_y = 0.8, color = "black")

volcano

