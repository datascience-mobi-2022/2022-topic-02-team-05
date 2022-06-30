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




#-----------------------------
#Mittelwert berechnen
#--------------------------------
mean.thca.tumor = apply(thca.tumor.va, 1, mean)
mean.thca.norm = apply(thca.norm.va, 1, mean)

#-----------------------------------------------------------------------------------
#Foldchange(FC) besrechnen!
#den log2FC_gene bekommen wir durch mean(condition1) - mean(condition2) 
#---------------------------------------------------------------------------------
log2fc.thca = mean.thca.norm - mean.thca.tumor


#--------------------------------------
#t-test paired durchführen: nicht auf den means, sondern auf den gecleanten daten
#--------------------------------------------
#as.vector(thca.norm.va[1,1:ncol(thca.norm.va)]), as.vector(thca.tumor.va[1,1:ncol(thca.tumor.va)]
#means in einen Datenframe packen:

#ttest.thca <- data.frame(expression = c(mean.thca.tumor, mean.thca.norm),
 #                        tissue.type = c(rep('tumor', length(mean.thca.tumor)), rep('normal', length(mean.thca.norm))))

#t.test(tumor ~ normal,
  #     data = ttest.thca,
   #    alternative = 'two.sided',
    #   paired = TRUE)
#wilcox.test(thca.norm.va[,1], thca.tumor.va[,1], paired = TRUE)

#for (mean in 1:length(mean.thca.norm)) {
#  x <- wilcox.test(mean.thca.norm[mean], mean.thca.tumor[mean], paired = TRUE)
#  print(x)
#}

#x = t.test(thca.norm.va[1,], thca.tumor.va[1,], alternative = 'two.sided')$p.value


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

#adjusted significanzlevel: ehemals: alpha = 0.05
alpha = 0.025
alpha.kor = alpha*bf

#gene, die sich signifikant (mit korrigiertem alpha) unterscheiden in tumor von normal tissue
fc.sig <- log2fc.thca[p.values < alpha.kor]
pv.sig <- p.values[p.values < alpha.kor]


#Volcanoplot
plot(fc.sig, -log(pv.sig))
plot(log2fc.thca, -log(p.values), 
     xlab='log(foldchange)',
     ylab = '-log10(Pvalues)')

data.thca = data.frame(log2fc.thca, p.values) #die 2 Vektoren f?r unseren Volcano PLot werden in einen df gepackt, damit daraus ein plot erstellt werden kann

data.thca = data.frame(log2fc.thca, p.values) #die 2 Vektoren f?r unseren Volcano PLot werden in einen df gepackt, damit daraus ein plot erstellt werden kann
with(data.thca, plot(data.thca$log2fc.thca, -log10(data.thca$p.values), main = "Volcano plot", xlab='log2(foldchange)',
                     ylab = '-log10(Pvalues)', pch = 20))

with(subset(data.thca, log2fc.thca>alpha.kor & p.values < alpha.kor), points(log2fc.thca, -log10(p.values), col="orange", pch = 20))
with(subset(data.thca, log2fc.thca< -alpha.kor & p.values < alpha.kor), points(log2fc.thca, -log10(p.values), col="green", pch = 20))
with(subset(data.thca, abs(log2fc.thca) < alpha.kor | p.values > alpha.kor), points(log2fc.thca, -log10(p.values), pch=19, col="gray"))

#-----------------------
#noch mehr preprocessing, die rownames werden zu den Ids
#-----------------------
#extrhiren aller gene die in den Expressionsdaten vorkommen
thca_genes = rownames(data.thca)

#dieser vetor enth?lt sowohl enseblm id als auch genenamen und muss daher gespalten werden
thca_genes = strsplit(thca_genes, split = '|', fixed = TRUE)

#speicher der ensembl ids ohne versionsnummer als eigenen vektor
thca_geneids = sapply(thca_genes, function(thca_genes){return(thca_genes[1])})
thca_geneids = strsplit(thca_geneids, split = '.', fixed = TRUE)
thca_geneids = sapply(thca_geneids, function(thca_geneids){return(thca_geneids[1])})

#speicher der genenames ohne versionsnummer als eigenen vektor
thca_genenames = sapply(thca_genes, function(thca_genes){return(thca_genes[2])})
thca_genenames = strsplit(thca_genenames, split = '.', fixed = TRUE)
thca_genenames = sapply(thca_genenames, function(thca_genenames){return(thca_genenames[1])})

thca_genes = cbind.data.frame(thca_geneids,thca_genenames)

rownames(data.thca) <- thca_geneids

#vergleichen von dem pathway wo welche geneid vorkommt, den pathway kann man jetzt beliebig ersetzen
random.pathway = match(genesets_ids$Apop_SURVIVAL.ensembl_gene_id, rownames(data.thca))
random.pathway.narm = na.omit(random.pathway)
data.thca[random.pathway.narm,3] = "blue"

with(data.thca, plot(data.thca$log2fc.thca, -log10(data.thca$p.values), main = "Volcano plot", xlab='log2(foldchange)',
                     ylab = '-log10(Pvalues)', pch = 20))
abline(h = -log10(alpha.kor))

with(subset(data.thca, data.thca$V3 == "blue" ), points(log2fc.thca, -log10(p.values), col="red", pch = 20))

#install.packages('calibrate')
#library(calibrate)


#with(subset(data.thca, p.values < alpha.kor/10000000000000000 & abs(log2fc.thca)>alpha.kor/1000000, textxy(log2fc.thca, -log10(p.values)), labs=rownames(data.thca), cex=.8))
#hiermit h?tten wir das noch beschriften k?nnen, hab ich aber nicht hinbekommen


#fc.pv.sig <- data.frame(fc.sig, pv.sig)

#Volcano_logpvalues = 
#  ggplot(fc.pv.sig, aes(x, y)) + 
#  geom_point(shape = 21, size = 1.5, fill = "red") +
#  ggtitle("Variance over Mean") + 
#  xlab("Mean") + 
#  ylab("Standardized Variance")
#geom_text(check_overlap = TRUE)

#annotate("text", x = -10:15, y = 5:10.5, label = RownamesTCGA)
#ggplot(Mea_clean, stand_var_clean, "p", pch = 19, main = "Variance over mean", xlab = "Mean", ylab = "standardized_Variance") -> descriptiveMeanVar 

#Punkte ab gewissem Schwellenwert benennen 
#if (y > 5) {
#  text(Mea_clean, stand_var_clean, labels = row.names(tcga_exp_cleaned), cex = 0.7, pos = 3)
#}
#Var_over_Mean_plot

#Histogramm der p.values
hist(p.values)
abline(v =alpha.kor)




