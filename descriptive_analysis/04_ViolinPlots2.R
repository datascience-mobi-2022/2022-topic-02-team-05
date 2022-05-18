#-----------------------------------
#f?r jedes gen p-ewrt berechnen vorher foldchange, log von foldchange xachse, y achse pwert, twosided ttest
#---------------------

tcga_tumor_norm = readRDS("~/GitHub/2022-topic-02-team-05/data/tcga_tumor_normal_datascience_proj_2022.RDS")
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

ttest.thca <- data.frame(expression = c(mean.thca.tumor, mean.thca.norm),
                         tissue.type = c(rep('tumor', length(mean.thca.tumor)), rep('normal', length(mean.thca.norm))))

t.test(tumor ~ normal,
       data = ttest.thca,
       alternative = 'two.sided',
       paired = TRUE)




