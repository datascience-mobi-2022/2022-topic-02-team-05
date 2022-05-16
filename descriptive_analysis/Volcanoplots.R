#-----------------------------------
#für jedes gen p-ewrt berechnen vorher foldchange, log von foldchange xachse, y achse pwert, twosided ttest
#---------------------

tcga_tumor_norm = readRDS("~/GitHub/2022-topic-02-team-05/data/tcga_tumor_normal_datascience_proj_2022.RDS")

#es muss noch die Verteilung der Gene überprüft werden, wir machen einen Shapiro-Wilk test

apply(luad.tumor, 1, function(x){
  ks.test(x, pnorm)
})
y = as.numeric(unique(luad.tumor[1,]))

x = as.numeric(unique(luad.norm[1,]))

ks.test(x, y)

#wenn der eraltene pWert größer gleich 0.05 ist, dann sind die Daten normalverteilt


#cleanen der Daten
gc() #gibt arbeitsspeciher frei der f?r die gro?en datenmengen gebraucht wird
tcga_tumor_norm = na.omit(tcga_tumor_norm)




#a R object containing, for 5 tumor types, the expression data of matched tumor and normal tissue
#Aufdröseln der einzelnen Daten
luad = tcga_tumor_norm[["LUAD"]]
luad.tumor = luad[["tumor"]]
luad.norm = luad[["normal"]]
luad.annot = luad[["clinical"]]




brca = tcga_tumor_norm[["BRCA"]]
brca.tumor = luad[["tumor"]]
brca.norm = luad[["normal"]]
brca.annot = luad[["clinical"]]

kirc = tcga_tumor_norm[["KIRC"]]
kirc.tumor = luad[["tumor"]]
kirc.norm = luad[["normal"]]
kirc.annot = luad[["clinical"]]

prad = tcga_tumor_norm[["PRAD"]]
prad.tumor = luad[["tumor"]]
prad.norm = luad[["normal"]]
prad.annot = luad[["clinical"]]

thca = tcga_tumor_norm[["THCA"]]
thca.tumor = luad[["tumor"]]
thca.norm = luad[["normal"]]
thca.annot = luad[["clinical"]]

#Berrechnen der Varianz
luad.tumor.var = apply(luad.tumor, 1, var)
luad.norm.var = apply(luad.norm, 1, var)

#PLotten als Histogramm
hist(log(luad.tumor.var), breaks = 50, probability = TRUE)
hist(log(luad.norm.var), breaks = 50, probability = TRUE)

#cutting der gene mit sehr niedriger exression d.h. log(var) < -1
#erstmal willk?rlich festgestzt 
#die Werte dürfen ja nur gelöscht werden wenn sie in beiden eine niedrige Varianz haben
if (luad.tumor.var < 0 & luad.norm.var < 0){

luad.tumor = luad.tumor[log(luad.tumor.var) > 0, ]
luad.norm = luad.norm[log(luad.norm.var) > 0, ]
}

#HIER FLIEGEN SEHR VIELE GENE RAUS! IST ES ÜBERHUAPT SINNVOLL VORHER NACH VARIANZ ZU FILTERN??



#berechnen des Mittelwertes für alle Gene
mean.luad.tumor = apply(luad.tumor, 1, mean)
mean.luad.norm = apply(luad.norm, 1, mean)

mean.brca.tumor = apply(brca.tumor, 1, mean)
mean.brca.norm = apply(brca.norm, 1, mean)

mean.kirc.tumor = apply(kirc.tumor, 1, mean)
mean.kirc.norm = apply(kirc.norm, 1, mean)

mean.prad.tumor = apply(prad.tumor, 1, mean)
mean.prad.norm = apply(prad.norm, 1, mean)

mean.thca.tumor = apply(thca.tumor, 1, mean)
mean.thca.norm = apply(thca.norm, 1, mean)

#berechnen des Foldchanges
thca.FC = mean.thca.norm/mean.thca.tumor
luad.FC = mean.luad.norm/mean.luad.tumor
brca.FC = mean.brca.norm/mean.brca.tumor
kirc.FC = mean.kirc.norm/mean.kirc.tumor
prad.FC = mean.prad.norm/mean.prad.tumor

#den log2FC_gene bekommen wir durch mean(condition1) - mean(condition2) 

log2fc.thca = mean.thca.norm - mean.thca.tumor

log2fc.luad = mean.luad.norm - mean.luad.tumor

log2fc.brca = mean.brca.norm - mean.brca.tumor

log2fc.kirc = mean.kirc.norm - mean.kirc.tumor

log2fc.prad = mean.prad.norm - mean.prad.tumor

#weil ich das sonst nicht hinbekomme vertausche ich hier Zeilen und SPalten des df und wandel die erhaltene Matrix wieder in einen df um

prad.tumor.t = t(prad.tumor)
prad.tumor.t = as.data.frame(prad.tumor.t)
prad.norm.t = t(prad.norm)
prad.norm.t = as.data.frame(prad.norm.t)

thca.tumor.t = t(thca.tumor)
thca.tumor.t = as.data.frame(thca.tumor.t)
thca.norm.t = t(thca.norm)
thca.norm.t = as.data.frame(thca.norm.t)

luad.tumor.t = t(luad.tumor)
luad.tumor.t = as.data.frame(luad.tumor.t)
luad.norm.t = t(luad.norm)
luad.norm.t = as.data.frame(luad.norm.t)

#wir machen jetz einen paired, two-tailed t-test um die p-Werte zu bestimmen

p.thca = mapply(function(x,y){t.test(x,y, alternative = "two.sided")$p.value}, thca.tumor,  thca.norm)

p.luad = mapply(function(x,y){t.test(x,y, alternative = "two.sided", paired = TRUE)$p.value}, luad.tumor.t,  luad.norm.t)

p.brca = mapply(function(x,y){t.test(x,y, alternative = "two.sided")$p.value}, thca.tumor,  thca.norm)

p.kirc = mapply(function(x,y){t.test(x,y, alternative = "two.sided")$p.value}, thca.tumor,  thca.norm)

p.prad = mapply(function(x,y){t.test(x,y, alternative = "two.sided", paired = TRUE)$p.value}, prad.tumor.t,  prad.norm.t)

#der p-value muss jetzt noch korrigiert werden, das kann mit der Bonferroni-Methode gemacht werden, weil unser Datensatz iwie so groß ist

pa.thca = p.adjust(p.thca, method = "bonferroni")
pa.luad = p.adjust(p.luad, method = "bonferroni")
pa.brca = p.adjust(p.brca, method = "bonferroni")
pa.kirc = p.adjust(p.kirc, method = "bonferroni")
pa.prad = p.adjust(p.prad, method = "bonferroni")

#daraus können wir dann den volcano plot erstellen, wenn man die low-variance Gene rauswirft sieht der sehr viel schöner aus, dann haben die Vektoren aber unterschiedliche Längen und der Plot bzw. die t-tests funktioniert nicht mehr...
#die Formatierung müssen wir uns evtl mal zusammen anschauen, sodass die nicht signifikanten pWerte grau werden und die signifikanten Gene markiert sind
vol.plot.prad = plot(log2fc.prad, -log10(pa.prad), pch = 20, main = "Volcano plot")

log10.prad = which(pa.prad < 0.05)
names(log10.prad)

points(log2fc.prad, -log10(pa.prad), pch=20, col=names(log10.prad))

ifelse (log2fc.prad > 1,
  points(log2fc.prad, -log10(pa.prad), pch=20, col="orange")
)
ifelse (log2fc.prad < -1){
  points(log2fc.prad, -log10(pa.prad), pch=20, col = "green")
}


vol.plot.prad + 
  geom_hline(yintercept = -log(0.05),
                           linetype = "dashed") +
  geom_vline(xintercept = c(log2(0.5), log2(2)),
                      linetype = "dashed")


vol.plot.luad = plot(log2fc.luad, -log10(pa.luad), pch = 20, main = "Volcano plot")

  