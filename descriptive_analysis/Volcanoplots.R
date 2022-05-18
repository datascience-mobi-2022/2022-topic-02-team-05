#-----------------------------------
#f?r jedes gen p-ewrt berechnen vorher foldchange, log von foldchange xachse, y achse pwert, twosided ttest
#---------------------

tcga_tumor_norm = readRDS("~/GitHub/2022-topic-02-team-05/data/tcga_tumor_normal_datascience_proj_2022.RDS")
thca.tumor <- tcga_tumor_norm$THCA$tumor
thca.norm <- tcga_tumor_norm$THCA$normal

#es muss noch die Verteilung der Gene ?berpr?ft werden, wir machen einen Shapiro-Wilk test -> Test auf normally distributed

apply(thca.tumor, 1, function(x){
  ks.test(x, pnorm)
})
y = as.numeric(unique(thca.tumor[1,]))

x = as.numeric(unique(thca.norm[1,]))

ks.test(x, y) # Kolmogorow-Smirnow-Test: Test auf Übereinstimmung zweier WK-Verteilungen

qqplot(thca.tumor,thca.norm)
qqplot(x,y)

#wenn der erhaltene pWert >= 0.05 ist, dann werden die Daten als normalverteilt angenommen

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
  
  

#a R object containing, for 5 tumor types, the expression data of matched tumor and normal tissue
#Aufdroeseln der einzelnen Daten
#luad = tcga_tumor_norm[["LUAD"]]
#luad.tumor = luad[["tumor"]]
#luad.norm = luad[["normal"]]
#luad.annot = luad[["clinical"]]

#brca = tcga_tumor_norm[["BRCA"]]
#brca.tumor = luad[["tumor"]]
#brca.norm = luad[["normal"]]
#brca.annot = luad[["clinical"]]

#kirc = tcga_tumor_norm[["KIRC"]]
#kirc.tumor = luad[["tumor"]]
#kirc.norm = luad[["normal"]]
#kirc.annot = luad[["clinical"]]

#prad = tcga_tumor_norm[["PRAD"]]
#prad.tumor = luad[["tumor"]]
#prad.norm = luad[["normal"]]
#prad.annot = luad[["clinical"]]

#thca = tcga_tumor_norm[["THCA"]]
#thca.tumor = luad[["tumor"]]
#thca.norm = luad[["normal"]]
#thca.annot = luad[["clinical"]]


#berechnen des Mittelwertes f?r alle Gene
mean.luad.tumor = apply(luad.tumor, 1, mean)
mean.luad.norm = apply(luad.norm, 1, mean)

mean.brca.tumor = apply(brca.tumor, 1, mean)
mean.brca.norm = apply(brca.norm, 1, mean)

mean.kirc.tumor = apply(kirc.tumor, 1, mean)
mean.kirc.norm = apply(kirc.norm, 1, mean)

mean.prad.tumor = apply(prad.tumor, 1, mean)
mean.prad.norm = apply(prad.norm, 1, mean)

#Mean der gene von allen thcapatienten
mean.thca.tumor = apply(thca.tumor.va, 1, mean)
mean.thca.norm = apply(thca.norm.va, 1, mean)

#-----------------------------------------------------------------------------------
#Foldchange(FC) besrechnen!
#den log2FC_gene bekommen wir durch mean(condition1) - mean(condition2) 
#---------------------------------------------------------------------------------

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

brca.tumor.t = t(brca.tumor)
brca.tumor.t = as.data.frame(brca.tumor.t)
brca.norm.t = t(brca.norm)
brca.norm.t = as.data.frame(brca.norm.t)

kirc.tumor.t = t(kirc.tumor)
kirc.tumor.t = as.data.frame(kirc.tumor.t)
kirc.norm.t = t(kirc.norm)
kirc.norm.t = as.data.frame(kirc.norm.t)

thca.tumor.t = t(thca.tumor.va)
thca.tumor.t = as.data.frame(thca.tumor.t)
thca.norm.t = t(thca.norm.va)
thca.norm.t = as.data.frame(thca.norm.t)

luad.tumor.t = t(luad.tumor)
luad.tumor.t = as.data.frame(luad.tumor.t)
luad.norm.t = t(luad.norm.va)
luad.norm.t = as.data.frame(luad.norm.t)

#wir machen jetz einen paired, two-tailed t-test um die p-Werte zu bestimmen

p.thca = mapply(function(x,y){t.test(x,y, alternative = "two.sided", paired = TRUE)$p.value}, thca.tumor.t,  thca.norm.t)

p.luad = mapply(function(x,y){t.test(x,y, alternative = "two.sided", paired = TRUE)$p.value}, luad.tumor.t,  luad.norm.t)

p.brca = mapply(function(x,y){t.test(x,y, alternative = "two.sided", paired = TRUE)$p.value}, thca.tumor.t,  thca.norm.t)

p.kirc = mapply(function(x,y){t.test(x,y, alternative = "two.sided", paired = TRUE)$p.value}, thca.tumor.t,  thca.norm.t)

p.prad = mapply(function(x,y){t.test(x,y, alternative = "two.sided", paired = TRUE)$p.value}, prad.tumor.t,  prad.norm.t)

#der p-value muss jetzt noch korrigiert werden, das kann mit der Bonferroni-Methode gemacht werden, weil unser Datensatz iwie so gro? ist

pa.thca = p.adjust(p.thca, method = "bonferroni")
pa.luad = p.adjust(p.luad, method = "bonferroni")
pa.brca = p.adjust(p.brca, method = "bonferroni")
pa.kirc = p.adjust(p.kirc, method = "bonferroni")
pa.prad = p.adjust(p.prad, method = "bonferroni")

#daraus k?nnen wir dann den volcano plot erstellen, wenn man die low-variance Gene rauswirft sieht der sehr viel sch?ner aus, dann haben die Vektoren aber unterschiedliche L?ngen und der Plot bzw. die t-tests funktioniert nicht mehr...
#die Formatierung m?ssen wir uns evtl mal zusammen anschauen, sodass die nicht signifikanten pWerte grau werden und die signifikanten Gene markiert sind


data.prad = data.frame(log2fc.prad, pa.prad) #die 2 Vektoren f?r unseren Volcano PLot werden in einen df gepackt, damit daraus ein plot erstellt werden kann

with(data.prad, plot(data.prad$log2fc.prad, -log10(data.prad$pa.prad), pch = 20, main = "Volcano plot"))
with(subset(data.prad, pa.prad < 0.05 ), points(log2fc.prad, -log10(pa.prad), pch=20, col="red"))
with(subset(data.prad, log2fc.prad>1), points(log2fc.prad, -log10(pa.prad), pch=20, col="orange"))
with(subset(data.prad, log2fc.prad< -1), points(log2fc.prad, -log10(pa.prad), pch=20, col="green"))

data.thca = data.frame(log2fc.thca, pa.thca) #die 2 Vektoren f?r unseren Volcano PLot werden in einen df gepackt, damit daraus ein plot erstellt werden kann
with(data.thca, plot(data.thca$log2fc.thca, -log10(data.thca$pa.thca), pch = 20, main = "Volcano plot"))
with(subset(data.thca, pa.thca < 0.05 ), points(log2fc.thca, -log10(pa.thca), pch=20, col="red"))
with(subset(data.thca, log2fc.thca>1), points(log2fc.thca, -log10(pa.thca), pch=20, col="orange"))
with(subset(data.thca, log2fc.thca< -1), points(log2fc.thca, -log10(pa.thca), pch=20, col="green"))

data.luad = data.frame(log2fc.luad, pa.luad) #die 2 Vektoren f?r unseren Volcano PLot werden in einen df gepackt, damit daraus ein plot erstellt werden kann
with(data.luad, plot(data.luad$log2fc.luad, -log10(data.luad$pa.luad), pch = 20, main = "Volcano plot"))
with(subset(data.luad, pa.luad < 0.05 ), points(log2fc.luad, -log10(pa.luad), pch=20, col="red"))
with(subset(data.luad, log2fc.luad>1), points(log2fc.luad, -log10(pa.luad), pch=20, col="orange"))
with(subset(data.luad, log2fc.luad< -1), points(log2fc.luad, -log10(pa.luad), pch=20, col="green"))

data.brca = data.frame(log2fc.brca, pa.brca) #die 2 Vektoren f?r unseren Volcano PLot werden in einen df gepackt, damit daraus ein plot erstellt werden kann
with(data.brca, plot(data.brca$log2fc.brca, -log10(data.brca$pa.brca), pch = 20, main = "Volcano plot"))
with(subset(data.brca, pa.brca < 0.05 ), points(log2fc.brca, -log10(pa.brca), pch=20, col="red"))
with(subset(data.brca, log2fc.brca>1), points(log2fc.brca, -log10(pa.brca), pch=20, col="orange"))
with(subset(data.brca, log2fc.brca< -1), points(log2fc.brca, -log10(pa.brca), pch=20, col="green"))

data.kirc = data.frame(log2fc.kirc, pa.kirc) #die 2 Vektoren f?r unseren Volcano PLot werden in einen df gepackt, damit daraus ein plot erstellt werden kann
with(data.kirc, plot(data.kirc$log2fc.kirc, -log10(data.kirc$pa.kirc), pch = 20, main = "Volcano plot"))
with(subset(data.kirc, pa.kirc < 0.05 ), points(log2fc.kirc, -log10(pa.kirc), pch=20, col="red"))
with(subset(data.kirc, log2fc.kirc>1), points(log2fc.kirc, -log10(pa.kirc), pch=20, col="orange"))
with(subset(data.kirc, log2fc.kirc< -1), points(log2fc.kirc, -log10(pa.kirc), pch=20, col="green"))


#hiermit werden die besonders signifikanten Gennamen neben die jeweiligen Punkte geschrieben, die SIgnifikanzniveaus m?ssen evtl nochmal ?berarbeitet werden
install.packages('calibrate')
library(calibrate)


with(subset(data.prad, pa.prad<.005 & abs(log2fc.prad)>5), textxy(log2fc.prad, -log10(pa.prad), labs=rownames(data), cex=.8))



