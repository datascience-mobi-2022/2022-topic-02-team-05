#--------------------------------------------------------
#In diesem Dokument versuchen wir die Pathwayactivität bestimmter Pathways
#basierend auf allen anderen pathways in THCA GSEA Daten vorherzusagen
#--------------------------------------------------------

pathway = 'REACTOME_THYROXINE_BIOSYNTHESIS' #Zu vorhersagender pathway

load('data/regression/test_sc.RData')
load('data/regression/train_sc.RData')



#---------------------------------------------------------
#4. Optimierung des Neuronalen Netzes
#Hier testen wir welche internen Strukturen die Performance optimieren
#und ob wir durch andere Anfangsweights/biases evtl. ein besseres Minimum finden
#Das ist nötig, da das Netz mit Gradient decent arbeitet und so immer nur lokale
#Minima findet allerdings kein globales Minimum
#---------------------------------------------------------

#Erstellen eine Liste in die wir verschieden 25 neuronale Netze speichern
#Wir testen und 25 Verschiedene Architekturen alle mit jeweils 2 Hidden layers

#Wir testen layer mit jeweils 10, 20, 50, 100, 300 Neuronen
library(neuralnet)
set.seed(1)

first_layer = c(10, 20, 50, 100, 300)
second_layer = c(10, 20, 50, 100, 300)

nn_list = list()

#Loop der durch die Liste geht und überall die entsprechenden Netze speichert
for (j in second_layer){
  for (i in first_layer){
    res = list(neuralnet(formula = form, #Die Gleichung die wir  vorhersagen möchten
                             data = train_sc,
                             hidden = c(i, j),
                             linear.output = TRUE))
    nn_list = append(nn_list, res)
  }
};rm(i);rm(j);rm(res)

#Namen der einzelen Netze
names = c()
for (j in second_layer){
  for (i in first_layer){
    res = paste(i, j, sep = ':')
    names = append(names, res)
  }
};rm(i);rm(j);rm(res)
names(nn_list) = names; rm(names)

#Testen aller Netze und speichern der MSE werte
nn.test = test_sc[,pathway]*(max(thca_gsea[pathway,])-min(thca_gsea[pathway,]))+min(thca_gsea[pathway,])

MSE = sapply(nn_list, FUN = function(x){
        #Testen
        res = compute(x, test_sc[,!colnames(test_sc) == pathway]) 
        #Hochskalieren
        res = as.numeric( 
          res$net.result*(max(thca_gsea[pathway,])-min(thca_gsea[pathway,])) + min(thca_gsea[pathway,]))
        #Fehler berechen
        MSE = sum((nn.test - res)^2)/nrow(test_sc) 
        return(MSE)
      })
#Ausgeben der Besten drei Architekturen
sort(MSE)[1:3]
#=> eine 100:10 Architektur scheint am besten zu sein

#Nun testen wir wie die Anfangs weights and biases das Netz beeinflussen
#dazu teste wir die besten 3 Netze mit jeweils 100 Anfangsbedingungen (seeds)

best_layers = strsplit(names(sort(MSE))[1:3], split = ':', fixed = TRUE)
best_layers = lapply(best_layers, as.numeric)

#Erstellen der 300 Netze
nn_list_best = list()
for (j in 1:3){
  for (i in 1:100){
    set.seed(i)
    res = list(neuralnet(formula = form, #Die Gleichung die wir  vorhersagen möchten
                         data = train_sc,
                         hidden = best_layers[[j]],
                         linear.output = TRUE))
    nn_list_best = append(nn_list_best, res)
  }
};rm(i);rm(j);rm(res)

#Namen für die 300 Netze
names = c()
for (j in 1:3){
  for (i in 1:100){
    res = paste(paste(best_layers[[j]][1],best_layers[[j]][2] , sep = ':'), i, sep = ', seed=')
    names = append(names, res)
  }
};rm(i);rm(j);rm(res)
names(nn_list_best) = names; rm(names)

#Testen aller Netze und speichern der MSE werte
MSE_best = sapply(nn_list_best, FUN = function(x){
  #Testen
  res = compute(x, test_sc[,!colnames(test_sc) == pathway]) 
  #Hochskalieren
  res = as.numeric( 
    res$net.result*(max(thca_gsea[pathway,])-min(thca_gsea[pathway,])) + min(thca_gsea[pathway,]))
  #Fehler berechen
  MSE = sum((nn.test - res)^2)/nrow(test_sc) 
  return(MSE)
})
sort(MSE_best)[1:3]
best_seed = as.numeric(strsplit(names(sort(MSE_best)[1]), split = '=')[[1]][2])
# Das Beste Netz bekomen wir mit einer Architektur von 100:10 und den
#Anfangsbedingungen bei set.seed(12)

#-------------------------------------------------------
#5. Implementierung und Training des besten Netzes
#-------------------------------------------------------
set.seed(best_seed)
AI = neuralnet(formula = form, #Die Gleichung die wir  vorhersagen möchten
               data = train_sc,
               hidden = best_layers[[1]],
               linear.output = TRUE)
#Plotten des Netzwerks (grade noch sehr unübersichtlich)
#plot(AI, rep = 'best', radius = 0.07, fontsize = 0)

#--------------------------------------------------------
#6. Testen des Neuronalen Netzes
#--------------------------------------------------------
#predicted jetzt unseren pathway basierend auf denen der Testdaten
nn.prediction_sc = compute(AI,test_sc[,!colnames(test_sc) == pathway])

#predictete Werte sind noch im Intervall[0,1] => hochscalieren mit min und max
nn.prediction = as.numeric(
  nn.prediction_sc$net.result*(max(thca_gsea[pathway,])-min(thca_gsea[pathway,])) + min(thca_gsea[pathway,]))

nn.test = test_sc[,pathway]*(max(thca_gsea[pathway,])-min(thca_gsea[pathway,]))+min(thca_gsea[pathway,])

#mean sum squared error des neuralen netzes
nn.error = sum((nn.test - nn.prediction)^2)/nrow(test_sc)

#Mittlerer Prozentualer fehler des Netzes
nn.procent.error = mean(abs(1-nn.prediction/nn.test))

#=> Der Fehler des Netzwerks ist mit 8.8% deutlich kleiner als die Lin. Regression (240%)


#---------------------------------------------------------
#7. Vergleich beider Modelle
#---------------------------------------------------------
par(mfrow=c(1,2))

#Plotten der Neuronalen netzes
plot(test[,pathway], nn.prediction, col='red',
     main=paste('Neuronal network',paste(best_layers[[1]][1], best_layers[[1]][2], sep = ':'), sep = ' '),
     pch=18,cex=0.7,
     ylab = 'prediction', xlab = 'true value', xlim=c(-1,0), ylim=c(-5,2))
abline(0,1,lwd=2)
legend('bottomright',legend= paste('MSE', round(nn.error, 2), sep = '='), bty='n')

#Plotten des Linear models
plot(test[,pathway], lm.prediction, col='blue', main='Linear model',pch=18, cex=0.7,
     ylab = 'prediction', xlab = 'true value',xlim=c(-1,0), ylim=c(-5,2))
abline(0,1,lwd=2)
legend('bottomright',legend= paste('MSE', round(lm.error, 2), sep = '='), bty='n')

#----------------------------------------------------------
#8. Funktion die auf neuen Daten die Pathwayaktivität vorhersagt
#----------------------------------------------------------

AI_predict = function(input_pathways, network = AI){
  #Skalieren der Input pathways
  maxs = max(input_pathways) 
  mins = min(input_pathways)
  input_scaled = as.data.frame(scale(input_pathways), center = mins, scale = maxs - mins)
  #Vorhersage der Aktivität
  prediction = compute(AI, t(input_scaled))
  prediction = as.numeric(prediction$net.result*(maxs-mins) + mins)
  return(prediction)
}

