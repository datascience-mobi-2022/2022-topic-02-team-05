#--------------------------------------------------------
#In diesem Dokument versuchen wir die Pathwayactivität bestimmter Pathways
#basierend auf allen anderen pathways in THCA Tumorgewebe vorherzusagen
#--------------------------------------------------------

#--------------------------------------------------------
#1. Präparation der Daten
#--------------------------------------------------------
load('data/thca_gsea.RData')

#Splitten der Daten in eine Test und Trainingsgruppe (75% für Training)
set.seed(1) #Für konsistente Indices da sample() random ist
index = sample(1:ncol(thca_gsea), round(0.75*ncol(thca_gsea)))
train = as.data.frame(t(thca_gsea[, index]))
test = as.data.frame(t(thca_gsea[, -index]))

#--------------------------------------------------------
#2. Linear Model der Trainingsdaten
#--------------------------------------------------------
#Definieren der Formel für die vorhersage
pathway = 'REACTOME_THYROXINE_BIOSYNTHESIS' #Zu vorhersagender pathway
form = paste(c(pathway,
               paste(rownames(thca_gsea)[!pathway == rownames(thca_gsea)], collapse = " + ")),
             collapse = ' ~ '
)
#Durchführen der Regression
linear.model = glm(formula = form, family = 'gaussian', data = train)

#Test der Regression
lm.prediction = predict(linear.model, test)
lm.error = sum((lm.prediction - test[,pathway])^2)/nrow(test) #berrechnung des meansquared error als Guete für das Modell

#--------------------------------------------------------
#3. Neuronales Netz
#--------------------------------------------------------
#Skalieren der Daten nach mit einer Min/max skalierung sodass sie im Intervall [0,1] liegen
maxs = apply(thca_gsea, 1, max) 
mins = apply(thca_gsea, 1, min)
thca_scaled = as.data.frame(scale(t(thca_gsea), center = mins, scale = maxs - mins))
train_sc = thca_scaled[index,]
test_sc = thca_scaled[-index,]

#Implementierung und Training des Netzes
library(neuralnet)
set.seed(1)
AI = neuralnet(formula = form, #Die Gleichung die wir  vorhersagen möchten
               data = train_sc,
               hidden = c(20, 10),
               linear.output = TRUE)
#Plotten des Netzwerks (grade noch sehr unübersichtlich)
#plot(AI, rep = 'best', radius = 0.07, fontsize = 0)

#--------------------------------------------------------
#4. Testen des Neuronalen Netzes
#--------------------------------------------------------
#predicted jetzt unseren pathway basierend auf denen der Testdaten
nn.prediction_sc = compute(AI,test_sc[,!colnames(test_sc) == pathway])

#predictete Werte sind noch im Intervall[0,1] => hochscalieren mit min und max
nn.prediction = as.numeric(
  nn.prediction_sc$net.result*(max(thca_gsea[pathway,])-min(thca_gsea[pathway,])) + min(thca_gsea[pathway,]))
                           
nn.test = test_sc[,pathway]*(max(thca_gsea[pathway,])-min(thca_gsea[pathway,]))+min(thca_gsea[pathway,])

#mean sum squared error des neuralen netzes
nn.error = sum((nn.test - nn.prediction)^2)/nrow(test_sc)

#=> Der Fehler des Netzwerks ist mit 0.13 deutlich kleiner als die Lin. Regression (4.42)

#---------------------------------------------------------
#Vergleich beider Modelle
#---------------------------------------------------------
par(mfrow=c(1,2))

#Plotten der Neuronalen netzes
plot(test[,pathway], nn.prediction, col='red', main='Neuronal network', pch=18,cex=0.7,
     ylab = 'prediction', xlab = 'true value', xlim=c(-1,0), ylim=c(-5,2))
abline(0,1,lwd=2)
legend('bottomright',legend= paste('MSE', round(nn.error, 2), sep = '='), bty='n')

#Plotten des Linear models
plot(test[,pathway], lm.prediction, col='blue', main='Linear model',pch=18, cex=0.7,
     ylab = 'prediction', xlab = 'true value',xlim=c(-1,0), ylim=c(-5,2))
abline(0,1,lwd=2)
legend('bottomright',legend= paste('MSE', round(lm.error, 2), sep = '='), bty='n')

