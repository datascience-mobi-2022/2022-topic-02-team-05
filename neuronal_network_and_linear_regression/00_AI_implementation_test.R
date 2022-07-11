#--------------------------------------------------------
#In diesem Dokument versuchen wir die Grundlagen eines deep leraning netzes
#nachzuvollziehen dazu benutzen wir einen Trainingsdaten satz
#(Boston dataset aus dem MASS) aus dem MASS package
# Quelle: https://www.r-bloggers.com/2015/09/fitting-a-neural-network-in-r-neuralnet-package/
#--------------------------------------------------------

#Laden der Trainingsdaten
set.seed(123)
library(MASS)
data <- Boston

#---------------------
#Vorarbeit 
#----------------------

#Checking for NAs
apply(data,2,function(x) sum(is.na(x)))
#=> keine NAs

#zufälliges Teilen der Daten in Trainings und Testdaten
index <- sample(1:nrow(data),round(0.75*nrow(data)))
train <- data[index,] #75% der Daten
test <- data[-index,] #25% der Daten

#linear model für die Varibale medv um nacher das netz zu vergleichen basierend auf allen Var
lm.fit <- glm(medv ~ . , data = train)
summary(lm.fit)

#Vorhersage der testdaten mit lin. modell
pr.lm <- predict(lm.fit, test)
MSE.lm <- sum((pr.lm - test$medv)^2)/nrow(test) #berrechnung des meansquared error als Guete für das Modell

#---------------------------
#Datenpräpartion für Deep leraning
#---------------------------
#Skalieren der Daten sodass die im [0,1] Intervall liegen (hier mit min, max scaling)
maxs <- apply(data, 2, max) 
mins <- apply(data, 2, min)
scaled <- as.data.frame(scale(data, center = mins, scale = maxs - mins))
#scaled trainig and test data
train_ <- scaled[index,]
test_ <- scaled[-index,]


#Anmerkungen für das Netz: Rule of thumb hidden layers haben 2/3 der neurons der vorherigen
#daher für uns mit 2 hidden layers einer  13:5:3:1 structur (ein output neuron)

library(neuralnet)
n <- names(train_)
f <- as.formula(paste("medv ~", paste(n[!n %in% "medv"], collapse = " + ")))
nn <- neuralnet(f, #Regressionsformel die wr vorhersagen wollen
                data=train_, #Trainingsdata
                hidden=c(5,3), #Anzahl und art der hidden layers
                linear.output=T) #TRUE für regrssionsmodelle FALSE für kategorische
plot(nn) #black sind weights und blue biases des Modells
#=> das neu.net ist jetzt traine

#--------------------------------
#Testen des Netzes
#--------------------------------
pr.nn <- compute(nn,test_[,1:13]) #predicted jetzt mdev basierend auf den ersten13 Var der testdaten
#predictete Werte sind noch im Intervall[0,1] => hochscalieren mit min und max
pr.nn_ <- pr.nn$net.result*(max(data$medv)-min(data$medv))+min(data$medv)
test.r <- (test_$medv)*(max(data$medv)-min(data$medv))+min(data$medv)
#mean sum squared error des neuralen netzes
MSE.nn <- sum((test.r - pr.nn_)^2)/nrow(test_)

print(paste(MSE.lm,MSE.nn)) #Neuronales netz ist also besser

#Visualisierung zwischen beiden Modellen
par(mfrow=c(1,2))
plot(test$medv,pr.nn_,col='red',main='Real vs predicted NN',pch=18,cex=0.7)
abline(0,1,lwd=2)
legend('bottomright',legend='NN',pch=18,col='red', bty='n')
plot(test$medv,pr.lm,col='blue',main='Real vs predicted lm',pch=18, cex=0.7)
abline(0,1,lwd=2)
legend('bottomright',legend='LM',pch=18,col='blue', bty='n', cex=.95)
#Auch hier sieht man dass das linear modell besser performed



