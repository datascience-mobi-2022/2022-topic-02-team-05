#--------------------------------------------
#In diesem Dokument werden wir eine lineare Regression in unseren THCA GSEA
#Daten durchführen um so eine Pathwayaktivität vorherzusagen
#--------------------------------------------

#Zu vorhersagender pathway
pathway = 'REACTOME_THYROXINE_BIOSYNTHESIS' 

load('data/regression/test_pca.RData')
load('data/regression/train_pca.RData')
load('data/regression/test.RData')
load('data/regression/train.RData')

#--------------------------------------------------------
#1. Correlation der Pathways überprüfen
#--------------------------------------------------------
#Correlation aller Pathways berrechnen
cor = cor(train[,!colnames(train) == pathway], method = 'spearman')

#Visualisierung der Heatmap mit pheatmap
library(pheatmap)
pheatmap(as.matrix(cor),
         breaks = seq(-1, 1, length.out = 201),
         color = colorRampPalette(c('blue4','white','red'),
                                  bias = 1,
                                  space = 'rgb',
                                  interpolate = 'linear'
         )(200),
         clustering_method = 'average', treeheight_row = 0, treeheight_col = 0,
         cellwidth = 0.7, cellheight = 0.7,
         show_colnames = FALSE,show_rownames = FALSE, border_color = 'lightcyan2',
         legend_breaks = c(-1, 0, 1),
         legend_labels = c('neg. correlation','no correlation','pos. correlation')
)

#Wir führen nun PCA durch um die Correlation zwischen den Pathways aufzuheben und
#so hofentlich die Regression zu verbessern

cor_pca = cor(train_pca[,!colnames(train_pca) == pathway], method = 'spearman')
pheatmap(as.matrix(cor_pca),
         breaks = seq(-1, 1, length.out = 201),
         color = colorRampPalette(c('blue4','white','red'),
                                  bias = 1,
                                  space = 'rgb',
                                  interpolate = 'linear'
         )(200),
         clustering_method = 'average', treeheight_row = 0, treeheight_col = 0,
         cellwidth = 2, cellheight = 2,
         show_colnames = FALSE,show_rownames = FALSE, border_color = 'lightcyan2',
         legend_breaks = c(-1, 0, 1),
         legend_labels = c('neg. correlation','no correlation','pos. correlation')
)

par(mfrow=c(1,2))
hist(cor)
hist(cor_pca)

#--------------------------------------------------------
#2. Lineare Regression ohne PCA
#--------------------------------------------------------
#Definieren der Formel für die vorhersage

form = paste(c(pathway,
               paste(colnames(train)[!pathway == colnames(train)], collapse = " + ")),
             collapse = ' ~ '
)

#Durchführen der Regression
linear.model = glm(formula = form, family = 'gaussian', data = train)

#Test der Regression
lm.prediction = predict(linear.model, test)
cor(train[, pathway], linear.model$y) #=> Cor zwischen Testwerten und Prediction ist wie erwartet 1

#berrechnung des meansquared error als Guete für das Modell
lm.MSE = sum((lm.prediction - test[,pathway])^2)/nrow(test) 

#Anschauen der Residuals
cor(train[, pathway], linear.model$residuals) #Korrelation ist sehr niedrig => gutes Modell
qqplot(qnorm(seq(0,1,0.01)), quantile(linear.model$residuals, seq(0,1,0.01))) #=> Residuals sind annähernd normalverteilt 

#--------------------------------------------------------
#3. Lineare Regression mit PCA
#--------------------------------------------------------
#Definieren der Formel für die vorhersage

form.pca = paste(c(pathway,
             paste(colnames(train_pca)[!pathway == colnames(train_pca)], collapse = " + ")),
             collapse = ' ~ '
)

#Durchführen der Regression
linear.model.pca = glm(formula = form.pca, family = 'gaussian', data = train_pca)

#Test der Regression
lm.pca.prediction = predict(linear.model.pca, test_pca)
cor(train_pca[, pathway], linear.model.pca$y) #=> Cor zwischen Testwerten und Prediction ist wie erwartet 1

#berrechnung des meansquared error als Guete für das Modell
lm.pca.MSE = sum((lm.pca.prediction - test_pca[,pathway])^2)/nrow(test_pca) 

#Anschauen der Residuals
cor(train_pca[, pathway], linear.model.pca$residuals) #Korrelation ist sehr niedrig => gutes Modell
qqplot(qnorm(seq(0,1,0.01)), quantile(linear.model.pca$residuals, seq(0,1,0.01))) #=> Residuals sind annähernd normalverteilt 

#--------------------------------------------------------
#4. Testen eines Nullmodells zum Vergleich
#--------------------------------------------------------
#Wir nehmen als Vorhersage den Mittleren Aktivitätswert des Pathways
null.prediction = mean(train[, pathway])
#Nun berrechnen wir auch hier den MSE
null.MSE = sum((rep(null.prediction, nrow(test)) - test[,pathway])^2)/nrow(test)

par(mfrow=c(1,3))

#Plotten der Neuronalen netzes
plot(test[,pathway], lm.prediction, col='red',
     main= 'LM',
     pch=18,cex=0.7,
     ylab = 'prediction', xlab = 'true value', xlim=c(-1,0), ylim=c(-5,2))
abline(0,1,lwd=2)
#legend('bottomright',legend= paste('MSE', round(nn.error, 2), sep = '='), bty='n')

#Plotten des Linear models
plot(test[,pathway], lm.pca.prediction, col='blue', main='Linear model PCA',pch=18, cex=0.7,
     ylab = 'prediction', xlab = 'true value',xlim=c(-1,0), ylim=c(-5,2))
abline(0,1,lwd=2)
#legend('bottomright',legend= paste('MSE', round(lm.error, 2), sep = '='), bty='n')

#Plotten des Linear models
plot(test[,pathway], rep(null.prediction, nrow(test)), col='green', main='null model',pch=18, cex=0.7,
     ylab = 'prediction', xlab = 'true value',xlim=c(-1,0), ylim=c(-5,2))
abline(0,1,lwd=2)
#legend('bottomright',legend= paste('MSE', round(lm.error, 2), sep = '='), bty='n')

