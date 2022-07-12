#---------------------------------------------------
#In diesem Dokument werden wir alle unsere Regressionsmodelle
#mittels F-tests vergleichen um herauszufinden ob das Neuronal network besser
#als die Lineare regression ist
#---------------------------------------------------

#Zu vorhersagender pathway
pathway = 'RODRIGUES_DCC_TARGETS_UP'

#Laden der Vorhersagen
load('data/regression/nn.prediction.RData')
load('data/regression/lm.prediction.RData')
#load('data/regression/lm.pca.prediction.RData')
load('data/regression/null.prediction.RData')

#Laden der MSEs
load('data/regression/nn.MSE.RData')
load('data/regression/lm.MSE.RData')
#load('data/regression/lm.pca.MSE.RData')
load('data/regression/null.MSE.RData')

#Laden der Testdaten
load('data/regression/test.RData')
test = test[,pathway]


#---------------------------------------------------
#1. Graphischer Vergleich
#Wir plotten echte Werte gegen die Vorhergesagten des Modells
#---------------------------------------------------

par(mfrow=c(1,3))

#Plotten der Neuronalen netzes
plot(test, nn.prediction, col='red',
     main= 'Neuronal network',
     pch=18,cex=1,
     ylab = 'prediction', xlab = 'true value', xlim=c(-1,1), ylim=c(-1,1))
abline(0,1,lwd=2)
legend('bottomright',legend= paste('MSE', round(nn.MSE, 3), sep = '='), bty='n')

#Plotten des null models
plot(test, null.prediction, col='green', main='Null model',pch=18, cex=1,
     ylab = 'prediction', xlab = 'true value',xlim=c(-1,1), ylim=c(-1,1))
abline(0,1,lwd=2)
legend('bottomright',legend= paste('MSE', round(null.MSE, 3), sep = '='), bty='n')

#Plotten des Linear models
plot(test, lm.prediction, col='blue',
     main= 'Linear model',
     pch=18,cex=1,
     ylab = 'prediction', xlab = 'true value', xlim=c(-1,1), ylim=c(-1,1))
abline(0,1,lwd=2)
legend('bottomright',legend= paste('MSE', round(lm.MSE, 3), sep = '='), bty='n')

# #Plotten des Linear models mit PCA
# plot(test, lm.pca.prediction, col='green', main='Linear model with PCA',pch=18, cex=1,
#      ylab = 'prediction', xlab = 'true value',xlim=c(-1,1), ylim=c(-1,1))
# abline(0,1,lwd=2)
# legend('bottomright',legend= paste('MSE', round(lm.pca.MSE, 3), sep = '='), bty='n')


#---------------------------------------------------
#2. Vergleich ob sich die Modelle signifikant unterscheiden
#Hierzu führen wir paarweise F-tests durch
#---------------------------------------------------
#F-test testen ob sich die Varianzen signifikant unterscheiden
#in diesem Fall betrechten wir die Varianzen der Residuals => je näher an null desto besser
predictions = matrix(c(nn.prediction, lm.prediction,  null.prediction), ncol = 3,
                     dimnames = list(seq(1,15,1), c('nn','lm','null')))
residuals = apply(predictions, 2, FUN = function(x){x-test})

#Da der F-test sehr anfällig gegen nicht normalverteilte Daten ist machen wir Shapiro-wilk tests
apply(residuals, 2, FUN = function (x){shapiro.test(x)$p.value})
#=> Alle Daten können außer das Nullmodell können als Normalvrteilt angenommen werden

#Durchführen der F-tests jedes gegen jedes Modell
ftest = apply(residuals, 2, FUN = function(x){
  apply(residuals, 2,  FUN = function(y){
    return(var.test(x, y, alternative = 'two.sided')$p.value)
  })
})
#ftest so lesen ist das 

#Visualisierung der Pvalues mit pheatmap
library(pheatmap)
pheatmap(ftest, main = 'Pvalue comparison of regression models',
         breaks = seq(0, 0.5, length.out = 21),
         color = colorRampPalette(c('red','black'),
                                  bias = 1,
                                  space = 'rgb',
                                  interpolate = 'linear'
         )(20),
         clustering_method = 'average', treeheight_row = 0, treeheight_col = 0,
         cellwidth = 40, cellheight = 40,
         show_colnames = TRUE,show_rownames = TRUE, border_color = NA,
         legend_breaks = c(0,0.05, 0.5),
         legend = FALSE,
         display_numbers = TRUE, number_color = 'white',
         legend_labels = c('significant','alpha', 'non significant')
)

