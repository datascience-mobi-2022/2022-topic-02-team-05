#--------------------------------------------
#In diesem Dokument werden wir eine lineare Regression in unseren THCA GSEA
#Daten durchführen um so eine Pathwayaktivität vorherzusagen
#--------------------------------------------

#Zu vorhersagender pathway
pathway = 'RODRIGUES_DCC_TARGETS_UP' 

# load('data/regression/test_pca.RData')
# load('data/regression/train_pca.RData')
load('data/regression/test.RData')
load('data/regression/train.RData')


#--------------------------------------------------------
#1. Correlation der Pathways überprüfen
#--------------------------------------------------------
path.values = train[, colnames(train) == pathway]

#Varianz cleaning der Daten (nur die 20% der Variantesten bleiben'
variance = apply(train[, !colnames(train) == pathway], 2, var)
train.topvar = train[, variance > quantile(variance, probs = 0.9)]; rm(variance)
cor.topvar = cor(train.topvar)
diag(cor.topvar) = NA

library(pheatmap)
pheatmap(as.matrix(cor.topvar),
         breaks = seq(-1, 1, length.out = 201),
         color = colorRampPalette(c('blue4','white','red'),
                                  bias = 1,
                                  space = 'rgb',
                                  interpolate = 'linear'
         )(200),
         clustering_method = 'average', treeheight_row = 0, treeheight_col = 0,
         show_colnames = FALSE,show_rownames = FALSE, border_color = 'lightcyan2',
         legend_breaks = c(-1, 0, 1),
         legend_labels = c('neg. correlation','no correlation','pos. correlation')
)

#Entfernen aller stark correlierten pathways
cor.topvar[upper.tri(cor.topvar)] = 0
diag(cor.topvar) = 0
cor_cleaned_train = train.topvar[, !apply(cor.topvar, 2, function(x) any(abs(x) > 0.5, na.rm = TRUE))]

rm(train.topvar, cor.topvar)

cor.cleaned = cor(cor_cleaned_train)
diag(cor.cleaned) = NA
pheatmap(as.matrix(cor.cleaned),
         breaks = seq(-1, 1, length.out = 201),
         color = colorRampPalette(c('blue4','white','red'),
                                  bias = 1,
                                  space = 'rgb',
                                  interpolate = 'linear'
         )(200),
         clustering_method = 'average', treeheight_row = 0, treeheight_col = 0,
         show_colnames = FALSE,show_rownames = FALSE, border_color = 'lightcyan2',
         legend_breaks = c(-1, 0, 1),
         legend_labels = c('neg. correlation','no correlation','pos. correlation')
)

rm(cor.cleaned)

#Hinzufügen unseres pathways zu den un korrelierten Daten
train_cleaned = cbind(path.values, cor_cleaned_train)
colnames(train_cleaned)[1] = pathway; rm(cor_cleaned_train, path.values)


#--------------------------------------------------------
#2. Lineare Regression ohne PCA
#--------------------------------------------------------

#Definieren der Formel für die vorhersage
form = formula(paste(c(pathway,
                       paste(colnames(train_cleaned)[!pathway == colnames(train_cleaned)], collapse = " + ")),
                     collapse = ' ~ '
))

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

#Speichern für Vergleich
save(lm.prediction, file = 'data/regression/lm.prediction.RData')
save(lm.MSE, file = 'data/regression/lm.MSE.RData')

#--------------------------------------------------------
#3. Lineare Regression mit nur signifikanten Pathways
#--------------------------------------------------------

#Verbessern des Modells durch entfernen der nicht signifikanten Pathways
lm.summary = summary(linear.model)$coefficients

#Behalten der signifikantesten 20% der Pathways die zum Modell beitragen
regression.path = names(which(lm.summary[,4] < quantile(lm.summary[,4], probs = 0.2)))
form.new = formula(paste(c(pathway,
                           paste(regression.path, collapse = " + ")),
                         collapse = ' ~ '
))

linear.model.new = lm(formula = form.new, data = train)

#Test der Regression
lm.new.prediction = predict(linear.model.new, test)
cor(train[, pathway], linear.model.new$fitted.values)
#=> Cor zwischen Testwerten und Prediction ist nicht annähernd 1

#berrechnung des meansquared error als Guete für das Modell
lm.new.MSE = sum((lm.new.prediction - test[,pathway])^2)/nrow(test) 

#Anschauen der Residuals
cor(train[, pathway], linear.model.new$residuals)
#Korrelation ist nicht besonders niedrig => schlechtes Modell
qqplot(qnorm(seq(0,1,0.01)), quantile(linear.model.new$residuals, seq(0,1,0.01)))
#=> Residuals sind annähernd normalverteilt 

#Speichern für Vergleich
save(lm.new.prediction, file = 'data/regression/lm.new.prediction.RData')
save(lm.new.MSE, file = 'data/regression/lm.new.MSE.RData')


# #--------------------------------------------------------
# #3. Lineare Regression mit PCA
# #--------------------------------------------------------
# #Definieren der Formel für die vorhersage
# 
# form.pca = paste(c(pathway,
#              paste(colnames(train_pca)[!pathway == colnames(train_pca)], collapse = " + ")),
#              collapse = ' ~ '
# )
# 
# #Durchführen der Regression
# linear.model.pca = glm(formula = form.pca, family = 'gaussian', data = train_pca)
# 
# #Test der Regression
# lm.pca.prediction = predict(linear.model.pca, test_pca)
# cor(train_pca[, pathway], linear.model.pca$y) #=> Cor zwischen Testwerten und Prediction ist wie erwartet 1
# 
# #berrechnung des meansquared error als Guete für das Modell
# lm.pca.MSE = sum((lm.pca.prediction - test_pca[,pathway])^2)/nrow(test_pca) 
# 
# #Anschauen der Residuals
# cor(train_pca[, pathway], linear.model.pca$residuals) #Korrelation ist sehr niedrig => gutes Modell
# qqplot(qnorm(seq(0,1,0.01)), quantile(linear.model.pca$residuals, seq(0,1,0.01))) #=> Residuals sind annähernd normalverteilt 
# 
# #Speichern für Vergleich
# save(lm.pca.prediction, file = 'data/regression/lm.pca.prediction.RData')
# save(lm.pca.MSE, file = 'data/regression/lm.pca.MSE.RData')


#--------------------------------------------------------
#4. Testen eines Nullmodells zum Vergleich
#--------------------------------------------------------

#Wir nehmen als Vorhersage den Mittleren Aktivitätswert des Pathways
null.prediction = rep(mean(train[, pathway]),nrow(test)) 
#Nun berrechnen wir auch hier den MSE
null.MSE = sum((null.prediction - test[,pathway])^2)/nrow(test)

#Speichern für Vergleich
save(null.prediction, file = 'data/regression/null.prediction.RData')
save(null.MSE, file = 'data/regression/null.MSE.RData')
