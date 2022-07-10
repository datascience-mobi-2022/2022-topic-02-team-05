#-----------------------------------------
#In diesem Dokument pärparieren wir unsere THCA GSEA DAten um damit lineare
#Regression druchzuführen und ein Neuronales Netz zu trainieren
#-----------------------------------------

load('data/thca_gsea.RData')

#--------------------------------------------------------
#0. Selektion eines geeigneten Pathways
#--------------------------------------------------------
#Eine Regression ist nur dann sinnvoll wenn der Pathway überhaupt eine ausreichend
#hohe varianz hat (sonst kann man auch ein null model mit Mittelwert benutzen)
#=> Wir suchen Pathways die sich signifikant in Tumor und Normalgewebe unterscheiden
# und zusätzlich hochvariant sind.

load('data/regression/thca_pval_gsva.RData') #Pvalues der Pathways
path.sig = sort(thca_pval_gsva)[1:25] #Gibt die top 25 signifikant unterschiedlichen pathways

path.var = apply(thca_gsea, 1, var)
path.topvar = path.var[path.var > quantile(path.var, probs = 0.8)] #top 20% der variantesten pathways

#Gibt uns die Schnittmenge aus Hochvarianten und Signifikant unterschiedlichen Pathways
path.selection = path.topvar[na.omit(match(names(path.sig), names(path.topvar)))]
path.selection; rm(thca_pval_gsva, path.sig, path.var, path.topvar)


pathway = 'RODRIGUES_DCC_TARGETS_UP' #Ein set von Oncogenen gefunden in Colonkarzinomen

#--------------------------------------------------------
#1. Split in Trainings und Testdaten
#--------------------------------------------------------


#Splitten der Daten in eine Test und Trainingsgruppe (75% für Training)
set.seed(1) #Für konsistente Indices da sample() random ist
index = sample(1:ncol(thca_gsea), round(0.75*ncol(thca_gsea)))
train = as.data.frame(t(thca_gsea[, index]))
test = as.data.frame(t(thca_gsea[, -index]))

save(train, file = 'data/regression/train.RData')
save(test, file = 'data/regression/test.RData')

# #--------------------------------------------------------
# #2.PCA transformation für unkorrelierte Werte für die Regression
# #--------------------------------------------------------
# 
# thca_pathway = thca_gsea[pathway, ] #Werte für unseren Pathway
# thca_pca_data = prcomp(t(thca_gsea[!rownames(thca_gsea) == pathway, ]))$x #PCA für alle Werte ohne pathway
# 
# train_pca = cbind.data.frame(thca_pathway[index], thca_pca_data[index, ])  
#   colnames(train_pca)[1] = pathway
# test_pca = cbind.data.frame(thca_pathway[-index], thca_pca_data[-index, ])
#   colnames(test_pca)[1] = pathway
#   
# save(train_pca, file = 'data/regression/train_pca.RData')
# save(test_pca, file = 'data/regression/test_pca.RData')

#--------------------------------------------------------
#3. Daten Präparation für Neuronales Netz
#--------------------------------------------------------
#Skalieren der Daten nach mit einer Min/max skalierung sodass sie im Intervall [0,1] liegen
maxs = apply(thca_gsea, 1, max) 
mins = apply(thca_gsea, 1, min)
thca_scaled = as.data.frame(scale(t(thca_gsea), center = mins, scale = maxs - mins))
train_sc = thca_scaled[index,]
test_sc = thca_scaled[-index,]

save(train_sc, file = 'data/regression/train_sc.RData')
save(test_sc, file = 'data/regression/test_sc.RData')

