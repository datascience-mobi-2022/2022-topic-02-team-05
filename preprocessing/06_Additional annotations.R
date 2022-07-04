#--------------------------
#In diesem Dokument fügen wir unseren Krebsannotationsdaten weitere hinzu,
#zB Art des Krebs
#--------------------------

tcga_anno = readRDS('~/GitHub/2022-topic-02-team-05/data/tcga_tumor_annotation.RDS')

#Data frame der Jedenkrebstyp mit der Foorm des Krebs verbindet
canc_form = c('Adenocarcinoma', #ACC
              'Transitional cell carcinoma', #BLCA
              'Carcinoma', #BRCA
              'Carcinoma', #CESC
              'Adenocarcinoma', #CHOL
              'Adenocarcinoma', #COAD
              'Leukemia', #dlbc
              'Carcinoma', #ESCA
              'Glioblastoma', #GBM
              'Squamous cell carcinoma', #HNSC
              'Adenocarcinoma', #KICH
              'Adenocarcinoma', #KIRC
              'Adenocarcinoma', #KIRP
              'Leukemia', #LAML
              'Glioblastoma', #LGG
              'Adenocarcinoma', #LIHC
              'Adenocarcinoma', #LUAD
              'Squamous cell carcinoma', #LUSC
              'Sarcoma', #MESO
              'Adenocarcinoma', #OV
              'Adenocarcinoma', #PAAP
              'Glioblastoma', #PCPG
              'Adenocarcinoma', #PRAD
              'Adenocarcinoma', #READ
              'Sarcoma', #SARC
              'Melanoma', #SKCM
              'Adenocarcinoma', #STAD
              'Carcinoma', #TGCT
              'Carcinoma', #THCA
              'Carcinoma', #THYM
              'Adenocarcinoma', #UCEC
              'Sarcoma', #UCS
              'Melanoma') #UVM
#Erstelen eines df mit Krebs abbkürzung und der passenden Krebsform
Cancer_form_anno = cbind.data.frame(names(table(tcga_anno$cancer_type_abbreviation)), canc_form)
colnames(Cancer_form_anno) = c('Cancer_Abb', 'Cancer_form')

#Erweitern der Annotation um eine Spalte mit der Krebsform
tcga_anno$cancer_form = sapply(tcga_anno$cancer_type_abbreviation, FUN = function(y){
  Cancer_form_anno$Cancer_form[y == Cancer_form_anno$Cancer_Abb]
})

save(tcga_anno, file = 'data/tcga_anno.RData')
save(Cancer_form_anno, file = 'data/Cancer_form_anno.RData')


#-------------------------------------------------
#Anschließend bennen wir die Histological grade für unsere THCA annotation um
#sodass nicht Thyroid Papillary Carcinoma davor steht
#-------------------------------------------------
load('data/thca_anno.RData')

#Umbennen der Histologischen stufe (entfernen der unnötigen vorsatzes THCA)
thca_anno$histological_type = sapply(thca_anno$histological_type, FUN = function(x){
  return(strsplit(x, split = '- ', fixed = TRUE)[[1]][2])})
#Alle NAs werden other genannt
thca_anno$histological_type[is.na(thca_anno$histological_type)] = 'other'

save(thca_anno, file = 'data/thca_anno.RData')

