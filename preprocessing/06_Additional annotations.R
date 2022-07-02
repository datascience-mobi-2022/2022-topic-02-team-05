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
              'Leucemia', #dlbc
              'Carcinoma', #ESCA
              'Glioblastoma', #GBM
              'Squamous cell carcinoma', #HNSC
              'Adenocarcinoma', #KICH
              'Adenocarcinoma', #KIRC
              'Adenocarcinoma', #KIRP
              'Leucemia', #LAML
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
