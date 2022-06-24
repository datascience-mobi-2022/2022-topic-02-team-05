#--------------------------------------
# In diesem Dkument nehemn wir die GSEA Daten und f?hren eine PCA sowie wie
# ein clustering der transformierten Daten durch so k?nnen wir hoffentlich verschiedenen
# unterschiedliche Pathway activt?ten zuordnen
#------------------------------------------
library(ggplot2)
library(metaplot)
library(gridExtra)
   

load('~/GitHub/2022-topic-02-team-05/data/tcga_gsva.RData')
tcga_anno = readRDS('~/GitHub/2022-topic-02-team-05/data/tcga_tumor_annotation.RDS')

#--------------------------------------------
#Durchf?hren der PCA
#---------------------------------------------
PCA = prcomp(t(tcga_gsva))

#Varianz die erkl?rt wird
PCs = 1:658; var_prop = PCA$sdev^2/sum(PCA$sdev^2); var_cum = cumsum(var_prop)
Var = cbind.data.frame(PCs, var_prop, var_cum)
v1 = ggplot(Var[1:10,], aes(PCs, var_prop)) + geom_point()+
  labs(y = "Proportion of Variance Explained")+
  scale_x_continuous(name = 'Principle Components', breaks = 1:10)

v2 = ggplot(Var[1:10,], aes(PCs, var_cum)) + geom_point()+
  labs(y = "Cummulative Proportion of Variance Explained")+
  scale_x_continuous(name = 'Principle Components', breaks = 1:10)
grid.arrange(grobs = list(v1,v2), ncol= 2)

#visualisierung der pathways die in einer PC stecken
pathways = as.data.frame(PCA$rotation)
pathways$names = row.names(pathways)

pathways = pathways[order(abs(pathways$PC1), decreasing = TRUE),]
t1 = ggplot(pathways[1:10,], aes(names, PC1)) + geom_col()+
  labs(x = 'Pathway', y= 'Contribution to PC1')+
  theme(axis.text = element_text(angle = 90,))

pathways = pathways[order(abs(pathways$PC2), decreasing = TRUE),]
t2 = ggplot(pathways[1:10,], aes(names, PC2)) + geom_col()+
  labs(x = 'Pathway', y= 'Contribution to PC2')+
  theme(axis.text = element_text(angle = 90,))

pathways = pathways[order(abs(pathways$PC3), decreasing = TRUE),]
t3 = ggplot(pathways[1:10,], aes(names, PC3)) + geom_col()+
  labs(x = 'Pathway', y= 'Contribution to PC3')+
  theme(axis.text = element_text(angle = 90,))
grid.arrange(grobs = list(t1, t2, t3), ncol = 3)


#plotten der daten gef?rbt nach tumortype
PCA_data = as.data.frame(PCA$x) 
PCA_data = PCA_data[order(row.names(PCA_data)),]

tcga_anno = tcga_anno[1:800,]
tcga_anno = tcga_anno[order(tcga_anno$sample),]

PCA_data$type = ifelse(tcga_anno$cancer_type_abbreviation == 'BRCA', 'BRCA',
                 ifelse(tcga_anno$cancer_type_abbreviation == 'LUAD', 'LUAD',
                  ifelse(tcga_anno$cancer_type_abbreviation == 'KIRC', 'KIRC',
                   ifelse(tcga_anno$cancer_type_abbreviation == 'PRAD', 'PRAD',
                    ifelse(tcga_anno$cancer_type_abbreviation == 'THCA', 'THCA','other')
                    )
                   )
                  )
                 )

p1 = ggplot(PCA_data, aes(PC1, PC2)) + geom_point(size = 1,aes(colour = type)) +
  scale_color_manual("Tumor type", values = c('BRCA'='red','LUAD'='blue','KIRC'='green','PRAD'='orange','THCA'='black','other'='grey'))+
  theme(panel.background = element_rect(fill = 'white'),
        panel.grid.major = element_line(colour = "darkgrey", size=0.25),
        panel.grid.minor = element_line(colour = "darkgrey", size=0.25),
        panel.border = element_rect(colour = 'black', fill = NA, size=0.25))
        #aspect.ratio = c(1,1))
  
p2 = ggplot(PCA_data, aes(PC1, PC3)) + geom_point(size = 1,aes(colour = type)) +
  scale_color_manual("Tumor type", values = c('BRCA'='red','LUAD'='blue','KIRC'='green','PRAD'='orange','THCA'='black','other'='grey'))+
  theme(legend.position="none",
        panel.background = element_rect(fill = 'white'),
        panel.grid.major = element_line(colour = "darkgrey", size=0.25),
        panel.grid.minor = element_line(colour = "darkgrey", size=0.25),
        panel.border = element_rect(colour = 'black', fill = NA, size=0.25))
        #aspect.ratio = c(1,1))

p3 = ggplot(PCA_data, aes(PC2, PC3)) + geom_point(size = 1,aes(colour = type)) +
  scale_color_manual("Tumor type", values = c('BRCA'='red','LUAD'='blue','KIRC'='green','PRAD'='orange','THCA'='black','other'='grey'))+
  theme(legend.position="none",
        panel.background = element_rect(fill = 'white'),
        panel.grid.major = element_line(colour = "darkgrey", size=0.25),
        panel.grid.minor = element_line(colour = "darkgrey", size=0.25),
        panel.border = element_rect(colour = 'black', fill = NA, size=0.25))
        #aspect.ratio = c(1,1))

# grid.arrange(grobs = list(p1, p2, p3),ncol = 3,  top = 'PCA of GSEA data',
#              widths = c(1,1,1.4))

grid.arrange(p1, arrangeGrob(p2,p3, ncol=2), heights=c(2.5/4, 2.5/4), ncol=1,top = 'PCA of GSEA data')

grid.arrange(arrangeGrob(v1, v2, ncol = 2),
             arrangeGrob(t1,t2,t3, ncol = 3),
             arrangeGrob(p1,p2,p3, ncol = 3),
             ncol=1)
