#----------------------------------------------
#In diesem Dokument wird die Varianz und der Mean für jedes Gen für jeden Tumortype geplotted
#----------------------------------------------

#Laden der bereits gecleanten Daten aus der großen tcga matrix 
load('data/tcga_exp_cleaned.RData')

#load ggplot2 
library(ggplot2)

#Variance of every gene
Var_clean <- apply(tcga_exp_cleaned, 1, var)
#mean of every gene 
Mea_clean <- apply(tcga_exp_cleaned, 1, mean)
x = as.vector(Mea_clean)

#Standardisieren der Varianz
stand_var_clean <- scale(Var_clean)
y = as.vector(stand_var_clean)

#data frame out of Mea_clean and stand_var_clean
as.vector(row.names(tcga_exp_cleaned)) -> names
data.frame(x, y, names) -> df_VoverM

#plotten der Werte mit Varianz über Mittelwert 
Var_over_Mean_plot = 
  ggplot(df_VoverM, aes(x, y, label = names)) + 
  geom_point(shape = 21, size = 1.5, fill = "red") +
  geom_text(data = subset(df_VoverM, y > 5), size = 2) + 
  ggtitle("Variance over Mean") + 
  xlab("Mean") + 
  ylab("Standardized Variance") 
 
  #annotate("text", x = -10:15, y = 5:10.5, label = RownamesTCGA)
  #ggplot(Mea_clean, stand_var_clean, "p", pch = 19, main = "Variance over mean", xlab = "Mean", ylab = "standardized_Variance") -> descriptiveMeanVar 
 
#Punkte ab gewissem Schwellenwert benennen 
if (y > 5) {
  text(Mea_clean, stand_var_clean, labels = row.names(tcga_exp_cleaned), cex = 0.7, pos = 3)
}
Var_over_Mean_plot


