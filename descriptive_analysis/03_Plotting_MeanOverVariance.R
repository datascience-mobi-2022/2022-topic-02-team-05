#----------------------------------------------
#In diesem Dokument wird die Varianz und der Mean für jedes Gen für jeden Tumortype geplotted
#----------------------------------------------

#Laden der bereits gecleanten Daten aus der großen tcga matrix - außerdem laden der ungecleanten Daten
load('data/tcga_exp_cleaned.RData')
#load('data/tcga_exp.RData')

#grobe Analyse großer Datensatz Matrix
#str(tcga_exp)
#head(tcga_exp)
#str(tcga_exp_cleaned)
#head(zcga_exp_cleaned)

#Berechnung der Varianz und des Means für alle Gene 
#Var_big <- var(tcga_exp)
#Mea_big <- mean(tcga_exp)

#Var_clean <- var(tcga_exp_cleaned)
Var_clean <- apply(tcga_exp_cleaned, 1, var)
#Mea_clean <- mean(tcga_exp_cleaned)
Mea_clean <- apply(tcga_exp_cleaned, 1, mean)

#Standardisieren der Varianz und logarithmieren des Mittelwerts
#log_mean_clean <- log(Mea_clean)
  #log_mean_clean <- apply(Mea_clean, log)
stand_var_clean <- scale(Var_clean)
#Func1 <- function(x){
  #(Var_clean - mean(Var_clean)) / var(Var_clean)
#}
#stand_var_clean <- apply(Var_clean, Func1)
                         
#plotten der Werte mit Varianz über Mittelwert 
plot(Mea_clean, stand_var_clean, "p", pch = 19, main = "Descriptive analysis", sub = "Variance over mean", xlab = "log_Mean", ylab = "standardized_Variance") -> descriptiveMeanVar 

#Punkte ab gewissem Schwellenwert benennen 
if (log_mean_clean > SCHWELLE & stand_var_big > SCHWELLE) {
  text(log_mean_clean, stand_var_big, labels = row.names(tcga_exp_cleaned), cex = 0.7, pos = 3)
}
