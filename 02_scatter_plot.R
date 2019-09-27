
library(dplyr)
library(survival)


#--- Set cox file which is matched to cancer type of dataframe ---#

data.cox = read.csv("C:/Users/Minhy/Desktop/DATA/skcm_cox.csv") 


#--- Filter genes by cox coefficient, top 100 or bottom 100 ---#

data.cor = result.analysis
colnames(data.cor) = c('TCGA.Name','correlation','cor.p.value')

data.final = merge(data.cor,data.cox,by='TCGA.Name',all=FALSE)


#--- Calculate R squared and generate scatter plot ---#

fit =  lm(correlation~Cox.coefficient,data.final)
plot(data.final$correlation,data.final$Cox.coefficient,
     xlab = '',ylab = '',
     pch = 20, col='steelblue')
abline(reg=fit) # draw scatter plot

stat = summary(fit)
stat$r.squared # print R squared value
length(patient.names) # number of patients in study
nrow(data.final) # number of genes used in study
nrow(data.final[data.final$Raw.p.value < 0.01,]) # number of genes with raw p < 0.01

#--- End of whole code ---#
