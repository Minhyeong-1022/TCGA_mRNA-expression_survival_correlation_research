
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

#--- End of whole code ---#
