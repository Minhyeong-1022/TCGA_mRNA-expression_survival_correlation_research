library(survival)


#--- Set cox file which is matched to cancer type of dataframe ---#

data.cox = read.csv("C:/Users/Minhy/Desktop/DATA/skcm_cox.csv") 


#--- Filter genes by cox coefficient, top 100 or bottom 100 ---#

data.cor = result.analysis
colnames(data.cor) = c('TCGA.Name','correlation','cor.p.value')

data.final = merge(data.cor,data.cox,by='TCGA.Name',all=FALSE)

top.results = tail(arrange(data.final,desc(Cox.coefficient)),n=100) # set arrange method for analysis (tail or head)
top.genes.index = match(top.results$TCGA.Name,gene.names)


#--- Generate kaplan meier plot with selected index ---#

selected.index = c()
for (i in top.genes.index) {
  if (class(rna.surv.list.kaplan[[i]]) == 'data.frame') {
    selected.index = append(selected.index,i)
  }
}

data.kaplan = data.frame() # select top genes 
for (i in selected.index) {
  time = rna.surv.list.kaplan[[i]]$clinical.survival
  event = rna.surv.list.kaplan[[i]]$clinical.status
  group = rna.surv.list.kaplan[[i]]$group.integrated
  column = as.data.frame(cbind(time,event,group))
  data.kaplan = as.data.frame(rbind(data.kaplan,column))
}

km = survfit(Surv(data.kaplan$time,data.kaplan$event) ~as.factor(data.kaplan$group)) # draw kaplan meier plot
plot(km,col=c(1,2))
lLab = c("mRNA High","mRNA Low")  ## legend labels
legend("topright",legend=lLab,col=1:2,lty=1:1,horiz=FALSE,,bty='n')


#--- Calculate R squared and generate scatter plot ---#

fit =  lm(correlation~Cox.coefficient,data.final)
plot(data.final$correlation,data.final$Cox.coefficient,
     xlab = '',ylab = '',
     pch = 20, col='steelblue')
abline(reg=fit) # draw scatter plot

stat = summary(fit)
stat$r.squared # print R squared value


#--- End of whole code ---#


