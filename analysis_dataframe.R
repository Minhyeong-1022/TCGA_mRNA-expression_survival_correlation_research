
data.cor = result.analysis
colnames(data.cor) = c('TCGA.Name','correlation','cor.p.value')

data.cox = read.csv("C:/Users/Minhy/Desktop/DATA/skcm_cox.csv") # set file which is matched to cancer type of dataframe

data.final = merge(data.cor,data.cox,by='TCGA.Name',all=FALSE)

top.results = tail(arrange(data.final,desc(Cox.coefficient)),n=100) # set arrange method for analysis (tail or head)
top.genes.index = match(top.results$TCGA.Name,gene.names)


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

title(main = 'Filtered by Cox - OV')

###################

plot(data.final$correlation,data.final$Cox.coefficient)
fit =  lm(correlation~Cox.coefficient,data.final)
stat = summary(fit)
stat




