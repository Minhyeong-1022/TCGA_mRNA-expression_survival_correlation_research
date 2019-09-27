#--- Generate kaplan meier plot with selected index ---#

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
lLab = c("mRNA High","mRNA Low")  ## legend labels
legend("topright",legend=lLab,col=1:2,lty=1:1,horiz=FALSE,,bty='n')

#--- End of whole code ---#
