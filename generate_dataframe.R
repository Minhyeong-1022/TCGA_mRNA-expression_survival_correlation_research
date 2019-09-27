if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("RTCGAToolbox", version = "3.8")
install.packages("stringi")
install.packages("rgl")
install.packages("vioplot")

library(plyr)
library(rgl)
library(vioplot)
library(stringi)
library(survival)
library(RTCGAToolbox)

cancer.type = "SKCM"  # Set Parameter

lastRunDate <- getFirehoseRunningDates()[1]
caller.v = getFirehoseData (dataset=cancer.type, runDate=lastRunDate,
                            gistic2_Date=getFirehoseAnalyzeDates(1),
                            Clinic=TRUE, RNASeq2GeneNorm=TRUE)
sheet.clinical <- getData(caller.v,"clinical")
sheet.rna <- getData(caller.v,"RNASeq2GeneNorm")
gene.names = rownames(sheet.rna)
sample.names = substr(colnames(sheet.rna),1,12)
patient.names = rownames(sheet.clinical)
patient.names =  toupper(patient.names)
substr(patient.names,5,5)<-'-'
substr(patient.names,8,8)<-'-'
death.index = which(is.na(sheet.clinical$days_to_last_followup)==TRUE)
clinical.survival = sheet.clinical$days_to_last_followup
clinical.survival[death.index] = sheet.clinical$days_to_death[death.index]
clinical.status = sheet.clinical$vital_status
clinical.data.frame = as.data.frame(cbind(patient.names,clinical.status,clinical.survival)) 

rna.surv.list = list() # to store of expression and survival each genes 
rna.surv.list.kaplan = list()
m.diff = c()
correlation = c()
p.value = c()

for (i in 1:length(gene.names)) {
  rna.surv.df = as.data.frame(cbind(sample.names,sheet.rna[i,]))
  rownames(rna.surv.df) = seq(1,nrow(rna.surv.df),1)
  colnames(rna.surv.df) = c('patient.names','rna.exp')
  rna.surv.list[[i]] = as.data.frame(merge(rna.surv.df,clinical.data.frame,by='patient.names',all.x=TRUE))
  rna.surv.list[[i]]$rna.exp = as.numeric(as.character(rna.surv.list[[i]]$rna.exp))
  rna.surv.list[[i]]$clinical.survival = as.numeric(as.character(rna.surv.list[[i]]$clinical.survival))
  
  mRNAexp = rna.surv.list[[i]]$rna.exp[rna.surv.list[[i]]$clinical.status == 1]
  survival_time = rna.surv.list[[i]]$clinical.survival[rna.surv.list[[i]]$clinical.status == 1]
  
  analysis.sheet = cbind(mRNAexp,survival_time)
  analysis.sheet = na.omit(analysis.sheet)
  
  correlation.test = cor.test(analysis.sheet[,1],analysis.sheet[,2])
  correlation[i] = correlation.test$estimate
  p.value[i] = correlation.test$p.value
}

for (i in 1:length(gene.names)) {
  rna.exp.median = median(rna.surv.list[[i]]$rna.exp)
  if (rna.exp.median > 1) {
    surv.in.high = rna.surv.list[[i]][rna.surv.list[[i]]$rna.exp > rna.exp.median,]
    surv.in.low = rna.surv.list[[i]][rna.surv.list[[i]]$rna.exp < rna.exp.median,]
    
    group.1 = rep('High',nrow(surv.in.high))
    group.2 = rep('Low',nrow(surv.in.low))
    group.integrated = c(group.1,group.2)
    surv.kaplan = data.frame(rbind(surv.in.high,surv.in.low))
    rna.surv.list.kaplan[[i]] = data.frame(surv.kaplan,group.integrated)
  } else rna.surv.list.kaplan[[i]] = NA
}

result.analysis = as.data.frame(cbind(gene.names,correlation,p.value))
#result.analysis$m.diff = as.numeric(as.character(result.analysis$m.diff))
result.analysis$correlation = as.numeric(as.character(result.analysis$correlation))
result.analysis$p.value = as.numeric(as.character(result.analysis$p.value))
colnames(result.analysis) = c('gene.names','correlation','p.value')
result.analysis = na.omit(result.analysis)

# --- End of whole code ---#

