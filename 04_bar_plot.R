library(ggplot2)

#--- Generate bar plot ---#

Table = read.csv("C:/Users/Minhy/Desktop/Research/cancer_data/Table_1.csv") 
Indicator = (Table$R_squared*Table$Selected_genes)/Table$Total_genes

Table = cbind(Table,Indicator)

ggplot(Table, aes(x=reorder(Cancer_type,Indicator), y=Indicator))+
geom_bar(fill='darkblue',stat='identity')+labs(x="",y="",title="Prognostic sensitivity by cancer types")+
coord_flip()+theme(axis.title=element_text(size=15),plot.title=element_text(size=15))

#--- End of whole code ---#
