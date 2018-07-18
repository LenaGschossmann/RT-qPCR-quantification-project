#MCMC.qpcr package

library(ggplot2)
library(gridExtra)
library(cowplot)
library(grid)
library(ggsignif)
library(scales)
library(RColorBrewer)
library(RDML)
library(reshape2)
library(MCMC.qpcr)

filepath = 'Y:/Lena/Bachelorthesis/Experiments/PCR/Mito_gene_expression/150618/csv files/'

dataCq = read.csv(paste(filepath, 'Threshold_Cq_values.csv', sep=''), sep=';', dec='.')

#######################################################################
namePlate = 'Plate 1'
Genes = as.character(unique(dataCq$Target))
replicates = c(1,2)
Animals = 6
IC = 1
Lines_order = c('HR', 'IR', 'LR') #Lines from top to bottom

###############################################################

Lines = c(Lines_order, Lines_order, 'IC', 'NTC')

dataCq$Line = c(rep(Lines[1:(Animals+IC)], each = length(replicates)*length(Genes)), rep(Lines[length(Lines)], length(Genes)))
dataCq$Replicate = rep(c(1:(nrow(dataCq)/2)), each=2)
  

#Remodel dataframe
dataNew = data.frame(matrix(ncol=(length(Genes)+2)))
colnames(dataNew) = c('sample', 'Line', Genes)

cntA = 1
cntB = 1

for(idx in 1:(Animals)){
  dataNew[(cntA:(cntA+length(replicates)-1)),1] = idx
  dataNew[(cntA:(cntA+length(replicates)-1)),2] = c(Lines[idx])
  
  for (iGene in Genes){
    dataNew[(cntA:(cntA+length(replicates)-1)), which(colnames(dataNew) == iGene)] = dataCq$Cq[(cntB:(cntB+length(replicates)-1))]
    cntB=cntB+length(replicates)
  }
  
  cntA=cntA+length(replicates)
}





# Efficiencies

efficiencies = data.frame(matrix(ncol=2, nrow=length(Genes)))
colnames(efficiencies) = c('Gene', 'Eff')

efficiencies$Gene = Genes
efficiencies$Eff = c(2.09, 1.99, 1.98, 2.03, 1.99, 2.11)


##################

qs = cq2counts(data=dataNew, effic=efficiencies, genecols=c(3:length(dataNew)), condcols=c(1:2))
# qs$Line=relevel(qs$Line, ref='IR')
naive = mcmc.qpcr(data=qs, fixed="Line")

s1 = HPDsummary(model=naive, data=qs)
s2 = HPDsummary(model=naive, data=qs, relative = TRUE)



# Save in complete data file
dataTot = dataNew

