
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

filepath = 'Y:/Lena/Bachelorthesis/Experiments/PCR/Mito_gene_expression/150618/csv files/Results_Mitogenes_plate1_LG_150618 - '

dataGenes = read.csv(paste(filepath, 'Quantification Amplification Results_SYBR.csv'), sep=';', dec=',')
dataGenes = dataGenes[,2:ncol(dataGenes)]

#well-content matrix
Qsum = read.csv(paste(filepath, 'Quantification Summary_0.csv'), sep=';')

#######################################################################
namePlate = 'Plate 1'
Genes = as.character(unique(Qsum$Target))
replicates = c(1,2)
Animals = 6
IC = 1
Lines_order = c('HR', 'IR', 'LR') #Lines from top to bottom
#######################################################################

CycNum = length(dataGenes$Cycle)

Lines = c(Lines_order, Lines_order, 'IC', 'NTC')
LinesMatrix = c(rep(Lines[1:(Animals+IC)], each = length(replicates)*length(Genes)*CycNum), rep(Lines[length(Lines)], CycNum*length(Genes)))

#New comprehensive data frame 
dataNew= melt(dataGenes, id=(c('Cycle')))

dataNew$Gene = rep(as.character(Qsum$Target), each=CycNum)
colnames(dataGenes) = c('Cycle', as.character(Qsum$Target))

dataNew$Replicate = c(rep(rep(replicates, each=CycNum), length(Genes)*(Animals+IC)), rep(1,length(Genes)*CycNum))
dataNew$Line= LinesMatrix
dataNew$Plate = namePlate


##Save & import dataNew
savingPath = 'Y:/Lena/Bachelorthesis/Experiments/PCR/Mito_gene_expression/analysis/raw_data/'
saveRDS(dataNew, file=paste(savingPath,'raw_data_', namePlate, '.rds', sep=''))

imported = readRDS('Y:/Lena/Bachelorthesis/Experiments/PCR/Mito_gene_expression/analysis/raw_data/raw_data_Plate 1.rds')


#######################  Graphics
show_only = 'GAPDH' #'all' if to show all
show_cycles = c(15:20) #'all' if to show all

if (show_only == 'all' & show_cycles == 'all'){
  matrix = c(1:length(dataNew$Gene))
} else if (show_only == 'all' & show_cycles != 'all') {
  matrix = which(dataNew$Cycle %in% show_cycles)
} else if (show_only != 'all' & show_cycles == 'all') {
  matrix = which(dataNew$Gene == show_only)
} else {
  matrix = which(dataNew$Gene == show_only & dataNew$Cycle %in% show_cycles)
}

ggplot(data= dataNew[matrix,], aes(x=dataNew$Cycle[matrix], y=dataNew$value[matrix], group=as.factor(dataNew$variable[matrix]), color=as.factor(dataNew$Line[matrix])))+
  geom_line(size = 1)+ scale_x_continuous(breaks=c(5,10,15,20,25,30,35,40))










