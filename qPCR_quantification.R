### qPCR data quantification with the comparative delta-delta-2G method

library(ggplot2)
library(gridExtra)
library(cowplot)
library(grid)
library(ggsignif)
library(scales)
library(RColorBrewer)
library(xlsx)

#### set some parameters and stuff
# setwd('C:/Users/Gschossmann/Dropbox/studies/Osnabrück/Universität/Bachelorarbeit_DroBo/experiments/Analysis/PCR')
# filepath = 'Y:/Lena/Bachelorthesis/Experiments/PCR/Mito_gene_expression/analysis'
#testing:
filepath = 'C:/Users/lena_/Dropbox/studies/Osnabrück/Universität/Bachelorarbeit_DroBo/experiments/Analysis/PCR/'

numRepl = 2 #number of technical replicates
HKG = c('HPRT', 'GAPDH') #number of housekeeping genes

# import data from CFX manager exported excel table that has been collected in 1 excel
# the imported data table shall have following columns: plate, date, threshold, animal, group, well, Gene, Cq value
data_tot = read.csv(paste(filepath, 'data_tot.csv', sep='/'), dec='.', sep=';')
colnames(data_tot) = c('Threshold', 'Plate', 'Well', 'Gene', 'Animal', 'Group', 'Cq')
data_tot = data_tot[which(data_tot$Animal != 'NTC'),]  #Kick out NTCs

# import the efficiencies
# effs = read.csv(paste(filepath, 'efficiency calculation/efficiency_values.csv', sep='/'))
#testing:
effs = read.csv(paste(filepath, 'efficiency_values.csv', sep='/'))

##### get some fix variables extracted
plates= unique(data_tot$Plate)
animals = unique(data_tot$Animal)
# animals_perPlate
genes = unique(data_tot$Gene)
GOI = genes[!(genes %in% HKG)]
genes = append(HKG, as.character(GOI)) #so that the HKG come as first elements
Group = c('HR', 'LR')
control = 'IR'
#####

dataAveraged = data.frame(matrix(nrow=(nrow(data_tot)/numRepl), ncol=(ncol(data_tot)+1)))
colnames(dataAveraged) = c(colnames(data_tot), 'CV')

outlier = data.frame(matrix(ncol=ncol(data_tot)))
colnames(outlier) = colnames(data_tot)

################################################### Check duplicate precision & combine replicates into one averaged value
cnt=1
tmp_warning_CV =''
iO = 1 #for counting outlier
#loop through plates, animals, genes, replicates
for (iPl in plates){
  
  for (iA in unique(data_tot$Animal[which(data_tot$Plate == iPl)])){
    
    for (iG in unique(data_tot$Gene[which(data_tot$Plate == iPl & data_tot$Animal == iA)])){
      tmp =0 #necessary to initialize tmp for append function
      tmpMatrix = which(data_tot$Plate == iPl & data_tot$Animal == iA & data_tot$Gene == iG)
      for (iR in 1:numRepl){
         tmp = append(data_tot$Cq[tmpMatrix][iR], tmp)
        }
        tmp = tmp[1:numRepl] #get rid of 0 in beginning
        #reduce replicates into one mean value and store in new dataframe
        # dataAveraged[cnt,1:(ncol(dataAveraged)-1)] = data_tot[tmpMatrix[1],]
        if((max(tmp)-min(tmp)) < 1){ #exclude outliers | criterion: difference between replicates
          dataAveraged$Threshold[cnt] = data_tot$Threshold[tmpMatrix[1]]
          dataAveraged$Plate[cnt] = data_tot$Plate[tmpMatrix[1]]
          dataAveraged$Well[cnt] = as.character(data_tot$Well[tmpMatrix[1]])
          dataAveraged$Gene[cnt] = as.character(data_tot$Gene[tmpMatrix[1]])
          dataAveraged$Group[cnt] = as.character(data_tot$Group[tmpMatrix[1]])
          dataAveraged$Animal[cnt] = as.character(data_tot$Animal[tmpMatrix[1]])
          dataAveraged$Cq[cnt] = mean(tmp)
          tmpCV = sd(tmp)/mean(tmp) 
          dataAveraged$CV[cnt] = tmpCV
        } else {
          for(iR in 1:numRepl){ # write outlier in dataframe
            outlier = rbind(outlier[], data_tot[tmpMatrix[iR],])
            iO = iO+1
          }
          # take first replicates value anyway (this should be commented out later)
          dataAveraged$Threshold[cnt] = data_tot$Threshold[tmpMatrix[1]]
          dataAveraged$Plate[cnt] = data_tot$Plate[tmpMatrix[1]]
          dataAveraged$Well[cnt] = as.character(data_tot$Well[tmpMatrix[1]])
          dataAveraged$Gene[cnt] = as.character(data_tot$Gene[tmpMatrix[1]])
          dataAveraged$Group[cnt] = as.character(data_tot$Group[tmpMatrix[1]])
          dataAveraged$Animal[cnt] = as.character(data_tot$Animal[tmpMatrix[1]])
          dataAveraged$Cq[cnt] = min(tmp)
          tmpCV = 0 
          dataAveraged$CV[cnt] = 0
        }
        cnt=cnt+1
    }
  }
}

################################################### IC-Normalization
# 1) Average IC of different plates
# 2) Correct averaged IC-GOI by averaged IC-HKG
# 3) Set the resulting differences for GOI as 100%
# 4) Correct each plates IC-GOI by the according IC-HKG
# 5) Relate the in 4) calculated differences to the averaged IC differences --> Correction factor per plate

# ! formula for correction takes differential efficiencies into account: 

#split dataAveraged into data and tmp_IC
tmp_IC = dataAveraged[dataAveraged$Animal == 'IC',]
data_Genes = dataAveraged[dataAveraged$Animal != 'IC',]
IC_corrFactors = data.frame(matrix(ncol=(length(GOI)+1), nrow=length(plates)))
colnames(IC_corrFactors) = append('Plate', as.character(GOI))
IC_corrFactors$Plate[] = 1:length(plates)

av_corr_IC = data.frame(matrix(ncol=2, nrow=length(GOI)))
colnames(av_corr_IC) = c('Gene', 'Av_Corr_Value')
av_corr_IC[,1] = GOI

# 1) & 2)
tmpMatrix_HKG = tmp_IC$Gene %in% HKG
av_IC_HKG = mean(tmp_IC$Cq[tmpMatrix_HKG]) # Averaged IC HKG means

for (iG in 1:length(GOI)){ #average all IC values per gene and then correct GOI by HKG
  tmpMatrix_GOI = tmp_IC$Gene %in% GOI[iG]
  av_corr_IC$Gene[iG] = GOI[iG]
  av_corr_IC$Av_Corr_Value[iG] = mean(tmp_IC$Cq[tmpMatrix_GOI]) - av_IC_HKG
}

# 4) & 5)
for (iPl in 1:length(plates)){
  tmp_avHKG = mean(tmp_IC$Cq[tmp_IC$Plate == plates[iPl] & tmpMatrix_HKG])
  for (iG in 1:length(GOI)){
    IC_corrFactors[iPl,(iG+1)] = (tmp_IC$Cq[tmp_IC$Plate == plates[iPl] & tmp_IC$Gene == GOI[iG]]-tmp_avHKG)/av_corr_IC$Av_Corr_Value[iG]
  }
}


tmpMatrix_HKG = data_Genes$Gene %in% HKG
data_corrected = data_Genes[!tmpMatrix_HKG,]

for (iPl in 1:length(plates)){
  for (iA in unique(data_tot$Animal[which(data_tot$Plate == iPl)])){
    tmp_avHKG = mean(data_Genes$Cq[ tmpMatrix_HKG & data_Genes$Plate == plates[iPl] & data_Genes$Animal == iA])
    for (iG in 1:length(GOI)){
      tmpCq = data_corrected$Cq[data_corrected$Plate == plates[iPl] & data_corrected$Animal == iA & data_corrected$Gene == GOI[iG]] - tmp_avHKG #correct by the animals housekeeping gene
      tmpCq = tmpCq * IC_corrFactors[iPl, (iG+1)] #Correct by IC correction factor
      data_corrected$IC_corr_Cq[data_corrected$Plate == plates[iPl] & data_corrected$Animal == iA & data_corrected$Gene == GOI[iG]] = tmpCq
    }
  }
}



################################################## Normalization of each LR, HR animal by IR-average
# split data_corrected into Control (IR) and Treatment (LR, HR)
tmp_Ctrl = data_corrected[data_corrected$Group == 'IR',]
data_corrected = data_corrected[!(data_corrected$Group == 'IR'),]

av_corr_Ctrl = data.frame(matrix(ncol=2, nrow=length(GOI)))
colnames(av_corr_Ctrl) = c('Gene', 'Ctrl_Av_corr_Cq')


for (iG in 1:length(GOI)){
  av_corr_Ctrl$Gene[iG] = as.character(GOI[iG])
  av_corr_Ctrl$Ctrl_Av_corr_Cq[iG] = mean(tmp_Ctrl$IC_corr_Cq[tmp_Ctrl$Gene == GOI[iG]])
}

for(iX in 1:nrow(data_corrected)){
  tmpGene = data_corrected$Gene[iX]
  data_corrected$Ctrl_corr_Cq[iX] = data_corrected$IC_corr_Cq[iX] - av_corr_Ctrl$Ctrl_Av_corr_Cq[av_corr_Ctrl$Gene == tmpGene]
}


################################################## Calculate Fold Change, Mean Difference & SEM

data_corrected$Fold_change[1:nrow(data_corrected)] = 2^(-(data_corrected$Ctrl_corr_Cq[]))

#Fold change summary
data_summary = data.frame(matrix(ncol=4, nrow=length(Group)*length(GOI))) #Mean and SEM per GOI per Group
colnames(data_summary) = c('Group', 'Gene', 'Mean', 'SEM')

cnt=1
for(iGr in 1:length(Group)){
  
  for(iG in 1:length(GOI)){
    tmpMatrix = data_corrected$Group == Group[iGr] & data_corrected$Gene == GOI[iG]
    data_summary$Group[cnt] = as.character(Group[iGr])
    data_summary$Gene[cnt] = as.character(GOI[iG])
    data_summary$Mean[cnt] = mean(data_corrected$Fold_change[tmpMatrix])
    data_summary$SEM[cnt] = sd(data_corrected$Fold_change[tmpMatrix])/sqrt(sum(tmpMatrix))
    cnt=cnt+1
  }
}

################################################## Plots

ggplot(data=data_corrected, aes(x=data_corrected$Gene, y=data_corrected$Fold_change, fill=data_corrected$Group))+
  geom_bar(stat='identity', position='dodge')



################################ Output

print(paste('CV was smaller than it should be in following replicates:', tmp_warning_CV))








