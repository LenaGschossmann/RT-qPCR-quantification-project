##### Stats business

library(car)

################################################## Mito gene expression

data_corrected$Rel_quantity=round(data_corrected$Rel_quantity, digits=3)

#Test for homogeneity of variance
for(iG in unique(data_corrected$Gene)){
  tmpData = data_corrected[data_corrected$Gene== iG,]
  print(iG)
  print(leveneTest(IC_corr_Cq ~ Group, data=tmpData))
}

##Kruskal Wallis
# sink('Kruskal-Wallis-Test.txt')
for (iG in 1:length(GOI)){
  varNames = GOI
  print(as.character(GOI[iG]))
  print(kruskal.test(data=data_corrected,data_corrected$Rel_quantity[data_corrected$Gene == as.character(GOI[iG])] ~ as.factor(data_corrected$Group[data_corrected$Gene == as.character(GOI[iG])])))
}
# sink()


#Man-Whitney-U Test
groupList = list(c('HR', 'IR'), c('IR', 'LR'), c('HR', 'LR'))
for (groups in groupList){
  print(paste(groups[1], 'vs.', groups[2]))
  for (iG in 1:length(GOI)){
    tmpTreat = data_corrected[data_corrected$Group %in% groups,]
    varNames = GOI
    print(as.character(GOI[iG]))
    print(wilcox.test(data=tmpTreat, tmpTreat$Rel_quantity[tmpTreat$Gene == as.character(GOI[iG])] ~ as.factor(tmpTreat$Group[tmpTreat$Gene == as.character(GOI[iG])])))
  }
}


for(iG in GOI){
  for(iTr in unique(data_corrected$Group)){
    print(sprintf('Gene: %s | Group; %s | n = %i',iG, iTr, nrow(data_corrected[data_corrected$Gene == iG & data_corrected$Group == iTr,])))
  }
}



###################################################### DeltaDelta Method
for (iG in 1:length(GOI)){
  varNames = GOI
  print(as.character(GOI[iG]))
  print(kruskal.test(data=data_corrected,data_corrected$DeltaDelta[data_corrected$Gene == as.character(GOI[iG])] ~ as.factor(data_corrected$Group[data_corrected$Gene == as.character(GOI[iG])])))
}

groupList = list(c('HR', 'IR'), c('IR', 'LR'), c('HR', 'LR'))
for (groups in groupList){
  print(paste(groups[1], 'vs.', groups[2]))
  for (iG in 1:length(GOI)){
    tmpTreat = data_corrected[data_corrected$Group %in% groups,]
    varNames = GOI
    print(as.character(GOI[iG]))
    print(wilcox.test(data=tmpTreat, tmpTreat$DeltaDelta[tmpTreat$Gene == as.character(GOI[iG])] ~ as.factor(tmpTreat$Group[tmpTreat$Gene == as.character(GOI[iG])])))
  }
}

##################################################### Apoptosis Ratio

print(kruskal.test(data=data_Ratios,data_Ratios$Bcl2_over_Bax ~ as.factor(data_Ratios$Group)))

for(iGr in Group){
  tmp_data_Ratios = data_Ratios[data_Ratios$Group == Control | data_Ratios$Group == iGr,]
  print(paste(iGr, 'vs.', Control))
  print(wilcox.test(data=tmp_data_Ratios, tmp_data_Ratios$Bcl2_over_Bax ~ as.factor(tmp_data_Ratios$Group)))
}

groupList = list(c('HR', 'IR'), c('IR', 'LR'), c('HR', 'LR'))
for (groups in groupList){
  print(paste(groups[1], 'vs.', groups[2]))
  print(wilcox.test(data=data_Ratios[data_Ratios$Group == groups[1] | data_Ratios$Group == groups[2],],
                  data_Ratios$Bcl2_over_Bax[data_Ratios$Group == groups[1] | data_Ratios$Group == groups[2]] ~
                    as.factor(data_Ratios$Group[data_Ratios$Group == groups[1] | data_Ratios$Group == groups[2]])))
}

##################################################### DNA damage

#mtDNA cn
statData=data_corrected[data_corrected$Gene== 'shMito',]
statData$Group = as.factor(statData$Group)
print('mtDNA copy number')
print(kruskal.test(data=statData, statData$mtCN ~ statData$Group))

#DNA damage
for(iG in GOI){
  statData=data_corrected[data_corrected$Gene== iG,]
  statData$Group = as.factor(statData$Group)
  print(paste('Lesion Frequency per 10kb', iG, sep=': '))
  print(kruskal.test(data=statData, statData$LesionFreq ~ statData$Group))
}



