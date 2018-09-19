##### stats business
library(car)
library(Hmisc)
library(openxlsx)

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
List_pVals = data.frame(matrix(ncol=2, nrow=length(groupList)))
colnames(List_pVals) = c('Groups', 'p.value')
for (iG in 1:length(GOI)){
  for (iGr in 1:length(groupList)){
    groups = unlist(groupList[iGr])
    print(paste(groups[1], 'vs.', groups[2]))
    tmpData = data_corrected[data_corrected$Group %in% groups,]
    print(as.character(GOI[iG]))
    print(wilcox.test(data=tmpData, tmpData$Rel_quantity[tmpData$Gene == as.character(GOI[iG])] ~ as.factor(tmpData$Group[tmpData$Gene == as.character(GOI[iG])])))
    tmp_pVal = wilcox.test(data=tmpData, tmpData$Rel_quantity[tmpData$Gene == as.character(GOI[iG])] ~ as.factor(tmpData$Group[tmpData$Gene == as.character(GOI[iG])]))
    List_pVals$Groups[iGr] = paste(groups[1], groups[2], sep='_')
    List_pVals$p.value[iGr] = tmp_pVal$p.value
  }
  n = length(groupList)
  List_pVals = List_pVals[order(List_pVals$p.value, decreasing = FALSE),]
  for(idx in 1:nrow(List_pVals)){
    List_pVals$p.value[idx] = List_pVals$p.value[idx] * n
    n=n-1
  }
  print('Bonferroni corrected p-values')
  print(List_pVals)
}


for(iG in GOI){
  for(iTr in unique(data_corrected$Group)){
    print(sprintf('Gene: %s | Group; %s | n = %i',iG, iTr, nrow(data_corrected[data_corrected$Gene == iG & data_corrected$Group == iTr,])))
  }
}

###################################################### DeltaDelta Method
#Kruskal-Wallis
for (iG in 1:length(GOI)){
  varNames = GOI
  print(as.character(GOI[iG]))
  print(kruskal.test(data=data_corrected,data_corrected$DeltaDelta[data_corrected$Gene == as.character(GOI[iG])] ~ as.factor(data_corrected$Group[data_corrected$Gene == as.character(GOI[iG])])))
}

#Man-Whitney-U Test
groupList = list(c('HR', 'IR'), c('IR', 'LR'), c('HR', 'LR'))
List_pVals = data.frame(matrix(ncol=2, nrow=length(groupList)))
colnames(List_pVals) = c('Groups', 'p.value')
for (iG in 1:length(GOI)){
  for (iGr in 1:length(groupList)){
    groups = unlist(groupList[iGr])
    print(paste(groups[1], 'vs.', groups[2]))
    tmpData = data_corrected[data_corrected$Group %in% groups,]
    print(as.character(GOI[iG]))
    print(wilcox.test(data=tmpData, tmpData$DeltaDelta[tmpData$Gene == as.character(GOI[iG])] ~ as.factor(tmpData$Group[tmpData$Gene == as.character(GOI[iG])])))
    tmp_pVal = wilcox.test(data=tmpData, tmpData$DeltaDelta[tmpData$Gene == as.character(GOI[iG])] ~ as.factor(tmpData$Group[tmpData$Gene == as.character(GOI[iG])]))
    List_pVals$Groups[iGr] = paste(groups[1], groups[2], sep='_')
    List_pVals$p.value[iGr] = tmp_pVal$p.value
  }
  n = length(groupList)
  List_pVals = List_pVals[order(List_pVals$p.value, decreasing = FALSE),]
  for(idx in 1:nrow(List_pVals)){
    List_pVals$p.value[idx] = List_pVals$p.value[idx] * n
    n=n-1
  }
  print('Bonferroni corrected p-values')
  print(List_pVals)
}


############ Method Comparison
statData = data_corrected_tot

for(iG in unique(statData$Gene)){
  for(iGr in unique(statData$Group)){
    tmpData = statData[statData$Gene == iG & statData$Group == iGr,]
    tmpD = rbind(tmpData,tmpData)
    tmpD$IDX[1:nrow(tmpData)] = 'GED'
    tmpD$ColOFinterest[1:nrow(tmpData)]= tmpData$Rel_quantity
    tmpD$IDX[(nrow(tmpData)+1):nrow(tmpD)] = 'DeltaDelta'
    tmpD$ColOFinterest[(nrow(tmpData)+1):nrow(tmpD)]= tmpData$DeltaDelta
    print(paste(iG, iGr, sep=' | '))
    print(wilcox.test(statData$Rel_quantity[statData$Gene == iG & statData$Group == iGr], statData$DeltaDelta[statData$Gene == iG & statData$Group == iGr]))
  }
}




##################################################### Stability HKG
for (iG in 1:length(HKG)){
  varNames = HKG
  print(as.character(HKG[iG]))
  # print(kruskal.test(data=data_corrected,data_corrected$IC_corr_Cq[data_corrected$Gene == as.character(HKG[iG])] ~ as.factor(data_corrected$Group[data_corrected$Gene == as.character(HKG[iG])])))
  print(kruskal.test(data=Efficiencies,Efficiencies$Efficiency[Efficiencies$Gene == as.character(HKG[iG])] ~ as.factor(Efficiencies$Group[Efficiencies$Gene == as.character(HKG[iG])])))
}



##################################################### Apoptosis Ratio
#Kruskal-Wallis-Test
print(kruskal.test(data=data_Ratios,data_Ratios$Bcl2_over_Bax ~ as.factor(data_Ratios$Group)))

for(iGr in unique(data_Ratios$Group)){
  tmpData = data_Ratios[data_Ratios$Group == iGr,]
  print(iGr)
  print(t.test(tmpData$Bcl2_over_Bax, mu =1))
}

#Man-Whitney-U Test
groupList = list(c('HR', 'IR'), c('IR', 'LR'), c('HR', 'LR'))
List_pVals = data.frame(matrix(ncol=2, nrow=length(groupList)))
colnames(List_pVals) = c('Groups', 'p.value')
for (iGr in 1:length(groupList)){
  groups = unlist(groupList[iGr])
  print(paste(groups[1], 'vs.', groups[2]))
  tmp_data_Ratios = data_Ratios[data_Ratios$Group %in% groups,]
  print(wilcox.test(data=tmp_data_Ratios, tmp_data_Ratios$Bcl2_over_Bax ~ as.factor(tmp_data_Ratios$Group)))
  tmp_pVal = wilcox.test(data=tmp_data_Ratios, tmp_data_Ratios$Bcl2_over_Bax ~ as.factor(tmp_data_Ratios$Group))
  List_pVals$Groups[iGr] = paste(groups[1], groups[2], sep='_')
  List_pVals$p.value[iGr] = tmp_pVal$p.value
}
n = length(groupList)
List_pVals = List_pVals[order(List_pVals$p.value, decreasing = FALSE),]
for(idx in 1:nrow(List_pVals)){
  List_pVals$p.value[idx] = List_pVals$p.value[idx] * n
  n=n-1
}
print('Bonferroni corrected p-values')
print(List_pVals)

for(iTr in unique(data_Ratios$Group)){
    print(sprintf('Group; %s | n = %i',iTr, nrow(data_Ratios[data_Ratios$Group == iTr,])))
}


##################################################### DNA damage

#Test for homogeneity of variance
for(iG in unique(data_corrected$Gene)){
  tmpData = data_corrected[data_corrected$Gene== iG,]
  print(iG)
  print(leveneTest(IC_corr_Cq ~ Group, data=tmpData))
}

###############mtDNA cn
statData=data_corrected[data_corrected$Gene== 'shMito',]
statData$Group = as.factor(statData$Group)
print('mtDNA copy number')
print(kruskal.test(data=statData, statData$mtCN ~ statData$Group))

groupList = list(c('HR', 'IR'), c('IR', 'LR'), c('HR', 'LR'))
List_pVals = data.frame(matrix(ncol=2, nrow=length(groupList)))
colnames(List_pVals) = c('Groups', 'p.value')
for (iGr in 1:length(groupList)){
  groups = unlist(groupList[iGr])
  print(paste(groups[1], 'vs.', groups[2]))
  tmp_statData = statData[statData$Group %in% groups,]
  print(wilcox.test(data=tmp_statData, tmp_statData$mtCN ~ as.factor(tmp_statData$Group)))
  tmp_pVal = wilcox.test(data=tmp_statData, tmp_statData$mtCN ~ as.factor(tmp_statData$Group))
  List_pVals$Groups[iGr] = paste(groups[1], groups[2], sep='_')
  List_pVals$p.value[iGr] = tmp_pVal$p.value
}
n = length(groupList)
List_pVals = List_pVals[order(List_pVals$p.value, decreasing = FALSE),]
for(idx in 1:nrow(List_pVals)){
  List_pVals$p.value[idx] = List_pVals$p.value[idx] * n
  n=n-1
}
print('Bonferroni corrected p-values')
print(List_pVals)

for(iTr in unique(statData$Group)){
  print(sprintf('Group; %s | n = %i',iTr, nrow(statData[statData$Group == iTr,])))
}


#################DNA damage
for(iG in GOI){
  statData=data_corrected[data_corrected$Gene== iG,]
  statData$Group = as.factor(statData$Group)
  print(paste('Lesion Frequency per 10kb', iG, sep=': '))
  print(kruskal.test(data=statData, statData$LesionFreq ~ statData$Group))
}

groupList = list(c('HR', 'IR'), c('IR', 'LR'), c('HR', 'LR'))
List_pVals = data.frame(matrix(ncol=2, nrow=length(groupList)))
colnames(List_pVals) = c('Groups', 'p.value')
for(iG in GOI){
  statData=data_corrected[data_corrected$Gene== iG,]
  statData$Group = as.factor(statData$Group)
  for (iGr in 1:length(groupList)){
    groups = unlist(groupList[iGr])
    print(paste(groups[1], 'vs.', groups[2]))
    tmp_statData = statData[statData$Group %in% groups,]
    print(wilcox.test(data=tmp_statData, tmp_statData$LesionFreq ~ as.factor(tmp_statData$Group)))
    tmp_pVal = wilcox.test(data=tmp_statData, tmp_statData$LesionFreq ~ as.factor(tmp_statData$Group))
    List_pVals$Groups[iGr] = paste(groups[1], groups[2], sep='_')
    List_pVals$p.value[iGr] = tmp_pVal$p.value
  }
  n = length(groupList)
  List_pVals = List_pVals[order(List_pVals$p.value, decreasing = FALSE),]
  for(idx in 1:nrow(List_pVals)){
    List_pVals$p.value[idx] = List_pVals$p.value[idx] * n
    n=n-1
  }
  print('Bonferroni corrected p-values')
  print(List_pVals)
  
  for(iTr in unique(statData$Group)){
    print(sprintf('Group; %s | n = %i',iTr, nrow(statData[statData$Group == iTr,])))
  }
}



# Lesion Frequency per gene
statData=data_corrected[!(is.na(data_corrected$LesionFreq)),]
print(wilcox.test(data=statData, statData$LesionFreq ~ statData$Gene))

for(iG in c('shMito', 'D-Loop', 'Cox3')){
  for(iTr in unique(data_corrected$Group)){
    print(sprintf('Gene: %s | Group; %s | n = %i',iG, iTr, nrow(data_corrected[data_corrected$Gene == iG & data_corrected$Group == iTr,])))
  }
}


#################################################### CORT
statData = cort

#Test for homogeneity of variance
print(leveneTest(CORT ~ Group, data=statData))

##Kruskal Wallis
print(kruskal.test(data=statData,statData$CORT ~ as.factor(statData$Group)))

#Man-Whitney-U Test
groupList = list(c('HR', 'IR'), c('IR', 'LR'), c('HR', 'LR'))
List_pVals = data.frame(matrix(ncol=2, nrow=length(groupList)))
colnames(List_pVals) = c('Groups', 'p.value')
  for (iGr in 1:length(groupList)){
    groups = unlist(groupList[iGr])
    print(paste(groups[1], 'vs.', groups[2]))
    tmpData = statData[statData$Group %in% groups,]
    tmpVal = wilcox.test(data=tmpData, tmpData$CORT ~ as.factor(tmpData$Group))
    print(tmpVal)
    List_pVals$Groups[iGr] = paste(groups[1], groups[2], sep='_')
    List_pVals$p.value[iGr] = tmpVal$p.value
  }
  n = length(groupList)
  List_pVals = List_pVals[order(List_pVals$p.value, decreasing = FALSE),]
  for(idx in 1:nrow(List_pVals)){
    List_pVals$p.value[idx] = List_pVals$p.value[idx] * n
    n=n-1
  }
  print('Bonferroni corrected p-values')
  print(List_pVals)


