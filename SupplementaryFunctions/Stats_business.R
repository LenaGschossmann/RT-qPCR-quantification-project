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


##################################################### DNA damage

#Test for homogeneity of variance
for(iG in unique(data_corrected$Gene)){
  tmpData = data_corrected[data_corrected$Gene== iG,]
  print(iG)
  print(leveneTest(IC_corr_Cq ~ Group, data=tmpData))
}

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

# Lesion Frequency per gene
statData=data_corrected[!(is.na(data_corrected$LesionFreq)),]
print(wilcox.test(data=statData, statData$LesionFreq ~ statData$Gene))


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


