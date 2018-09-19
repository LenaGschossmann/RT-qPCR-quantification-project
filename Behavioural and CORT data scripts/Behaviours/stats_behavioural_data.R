
Sys.setenv(JAVA_HOME= 'C:\\Users\\lena_\\Downloads\\Matlab 2017a\\_temp_matlab_R2017a_win64\\sys\\java\\jre\\win64\\jre')
library(rJava)
library(xlsx)
library(car)

##### Stat Tests

# library(Rmisc)

setwd('C:/Users/lena_/Dropbox/studies/Osnabrück/Universität/Bachelorarbeit_DroBo/experiments/Analysis Code/Behaviours')

#for SPT
# statData = Subgroup_anhedonic[,c(1,2,5,8,11,12,13,14,15)]

statData = raw_cleaned
colnames(statData)[2]= 'Line'
# statData=statData[statData$Mouse_Subline==1,]
# statData=statData[statData$Mouse_Subline==2,]
statData$Line[statData$Line == 1] ='HR'
statData$Line[statData$Line == 2] ='IR'
statData$Line[statData$Line == 3] ='LR'

#Test for homogeneity of variance
for(idx in 3:ncol(statData)){
  print(colnames(statData[idx]))
  print(leveneTest(statData[,idx] ~ statData$Line, data=statData))
}

# Test for normaldistribution of data
for(idx in 3:ncol(statData)){
  print(colnames(statData[idx]))
  print(shapiro.test(statData[,idx]))
}

#Kruskal Wallis
# sink('Kruskal-Wallis-Test.txt')
for (idx in 3:ncol(statData)){
  statData[,idx] = as.factor(statData[,idx])
  varNames = variable.names(statData)
  print(varNames[idx])
  print(kruskal.test(data=statData, statData[,idx] ~ as.factor(statData$Line)))
}
# sink()

#Man-Whitney-U Test
# sink('Man-Whitney-U-Test.txt')
groupList = list(c('HR', 'IR'), c('IR', 'LR'), c('HR', 'LR'))
# groupList = list(c(1, 2), c(2, 3), c(1, 3))
List_pVals = data.frame(matrix(ncol=2, nrow=length(groupList)))
colnames(List_pVals) = c('Groups', 'p.value')

for (idx in 3:length(statData)){
  varNames = variable.names(statData)
  for(iGr in 1:length(groupList)){
    groups = unlist(groupList[iGr])
    tmpData= statData[statData$Line %in% groups,]
    # tmpData$Line = as.factor(tmpData$Line)
    print(paste(varNames[idx], ':', groups[1], 'vs', groups[2], sep=' '))
    print(wilcox.test(data=tmpData, tmpData[,idx] ~ tmpData$Line))
    tmp_pVal = wilcox.test(data=tmpData, tmpData[,idx] ~ tmpData$Line)
    List_pVals$Groups[iGr] = paste(groups[1], groups[2], sep='_')
    List_pVals$p.value[iGr] = tmp_pVal$p.value
  }
  #+ Bonferroni correction: most sig. p-value*number of groups | sec. most significant p-values * (number of groups - 1) etc.
  n = length(groupList)
  List_pVals = List_pVals[order(List_pVals$p.value, decreasing = FALSE),]
  
  for(idx in 1:nrow(List_pVals)){
    List_pVals$p.value[idx] = List_pVals$p.value[idx] * n
    n=n-1
  }
  print('Bonferroni corrected p-values')
  print(List_pVals)
}
# sink()

####### SPT against 65%
for(iGr in unique(statData$Line)){
  tmpData = statData[statData$Line == iGr,]
  print(iGr)
  print(t.test(tmpData$per_AV, mu =65))
}

###### SPT test subgroups
statData = Subgroup_anhedonic
statData2 = Subgroup_hedonic
colnames(statData)[2]= 'Line'
statData$Line[statData$Line == 1] ='HR'
statData$Line[statData$Line == 2] ='IR'
statData$Line[statData$Line == 3] ='LR'
colnames(statData2) = colnames(statData)
statData2$Line[statData2$Line == 1] ='HR'
statData2$Line[statData2$Line == 2] ='IR'
statData2$Line[statData2$Line == 3] ='LR'

for(iGr in unique(statData$Line)){
  tmpData = statData[statData$Line == iGr,]
  tmpData2 = statData2[statData2$Line == iGr,]
  print(iGr)
  print(t.test(tmpData$per_AV, tmpData2$per_AV))
  # print(t.test(tmpData$per_AV, mu =65))
}

##### Av rel intake within groups between habituation, 2.5% forced and SPT
# 2.5% rel total intake vs rel av total intake during SPT
for(iGr in unique(statData$Line)){
  tmp_statData = statData[statData$Line == iGr,]
  tmpD = rbind(tmp_statData, tmp_statData)
  tmpD$IDX[1:nrow(tmp_statData)] = '2.5_Sucrose'
  tmpD$ColOFinterest[1:nrow(tmp_statData)] = tmp_statData$Rel_Suc
  tmpD$IDX[(nrow(tmp_statData)+1):nrow(tmpD)] = 'Av.Rel.Suc'
  tmpD$ColOFinterest[(nrow(tmp_statData)+1):nrow(tmpD)] = tmp_statData$Av_rel_TOT
  print(iGr)
  print(wilcox.test(data=tmpD,ColOFinterest ~ IDX))
}


# 2.5% rel total intake vs rel av total intake during habituation
for(iGr in unique(statData$Line)){
  tmp_statData = statData[statData$Line == iGr,]
  tmpD = rbind(tmp_statData, tmp_statData)
  tmpD$IDX[1:nrow(tmp_statData)] = '2.5_Sucrose'
  tmpD$ColOFinterest[1:nrow(tmp_statData)] = tmp_statData$Rel_Suc
  tmpD$IDX[(nrow(tmp_statData)+1):nrow(tmpD)] = 'Av_Total_Habituation_without_Suc'
  tmpD$ColOFinterest[(nrow(tmp_statData)+1):nrow(tmpD)] = tmp_statData$Av_Rel_H_noSuc
  print(iGr)
  print(wilcox.test(data=tmpD,ColOFinterest ~ IDX))
}

# SPT rel total intake vs rel av total intake during habituation
for(iGr in unique(statData$Line)){
  tmp_statData = statData[statData$Line == iGr,]
  tmpD = rbind(tmp_statData, tmp_statData)
  tmpD$IDX[1:nrow(tmp_statData)] = 'AvRelTotalSucrose'
  tmpD$ColOFinterest[1:nrow(tmp_statData)] = tmp_statData$Av_rel_TOT
  tmpD$IDX[(nrow(tmp_statData)+1):nrow(tmpD)] = 'Av_Total_Habituation_without_Suc'
  tmpD$ColOFinterest[(nrow(tmp_statData)+1):nrow(tmpD)] = tmp_statData$Av_Rel_H_noSuc
  print(iGr)
  print(wilcox.test(data=tmpD,ColOFinterest ~ IDX))
}


############################################# Print n per group

for(iTr in unique(statData$Line)){
  print(sprintf('Group; %s | n = %i',iTr, nrow(statData[statData$Line == iTr,])))
}
