
Sys.setenv(JAVA_HOME= 'C:\\Users\\lena_\\Downloads\\Matlab 2017a\\_temp_matlab_R2017a_win64\\sys\\java\\jre\\win64\\jre')
library(rJava)
library(xlsx)

# source('C:/Users/Gschossmann/Dropbox/studies/Osnabrück/Universität/Bachelorarbeit_DroBo/experiments/Analysis Code/PCR/qPCR_GeneExpressionAssay/getCORT.R')
source('C:/Users/lena_/Dropbox/studies/Osnabrück/Universität/Bachelorarbeit_DroBo/experiments/Analysis Code/CORT/getCORT.R')

fileX = 'Cohort VR/Data Males VR Project 1-Metabolism.csv'

setwd('C:/Users/lena_/Dropbox/studies/Osnabrück/Universität/Bachelorarbeit_DroBo/experiments/Analysis Code/Behaviours')
filepath = paste('C:/Users/lena_/Dropbox/studies/Osnabrück/Universität/Bachelorarbeit_DroBo/experiments/Behavioural tests', fileX, sep='/')

raw_data = read.csv(filepath)
raw_data = raw_data[!(is.na(raw_data$ID)),]
mouseLines = c('HR', 'IR', 'LR')

####### Add CORT data
fileX = 'Cohort VR/Data Males Project 1-CORT.csv'
filepathCORT = paste('C:/Users/lena_/Dropbox/studies/Osnabrück/Universität/Bachelorarbeit_DroBo/experiments/Behavioural tests', fileX, sep='/')
cort = getCORT(filepathCORT)
for(iR in 1:nrow(raw_data)){
  raw_data$CORT[iR] = cort$CORT[cort$ID == raw_data$ID[iR]]
}

descriptStats = data.frame(matrix(ncol = length(raw_data), nrow= 3*5))
names = variable.names(raw_data[1:length(raw_data)])
names[1] = 'Line'
names[2] = 'DescriptiveParam'
colnames(descriptStats) = names

cleanedStats = data.frame(matrix(ncol = length(raw_data), nrow= 3*5))
names = variable.names(raw_data[1:length(raw_data)])
names[1] = 'Line'
names[2] = 'DescriptiveParam'
colnames(cleanedStats) = names


#simles descriptive stats table
count = c(1:5)
for (idx in 1:3){
  tmp_data = raw_data[which(raw_data$Line == mouseLines[idx]),]
  #mean
  descriptStats[count[1],3:length(descriptStats)] = apply(tmp_data[,3:length(tmp_data)], 2, function(x) round(mean(x, na.rm = TRUE), digits=3))
  #standard deviation
  descriptStats[count[2],3:length(descriptStats)] = apply(tmp_data[,3:length(tmp_data)], 2, function(x) round(sd(x, na.rm = TRUE), digits=3))
  #SEM
  descriptStats[count[3],3:length(descriptStats)] = apply(tmp_data[,3:length(tmp_data)], 2, function(x) round(sd(x, na.rm = TRUE)/sqrt(sum(!is.na(x))), digits=3))
  #AV+sd
  descriptStats[count[4],3:length(descriptStats)] = apply(tmp_data[,3:length(tmp_data)], 2, function(x) round((mean(x, na.rm = TRUE)+ (2*sd(x, na.rm = TRUE))), digits=3))
  #AV-sd
  descriptStats[count[5],3:length(descriptStats)] = apply(tmp_data[,3:length(tmp_data)], 2, function(x) round((mean(x, na.rm = TRUE)-(2*sd(x, na.rm = TRUE))), digits=3))
  
  
  descriptStats[count[1:5], 2] = c('Average', 'Standard_Deviation', 'SEM', 'AV+sd', 'AV-sd')
  descriptStats[count[1:5],1] = rep(mouseLines[idx], 5)
  count = c((count[5]+1):(count[5]+5))
}


#Outlier detection after Dixon
#specify variables that are essential -> when outlier herein, animal should be removed (completely/from test)
exclusion_criteria = c('AV_BW', 'CORT')
FST_exclusion = c('Dur_floating')
OF_exclusion = c('thigmotaxis')

FST = 18:22 #Variable Indices corresponding to FST
OF = 23:31 #Variable Indices corresponding to OF

outlierTable = raw_data
outlierTable[,3:length(outlierTable)]=0

raw_cleaned = raw_data

elim = 0
elimFST = 0
elimOF = 0

#outer loop -> animals
for (idxA in 1:nrow(raw_data)){
  
  LineID = raw_data$Line[idxA]
  tmp_stats_AVplus = descriptStats[descriptStats$Line == LineID & descriptStats$DescriptiveParam == 'AV+sd',]
  tmp_stats_AVminus = descriptStats[descriptStats$Line == LineID & descriptStats$DescriptiveParam == 'AV-sd',]
  
  
  #inner loop -> variables
  for (idxV in 3:(length(raw_data))){
    isOutlier = 0
    #if Value is NA -> ignore
    if(is.na(raw_data[idxA, idxV]) == 1){
      break
    }
    #check for each variable if animal shows values outside of mean+/- 2*sd
    if(raw_data[idxA, idxV] > tmp_stats_AVplus[idxV] | raw_data[idxA, idxV] < tmp_stats_AVminus[idxV]){
      raw_cleaned[idxA, idxV] = NA
      outlierTable[idxA, idxV] = 1
      isOutlier = 1
    }
    #check if animal needs to be removed completely
    if(variable.names(raw_data)[idxV] %in% exclusion_criteria & isOutlier ==1){
      # raw_cleaned[idxA:(nrow(raw_cleaned)-1),] = raw_cleaned[(idxA+1):nrow(raw_cleaned), ]
      # raw_cleaned = raw_cleaned[1:(nrow(raw_cleaned)-1),]
      raw_cleaned[idxA,] = NA
      elim = elim+1
      break
    }
    #check if animal needs to be removed from FST
    if(variable.names(raw_data)[idxV] %in% FST_exclusion  & isOutlier ==1){
      raw_cleaned[idxA, FST] = NA
      elimFST = elimFST+1
    }
    #check if animal needs to be removed from OF
    if(variable.names(raw_data)[idxV] %in% OF_exclusion  & isOutlier ==1){
      raw_cleaned[idxA, OF] = NA
      elimOF = elimOF+1
    }
  }
}

raw_cleaned = raw_cleaned[!is.na(raw_cleaned$ID),]
# # Export raw_cleaned as excel
# write.xlsx(raw_cleaned, 'cleaned_data.xlsx')

# nHR = sum(raw_cleaned$Line == 'HR')
# nIR = sum(raw_cleaned$Line == 'IR')
# nLR = sum(raw_cleaned$Line == 'LR')


#### Run descriptive Stats again

count = c(1:5)
for (idx in 1:3){
  tmp_data = raw_cleaned[which(raw_cleaned$Line == mouseLines[idx]),]
  #mean
  cleanedStats[count[1],3:length(descriptStats)] = apply(tmp_data[,3:length(tmp_data)], 2, function(x) round(mean(x, na.rm = TRUE), digits=3))
  #standard deviation
  cleanedStats[count[2],3:length(descriptStats)] = apply(tmp_data[,3:length(tmp_data)], 2, function(x) round(sd(x, na.rm = TRUE), digits=3))
  #SEM
  cleanedStats[count[3],3:length(descriptStats)] = apply(tmp_data[,3:length(tmp_data)], 2, function(x) round(sd(x, na.rm = TRUE)/sqrt(sum(!is.na(x))), digits=3))
  #AV+sd
  cleanedStats[count[4],3:length(descriptStats)] = apply(tmp_data[,3:length(tmp_data)], 2, function(x) round((mean(x, na.rm = TRUE)+ 2*sd(x, na.rm = TRUE)), digits=3))
  #AV-sd
  cleanedStats[count[5],3:length(descriptStats)] = apply(tmp_data[,3:length(tmp_data)], 2, function(x) round((mean(x, na.rm = TRUE)-2*sd(x, na.rm = TRUE)), digits=3))
  
  
  cleanedStats[count[1:5], 2] = c('Average', 'Standard_Deviation', 'SEM', 'AV+sd', 'AV-sd')
  cleanedStats[count[1:5],1] = rep(mouseLines[idx], 5)
  count = c((count[5]+1):(count[5]+5))
}

# Create dataframe without NA for each test

# FST_cleaned = raw_cleaned[!is.na(raw_cleaned$Dur_struggling), 18:22]
#   data.frame(matrix(ncol=((18:22)+2), nrow = length(raw_cleaned$Dur_struggling[!is.na(raw_cleaned$Dur_struggling)]))

# Export Means+SEM table
write.xlsx(cleanedStats, paste('C:/Users/lena_/Dropbox/studies/Osnabrück/Universität/Bachelorarbeit_DroBo/experiments/Behavioural tests/Cohort VR/analysis files',
                               'Cohort_VR_Cleaned_DescriptiveStats.xlsx',sep='/'))

