
Sys.setenv(JAVA_HOME= 'C:\\Users\\lena_\\Downloads\\Matlab 2017a\\_temp_matlab_R2017a_win64\\sys\\java\\jre\\win64\\jre')
library(rJava)
library(xlsx)

##### Cleaning Behavioural data

setwd('C:/Users/lena_/Dropbox/studies/Osnabrück/Universität/Bachelorarbeit_DroBo/experiments/Analysis Code/Behaviours')
filepath = 'C:/Users/lena_/Dropbox/studies/Osnabrück/Universität/Bachelorarbeit_DroBo/experiments/Behavioural tests'

############## BW, Food, Water + Tissue
#BW, food, water
fileX = 'BW_Food_Water_cleaned.csv'
raw_data = read.csv(paste(filepath, fileX, sep='/'))
raw_data = raw_data[!(is.na(raw_data$Mouse_Line)),]
raw_data = raw_data[,which(variable.names(raw_data) %in% c('Mouse_original_ID', 'Mouse_Line', 'Mouse_SuBline', 'AV_BW', 'AV_rel_food', 'AV_rel_water', 'AV_food', 'AV_water'))]
exclusion_criteria = c('AV_BW', 'AV_rel_food', 'AV_rel_water')

#Tissue
fileX1 = 'Tissue_CohortI-IV.csv'
raw_tissue = read.csv(paste(filepath, fileX1, sep='/'))
raw_tissue = raw_tissue[!(is.na(raw_tissue$Line)),]
raw_tissue = raw_tissue[,which(variable.names(raw_tissue) %in% c('Mouse_original_ID', 'Rel_eWAT', 'Rel_muscle', 'Lean_body.mass'))]

# join dataframes
for(iR in 1:nrow(raw_data)){
  if(raw_data$Mouse_original_ID[iR] %in% unique(raw_tissue$Mouse_original_ID)){
    raw_data$Rel_muscle[iR] = raw_tissue$Rel_muscle[raw_tissue$Mouse_original_ID == raw_data$Mouse_original_ID[iR]]
    raw_data$Lean_body.mass[iR] = raw_tissue$Lean_body.mass[raw_tissue$Mouse_original_ID == raw_data$Mouse_original_ID[iR]]
    raw_data$Rel_eWAT[iR] = raw_tissue$Rel_eWAT[raw_tissue$Mouse_original_ID == raw_data$Mouse_original_ID[iR]]
  }else{}
}

############## SPT
#save parameters from BW, food, water
raw_BW_F_W = raw_data

fileX = 'SPT_cleaned.csv'
raw_data = read.csv(paste(filepath, fileX, sep='/'))

raw_data = raw_data[!(is.na(raw_data$AV__BW)),]
raw_data=raw_data[,2:ncol(raw_data)]

exclusion_criteria = c('per_AV', 'D1_TOT', 'D2_TOT', 'D3_TOT')

raw_data$Av_Rel_H_noSuc = mean(c(raw_data$D1_TOT,raw_data$D2_TOT, raw_data$D3_TOT)) / raw_data$AV__BW

#### Import FST data
fileX = 'FST_CohortI-IV.csv'
raw_FST = read.csv(paste(filepath, fileX, sep='/'))

raw_FST = raw_FST[!(is.na(raw_FST$Mouse_Line)),]


##### join data frames
# from BW - Food - Water
for(iR in 1:nrow(raw_data)){
  if(raw_data$Mouse_original_ID[iR] %in% unique(raw_BW_F_W$Mouse_original_ID)){
    raw_data$Dur_floating[iR] = raw_BW_F_W$Dur_floating[raw_BW_F_W$Mouse_original_ID == raw_data$Mouse_original_ID[iR]]
  }else{}
}

# from FST
for(iR in 1:nrow(raw_data)){
  if(raw_data$Mouse_original_ID[iR] %in% unique(raw_FST$Mouse_original_ID)){
    raw_data$Dur_swimming[iR] = raw_FST$Dur_swimming[raw_FST$Mouse_original_ID == raw_data$Mouse_original_ID[iR]]
  }else{}
}
  

################################################### proceed here for SPT & BW, FOod, Water
mouseLines = c(1,2,3)

descriptStats = data.frame(matrix(ncol =ncol(raw_data), nrow= 3*5))
names1 = variable.names(raw_data)
names1[2] = 'DescriptiveParam'
names1[1] = 'Line'
colnames(descriptStats) = names1

cleanedStats = data.frame(matrix(ncol=ncol(raw_data), nrow= 3*5))
colnames(cleanedStats) = names1

#descriptive stats table
count = c(1:5)
for (idx in 1:3){
  tmp_data = raw_data[raw_data$Mouse_Line == mouseLines[idx],4:ncol(raw_data)]
  #mean
  descriptStats[count[1],4:length(descriptStats)] = apply(tmp_data[,1:length(tmp_data)], 2, function(x) round(mean(x, na.rm = TRUE), digits=3))
  #standard deviation
  descriptStats[count[2],4:length(descriptStats)] = apply(tmp_data[,1:length(tmp_data)], 2, function(x) round(sd(x, na.rm = TRUE), digits=3))
  #SEM
  descriptStats[count[3],4:length(descriptStats)] = apply(tmp_data[,1:length(tmp_data)], 2, function(x) round(sd(x, na.rm = TRUE)/sqrt(sum(!is.na(x))), digits=3))
  #AV+sd
  descriptStats[count[4],4:length(descriptStats)] = apply(tmp_data[,1:length(tmp_data)], 2, function(x) round((mean(x, na.rm = TRUE)+ (2*sd(x, na.rm = TRUE))), digits=3))
  #AV-sd
  descriptStats[count[5],4:length(descriptStats)] = apply(tmp_data[,1:length(tmp_data)], 2, function(x) round((mean(x, na.rm = TRUE)-(2*sd(x, na.rm = TRUE))), digits=3))
  
  descriptStats[count[1:5],1] = rep(mouseLines[idx], 5)
  descriptStats[count[1:5], 2] = c('Average', 'Standard_Deviation', 'SEM', 'AV+sd', 'AV-sd')
  count = c((count[5]+1):(count[5]+5))
}

raw_cleaned = raw_data
raw_cleaned = raw_cleaned[!(is.na(raw_cleaned$Mouse_Line)),]
# #Outlier detection after Dixon
# #specify variables that are essential -> when outlier herein, animal should be removed (completely/from test)
# outlierTable = raw_data
# outlierTable[,3:length(outlierTable)]=0
# 
# 
# elim = 0
# 
# #outer loop -> animals
# for (idxA in 1:nrow(raw_data)){
#   LineID = raw_data$Mouse_Line[idxA]
#   tmp_stats_AVplus = descriptStats[descriptStats$Line == LineID & descriptStats$DescriptiveParam == 'AV+sd',]
#   tmp_stats_AVminus = descriptStats[descriptStats$Line == LineID & descriptStats$DescriptiveParam == 'AV-sd',]
#   #inner loop -> variables
#   for (idxV in 3:(ncol(raw_data))){
#     isOutlier = 0
#     #if Value is NA -> ignore
#     if(is.na(raw_data[idxA, idxV]) == 1){
#       break
#     }
#     #check for each variable if animal shows values outside of mean+/- 2*sd
#     if(raw_data[idxA, idxV] > tmp_stats_AVplus[idxV] | raw_data[idxA, idxV] < tmp_stats_AVminus[idxV]){
#       outlierTable[idxA, idxV] = 1
#       isOutlier = 1
#     }
#     #check if animal needs to be removed completely
#     if(variable.names(raw_data)[idxV] %in% exclusion_criteria & isOutlier ==1){
#       raw_cleaned[idxA:(nrow(raw_cleaned)-1),] = raw_cleaned[(idxA+1):nrow(raw_cleaned), ]
#       raw_cleaned = raw_cleaned[1:(nrow(raw_cleaned)-1),]
#       # raw_cleaned[idxA,] = NA
#       elim = elim+1
#       break
#     }
#   }
# }
# 
# raw_cleaned = na.omit(raw_cleaned)
# # # Export raw_cleaned as excel
# # write.xlsx(raw_cleaned, 'cleaned_data.xlsx')
# 
# # nHR = sum(raw_cleaned$Line == 'HR')
# # nIR = sum(raw_cleaned$Line == 'IR')
# # nLR = sum(raw_cleaned$Line == 'LR')
# 
# 
# #### Run descriptive Stats again
# 
# count = c(1:5)
# for (idx in 1:3){
#   tmp_data = raw_cleaned[raw_cleaned$Mouse_Line == mouseLines[idx],3:ncol(raw_data)]
#   #mean
#   cleanedStats[count[1],3:length(descriptStats)] = apply(tmp_data[,1:length(tmp_data)], 2, function(x) round(mean(x, na.rm = TRUE), digits=3))
#   #standard deviation
#   cleanedStats[count[2],3:length(descriptStats)] = apply(tmp_data[,1:length(tmp_data)], 2, function(x) round(sd(x, na.rm = TRUE), digits=3))
#   #SEM
#   cleanedStats[count[3],3:length(descriptStats)] = apply(tmp_data[,1:length(tmp_data)], 2, function(x) round(sd(x, na.rm = TRUE)/sqrt(sum(!is.na(x))), digits=3))
#   #AV+sd
#   cleanedStats[count[4],3:length(descriptStats)] = apply(tmp_data[,1:length(tmp_data)], 2, function(x) round((mean(x, na.rm = TRUE)+ 2*sd(x, na.rm = TRUE)), digits=3))
#   #AV-sd
#   cleanedStats[count[5],3:length(descriptStats)] = apply(tmp_data[,1:length(tmp_data)], 2, function(x) round((mean(x, na.rm = TRUE)-2*sd(x, na.rm = TRUE)), digits=3))
#   
#   
#   cleanedStats[count[1:5], 2] = c('Average', 'Standard_Deviation', 'SEM', 'AV+sd', 'AV-sd')
#   cleanedStats[count[1:5],1] = rep(mouseLines[idx], 5)
#   count = c((count[5]+1):(count[5]+5))
# }


######################## SPT

# # check if fluid intake during SPT differs by large from habituation phase (using DIxons test)
# raw_cleaned = raw_cleaned[,c(1:15,18,21,24,27)]
# 
# for(idx in 1:nrow(raw_cleaned)){
#   raw_cleaned$AV_H_TOT[idx] = mean(as.numeric(c(raw_cleaned$H4_TOT[idx],raw_cleaned$H3_TOT[idx], raw_cleaned$H2_TOT[idx], raw_cleaned$H1_TOT[idx] )))
#   raw_cleaned$mean_plus_2sd[idx] = raw_cleaned$AV_H_TOT[idx] + 3 * sd(as.numeric(c(raw_cleaned$H4_TOT[idx],raw_cleaned$H3_TOT[idx], raw_cleaned$H2_TOT[idx], raw_cleaned$H1_TOT[idx] )))
#   raw_cleaned$mean_minus_2sd[idx] = raw_cleaned$AV_H_TOT[idx] - 3 * sd(as.numeric(c(raw_cleaned$H4_TOT[idx],raw_cleaned$H3_TOT[idx], raw_cleaned$H2_TOT[idx], raw_cleaned$H1_TOT[idx] )))
#   
#   if(raw_cleaned$D1_TOT[idx] < raw_cleaned$mean_minus_2sd[idx] | raw_cleaned$D1_TOT[idx] > raw_cleaned$mean_plus_2sd[idx]){
#     raw_cleaned$D1_TOT[idx] = NA
#     raw_cleaned$per_D1[idx] = NA
#   }
#   if(raw_cleaned$D2_TOT[idx] < raw_cleaned$mean_minus_2sd[idx] | raw_cleaned$D2_TOT[idx] > raw_cleaned$mean_plus_2sd[idx]){
#     raw_cleaned$D2_TOT[idx] = NA
#     raw_cleaned$per_D2[idx] = NA
#   }
#   if(raw_cleaned$D3_TOT[idx] < raw_cleaned$mean_minus_2sd[idx] | raw_cleaned$D3_TOT[idx] > raw_cleaned$mean_plus_2sd[idx]){
#     raw_cleaned$D3_TOT[idx] = NA
#     raw_cleaned$per_D3[idx] = NA
#   }
# }
# 
# # mean total fluid intake during SPT
# for(idx in 1: nrow(raw_cleaned)){
#   raw_cleaned$AV_TOT[idx] = mean(c(raw_cleaned$D1_TOT[idx], raw_cleaned$D2_TOT[idx], raw_cleaned$D3_TOT[idx]))
# }
# 
# # mean total sucrose intake during SPT
# for(idx in 1: nrow(raw_cleaned)){
#   raw_cleaned$Suc_TOT[idx] = mean(c(raw_cleaned$S1[idx], raw_cleaned$S2[idx], raw_cleaned$S3[idx]))
# }
# 
# raw_cleaned = na.omit(raw_cleaned)



#### SPT: check if there are two distinct subgroups regarding anhedonia liking
raw_cleaned$PrefGroup[raw_cleaned$per_AV >= 65] = 1
raw_cleaned$PrefGroup[raw_cleaned$per_AV < 65] = 0

Subgroup_hedonic = raw_cleaned[raw_cleaned$PrefGroup == 1 & !(is.na(raw_cleaned$Mouse_original_ID)),]
Subgroup_anhedonic = raw_cleaned[raw_cleaned$PrefGroup == 0  & !(is.na(raw_cleaned$Mouse_original_ID)),]





####################################################
# Export Means+SEM table
write.xlsx(descriptStats, paste('C:/Users/lena_/Dropbox/studies/Osnabrück/Universität/Bachelorarbeit_DroBo/experiments/Behavioural tests/CohortsI-IV/analysis files',
                               'CohortsI-IV_SPT_Cleaned_DescriptiveStats.xlsx',sep='/'))





