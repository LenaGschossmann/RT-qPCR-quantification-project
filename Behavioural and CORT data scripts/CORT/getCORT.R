### CORT analysis
# fileX = 'Cohort VR/Data Males Project 1-CORT.csv'
# filepathCORT = paste('C:/Users/lena_/Dropbox/studies/Osnabrück/Universität/Bachelorarbeit_DroBo/experiments/Behavioural tests', fileX, sep='/')

getCORT = function(filepathCORT){
  raw_CORT= read.csv(filepathCORT)
  
  colnames(raw_CORT)[ncol(raw_CORT)] = 'CORT'
  colnames(raw_CORT)[which(colnames(raw_CORT) == 'Line')] = 'Group'
  
  raw_CORT$Group[raw_CORT$Group == 1] = 'HR'
  raw_CORT$Group[raw_CORT$Group == 2] = 'IR'
  raw_CORT$Group[raw_CORT$Group == 3] = 'LR'
  mouseGroups = c('HR', 'IR', 'LR')
  
  descriptStats = data.frame(matrix(ncol = 4, nrow= 3*5))
  colnames(descriptStats) = c('Group', 'SubGroup', 'DescriptiveParam', 'Value')
  
  cleanedStatsCORT = data.frame(matrix(ncol = 4, nrow= 3*5))
  colnames(cleanedStatsCORT) = c('Group', 'SubGroup', 'DescriptiveParam', 'Value')
  
  #simles descriptive stats table
  count = c(1:5)
  for (idx in 1:length(mouseGroups)){
    tmp_data = raw_CORT[which(raw_CORT$Group == mouseGroups[idx]),]
    #mean
    descriptStats$Value[count[1]] = round(mean(tmp_data$CORT,na.rm=TRUE), digits=3)
    #standard deviation
    descriptStats$Value[count[2]] = round(sd(tmp_data$CORT, na.rm = TRUE), digits=3)
    #SEM
    descriptStats$Value[count[3]] = round(sd(tmp_data$CORT, na.rm = TRUE)/sqrt(sum(!is.na(tmp_data$CORT))), digits=3)
    #AV+sd
    descriptStats$Value[count[4]] = round((mean(tmp_data$CORT, na.rm = TRUE)+ (2*sd(tmp_data$CORT, na.rm = TRUE))), digits=3)
    #AV-sd
    descriptStats$Value[count[5]] = round((mean(tmp_data$CORT, na.rm = TRUE)-(2*sd(tmp_data$CORT, na.rm = TRUE))), digits=3)
    
    descriptStats[count[1:5], 3] = c('Average', 'Standard_Deviation', 'SEM', 'AV+sd', 'AV-sd')
    descriptStats[count[1:5],1] = rep(mouseGroups[idx], 5)
    count = c((count[5]+1):(count[5]+5))
  }
  
  
  #Outlier detection after Dixon
  for (idxA in 1:nrow(raw_CORT)){
    GroupID = raw_CORT$Group[idxA]
    tmp_stats_AVplus = descriptStats$Value[descriptStats$Group == GroupID & descriptStats$DescriptiveParam == 'AV+sd']
    tmp_stats_AVminus = descriptStats$Value[descriptStats$Group == GroupID & descriptStats$DescriptiveParam == 'AV-sd']
    
    if(is.na(raw_CORT$CORT[idxA]) == 1){
      break
    }
    #check for each variable if animal shows values outside of mean+/- 2*sd
    if(raw_CORT$CORT[idxA] > tmp_stats_AVplus | raw_CORT$CORT[idxA] < tmp_stats_AVminus){
      raw_CORT$CORT[idxA]=NA
    }
  }
  
  cleaned_CORT = raw_CORT[!(is.na(raw_CORT$CORT)),]
  # # Export raw_cleaned as excel
  # write.xlsx(raw_cleaned, 'cleaned_CORT.xlsx')
  # 
  # nHR = sum(raw_cleaned$Group == 'HR')
  # nIR = sum(raw_cleaned$Group == 'IR')
  # nLR = sum(raw_cleaned$Group == 'LR')
  
  
  #### Run descriptive Stats again
  count = c(1:5)
  for (idx in 1:length(mouseGroups)){
    tmp_data = raw_CORT[which(raw_CORT$Group == mouseGroups[idx]),]
    #mean
    cleanedStatsCORT$Value[count[1]] = round(mean(tmp_data$CORT,na.rm=TRUE), digits=3)
    #standard deviation
    cleanedStatsCORT$Value[count[2]] = round(sd(tmp_data$CORT, na.rm = TRUE), digits=3)
    #SEM
    cleanedStatsCORT$Value[count[3]] = round(sd(tmp_data$CORT, na.rm = TRUE)/sqrt(sum(!is.na(tmp_data$CORT))), digits=3)
    #AV+sd
    cleanedStatsCORT$Value[count[4]] = round((mean(tmp_data$CORT, na.rm = TRUE)+ (2*sd(tmp_data$CORT, na.rm = TRUE))), digits=3)
    #AV-sd
    cleanedStatsCORT$Value[count[5]] = round((mean(tmp_data$CORT, na.rm = TRUE)-(2*sd(tmp_data$CORT, na.rm = TRUE))), digits=3)
    
    
    cleanedStatsCORT[count[1:5], 2] = c('Average', 'Standard_Deviation', 'SEM', 'AV+sd', 'AV-sd')
    cleanedStatsCORT[count[1:5],1] = rep(mouseGroups[idx], 5)
    count = c((count[5]+1):(count[5]+5))
  }

  return(cleaned_CORT)
}
