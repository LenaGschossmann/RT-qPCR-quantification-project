
#need to run qPCR_quantification.R first!!

useDelta = 1 #use results obtained by the DeltaDelta method
ini_data = data_corrected

if (useDelta == 1){
  savename1 = 'apoptosis_ratios_DeltaDelta.xlsx'
  savename2 = 'apoptosis_summary_DeltaDelta.xlsx'
  data_Ratios = data.frame(matrix(ncol = 4, nrow=nrow(ini_data)))
  colnames(data_Ratios) = c('Plate', 'Animal', 'Group', 'Bcl2_over_Bax')
  tmpGenes = c('Bax', 'Bcl2')
  cntR = 0
  for(iA in unique(ini_data$Animal)){
    if(sum(tmpGenes %in% ini_data$Gene[ini_data$Animal == iA]) == 2 ){
      cntR = cntR+1
      tmpBax = ini_data$DeltaDelta[ini_data$Animal == iA & ini_data$Gene == 'Bax']
      tmpBcl2 = ini_data$DeltaDelta[ini_data$Animal == iA & ini_data$Gene == 'Bcl2']
      data_Ratios$Animal[cntR] = iA
      data_Ratios$Plate[cntR] = ini_data$Plate[ini_data$Animal == iA & ini_data$Gene == 'Bax']
      data_Ratios$Group[cntR] = ini_data$Group[ini_data$Animal == iA & ini_data$Gene == 'Bax']
      data_Ratios$Bcl2_over_Bax[cntR] = tmpBcl2 / tmpBax
    }
  }
}else{
  savename1 = 'apoptosis_ratios.xlsx'
  savename2 = 'apoptosis_summary.xlsx'
  data_Ratios = data.frame(matrix(ncol = 4, nrow=nrow(ini_data)))
  colnames(data_Ratios) = c('Plate', 'Animal', 'Group', 'Bcl2_over_Bax')
  tmpGenes = c('Bax', 'Bcl2')
  cntR = 0
  for(iA in unique(ini_data$Animal)){
    if(sum(tmpGenes %in% ini_data$Gene[ini_data$Animal == iA]) == 2 ){
      cntR = cntR+1
      tmpBax = ini_data$Rel_quantity[ini_data$Animal == iA & ini_data$Gene == 'Bax']
      tmpBcl2 = ini_data$Rel_quantity[ini_data$Animal == iA & ini_data$Gene == 'Bcl2']
      data_Ratios$Animal[cntR] = iA
      data_Ratios$Plate[cntR] = ini_data$Plate[ini_data$Animal == iA & ini_data$Gene == 'Bax']
      data_Ratios$Group[cntR] = ini_data$Group[ini_data$Animal == iA & ini_data$Gene == 'Bax']
      data_Ratios$Bcl2_over_Bax[cntR] = tmpBcl2 / tmpBax
    }
  }
}

data_Ratios = na.omit(data_Ratios)


#summary
data_summary = data.frame(matrix(ncol=4, nrow=3))
colnames(data_summary)=c('Group','Mean','SEM', 'SD')
cnt=1
for(iT in unique(data_corrected$Group)){
    tmpMatrix = data_Ratios$Group == iT
    data_summary$Group[cnt] = iT
    data_summary$Mean[cnt] = mean(data_Ratios$Bcl2_over_Bax[tmpMatrix])
    data_summary$SD[cnt] = sd(data_Ratios$Bcl2_over_Bax[tmpMatrix])
    data_summary$SEM[cnt] = sd(data_Ratios$Bcl2_over_Bax[tmpMatrix])/sqrt(sum(tmpMatrix))
    cnt=cnt+1
}



################################################### Output
#Save Corrected data
write.xlsx(data_Ratios, paste(savepath,paste(part, savename1, sep='_'),sep='/'))

write.xlsx(data_summary, paste(savepath,paste(part, savename2, sep='_'),sep='/'))

# write.csv(data_summaryDeltaDelta, paste(savepath,paste(part, 'data_summaryDeltaDelta.csv', sep='_'),sep='/'))



