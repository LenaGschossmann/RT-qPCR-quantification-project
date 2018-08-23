
#need to run qPCR_quantification.R first!!

ini_data = data_corrected

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

data_Ratios = na.omit(data_Ratios)





