deltadeltaCq = function(tmp_data_delta, GOI, HKG, Control, Group){
# 2^-[(Cq_GOI - Cq_HKG)_treat - (Cq_GOI - Cq_HKG)_ctrl]
  
  
  ##################################################
  
  for(iPl in plates){
    for (iA in unique(tmp_data_delta$Animal[tmp_data_delta$Plate == iPl])){
      tmpGenes = unique(tmp_data_delta$Gene[tmp_data_delta$Plate == iPl & tmp_data_delta$Animal == iA])
      tmpGenes = tmpGenes[!(tmpGenes %in% HKG)]
      #calculate mean of Cq_HKG for all HKG
      tmp_Av_HKG = mean(tmp_data_delta$IC_corr_Cq[tmp_data_delta$Gene %in% HKG & tmp_data_delta$Animal == iA & tmp_data_delta$Plate == iPl])
      
      #calculate 2^Cq_GOI for all GOI
      for (iG in tmpGenes){
        tmp_delta1 = tmp_data_delta$IC_corr_Cq[tmp_data_delta$Animal == iA & tmp_data_delta$Gene == iG] - tmp_Av_HKG
        tmp_data_delta$DeltaDelta[tmp_data_delta$Animal == iA & tmp_data_delta$Gene == iG & tmp_data_delta$Plate == iPl] = tmp_delta1
      }
    }
  }
  
  tmp_data_delta_all = tmp_data_delta
  tmp_data_delta = tmp_data_delta[!(tmp_data_delta$Gene %in% HKG),]
  
  ################################################## 2. Normalization (to Control Group)
  #N = N_treatment / N_ctrl |  K: would going chancel out here
  # for derivation see: Scheefe et al. (2006): Quantitative real-time RT-PCR data analysis: current concepts and the novel "gene expression's CT difference" formula.
  # split tmp_data_delta into Control (IR) and Treatment (LR, HR)
  data_Ctrl = tmp_data_delta[tmp_data_delta$Group == Control,]
  # tmp_data_delta = tmp_data_delta[!(tmp_data_delta$Group == 'IR'),] #use this line if you want just information about treatment groups in the end
  
  tmp_Av_Ctrl = data.frame(matrix(ncol=2, nrow=length(GOI)))
  colnames(tmp_Av_Ctrl) = c('Gene', 'Ctrl_Av_N')
  
  #Average all HKG-normalized N values for Control Group for each GOI
  for (iG in 1:length(GOI)){
    tmp_Av_Ctrl$Gene[iG] = GOI[iG]
    tmp_Av_Ctrl$Ctrl_Av_N[iG] = mean(data_Ctrl$DeltaDelta[data_Ctrl$Gene == GOI[iG]])
  }
  
  #Relative quantity of GOI of Treatment Groups compared to Control Group
  for(iX in 1:nrow(tmp_data_delta)){
    tmpGene = tmp_data_delta$Gene[iX]
    tmp_data_delta$DeltaDelta[iX] = 2^-(tmp_data_delta$DeltaDelta[iX] - (tmp_Av_Ctrl$Ctrl_Av_N[tmp_Av_Ctrl$Gene == tmpGene]))
  }
  
  
  return(tmp_data_delta$DeltaDelta) 
}