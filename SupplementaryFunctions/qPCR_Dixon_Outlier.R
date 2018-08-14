
qPCR_Dixon_Outlier = function(dataAveraged, IC_corrected, genes, Treatments){
#Treatments: treatment and control groups
#genes: all genes, with RefGenes listed first and then GOI
#set IC_corrected to '1' if the IC corrected Cq values should be used, otherwise '0'

  if(IC_corrected == 1){
    dataAveraged$Cq = dataAveraged$IC_corr_Cq #overscribe IC_corr values into Cq column for further calculations
  } 
  
  #calculate mean and sd
  dataAv_prelimStats = data.frame(matrix(ncol=4, nrow=length(genes)*length(Treatments)))
  colnames(dataAv_prelimStats) = c('Treatment', 'Gene', 'Mean', 'SD')
  cnt=1
  for(iG in genes){
    for(iT in Treatments){
      dataAv_prelimStats$Gene[cnt] = iG
      dataAv_prelimStats$Treatment[cnt] = iT
      dataAv_prelimStats$Mean[cnt] = mean(dataAveraged$Cq[dataAveraged$Gene == iG & dataAveraged$Group == iT & is.na(dataAveraged$Cq) == FALSE])
      dataAv_prelimStats$SD[cnt] = sd(dataAveraged$Cq[dataAveraged$Gene == iG & dataAveraged$Group == iT & is.na(dataAveraged$Cq) == FALSE])
      cnt = cnt+1
    }
  }
  
  
  outlier = data.frame(matrix(ncol=7, nrow=nrow(data_tot)))
  colnames(outlier) = c('Plate', 'Animal', 'Group', 'Gene', 'Well', 'Cq', 'Included')
  iO=1
  for(iT in unique(dataAveraged$Group[!(dataAveraged$Group == '') & is.na(dataAveraged$Cq) == FALSE])){
    for(iPl in plates){
      for(iA in unique(dataAveraged$Animal[dataAveraged$Plate == iPl & dataAveraged$Group == iT & is.na(dataAveraged$Cq) == FALSE])){
        excludeAnimal = 0
        tmpGOI = unique(dataAveraged$Gene[dataAveraged$Plate == iPl & dataAveraged$Animal == iA & !(dataAveraged$Gene %in% HKG) & is.na(dataAveraged$Cq) == FALSE])
        for(iHKG in HKG){
          tmpMean = dataAv_prelimStats$Mean[dataAv_prelimStats$Gene == iHKG & dataAv_prelimStats$Treatment == iT]
          tmpSD = dataAv_prelimStats$SD[dataAv_prelimStats$Gene == iHKG & dataAv_prelimStats$Treatment == iT]
          if(abs(dataAveraged$Cq[dataAveraged$Animal == iA & dataAveraged$Gene == iHKG & dataAveraged$Plate == iPl & is.na(dataAveraged$Cq) == FALSE]-tmpMean) >= 2*tmpSD){ #Apply Dixons outlier criterion on averaged replicates
            excludeAnimal = 1
          }else{}
        }
        if(excludeAnimal == 1){
          tmpOut = dataAveraged$Animal == iA & dataAveraged$Plate == iPl & is.na(dataAveraged$Cq) == FALSE
          iOrange = c(iO:(iO + sum(tmpOut) -1))
          outlier$Plate[iOrange] = dataAveraged$Plate[tmpOut]
          outlier$Animal[iOrange] = as.character(dataAveraged$Animal[tmpOut])
          outlier$Group[iOrange] = as.character(dataAveraged$Group[tmpOut])
          outlier$Gene[iOrange] = as.character(dataAveraged$Gene[tmpOut])
          outlier$Cq[iOrange] = dataAveraged$Cq[tmpOut]
          outlier$Well[iOrange] = 'Averaged'
          outlier$Included[iOrange]= FALSE
          dataAveraged[tmpOut,] = NA
          iO = 1+iOrange[length(iOrange)]
        }else{
          for(iG in tmpGOI){
            tmpMean = dataAv_prelimStats$Mean[dataAv_prelimStats$Gene == iG & dataAv_prelimStats$Treatment == iT]
            tmpSD = dataAv_prelimStats$SD[dataAv_prelimStats$Gene == iG & dataAv_prelimStats$Treatment == iT]
            if(abs(dataAveraged$Cq[dataAveraged$Animal == iA & dataAveraged$Gene == iG & dataAveraged$Plate == iPl & is.na(dataAveraged$Cq) == FALSE]-tmpMean) >= 2*tmpSD){ #Apply Dixons outlier criterion on averaged replicates
              tmpOut = dataAveraged$Animal == iA & dataAveraged$Plate == iPl & dataAveraged$Gene == iG & is.na(dataAveraged$Cq) == FALSE
              outlier$Plate[iO] = dataAveraged$Plate[tmpOut]
              outlier$Animal[iO] = as.character(dataAveraged$Animal[tmpOut])
              outlier$Group[iO] = as.character(dataAveraged$Group[tmpOut])
              outlier$Gene[iO] = as.character(dataAveraged$Gene[tmpOut])
              outlier$Cq[iO] = dataAveraged$Cq[tmpOut]
              outlier$Well[iO] = 'Averaged'
              outlier$Included[iO]= FALSE
              dataAveraged[tmpOut,] = NA
              iO = iO+1
            }else{}
            
          }
        }
      }
    }
  }
  
  dataAveraged = na.omit(dataAveraged)
  outlier = na.omit(outlier)
  
  returnList = list(dataAveraged, outlier, dataAv_prelimStats)
  
 return(returnList) 
}