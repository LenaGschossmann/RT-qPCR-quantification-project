
qPCR_replAveraging = function(data_tot){
#data_tot mus contain following columns: 'Threshold', 'Plate', 'Well', 'Gene', 'Animal', 'Group', 'Cq'
#with replicated listed in different rows
  
  dataAveraged = data.frame(matrix(nrow=(nrow(data_tot)/numRepl), ncol=(ncol(data_tot)+1)))
  colnames(dataAveraged) = c(colnames(data_tot), 'CV')
  
  exclusions = data.frame(matrix(ncol=7, nrow=nrow(data_tot)))
  colnames(exclusions) = c('Plate', 'Animal', 'Group', 'Gene', 'Well', 'Cq', 'Included')
  
  cnt=1
  tmp_warning_CV =''
  iEx = 1 #for counting outlier
  #loop through plates, animals, genes, replicates
  for (iPl in plates){
    ExcludeRest = 0
    for (iA in unique(data_tot$Animal[which(data_tot$Plate == iPl)])){
      
      for (iG in unique(data_tot$Gene[which(data_tot$Plate == iPl & data_tot$Animal == iA)])){
        
        if(iG %in% exclusions$Gene[exclusions$Animal == iA & exclusions$Plate == iPl]){
        }else{
          tmp =0 #necessary to initialize tmp for append function
          tmpMatrix = which(data_tot$Plate == iPl & data_tot$Animal == iA & data_tot$Gene == iG)
          for (iR in 1:numRepl){
            tmp = append(data_tot$Cq[tmpMatrix][iR], tmp)
          }
          tmp = tmp[1:numRepl] #get rid of 0 in beginning
          #reduce replicates into one mean value and store in new dataframe
          # dataAveraged[cnt,1:(ncol(dataAveraged)-1)] = data_tot[tmpMatrix[1],]
          if((max(tmp)-min(tmp)) < acceptable_diff){ #exclude dirty replicates | criterion: difference between replicates
            dataAveraged$Threshold[cnt] = data_tot$Threshold[tmpMatrix[1]]
            dataAveraged$Plate[cnt] = data_tot$Plate[tmpMatrix[1]]
            dataAveraged$Well[cnt] = as.character(data_tot$Well[tmpMatrix[1]])
            dataAveraged$Gene[cnt] = as.character(data_tot$Gene[tmpMatrix[1]])
            dataAveraged$Group[cnt] = as.character(data_tot$Group[tmpMatrix[1]])
            dataAveraged$Animal[cnt] = as.character(data_tot$Animal[tmpMatrix[1]])
            dataAveraged$Cq[cnt] = mean(tmp)
            tmpCV = sd(tmp)/mean(tmp) 
            dataAveraged$CV[cnt] = tmpCV
            cnt=cnt+1
            
          } else {
            
            for(iR in 1:numRepl){ # write outlier in dataframe
              exclusions$Plate[iEx] = data_tot$Plate[tmpMatrix[iR]]
              exclusions$Animal[iEx] = as.character(data_tot$Animal[tmpMatrix[iR]])
              exclusions$Group[iEx] = as.character(data_tot$Group[tmpMatrix[iR]])
              exclusions$Gene[iEx] = as.character(data_tot$Gene[tmpMatrix[iR]])
              exclusions$Well[iEx] = as.character(data_tot$Well[tmpMatrix[iR]])
              exclusions$Cq[iEx] = data_tot$Cq[tmpMatrix[iR]]
              exclusions$Included[iEx]= FALSE
              iEx = iEx+1
            }
            if(iA == 'IC'){
              tmpAverage = mean(data_tot$Cq[data_tot$Gene == iG & data_tot$Plate != iPl & data_tot$Animal == 'IC'])
            }else{
              tmpAverage = mean(data_tot$Cq[data_tot$Gene == iG & data_tot$Animal != iA & data_tot$Animal != 'IC'])
            }
            print(sprintf("Plate: %d | Animal: %s | Gene: %s", iPl, iA, iG))
            print(paste('Value:',data_tot$Cq[tmpMatrix]))
            print(paste('Difference:', round(max(tmp)-min(tmp), digits=4)))
            print(paste('Average:', tmpAverage))
            ui = readline(prompt='Type Number of Replicate (1. -> <1>) you want to keep, <0> if you dont want to keep any or <a> if you want to keep all:')
            input_ok = 0
            while(input_ok == 0){
              if(ui == '0'){
                input_ok = 1
                if(iG %in% HKG){
                  ExcludeRest = 1
                }
              } else if(ui == 'a'){
                dataAveraged$Threshold[cnt] = data_tot$Threshold[tmpMatrix[1]]
                dataAveraged$Plate[cnt] = data_tot$Plate[tmpMatrix[1]]
                dataAveraged$Well[cnt] = as.character(data_tot$Well[tmpMatrix[1]])
                dataAveraged$Gene[cnt] = as.character(data_tot$Gene[tmpMatrix[1]])
                dataAveraged$Group[cnt] = as.character(data_tot$Group[tmpMatrix[1]])
                dataAveraged$Animal[cnt] = as.character(data_tot$Animal[tmpMatrix[1]])
                dataAveraged$Cq[cnt] = mean(tmp)
                tmpCV = sd(tmp)/mean(tmp) 
                dataAveraged$CV[cnt] = tmpCV
                exclusions$Included[(iEx-numRepl):(iEx-1)] = TRUE
                cnt=cnt+1
                input_ok = 1
                
              }else if(ui %in% 1:numRepl){
                # take indicated replicate value anyway
                ui=as.numeric(ui)
                dataAveraged$Threshold[cnt] = data_tot$Threshold[tmpMatrix[ui]]
                dataAveraged$Plate[cnt] = data_tot$Plate[tmpMatrix[ui]]
                dataAveraged$Well[cnt] = as.character(data_tot$Well[tmpMatrix[ui]])
                dataAveraged$Gene[cnt] = as.character(data_tot$Gene[tmpMatrix[ui]])
                dataAveraged$Group[cnt] = as.character(data_tot$Group[tmpMatrix[ui]])
                dataAveraged$Animal[cnt] = as.character(data_tot$Animal[tmpMatrix[ui]])
                dataAveraged$Cq[cnt] = data_tot$Cq[tmpMatrix][ui]
                dataAveraged$CV[cnt] = 0
                exclusions$Included[(iEx-1-numRepl+ui)] = TRUE
                cnt=cnt+1
                input_ok = 1
              }else{
                ui = readline(prompt='Please make a valid entry:')
              }
            }
          }
        }
      }
      
      if(ExcludeRest == 1){#if excluded Gene is HKG, exclude rest of genes of animal
        dataAveraged[dataAveraged$Animal == iA & dataAveraged$Plate == iPl & is.na(dataAveraged$Cq) == FALSE,] = NA
      }
      
    }
  }
  exclusions= na.omit(exclusions)
  returnList = list(dataAveraged, exclusions)
  
  return(returnList)
}