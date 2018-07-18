# for literature see also: Rao et al (2013): An improvement of the 2^(-delta delta CT) method for quantitative real-time polymerase chain reaction data analysis

qpcR_Efficiencies = function(dataThresholdCq, f_effs, exclusions, numRepl, weirdEffs, methodEff){
  ##as dataThresholdCq give data_tot to the function
  ##order of files in f_effs needs to correspond to plate numbering
  Efficiencies = data.frame(matrix(ncol=4))
  colnames(Efficiencies) = c('Animal', 'Gene', 'Slope', 'Efficiency')
  cntE = 1
  tmpPlates = sort(unique(dataThresholdCq$Plate))
  for(iF in 1:length(tmpPlates)){
    tmpData = dataThresholdCq[dataThresholdCq$Plate == tmpPlates[iF] & dataThresholdCq$Animal != 'IC',]
    ##single cycle values
    tmpRaw = read.csv(paste(f_effs[tmpPlates[iF]], 'Quantification Amplification Results_SYBR.csv'), dec=',', sep=';', stringsAsFactors = FALSE)
    ##rename columns
    tmpRaw = tmpRaw[,2:ncol(tmpRaw)]
    
    tmpRaw=tmpRaw[,colnames(tmpRaw) %in% c('Cycle',as.character(unique(tmpData$Well)))] #kick out wells with NTCs
    # tmpRawSum = dataThresholdCq[dataAveraged$Plate == tmpPlates[iF],]
    
    tmpAnimals = as.character(unique(tmpData$Animal[tmpData$Plate == tmpPlates[iF]]))
    # tmpAnimals = tmpAnimals[tmpAnimals != 'IC']
    
    ##Create a MaterIndex where you save what yu find in each well
    Master = data.frame(matrix(ncol=3))
    colnames(Master)=c('Animal', 'Gene', 'Well')
    cntA=1
    for(iA in tmpAnimals){
      tmpGenes = tmpData$Gene[tmpData$Animal == iA]
      for(icount in 1:length(tmpGenes)){
        Master = rbind(Master, vector(length=3))
        Master$Animal[cntA] = iA
        Master$Gene[cntA] = as.character(tmpData$Gene[tmpData$Animal == iA][icount])
        Master$Well[cntA] = as.character(tmpData$Well[tmpData$Animal == iA][icount])
        cntA=cntA+1
      }
    }
    Master = Master[1:(nrow(Master)-1),]
    
    tmpSubRaw = data.frame(matrix(ncol=(numRepl+1), nrow=length(tmpRaw$Cycle))) #Intensities of just one sample + replicates
    # colnames(tmpSubRaw) = c('Cycle', 'I')
    tmpSubRaw[,1] = tmpRaw$Cycle
    
    for(iA in tmpAnimals){
      tmpGenes = unique(Master$Gene[Master$Animal == iA])
      for(iG in tmpGenes){
        include=1
        cntTmp=1
        tmpVals = data.frame(matrix(ncol=(numRepl+1)))
        ##Check if outlier and if any of the replicates is used
        if(iA %in% exclusions$Animal & iG %in% exclusions$Gene[exclusions$Animal == iA & exclusions$Plate == tmpPlates[iF]] &
           sum(exclusions$Included[exclusions$Animal == iA & exclusions$Gene == iG])==0){ 
        ##in case all replicates are excluded
          include = 0
        }else if(iA %in% exclusions$Animal &
                 iG %in% exclusions$Gene[exclusions$Animal == iA & exclusions$Plate == tmpPlates[iF]] &
                 sum(exclusions$Included[exclusions$Animal == iA & exclusions$Gene == iG])==1){ 
        ##in case only one replicate is used
          replicate = which(exclusions$Included[exclusions$Animal == iA & exclusions$Gene == iG & exclusions$Plate == tmpPlates[iF]]==TRUE)
          tmpWell = exclusions$Well[exclusions$Animal == iA & exclusions$Gene == iG & exclusions$Plate == tmpPlates[iF] & exclusions$Included == TRUE]
          # tmpThreshold = round(tmpData$Cq[tmpData$Animal == iA & tmpData$Well == tmpWell], digits=0)
          tmpSubRaw[,2] = tmpRaw[,colnames(tmpRaw) == tmpWell]
          ##make use of qpcR package and fit 5-praameter model on amplification curve and get efficiency (calculated as F_n/F_n-1 with F being a point within the exponential phase)
          obj = pcrfit(tmpSubRaw, 1, 2)
          eff=efficiency(obj, plot=FALSE, type=methodEff)
          finalEff= eff$eff[[1]]
        }else{
          ##in case all replicates are used
          ##make use of qpcR package and fit 5-parameter model on signals of amplification curves of all replicates and get efficiency (calculated as F_n/F_n-1 with F being a point within the exponential phase)
          # for(iR in 1:numRepl){
          #   tmpWell = tmpData$Well[tmpData$Animal == iA & tmpData$Gene == iG][iR]
          #   tmpSubRaw[,(iR+1)] = tmpRaw[,colnames(tmpRaw) == tmpWell]
          # }
          # obj = pcrfit(tmpSubRaw, 1, (2:(ncol(tmpSubRaw))))
          # tmpEff=efficiency(obj, plot=TRUE, type='cpD2')
          # finalEff = tmpEff$eff
           
          ##In case you want to fit model to each replicates curve first and then average efficiencies
          tmpEff=vector(length=numRepl)
          for(iR in 1:numRepl){
            tmpWell = tmpData$Well[tmpData$Animal == iA & tmpData$Gene == iG][iR]
            tmpSubRaw[,(iR+1)] = tmpRaw[,colnames(tmpRaw) == tmpWell]
            #make use of qpcR package and fit 5-praameter model on amplification curve and get efficiency (calculated as F_n/F_n-1 with F being a point within the exponential phase)
            obj = pcrfit(tmpSubRaw, 1, (iR+1))
            eff = efficiency(obj, plot=FALSE, type=methodEff)
            tmpEff[iR]= eff$eff[[1]]
          }
          finalEff=mean(tmpEff)
          if(weirdEffs == 1){
            if(finalEff < 1.9 | finalEff > 2.1){
              take = tmpEff[which(tmpEff < 2.1 & tmpEff>1.9)]
              if(sum(take) == FALSE){
                finalEff = FALSE
              }else{
                finalEff = mean(tmpEff[take])
              }
            }
          }
        }
        Efficiencies = rbind(Efficiencies[], vector(length=ncol(Efficiencies)))
        Efficiencies$Plate[cntE] = tmpPlates[iF]
        Efficiencies$Animal[cntE] = iA
        Efficiencies$Gene[cntE] = iG
        Efficiencies$Efficiency[cntE] = finalEff
        # Efficiencies$Slope[cntE] = mean(tmpSlope)
        # Efficiencies$Efficiency[cntE] = exp(mean(tmpSlope))
        cntE = cntE+1
        # ggplot(data=tmpVals, aes(x=Cycle, y=I))+geom_point()+geom_smooth(method=lm)
        
      } #closes Gene loop
    } #closes Animal loop
    
  }
  Efficiencies = Efficiencies[Efficiencies$Animal != FALSE,] # to delete any empty rows at end
  
  return(Efficiencies)
  }
  