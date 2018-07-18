# for literature see also: Rao et al (2013): An improvement of the 2^(-delta delta CT) method for quantitative real-time polymerase chain reaction data analysis
#  take raw data (Intensity per cycle) of all samples (except the ones excluded in 1))  <--- linReg_Efficiencies.R
#     -> find Window of Linearity for each curve in separate:
#         - take differences between consecutive cycles, starting with 2 cycles before the Cq (eliminates background difference)
#         - compare differences: is next one greater /equal than current?
#         - if YES: inclued cycle that corresponds to n (from which n-1 has been subtracted)
#         - if NO: take current delta plus consecutive delta (so cycle points n and n+1)
#     ! equation to describe original Intensity curve: I(k) = y_background + F * X_initial * E^k | k= cycle, F = conversion factor: amount of product to intensity
#     ! as cycle differences are used: I = I(k) - I(k-1) --> F * X_initial * E^k * (E-1)
#     ! log-transform equation: log(I(k)-I(k-1)) = log[F * X_initial * E^k] + k* logE
#     -> take log of the selected delta values and fit regression model
#     ! --> U(k) = beta_0 + k*beta_1 --> beta_1 = slope
#     ! --> E = beta_1  -> slope represents 
#     -> Take mean of slope of replicates


linReg_Efficiencies = function(dataThresholdCq, dataAveraged, f_effs, exclusions, numRepl){
  #as dataThresholdCq give data_tot to the function
  #order of files in f_effs needs to correspond to plate numbering
  Efficiencies = data.frame(matrix(ncol=4))
  colnames(Efficiencies) = c('Animal', 'Gene', 'Slope', 'Efficiency')
  cntE = 1
  tmpPlates = sort(unique(dataAveraged$Plate))
  for(iF in 1:length(tmpPlates)){
    tmpData = dataThresholdCq[dataThresholdCq$Plate == tmpPlates[iF],]
    #single cycle values
    tmpRaw = read.csv(paste(f_effs[tmpPlates[iF]], 'Quantification Amplification Results_SYBR.csv'), dec=',', sep=';', stringsAsFactors = FALSE)
    #rename columns
    tmpRaw = tmpRaw[,2:ncol(tmpRaw)]
    
    tmpRaw=tmpRaw[,colnames(tmpRaw) %in% c('Cycle',as.character(unique(tmpData$Well)))]
    tmpRawSum = dataThresholdCq[dataThresholdCq$Plate == tmpPlates[iF],]
    
    
    tmpAnimals = as.character(unique(dataAveraged$Animal[dataAveraged$Plate == plates[iF]]))
    tmpAnimals = tmpAnimals[tmpAnimals != 'IC']
    tmpGenes = as.character(unique(dataAveraged$Gene[dataAveraged$Plate == plates[iF]]))
    
    Master = data.frame(matrix(ncol=3, nrow=length(tmpAnimals)*length(tmpGenes)*numRepl))
    colnames(Master)=c('Animal', 'Gene', 'Well')
    cntA=1
    for(iA in tmpAnimals){
      cntARange = cntA:(cntA +((length(tmpGenes)*numRepl))-1)
      Master$Animal[cntARange] = iA
      Master$Gene[cntARange] = as.character(tmpRawSum$Gene[tmpRawSum$Animal == iA])
      Master$Well[cntARange] = as.character(tmpRawSum$Well[tmpRawSum$Animal == iA])
      cntA = cntA + (length(tmpGenes)*numRepl)
    }
    
    tmpSubRaw = data.frame(matrix(ncol=2, nrow=length(tmpRaw$Cycle)))
    colnames(tmpSubRaw) = c('Cycle', 'I')
    tmpSubRaw$Cycle = tmpRaw$Cycle
    
    for(iA in tmpAnimals){
      for(iG in tmpGenes){
        
        include=1
        cntTmp=1
        #Check if outlier and if any of the replicates is used
        if(iA %in% exclusions$Animal[] & iG %in% exclusions$Gene[exclusions$Animal == iA] & sum(exclusions$Included[exclusions$Animal == iA & exclusions$Gene == iG])==0){ #in case all replicates are excluded
          include = 0
        }else if(iA %in% exclusions$Animal[exclusions$Animal == iA] & iG %in% exclusions$Gene[] & sum(exclusions$Included[exclusions$Animal == iA & exclusions$Gene == iG])==1){ #in case only one replicate is used
          tmpVals = data.frame(matrix(ncol=2))
          colnames(tmpVals) = c('Cycle', 'I')
          replicate = which(exclusions$Included[exclusions$Animal == iA & exclusions$Gene == iG]==TRUE)
          tmpWell = tmpData$Well[tmpData$Animal == iA & tmpData$Gene == iG][replicate]
          tmpThreshold = round(tmpData$Cq[tmpData$Animal == iA & tmpData$Well == tmpWell], digits=0)
          tmpSubRaw$I = tmpRaw[,colnames(tmpRaw) == Master$Well[Master$Animal == iA & Master$Gene == iG][replicate]]
          
          tmpSlope = 1
          
          #Find consecutive cycles with steepest slope (Max in sec derivativ) & background fluorescence correction
          samplePoints = c(tmpThreshold:(tmpThreshold+4)) #linear regression depending on threshold cycle (rounded) and 3 sequential cycles
          
          delta = tmpSubRaw$I[samplePoints[1]] - tmpSubRaw$I[(samplePoints[1]-1)]
          for(cntC in samplePoints){
            delta_new = tmpSubRaw$I[cntC+1] - tmpSubRaw$I[cntC]
            if(delta_new >= delta){
              tmpVals = rbind(tmpVals, c(tmpSubRaw$Cycle[cntC+1], delta_new))
            }else{}
            delta = delta_new
          }
          tmpVals = na.omit(tmpVals)
          if(nrow(tmpVals) > 1){
            #Linear regression fitted on values of every replicates separately
            tmpVals$I[]=log10(tmpVals$I[])
            tmpLM = lm(formula=I~Cycle, data=tmpVals)
            tmpSlope[iR] = tmpLM$coefficients[2]
          }else{#in case there is just one delta thats greater in slope than take that value plus the consecutive delta two have at least 2 values for regression fitting
            tmpVals = rbind(tmpVals, vector(length=ncol(tmpVals)))
            tmpVals$Cycle[2]=(tmpVals$Cycle[1]+1)
            tmpVals$I[2] = tmpSubRaw$I[tmpSubRaw$Cycle == tmpVals$Cycle[2]]-tmpSubRaw$I[tmpSubRaw$Cycle == (tmpVals$Cycle[2]-1)]
            tmpVals$I[]=log10(tmpVals$I[])
            tmpLM = lm(formula=I~Cycle, data=tmpVals)
            tmpSlope[iR] = tmpLM$coefficients[2]
          }
          
          # delta = diff(as.numeric(tmpSubRaw$I[tmpSubRaw$Cycle %in% as.integer(samplePoints)]))
          # delta_max = delta[max(delta)]
          # delta_max_Cyc = 1
          # cntRange = cntTmp:cntTmp+length(delta_max)
          # 
          # for(iRowBind in 1:(length(cntRange))){
          #   tmpVals = rbind(tmpVals[], vector(length=ncol(tmpVals)))
          # }
          # 
          # tmpVals$Animal[cntRange] = iA
          # tmpVals$Gene[cntRange] = iG
          # tmpVals$Cycle[cntRange]= samplePoints[1:(length(samplePoints)-1)]
          # tmpVals$I[cntRange]= tmpI
          
        }else{
          for(iR in 1:numRepl){
            tmpVals = data.frame(matrix(ncol=2))
            colnames(tmpVals) = c('Cycle', 'I')
            tmpWell = tmpData$Well[tmpData$Animal == iA & tmpData$Gene == iG][iR]
            tmpThreshold = round(tmpData$Cq[tmpData$Animal == iA & tmpData$Well == tmpWell], digits=0)
            # samplePoints = c(tmpThreshold-1, tmpThreshold, tmpThreshold+1, tmpThreshold+2, tmpThreshold+3) #linear regression depending on threshold cycle (rounded) and 3 sequential cycles
            # cntRange = cntTmp:(cntTmp+length(samplePoints)-2)
            
            tmpSubRaw$I = tmpRaw[,which(colnames(tmpRaw) == Master$Well[Master$Animal == iA & Master$Gene == iG][iR])]
            
            tmpSlope = vector(length=numRepl)
            
            #Find consecutive cycles with steepest slope (Max in sec derivativ) & background fluorescence correction
            samplePoints = c((tmpThreshold-1):(tmpThreshold+4)) #find window of linearity in area around Cq 
            delta = tmpSubRaw$I[samplePoints[1]] - tmpSubRaw$I[(samplePoints[1]-1)]
            for(cntC in samplePoints){
              delta_new = tmpSubRaw$I[cntC+1] - tmpSubRaw$I[cntC]
              if(delta_new >= delta){
                tmpVals = rbind(tmpVals, c(tmpSubRaw$Cycle[cntC+1], delta_new))
              }else{}
              delta = delta_new
            }
            tmpVals = na.omit(tmpVals)
            # if(nrow(tmpVals) == 0){
            #   print(sprintf('!!! Warning! Something is wrong with well: % s on plate: %i. Please delete manually from raw data!', tmpWell, plates[iF]))
            # }else
            if(nrow(tmpVals) > 1){
              #Linear regression fitted on values of every replicates separately
              tmpVals$I[]=log10(tmpVals$I[])
              tmpLM = lm(formula=I~Cycle, data=tmpVals)
              tmpSlope[iR] = tmpLM$coefficients[2]
            }else{ #in case there is just one delta thats greater in slope than take that value plus the consecutive delta two have at least 2 values for regression fitting
              tmpVals = rbind(tmpVals, vector(length=ncol(tmpVals)))
              tmpVals$Cycle[2]=(tmpVals$Cycle[1]+1)
              tmpVals$I[2] = tmpSubRaw$I[tmpSubRaw$Cycle == tmpVals$Cycle[2]]-tmpSubRaw$I[tmpSubRaw$Cycle == (tmpVals$Cycle[2]-1)]
              tmpVals$I[]=log10(tmpVals$I[])
              tmpLM = lm(formula=I~Cycle, data=tmpVals)
              tmpSlope[iR] = tmpLM$coefficients[2]
            }
            # #for background fluorescence correction
            # tmpI = diff(as.numeric(tmpSubRaw$I[tmpSubRaw$Cycle %in% as.integer(samplePoints)]))
            # 
            # for(iRowBind in 1:(length(cntRange))){
            #   tmpVals = rbind(tmpVals[], vector(length=ncol(tmpVals)))
            # }
            # tmpVals$Animal[cntRange] = iA
            # tmpVals$Gene[cntRange] = iG
            # tmpVals$Cycle[cntRange]= samplePoints[1:(length(samplePoints)-1)]
            # tmpVals$I[cntRange]= tmpI
            # 
            # cntTmp = cntTmp+length(samplePoints)-1
          }
        }
        
        # #Linear regression fitted on values from all replicates of a sample
        if(include == 1){
          #Calculate and Save Efficiency
          Efficiencies = rbind(Efficiencies[], vector(length=ncol(Efficiencies)))
          Efficiencies$Plate[cntE] = tmpPlates[iF]
          Efficiencies$Animal[cntE] = iA
          Efficiencies$Gene[cntE] = iG
          Efficiencies$Slope[cntE] = mean(tmpSlope)
          Efficiencies$Efficiency[cntE] = exp(mean(tmpSlope))
          cntE = cntE+1
          
        #   # tmpVals = tmpVals[1:(nrow(tmpVals)-1),]
        #   tmpVals$I[]=log10(tmpVals$I[])
        #   tmpLM = lm(formula=I~Cycle, data=tmpVals)
        # 
        #   #Calculate and Save Efficiency
        #   Efficiencies = rbind(Efficiencies[], vector(length=ncol(Efficiencies)))
        #   Efficiencies$Plate[cntE] = tmpPlates[iF]
        #   Efficiencies$Animal[cntE] = iA
        #   Efficiencies$Gene[cntE] = iG
        #   Efficiencies$Slope[cntE] = tmpLM$coefficients[2]
        #   Efficiencies$Efficiency[cntE] = exp(tmpLM$coefficients[2])
        #   cntE = cntE+1
        }else{}
        
        ggplot(data=tmpVals, aes(x=Cycle, y=I))+geom_point()+geom_smooth(method=lm)
        
      } #closes Gene loop
    } #closes Animal loop
    
  }
  Efficiencies = Efficiencies[Efficiencies$Animal != FALSE,] # to delete any empty rows at end
  
  return(Efficiencies)
  }
  