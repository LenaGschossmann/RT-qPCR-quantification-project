  ## qPCR data quantification
  # based on:
  # Rao et al (2013): An improvement of the 2^(-delta delta CT) method for quantitative real-time polymerase chain reaction data analysis
  # Scheefe et al. (2006): Quantitative real-time RT-PCR data analysis: current concepts and the novel "gene expression's CT difference" formula.
  # Shain EB & Clemens (2008): A new method for robust quantitative and qualitative analysis of real-time PCR. JM.Nucleic Acids Research, 36, e91.
  
  ## Analysis-Procedure
  # 1) check similarity of Cq of replicates and calculate CV, exclude outliers
  # 2) qpcR_Efficiencies:
  #     - import raw data of each amplification cycle and fit multi-parametrical model
  #       ((weighted) nonlinear least-squares (Levenberg-Marquardt) fitting) on combined replicates
  #     --> efficiency is obtained from F_n / F_(n-1) where F is raw fluorescence and n the cycle number
  #     For literature, see: A new method for robust quantitative and qualitative analysis of real-time PCR. Shain EB & Clemens JM.Nucleic Acids Research (2008), 36, e91.
  # 3) Calculate preliminary statistics (Mean, SD) per Treatment group, per gene
  # 4) Exclude outliers according to Dixon's criterion: (sample-mean)>=2*sd
  # 5) Average all Cq values of the Internal Control (IC) and set them as 100%
  #     -> calculate correction factor for IC genes of each plate according to the 100% af the averaged value
  #     -> correct each gene on each plate by the respective genes IC_correction factor for the respective plate
  # 6) Relatively quantify the expression of each of the genes of interest (GOI) in relation to the average of the housekeeping genes (HKG)
  #     (herefore the above estimated efficiencies are used)
  #     ! N = K*(Eff_ref)^Cq_ref/(Eff_sample)^Cq_sample (for Efficiency as between 1 and 2)
  #     for the HKG first [(1+Eff_ref)^Cq_ref] is calculated and then averaged
  # 7) Relatively quantify gene expression in treatment group compared to control group
  #     ! N = N_treatment / N_ctrl
  #     -> For each GOI, average the N from step 4) across all Control group animals
  #     -> For each GOI, relate N from step 4) of each each Treatment group animal to the respective control group averaged N
  # 8) For each GOI, average N across treatment group animals
  # 9) Again Outlier correction?
  
  # Sys.setenv(JAVA_HOME= 'C:\\Users\\lena_\\Downloads\\Matlab 2017a\\_temp_matlab_R2017a_win64\\sys\\java\\jre\\win64\\jre')
  # library(rJava)
  library(openxlsx)
  library(qpcR)
  
  options(decimals=5)
  
  #### set some parameters and stuff
  # setwd('C:/Users/Gschossmann/Dropbox/studies/Osnabrück/Universität/Bachelorarbeit_DroBo/experiments/Analysis code/PCR/')
  # filepath = 'Y:/Lena/Bachelorthesis/Experiments/PCR/qPCR_DNA_damage/analysis'
  
  #testing:
  # setwd('C:/Users/Gschossmann/Dropbox/studies/Osnabrück/Universität/Bachelorarbeit_DroBo/experiments/Analysis Code/PCR')
  filepath = 'C:/Users/lena_/Dropbox/studies/Osnabrück/Universität/Bachelorarbeit_DroBo/experiments/PCR/qPCR_DNA_damage/analysis'
  
  ##### Input
  numRepl = 2 #number of technical replicates
  HKG = c('GAPDH', 'shMito') #number of housekeeping genes
  
  acceptable_diff = 0.5 #between replicates
  
  exclude_weirdEffs = 1 # set to one if outlier >/< mean+/- 2sd shall be eliminated
  effsOuttaRange = 0  #set to 1 when samples with efficiencies <1.9 or >2.1 shall be excluded
  methodEff = 'cpD2' #determines from which point the efficiency is calculated (for clarification see: qpcR documentation)
  # cpD2: max of sec. derivative | cpD1: max of first derivative | maxE: max of efficiency curve | expR: from exponential region=cpD2-(cpD1-cpD2)
  excludeAnimal = c('17')
  IC_correction = 1 #set to 1 if there is an Internal Control on every plate that should be used for correction of interplate variability; else set 0
  
  Group = c('HR', 'LR')
  Control = 'IR'

  #####
  
  #Filepaths for raw Cq data:
  #testing
  f1='C:/Users/lena_/Dropbox/studies/Osnabrück/Universität/Bachelorarbeit_DroBo/experiments/PCR/qPCR_DNA_damage/200718/csv files 1/Results_Hc_short_long_mito_1_LG_200718 - '   #= plate1
  f2='C:/Users/lena_/Dropbox/studies/Osnabrück/Universität/Bachelorarbeit_DroBo/experiments/PCR/qPCR_DNA_damage/200718/csv files 2/Results_Hc_short_long_mito_2_LG_200718 - '   #= plate2
  f3='C:/Users/lena_/Dropbox/studies/Osnabrück/Universität/Bachelorarbeit_DroBo/experiments/PCR/qPCR_DNA_damage/260718/csv files 1/Results_Hc_short_long_mito_3_LG_260718 - '   #= plate3...
  f4='C:/Users/lena_/Dropbox/studies/Osnabrück/Universität/Bachelorarbeit_DroBo/experiments/PCR/qPCR_DNA_damage/260718/csv files 2/Results_Hc_short_long_mito_4_LG_260718 - '
  f5='C:/Users/lena_/Dropbox/studies/Osnabrück/Universität/Bachelorarbeit_DroBo/experiments/PCR/qPCR_DNA_damage/270718/csv files/Results_Hc_short_long_mito_5_LG_270718 - '
  
  
  # f1='Y:/Lena/Bachelorthesis/Experiments/PCR/qPCR_DNA_damage/200718/csv files 1/Results_Hc_short_long_mito_1_LG_200718 - '   #= plate1
  # f2='Y:/Lena/Bachelorthesis/Experiments/PCR/qPCR_DNA_damage/200718/csv files 2/Results_Hc_short_long_mito_2_LG_200718 - '   #= plate2
  # f3='Y:/Lena/Bachelorthesis/Experiments/PCR/qPCR_DNA_damage/260718/csv files 1/Results_Hc_short_long_mito_3_LG_260718 - '   #= plate3...
  # f4='Y:/Lena/Bachelorthesis/Experiments/PCR/qPCR_DNA_damage/260718/csv files 2/Results_Hc_short_long_mito_4_LG_260718 - '
  # f5='Y:/Lena/Bachelorthesis/Experiments/PCR/qPCR_DNA_damage/270718/csv files/Results_Hc_short_long_mito_5_LG_270718 - '
  
  f_effs = c(f1, f2, f3, f4, f5) #idx has to represent plate number
  
  #####
  # savepath = 'Y:/Lena/Bachelorthesis/Experiments/PCR/qPCR_DNA_damage/analysis'
  savepath = 'C:/Users/lena_/Dropbox/studies/Osnabrück/Universität/Bachelorarbeit_DroBo/experiments/PCR/qPCR_DNA_damage/analysis/R_analysis'
  
  #####
  
  # source('C:/Users/Gschossmann/Dropbox/studies/Osnabrück/Universität/Bachelorarbeit_DroBo/experiments/Analysis Code/PCR/qPCR_Efficiency/qpcR_Efficiencies.R')
  # source('C:/Users/Gschossmann/Dropbox/studies/Osnabrück/Universität/Bachelorarbeit_DroBo/experiments/Analysis Code/PCR/SupplementaryFunctions/qPCR_replAveraging.R')
  # source('C:/Users/Gschossmann/Dropbox/studies/Osnabrück/Universität/Bachelorarbeit_DroBo/experiments/Analysis Code/PCR/SupplementaryFunctions/qPCR_Dixon_Outlier.R')
  # source('C:/Users/Gschossmann/Dropbox/studies/Osnabrück/Universität/Bachelorarbeit_DroBo/experiments/Analysis Code/PCR/qPCR_GeneExpressionAssay/deltadeltaCq.R')
  
  source('C:/Users/lena_/Dropbox/studies/Osnabrück/Universität/Bachelorarbeit_DroBo/experiments/Analysis Code/PCR/qPCR_Efficiency/qpcR_Efficiencies.R')
  source('C:/Users/lena_/Dropbox/studies/Osnabrück/Universität/Bachelorarbeit_DroBo/experiments/Analysis Code/PCR/SupplementaryFunctions/qPCR_replAveraging.R')
  source('C:/Users/lena_/Dropbox/studies/Osnabrück/Universität/Bachelorarbeit_DroBo/experiments/Analysis Code/PCR/SupplementaryFunctions/qPCR_Dixon_Outlier.R')
  source('C:/Users/lena_/Dropbox/studies/Osnabrück/Universität/Bachelorarbeit_DroBo/experiments/Analysis Code/PCR/qPCR_GeneExpressionAssay/deltadeltaCq.R')
  
  ####################################################################################
  
  
  # import data from CFX manager exported excel table that has been collected in 1 excel
  # the imported data table shall have following columns: plate, date, threshold, animal, group, well, Gene, Cq value
  data_tot = read.csv(paste(filepath, 'data_tot.csv', sep='/'), dec='.', sep=',', stringsAsFactors = FALSE)
  colnames(data_tot) = c('Threshold', 'Plate', 'Well', 'Gene', 'Animal', 'Group', 'Cq')
  data_tot = data_tot[which(data_tot$Animal != 'NTC'),]  #Kick out NTCs
  data_tot = data_tot[!(is.na(data_tot$Cq)),]
  # data_tot = data_tot[!(data_tot$Animal %in% excludeAnimal),] #Kick out animals that shall be excluded
  
  
  ##### get some fix variables extracted
  plates= sort(unique(data_tot$Plate))
  animals = as.character(unique(data_tot$Animal))
  # animals_perPlate
  genes = as.character(unique(data_tot$Gene))
  GOI = genes[!(genes %in% HKG)]
  genes = append(HKG, as.character(GOI)) #so that the HKG come as first elements
  Treatments = c(Group, Control)
  geneSize = data.frame(matrix(ncol=2, nrow=length(GOI)))
  colnames(geneSize)=c('Gene', 'Size')
  geneSize[1,]=c('D-Loop', 11.19)
  geneSize[2,] = c('Cox3', 11.29)
  #####
  
  ## Substitute Well names (A01 -> A1) for compatibility reasons
  withZeros = c('01','02','03','04','05','06','07','08','09')
  noZeros=c('1','2','3','4','5','6','7','8','9')
  for(iNumber in 1:length(withZeros)){
    data_tot$Well[grep(withZeros[iNumber], data_tot$Well[])] = gsub(withZeros[iNumber], noZeros[iNumber], data_tot$Well[grep(withZeros[iNumber], data_tot$Well[])])
  }
  ####
  
  ################################################### Check duplicate precision & combine replicates into one averaged value
  
  AvList = qPCR_replAveraging(data_tot)
  
  dataAveraged = AvList[[1]]
  exclusions = AvList[[2]]
  
  ################################################### IC-Normalization to exclude effects of plate that do not have to do with the individual efficiencies but global factors
  # 1) Average IC of different plates
  # 2) Set  as 100%
  # 3) Relate the in 4) calculated differences to the averaged IC differences --> Correction factor per plate
  
  #split dataAveraged into data and data_IC
  data_IC = dataAveraged[dataAveraged$Animal == 'IC',]
  data_Genes = dataAveraged[dataAveraged$Animal != 'IC',]
  data_corrected = data_Genes
  
  if(IC_correction == 1){
    IC_corrFactors = data.frame(matrix(ncol=3, nrow=nrow(dataAveraged)))
    colnames(IC_corrFactors) = c('Plate', 'Gene', 'CorrFac')
    
    av_corr_IC = data.frame(matrix(ncol=2, nrow=length(genes)))
    colnames(av_corr_IC) = c('Gene', 'Av_Corr_Value')
    av_corr_IC[,1] = genes
    
    # 1)
    for (iG in 1:length(genes)){ #average all IC values per gene
      tmpMatrix_genes = data_IC$Gene %in% genes[iG]
      av_corr_IC$Gene[iG] = genes[iG]
      av_corr_IC$Av_Corr_Value[iG] = mean(data_IC$Cq[tmpMatrix_genes])
    }
    
    # 3)
    cnt = 1
    for (iPl in plates){
      for (iG in unique(data_IC$Gene[data_IC$Plate == iPl])){
        IC_corrFactors$Plate[cnt] = iPl
        IC_corrFactors$Gene[cnt] = iG
        IC_corrFactors$CorrFac[cnt] = data_IC$Cq[data_IC$Plate == iPl & data_IC$Gene == iG]/av_corr_IC$Av_Corr_Value[av_corr_IC$Gene == iG]
        cnt=cnt+1
      }
    }
    
    IC_corrFactors = na.omit(IC_corrFactors)
    #correct all animals by respective IC_correction-factor
    for (iPl in plates){
      tmpGenes = data_corrected$Gene[data_corrected$Plate == iPl]
      for (iG in 1:length(tmpGenes)){
        tmpCq = data_corrected$Cq[data_corrected$Plate == iPl][iG] * IC_corrFactors$CorrFac[IC_corrFactors$Plate == iPl & IC_corrFactors$Gene == tmpGenes[iG]] 
        data_corrected$IC_corr_Cq[data_corrected$Plate == iPl][iG] = tmpCq
      }
    }
  }else{}
  
  
  #################################################### Exclude outlier after Dixon
  
  DixonList = qPCR_Dixon_Outlier(data_corrected, 1, genes, Treatments)
  data_corrected = DixonList[[1]]
  outlier = DixonList[[2]]
  dataAv_prelimStats = DixonList[[3]]
  exclusions = rbind(exclusions, outlier)
  
  ################################################### Efficiencies
  
  Efficiencies= qpcR_Efficiencies(data_tot, f_effs, exclusions, numRepl, exclude_weirdEffs, effsOuttaRange,methodEff)
  
  #exclude animals with weird efficiencies, if HKG -> exclude all genes of that animal on plate
  if(exclude_weirdEffs == 1){
    tmpExclude = Efficiencies[Efficiencies$Outlier == 1,]
    # check if sample is already kicked out by Dixon correction
    for(iExcl in 1:nrow(tmpExclude)){
      if(sum(outlier$Plate %in% tmpExclude$Plate[iExcl] &
             outlier$Animal %in% tmpExclude$Animal[iExcl] &
             outlier$Gene %in% tmpExclude$Gene[iExcl]) > 0){
        tmpExclude[iExcl,] = NA
      }
    }
    tmpExclude = na.omit(tmpExclude)
    if(sum(tmpExclude$Gene %in% HKG) > 0){
      addEx = which(tmpExclude$Gene %in% HKG)
      for(iAdd in addEx){
        tmpAddexclude = data_corrected[data_corrected$Animal == tmpExclude$Animal[iAdd] &
                                       data_corrected$Plate == tmpExclude$Plate[iAdd] &
                                       data_corrected$Gene != tmpExclude$Gene[iAdd],]
        for(iAdd2 in 1:nrow(tmpAddexclude)){
          tmpExclude = rbind(tmpExclude, vector(length=ncol(tmpExclude)))
          tmpExclude$Animal[nrow(tmpExclude)] = tmpAddexclude$Animal[iAdd2]
          tmpExclude$Gene[nrow(tmpExclude)] = tmpAddexclude$Gene[iAdd2]
          tmpExclude$Plate[nrow(tmpExclude)] = tmpAddexclude$Plate[iAdd2]
        }
      }
    }
    
    for(iEx in 1:nrow(tmpExclude)){
      excludeMatrix = data_corrected$Animal == tmpExclude$Animal[iEx] & data_corrected$Gene == tmpExclude$Gene[iEx] &
        data_corrected$Plate == tmpExclude$Plate[iEx]
      if(sum(excludeMatrix > 0)){
        exclusions = rbind(exclusions, vector(length=ncol(exclusions)))
        exclusions$Plate[nrow(exclusions)] = data_corrected$Plate[excludeMatrix]
        exclusions$Well[nrow(exclusions)] = 'Averaged'
        exclusions$Animal[nrow(exclusions)] = data_corrected$Animal[excludeMatrix]
        exclusions$Group[nrow(exclusions)] = data_corrected$Group[excludeMatrix]
        exclusions$Gene[nrow(exclusions)] = data_corrected$Gene[excludeMatrix]
        exclusions$Cq[nrow(exclusions)] = data_corrected$Cq[excludeMatrix]
        exclusions$Included[nrow(exclusions)] = FALSE
        data_corrected[excludeMatrix,] = NA
        data_corrected = data_corrected[!is.na(data_corrected$Plate),]
      }
    }
    exclusions = exclusions[!(exclusions$Plate == 0),]
  }
  
  
  ################################################# Calculations
  
  for(iRow in 1:nrow(data_corrected)){
    #calculate (E)^Cq
    tmpGene = data_corrected$Gene[iRow]
    tmpAnimal = data_corrected$Animal[iRow]
    tmpPlate = data_corrected$Plate[iRow]
    tmp_eff = Efficiencies$Efficiency[Efficiencies$Gene == tmpGene & Efficiencies$Animal == tmpAnimal & Efficiencies$Plate == tmpPlate]
    data_corrected$Quant[iRow] = (tmp_eff)^(data_corrected$IC_corr_Cq[data_corrected$Animal == tmpAnimal & data_corrected$Gene == tmpGene])
  }
  
  
  
  ################################################## mt Copy number
  #####1. Normalization (to HKGs)
  #N = K*(1+Eff_ref)^Cq_ref/(1+Eff_sample)^Cq_sample  |  K: is going to chancel out in the second normalization step later, so wee dont regard it here
  # N here is the ratio of initial amount of target gene over initial amount of reference gene
  
  for(iPl in plates){
    for (iA in unique(data_corrected$Animal[data_corrected$Plate == iPl])){
      tmp_mtCN = data_corrected$Quant[data_corrected$Animal == iA & data_corrected$Gene == 'GAPDH' & data_corrected$Plate == iPl] /
        data_corrected$Quant[data_corrected$Animal == iA & data_corrected$Gene == 'shMito' & data_corrected$Plate == iPl]
      data_corrected$mtCN[data_corrected$Animal == iA & data_corrected$Gene == 'shMito' & data_corrected$Plate == iPl] = tmp_mtCN
    }
  }
  
  # data_corrected_all = data_corrected
  # data_corrected = data_corrected[!(data_corrected$Gene %in% HKG),]
  
  ##### 2. Normalization (to Control Group)
  #N = N_treatment / N_ctrl |  K: would going chancel out here
  # for derivation see: Scheefe et al. (2006): Quantitative real-time RT-PCR data analysis: current concepts and the novel "gene expression's CT difference" formula.
  # split data_corrected into Control (IR) and Treatment (LR, HR)
  data_Ctrl = data_corrected[data_corrected$Group == Control,]
  #Average all HKG-normalized N values for Control Group for shMito
  tmp_Av_Ctrl = mean(data_Ctrl$Quant[data_Ctrl$Gene == 'shMito'])
  
  #Relative quantity of shMito of Treatment Groups compared to Control Group
  for(iX in which(data_corrected$Gene=='shMito')){
    data_corrected$mtCN[iX] = data_corrected$Quant[iX] / tmp_Av_Ctrl
  }
  
  Av_shMito = mean(data_corrected$mtCN[which(data_corrected$Gene=='shMito')])
  
  ################################################## Calculate Lesion frequencies
  # following 
  #1)  E_Cq(D-Loop) / (E_Cq(shMito)/Av_E_Cq_shMito)
  #2) IR correction
  
  for(iPl in plates){
    for(iA in unique(data_corrected$Animal[data_corrected$Plate == iPl])){
      tmp_shMito = data_corrected$Quant[data_corrected$Plate == iPl & data_corrected$Animal == iA & data_corrected$Gene == 'shMito']
      for(iG in GOI){
        tmpVal = data_corrected$Quant[data_corrected$Plate == iPl & data_corrected$Animal == iA & data_corrected$Gene == iG] / (tmp_shMito /tmp_shMito)
        data_corrected$LesionFreq[data_corrected$Plate == iPl & data_corrected$Animal == iA & data_corrected$Gene == iG] = tmpVal
      }
    }
  }
  
  #take IR average and divide the other genes value by it
  for(iG in 1:length(GOI)){
    tmp_Av_Ctrl = mean(data_corrected$LesionFreq[data_corrected$Gene == GOI[iG] & data_corrected$Group == Control])
    for(iRow in which(data_corrected$Gene == GOI[iG])){
      data_corrected$LesionFreq[iRow] = data_corrected$LesionFreq[iRow] / tmp_Av_Ctrl
      data_corrected$LesionFreq[iRow] = (-log(data_corrected$LesionFreq[iRow]))*(10/as.numeric(geneSize$Size[geneSize$Gene == GOI[iG]]))
    }
  }
  
  
  
  
  ################################################## Calculate Mean Difference & SEM
  
  #summary
  data_summary = data.frame(matrix(ncol=4, nrow=3*length(GOI))) #Mean and SEM per GOI per Group
  colnames(data_summary) = c('Group', 'Gene', 'Mean', 'SEM')
  
  cnt=1
  for(iT in 1:length(unique(data_corrected_Treat$Group))){
    
    for(iG in 1:length(GOI)){
      tmpMatrix = data_corrected_Treat$Group == unique(data_corrected_Treat$Group)[iT] & data_corrected_Treat$Gene == GOI[iG]
      data_summary$Group[cnt] = unique(data_corrected_Treat$Group)[iT]
      data_summary$Gene[cnt] = as.character(GOI[iG])
      data_summary$Mean[cnt] = mean(data_corrected_Treat$mtCN[tmpMatrix])
      data_summary$SD[cnt] = sd(data_corrected_Treat$mtCN[tmpMatrix])
      data_summary$SEM[cnt] = sd(data_corrected_Treat$mtCN[tmpMatrix])/sqrt(sum(tmpMatrix))
      cnt=cnt+1
    }
  }
  
  
  
  ################################ Output
  
  #Total data
  write.xlsx(data_tot, paste(savepath,'data_total.xlsx',sep='/'))
  
  #final data, corrected to IR (only treatment groups)
  write.xlsx(data_corrected_Treat, paste(savepath,'data_corrected_Norm_to_IR.xlsx',sep='/'))
  
  #exclusions
  write.xlsx(exclusions, paste(savepath,'excluded samples.xlsx',sep='/'))
  
  #IC correction factors
  write.xlsx(IC_corrFactors, paste(savepath,'IC_correctino_factors.xlsx',sep='/'))
  
  #Efficiencies
  write.xlsx(Efficiencies, paste(savepath,'Efficiencies.xlsx',sep='/'))
  
  #final descriptive summary
  write.xlsx(data_summary, paste(savepath,'data_corrected_summary.xlsx',sep='/'))
  
  
  
  
