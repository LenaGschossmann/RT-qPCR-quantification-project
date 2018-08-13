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

library(ggplot2)
library(gridExtra)
library(cowplot)
library(grid)
library(ggsignif)
library(scales)
library(RColorBrewer)
# Sys.setenv(JAVA_HOME= 'C:\\Users\\lena_\\Downloads\\Matlab 2017a\\_temp_matlab_R2017a_win64\\sys\\java\\jre\\win64\\jre')
# library(rJava)
library(openxlsx)
library(qpcR)

#### set some parameters and stuff
setwd('C:/Users/Gschossmann/Dropbox/studies/Osnabrück/Universität/Bachelorarbeit_DroBo/experiments/Analysis code/PCR/')
filepath = 'Y:/Lena/Bachelorthesis/Experiments/PCR/qPCR_DNA_damage/analysis'

#testing:
# setwd('C:/Users/Gschossmann/Dropbox/studies/Osnabrück/Universität/Bachelorarbeit_DroBo/experiments/Analysis/PCR')
# filepath = 'C:/Users/lena_/Dropbox/studies/Osnabrück/Universität/Bachelorarbeit_DroBo/experiments/Analysis/PCR/'

##### Input
numRepl = 2 #number of technical replicates
HKG = c('GAPDH') #number of housekeeping genes

acceptable_diff = 0.5 #between replicates

exclude_weirdEffs = 0 #set to 1 when samples with efficiencies <1.9 or >2.1 shall be excluded
methodEff = 'cpD2' #determines from which point the efficiency is calculated (for clarification see: qpcR documentation)
# cpD2: max of sec. derivative | cpD1: max of first derivative | maxE: max of efficiency curve | expR: from exponential region=cpD2-(cpD1-cpD2)
excludeAnimal = c('0')
IC_correction = 1 #set to 1 if there is an Internal Control on every plate that should be used for correction of interplate variability; else set 0

Group = c('HR', 'LR')
Control = 'IR'
#####

#Filepaths for raw Cq data:
# #testing
# f1='C:/Users/lena_/Dropbox/studies/Osnabrück/Universität/Bachelorarbeit_DroBo/experiments/Analysis/PCR/Results_Mitogenes_plate1_LG_150618 - '
f1='Y:/Lena/Bachelorthesis/Experiments/PCR/qPCR_DNA_damage/200718/csv files 1/Results_Hc_short_long_mito_1_LG_200718 - '   #= plate1
f2='Y:/Lena/Bachelorthesis/Experiments/PCR/qPCR_DNA_damage/200718/csv files 2/Results_Hc_short_long_mito_2_LG_200718 - '   #= plate2
f3='Y:/Lena/Bachelorthesis/Experiments/PCR/qPCR_DNA_damage/260718/csv files 1/Results_Hc_short_long_mito_3_LG_260718 - '   #= plate3...
f4='Y:/Lena/Bachelorthesis/Experiments/PCR/qPCR_DNA_damage/260718/csv files 2/Results_Hc_short_long_mito_4_LG_260718 - '
f5='Y:/Lena/Bachelorthesis/Experiments/PCR/qPCR_DNA_damage/270718/csv files/Results_Hc_short_long_mito_5_LG_270718 - '

f_effs = c(f1, f2, f3, f4, f5) #idx has to represent plate number

#####
savepath = 'Y:/Lena/Bachelorthesis/Experiments/PCR/qPCR_DNA_damage/analysis'
#####

source('C:/Users/Gschossmann/Dropbox/studies/Osnabrück/Universität/Bachelorarbeit_DroBo/experiments/Analysis code/PCR/qPCR_Efficiency/qpcR_Efficiencies.R')
    
####################################################################################


# import data from CFX manager exported excel table that has been collected in 1 excel
# the imported data table shall have following columns: plate, date, threshold, animal, group, well, Gene, Cq value
data_tot = read.csv(paste(filepath, 'data_tot.csv', sep='/'), dec='.', sep=';', stringsAsFactors = FALSE)
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
#####

dataAveraged = data.frame(matrix(nrow=(nrow(data_tot)/numRepl), ncol=(ncol(data_tot)+1)))
colnames(dataAveraged) = c(colnames(data_tot), 'CV')

exclusions = data.frame(matrix(ncol=7, nrow=nrow(data_tot)))
colnames(exclusions) = c('Plate', 'Animal', 'Group', 'Gene', 'Well', 'Cq', 'Included')

####
## Substitute Well names (A01 -> A1) for compatibility reasons
withZeros = c('01','02','03','04','05','06','07','08','09')
noZeros=c('1','2','3','4','5','6','7','8','9')
for(iNumber in 1:length(withZeros)){
  data_tot$Well[grep(withZeros[iNumber], data_tot$Well[])] = gsub(withZeros[iNumber], noZeros[iNumber], data_tot$Well[grep(withZeros[iNumber], data_tot$Well[])])
}
####

################################################### Check duplicate precision & combine replicates into one averaged value

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


# Exclude outlier after Dixon
outlier = data.frame(matrix(ncol=ncol(exclusions), nrow=nrow(data_tot)))
colnames(outlier) = colnames(exclusions)
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
exclusions= na.omit(exclusions)

exclusions = rbind(exclusions, outlier)
  

################################################### Efficiencies

# Efficiencies =data.frame()
# Efficiencies=linReg_Efficiencies(data_tot, dataAveraged, f_effs, exclusions, numRepl)
# 
# Av_effs = data.frame(matrix(ncol=2, nrow=length(genes)))
# colnames(Av_effs) = c('Gene', 'Av_Efficiency')
# for(iG in 1: length(genes)){
#   Av_effs$Gene[iG] = genes[iG]
#   Av_effs$Av_Efficiency[iG] = mean(Efficiencies$Efficiency[Efficiencies$Gene == genes[iG]])
# }

Efficiencies= qpcR_Efficiencies(data_tot, f_effs, exclusions, numRepl, exclude_weirdEffs, methodEff)

#exclude animals with weird efficiencies, if HKG -> exclude all genes of that animal on plate
if(exclude_weirdEffs == 1){
  tmpExclude = Efficiencies[Efficiencies$Efficiency == 0,]
  if(sum(tmpExclude$Gene %in% HKG) > 0){
    addEx = which(tmpExclude$Gene %in% HKG)
    for(iAdd in addEx){
      tmpAddexclude = dataAveraged[dataAveraged$Animal == tmpExclude$Animal[iAdd] &
                                     dataAveraged$Plate == tmpExclude$Plate[iAdd] &
                                     dataAveraged$Gene != tmpExclude$Gene[iAdd],]
      for(iAdd2 in 1:nrow(tmpAddexclude)){
        tmpExclude = rbind(tmpExclude, vector(length=ncol(tmpExclude)))
        tmpExclude$Animal[nrow(tmpExclude)] = tmpAddexclude$Animal[iAdd2]
        tmpExclude$Gene[nrow(tmpExclude)] = tmpAddexclude$Gene[iAdd2]
        tmpExclude$Plate[nrow(tmpExclude)] = tmpAddexclude$Plate[iAdd2]
      }
    }
  }
  
  for(iEx in 1:nrow(tmpExclude)){
    excludeMatrix = dataAveraged$Animal == tmpExclude$Animal[iEx] & dataAveraged$Gene == tmpExclude$Gene[iEx] &
      dataAveraged$Plate == tmpExclude$Plate[iEx]
    if(sum(excludeMatrix > 0)){
      exclusions = rbind(exclusions, vector(length=ncol(exclusions)))
      exclusions$Plate[nrow(exclusions)] = dataAveraged$Plate[excludeMatrix]
      exclusions$Well[nrow(exclusions)] = 'Averaged'
      exclusions$Animal[nrow(exclusions)] = dataAveraged$Animal[excludeMatrix]
      exclusions$Group[nrow(exclusions)] = dataAveraged$Group[excludeMatrix]
      exclusions$Gene[nrow(exclusions)] = dataAveraged$Gene[excludeMatrix]
      exclusions$Cq[nrow(exclusions)] = dataAveraged$Cq[excludeMatrix]
      exclusions$Included[nrow(exclusions)] = FALSE
      dataAveraged[excludeMatrix,] = NA
      dataAveraged = na.omit(dataAveraged)
    }
  }
  exclusions = exclusions[!(exclusions$Plate == 0),]
}

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
  
  #correct all animals by respective IC_correction-factor
  for (iPl in plates){
    tmpGenes = data_corrected$Gene[data_corrected$Plate == iPl]
    for (iG in 1:length(tmpGenes)){
      tmpCq = data_corrected$Cq[data_corrected$Plate == iPl][iG] * IC_corrFactors$CorrFac[IC_corrFactors$Plate == iPl & IC_corrFactors$Gene == tmpGenes[iG]] 
      data_corrected$IC_corr_Cq[data_corrected$Plate == iPl][iG] = tmpCq
    }
  }
}else{}


################################################## 1. Normalization (to HKGs)
#N = K*(1+Eff_ref)^Cq_ref/(1+Eff_sample)^Cq_sample  |  K: is going to chancel out in the second normalization step later, so wee dont regard it here
# N here is the ratio of initial amount of target gene over initial amount of reference gene
# tmpMatrix_HKG = data_corrected$Gene %in% HKG

for(iPl in plates){
  for (iA in unique(data_corrected$Animal[data_corrected$Plate == iPl])){
    tmpGenes = unique(data_corrected$Gene[data_corrected$Plate == iPl & data_corrected$Animal == iA])
    tmpGenes = tmpGenes[!(tmpGenes %in% HKG)]
    #calculate mean of (E)^Cq_HKG for all HKG
    tmp_N_HKG = vector(length=length(HKG))
   for(iHKG in 1:length(HKG)){
     tmp_eff_HKG = Efficiencies$Efficiency[Efficiencies$Gene == HKG[iHKG] & Efficiencies$Animal == iA & Efficiencies$Plate == iPl]
     tmp_N_HKG[iHKG] = (tmp_eff_HKG)^(data_corrected$IC_corr_Cq[data_corrected$Gene == HKG[iHKG] & data_corrected$Animal == iA & data_corrected$Plate == iPl])
   }
    tmp_Av_Cq_HKG = mean(tmp_N_HKG)
  
    #calculate (E)^Cq_GOI for all GOI
    for (iG in tmpGenes){
      tmp_eff_GOI = Efficiencies$Efficiency[Efficiencies$Gene == iG & Efficiencies$Animal == iA & Efficiencies$Plate == iPl]
      tmpX_GOI = (tmp_eff_GOI)^(data_corrected$IC_corr_Cq[data_corrected$Animal == iA & data_corrected$Gene == iG])
      data_corrected$HKG_norm_N[data_corrected$Animal == iA & data_corrected$Gene == iG & data_corrected$Plate == iPl] = tmp_Av_Cq_HKG/tmpX_GOI
    }
  }
}

data_corrected = data_corrected[!(data_corrected$Gene %in% HKG),]

################################################## 2. Normalization (to Control Group)
#N = N_treatment / N_ctrl |  K: would going chancel out here
# for derivation see: Scheefe et al. (2006): Quantitative real-time RT-PCR data analysis: current concepts and the novel "gene expression's CT difference" formula.
# split data_corrected into Control (IR) and Treatment (LR, HR)
tmp_Ctrl = data_corrected[data_corrected$Group == Control,]
data_corrected_Treat = data_corrected
# data_corrected_Treat = data_corrected[!(data_corrected$Group == 'IR'),] #use this line if you want just information about treatment groups in the end

tmp_Av_Ctrl = data.frame(matrix(ncol=2, nrow=length(GOI)))
colnames(tmp_Av_Ctrl) = c('Gene', 'Ctrl_Av_N')

#Average all HKG-normalized N values for Control Group for each GOI
for (iG in 1:length(GOI)){
    tmp_Av_Ctrl$Gene[iG] = GOI[iG]
    tmp_Av_Ctrl$Ctrl_Av_N[iG] = mean(tmp_Ctrl$HKG_norm_N[tmp_Ctrl$Gene == GOI[iG]])
}

#Relative quantity of GOI of Treatment Groups compared to Control Group
for(iX in 1:nrow(data_corrected_Treat)){
  tmpGene = data_corrected_Treat$Gene[iX]
  data_corrected_Treat$Rel_quantity[iX] = data_corrected_Treat$HKG_norm_N[iX] / (tmp_Av_Ctrl$Ctrl_Av_N[tmp_Av_Ctrl$Gene == tmpGene])
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
    data_summary$Mean[cnt] = mean(data_corrected_Treat$Rel_quantity[tmpMatrix])
    data_summary$SD[cnt] = sd(data_corrected_Treat$Rel_quantity[tmpMatrix])
    data_summary$SEM[cnt] = sd(data_corrected_Treat$Rel_quantity[tmpMatrix])/sqrt(sum(tmpMatrix))
    cnt=cnt+1
  }
}

# ##Kick out outliers after Dixon's criterion (deviation from mean >2sd)
# data_corrected_Treat_backup = data_corrected_Treat
# for(idx in 1:nrow(data_corrected_Treat)){
#   tmpSD = data_corrected_Treat$Rel_quantity[idx] - data_summary$Mean[data_summary$Gene == data_corrected_Treat$Gene[idx] &
#                                                                        data_summary$Group == data_corrected_Treat$Group[idx]]
#   SDx2 = 2*data_summary$SD[data_summary$Gene == data_corrected_Treat$Gene[idx] & data_summary$Group == data_corrected_Treat$Group[idx]]
#   if(tmpSD >= SDx2){
#     data_corrected_Treat[idx,] = NA
#   }
# }
# data_corrected_Treat=na.omit(data_corrected_Treat)


################################################## Plots

#specify colors
colHR = rgb(255, 0, 0, 255, names = 'HR', max=255)
colIR = rgb(5,190,120, 255, names= 'IR', max=255)
colLR = rgb(87,87,249, 255, names= 'LR', max=255)
cols = c(colHR, colIR, colLR)

ggplot(data=data_corrected_Treat, aes(x=data_corrected_Treat$Gene, y=data_corrected_Treat$Rel_quantity, fill=data_corrected_Treat$Group))+
  stat_boxplot(position=position_dodge(.95))+
  xlab('Genes') + ylab('Relative gene expression (normalized to IR) [AU]]') + theme(text=element_text(size = 15))+
  theme(axis.title.y=element_text(margin = margin(t=0, r=10, b=0, l=0))) +
  theme(axis.title.x=element_text(margin = margin(t=5, r=0, b=0, l=0)))+
  theme(legend.direction = 'vertical', legend.justification='center')+
  scale_fill_manual(values = cols, name='Line')
  # geom_signif(comparison=list(c('HR', 'LR')), annotations = '***', y_position=Y_sign1, tip_length = 0, vjust=0.1)+
  # geom_signif(comparison=list(c('HR', 'IR')), annotations = 'NS.', y_position=Y_sign2, tip_length = 0, vjust=-0.1)+
  # geom_signif(comparison=list(c('IR', 'LR')), annotations = '***', y_position=Y_sign3, tip_length = 0, vjust=0.1)

# ggplot(data=data_corrected_Treat, aes(x=data_corrected_Treat$Gene, y=data_corrected_Treat$Rel_quantity, fill=data_corrected_Treat$Group))+
#   geom_boxplot(aes(ymin = min(data_corrected_Treat$Rel_quantity),
#                    lower = quantile(data_corrected_Treat$Rel_quantity, 0.25),
#                    middle = median(data_corrected_Treat$Rel_quantity),
#                    upper = quantile(data_corrected_Treat$Rel_quantity, 0.75),
#                    ymax = max(data_corrected_Treat$Rel_quantity)),position=position_dodge(.95))
# #   
# 
ggplot(data=data_summary, aes(x=data_summary$Gene, y=data_summary$Mean, fill=data_summary$Group))+
  geom_bar(stat='identity', position=position_dodge(.95))+
  geom_errorbar(aes(ymin=data_summary$Mean-data_summary$SEM, ymax=data_summary$Mean+data_summary$SEM),
                width=.8, stat='identity', position=position_dodge(.95))+
  scale_fill_manual(values = cols, name='Line')




################################################## Stats

data_corrected_Treat$Rel_quantity[]=round(data_corrected_Treat$Rel_quantity[], digits=3)

##Kruskal Wallis
# sink('Kruskal-Wallis-Test.txt')
for (iG in 1:length(GOI)){
  varNames = GOI
  print(as.character(GOI[iG]))
  print(kruskal.test(data=data_corrected_Treat,data_corrected_Treat$Rel_quantity[data_corrected_Treat$Gene == as.character(GOI[iG])] ~ as.factor(data_corrected_Treat$Group[data_corrected_Treat$Gene == as.character(GOI[iG])])))
  }
# sink()

#Man-Whitney-U Test
# sink('Man-Whitney-U-Test.txt')
for(iGr in Group){
  for (iG in 1:length(GOI)){
    tmpTreat = data_corrected_Treat[data_corrected_Treat$Group == Control | data_corrected_Treat$Group == iGr,]
    varNames = GOI
    print(as.character(GOI[iG]))
    print(paste(iGr, 'vs.', Control))
    print(wilcox.test(data=tmpTreat, tmpTreat$Rel_quantity[tmpTreat$Gene == as.character(GOI[iG])] ~ as.factor(tmpTreat$Group[tmpTreat$Gene == as.character(GOI[iG])])))
  }
}
# sink()

#Man-Whitney-U Test
# sink('Man-Whitney-U-Test.txt')
print(paste(Group[1], 'vs.', Group[2]))
for (iG in 1:length(GOI)){
  tmpTreat = data_corrected_Treat[data_corrected_Treat$Group != Control,]
  varNames = GOI
  print(as.character(GOI[iG]))
  print(wilcox.test(data=tmpTreat, tmpTreat$Rel_quantity[tmpTreat$Gene == as.character(GOI[iG])] ~ as.factor(tmpTreat$Group[tmpTreat$Gene == as.character(GOI[iG])])))
}
# sink()

for(iG in GOI){
  for(iTr in unique(data_corrected_Treat$Group)){
    print(sprintf('Gene: %s | Group; %s | n = %i',iG, iTr, nrow(data_corrected_Treat[data_corrected_Treat$Gene == iG & data_corrected_Treat$Group == iTr,])))
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



################################ Further plots

##Stability of RefG
ggplot(data=dataAveraged[dataAveraged$Gene %in% HKG,],
       aes(x=(1:nrow(dataAveraged[dataAveraged$Gene %in% HKG,])),
           y=dataAveraged$Cq[dataAveraged$Gene %in% HKG], group=dataAveraged$Gene[dataAveraged$Gene %in% HKG]))+
  geom_line()

