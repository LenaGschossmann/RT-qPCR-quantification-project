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
# 3) Average all Cq values of the Internal Control (IC) and set them as 100%
#     -> calculate correction factor for IC genes of each plate according to the 100% af the averaged value
#     -> correct each gene on each plate by the respective genes IC_correction factor for the respective plate
# 4) Calculate preliminary statistics (Mean, SD) per Treatment group, per gene
# 5) Exclude outliers according to Dixon's criterion: (sample-mean)>=2*sd
# 6) Relatively quantify the expression of each of the genes of interest (GOI) in relation to the average of the housekeeping genes (HKG)
#     (herefore the above estimated efficiencies are used)
#     ! N = K*(1+Eff_ref)^Cq_ref/(1+Eff_sample)^Cq_sample (for Efficiency as between 1 and 2)
#     for the HKG first [(1+Eff_ref)^Cq_ref] is calculated and then averaged
# 7) Relatively quantify gene expression in treatment group compared to control group
#     ! N = N_treatment / N_ctrl
#     -> For each GOI, average the N from step 6) across all Control group animals
#     -> For each GOI, relate N from step 6) of each each Treatment group animal to the respective control group averaged N
# 8) For each GOI, average N across treatment group animals
# 9) Again Outlier correction?

library(openxlsx)
library(qpcR)

options(decimals=5)

#### set some parameters and stuff
# setwd('C:/Users/Gschossmann/Dropbox/studies/Osnabrück/Universität/Bachelorarbeit_DroBo/experiments/Analysis Code/PCR/')
# filepath = 'Y:/Lena/Bachelorthesis/Experiments/PCR/qPCR_DNA_damage/analysis'

setwd('C:/Users/lena_/Dropbox/studies/Osnabrück/Universität/Bachelorarbeit_DroBo/experiments/Analysis Code/PCR/')
filepath = 'C:/Users/lena_/Dropbox/studies/Osnabrück/Universität/Bachelorarbeit_DroBo/experiments/PCR/Mito_gene_expression/analysis'

part = 'part3'

##### Input
numRepl = 2 #number of technical replicates
HKG = c('HPRT', 'GAPDH') #number of housekeeping genes

acceptable_diff = 0.5 #between replicates

exclude_weirdEffs = 1 # set to one if outlier >/< mean+/- 2sd shall be eliminated
effsOuttaRange = 1  #set to 1 when samples with efficiencies <1.9 or >2.1 shall be excluded
methodEff = 'cpD2' #determines from which point the efficiency is calculated (for clarification see: qpcR documentation)
# cpD2: max of sec. derivative | cpD1: max of first derivative | maxE: max of efficiency curve | expR: from exponential region=cpD2-(cpD1-cpD2)

effUpperBond = 0 #set to 1 if all efficiencies > 2 should be set to 2 as they are artefacts
effTakeMean = 0 #set to 1 if the mean of the gene's efficiency should be used

excludeAnimal = c('17')
IC_correction = 1 #set to 1 if there is an Internal Control on every plate that should be used for correction of interplate variability; else set 0

Group = c('HR', 'LR')
Control = 'IR'
#####

#Filepaths for raw Cq data:
f1='C:/Users/lena_/Dropbox/studies/Osnabrück/Universität/Bachelorarbeit_DroBo/experiments/PCR/Mito_gene_expression/150618/csv files/Results_Mitogenes_plate1_LG_150618 - '   #= plate1
f2='C:/Users/lena_/Dropbox/studies/Osnabrück/Universität/Bachelorarbeit_DroBo/experiments/PCR/Mito_gene_expression/180618/csv files/Results_Mitogenes_plate2_LG_180618 - '   #= plate2
f3='C:/Users/lena_/Dropbox/studies/Osnabrück/Universität/Bachelorarbeit_DroBo/experiments/PCR/Mito_gene_expression/190618/csv files 1/Results_Mitogenes_plate3_LG_190618 - '   #= plate3...
f4='C:/Users/lena_/Dropbox/studies/Osnabrück/Universität/Bachelorarbeit_DroBo/experiments/PCR/Mito_gene_expression/190618/csv files 2/Results_Mitogenes_plate4_LG_190618 - '
f5='C:/Users/lena_/Dropbox/studies/Osnabrück/Universität/Bachelorarbeit_DroBo/experiments/PCR/Mito_gene_expression/200618/csv files/Results_Mitogenes_plate5_LG_200618 - '
f6='C:/Users/lena_/Dropbox/studies/Osnabrück/Universität/Bachelorarbeit_DroBo/experiments/PCR/Mito_gene_expression/050718/csv files 1/Results_Mitogenes_plate6_LG_050718 - '
f7='C:/Users/lena_/Dropbox/studies/Osnabrück/Universität/Bachelorarbeit_DroBo/experiments/PCR/Mito_gene_expression/050718/csv files 2/Results_Mitogenes_plate7_LG_050718 - '
f8= 'C:/Users/lena_/Dropbox/studies/Osnabrück/Universität/Bachelorarbeit_DroBo/experiments/PCR/Mito_gene_expression/060718/csv files 1/Results_Mitogenes_plate8_LG_060718 - '
f9='C:/Users/lena_/Dropbox/studies/Osnabrück/Universität/Bachelorarbeit_DroBo/experiments/PCR/Mito_gene_expression/060718/csv files 2/Results_Mitogenes_plate9_LG_060718 - '
f10='C:/Users/lena_/Dropbox/studies/Osnabrück/Universität/Bachelorarbeit_DroBo/experiments/PCR/Mito_gene_expression/090718/csv files 1/Results_Mitogenes_plate10_LG_090718 - '
f11='C:/Users/lena_/Dropbox/studies/Osnabrück/Universität/Bachelorarbeit_DroBo/experiments/PCR/Mito_gene_expression/100718/csv files/Results_Mitogenes_plate11_LG_100718 - '
f12='C:/Users/lena_/Dropbox/studies/Osnabrück/Universität/Bachelorarbeit_DroBo/experiments/PCR/Mito_gene_expression/110718/csv files/Results_Mitogenes_plate12_LG_110718 - '
f13='C:/Users/lena_/Dropbox/studies/Osnabrück/Universität/Bachelorarbeit_DroBo/experiments/PCR/Mito_gene_expression/160718/csv files 1/Results_Mitogenes_plate13_LG_160718 - '
f14='C:/Users/lena_/Dropbox/studies/Osnabrück/Universität/Bachelorarbeit_DroBo/experiments/PCR/Mito_gene_expression/160718/csv files 2/Results_Mitogenes_plate14_LG_160718 - '


# f1='Y:/Lena/Bachelorthesis/Experiments/PCR/Mito_gene_expression/150618/csv files/Results_Mitogenes_plate1_LG_150618 - '   #= plate1
# f2='Y:/Lena/Bachelorthesis/Experiments/PCR/Mito_gene_expression/180618/csv files/Results_Mitogenes_plate2_LG_180618 - '   #= plate2
# f3='Y:/Lena/Bachelorthesis/Experiments/PCR/Mito_gene_expression/190618/csv files 1/Results_Mitogenes_plate3_LG_190618 - '   #= plate3...
# f4='Y:/Lena/Bachelorthesis/Experiments/PCR/Mito_gene_expression/190618/csv files 2/Results_Mitogenes_plate4_LG_190618 - '
# f5='Y:/Lena/Bachelorthesis/Experiments/PCR/Mito_gene_expression/200618/csv files/Results_Mitogenes_plate5_LG_200618 - '
# f6='Y:/Lena/Bachelorthesis/Experiments/PCR/Mito_gene_expression/050718/csv files 1/Results_Mitogenes_plate6_LG_050718 - '
# f7='Y:/Lena/Bachelorthesis/Experiments/PCR/Mito_gene_expression/050718/csv files 2/Results_Mitogenes_plate7_LG_050718 - '
# f8= 'Y:/Lena/Bachelorthesis/Experiments/PCR/Mito_gene_expression/060718/csv files 1/Results_Mitogenes_plate8_LG_060718 - '
# f9='Y:/Lena/Bachelorthesis/Experiments/PCR/Mito_gene_expression/060718/csv files 2/Results_Mitogenes_plate9_LG_060718 - '
# f10='Y:/Lena/Bachelorthesis/Experiments/PCR/Mito_gene_expression/090718/csv files 1/Results_Mitogenes_plate10_LG_090718 - '
# f11='Y:/Lena/Bachelorthesis/Experiments/PCR/Mito_gene_expression/100718/csv files/Results_Mitogenes_plate11_LG_100718 - '
# f12='Y:/Lena/Bachelorthesis/Experiments/PCR/Mito_gene_expression/110718/csv files/Results_Mitogenes_plate12_LG_110718 - '
# f13='Y:/Lena/Bachelorthesis/Experiments/PCR/Mito_gene_expression/160718/csv files 1/Results_Mitogenes_plate13_LG_160718 - '
# f14='Y:/Lena/Bachelorthesis/Experiments/PCR/Mito_gene_expression/160718/csv files 2/Results_Mitogenes_plate14_LG_160718 - '
f_effs = c(f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, f14) #idx has to represent plate number

#####
# savepath = 'Y:/Lena/Bachelorthesis/Experiments/PCR/Mito_gene_expression/analysis/analysed files'
savepath = 'C:/Users/lena_/Dropbox/studies/Osnabrück/Universität/Bachelorarbeit_DroBo/experiments/PCR/Mito_gene_expression/analysis/R_analysis'
#####

# source('C:/Users/Gschossmann/Dropbox/studies/Osnabrück/Universität/Bachelorarbeit_DroBo/experiments/Analysis Code/PCR/qPCR_Efficiency/qpcR_Efficiencies.R')
# source('C:/Users/Gschossmann/Dropbox/studies/Osnabrück/Universität/Bachelorarbeit_DroBo/experiments/Analysis Code/PCR/SupplementaryFunctions/qPCR_replAveraging.R')
# source('C:/Users/Gschossmann/Dropbox/studies/Osnabrück/Universität/Bachelorarbeit_DroBo/experiments/Analysis Code/PCR/SupplementaryFunctions/qPCR_Dixon_Outlier.R')
# source('C:/Users/Gschossmann/Dropbox/studies/Osnabrück/Universität/Bachelorarbeit_DroBo/experiments/Analysis Code/PCR/SupplementatyFunctions/deltadeltaCq.R')

source('C:/Users/lena_/Dropbox/studies/Osnabrück/Universität/Bachelorarbeit_DroBo/experiments/Analysis Code/PCR/qPCR_Efficiency/qpcR_Efficiencies.R')
source('C:/Users/lena_/Dropbox/studies/Osnabrück/Universität/Bachelorarbeit_DroBo/experiments/Analysis Code/PCR/SupplementaryFunctions/qPCR_replAveraging.R')
source('C:/Users/lena_/Dropbox/studies/Osnabrück/Universität/Bachelorarbeit_DroBo/experiments/Analysis Code/PCR/SupplementaryFunctions/qPCR_Dixon_Outlier.R')
source('C:/Users/lena_/Dropbox/studies/Osnabrück/Universität/Bachelorarbeit_DroBo/experiments/Analysis Code/PCR/SupplementaryFunctions/deltadeltaCq.R')


####################################################################################

# import data from CFX manager exported excel table that has been collected in 1 excel
# the imported data table shall have following columns: plate, date, threshold, animal, group, well, Gene, Cq value
data_tot = read.csv(paste(filepath, paste(part,'csv', sep='.'), sep='/'), dec='.', sep=';', stringsAsFactors = FALSE)
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

Efficiencies= qpcR_Efficiencies(data_tot, f_effs, exclusions, numRepl, exclude_weirdEffs, effsOuttaRange, methodEff)

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
      data_corrected = na.omit(data_corrected)
    }
  }
  exclusions = exclusions[!(exclusions$Plate == 0),]
}


########## use mean of efficiency
if(effTakeMean == 1){
  for(iPl in plates){
    tmpGenes = unique(Efficiencies$Gene[Efficiencies$Plate == iPl])
    for(iG in tmpGenes){
      tmpMean = mean(Efficiencies$Efficiency[Efficiencies$Gene == iG & Efficiencies$Plate == iPl])
      Efficiencies$Efficiency[Efficiencies$Gene == iG & Efficiencies$Plate == iPl]= tmpMean
    }
  }
  
}


######### set all efficiencies > 2 to 2
if(effUpperBond == 1){
   Efficiencies$Efficiency[Efficiencies$Efficiency > 2] = 2
}



################################################## 1. Normalization (to HKGs)
#N = K*(1+Eff_ref)^Cq_ref/(1+Eff_sample)^Cq_sample  |  K: is going to chancel out in the second normalization step later, so wee dont regard it here | the '1' in front of the 'E' is not needed since efficiency already between 1 and 2
# N here is the ratio of initial amount of target gene over initial amount of reference gene

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
    tmp_Av_HKG = mean(tmp_N_HKG)
  
    #calculate (E)^Cq_GOI for all GOI
    for (iG in tmpGenes){
      tmp_eff_GOI = Efficiencies$Efficiency[Efficiencies$Gene == iG & Efficiencies$Animal == iA & Efficiencies$Plate == iPl]
      tmpX_GOI = (tmp_eff_GOI)^(data_corrected$IC_corr_Cq[data_corrected$Animal == iA & data_corrected$Gene == iG])
      data_corrected$HKG_norm_N[data_corrected$Animal == iA & data_corrected$Gene == iG & data_corrected$Plate == iPl] = tmp_Av_HKG/tmpX_GOI
    }
    # Calculate RefG stability
    for (iG in HKG){
      tmp_eff_HKG = Efficiencies$Efficiency[Efficiencies$Gene == iG & Efficiencies$Animal == iA & Efficiencies$Plate == iPl]
      tmpX_HKG = (tmp_eff_HKG)^(data_corrected$IC_corr_Cq[data_corrected$Animal == iA & data_corrected$Gene == iG])
      data_corrected$HKG_norm_N[data_corrected$Animal == iA & data_corrected$Gene == iG & data_corrected$Plate == iPl] = tmpX_HKG
    }
    
  }
}




data_corrected_all = data_corrected
data_corrected = data_corrected[!(data_corrected$Gene %in% HKG),]

################################################## 2. Normalization (to Control Group)
#N = N_treatment / N_ctrl |  K: would going chancel out here
# for derivation see: Scheefe et al. (2006): Quantitative real-time RT-PCR data analysis: current concepts and the novel "gene expression's CT difference" formula.
# split data_corrected into Control (IR) and Treatment (LR, HR)
data_Ctrl = data_corrected[data_corrected$Group == Control,]
# data_corrected = data_corrected[!(data_corrected$Group == 'IR'),] #use this line if you want just information about treatment groups in the end

tmp_Av_Ctrl = data.frame(matrix(ncol=2, nrow=length(GOI)))
colnames(tmp_Av_Ctrl) = c('Gene', 'Ctrl_Av_N')

#Average all HKG-normalized N values for Control Group for each GOI
for (iG in 1:length(GOI)){
    tmp_Av_Ctrl$Gene[iG] = GOI[iG]
    tmp_Av_Ctrl$Ctrl_Av_N[iG] = mean(data_Ctrl$HKG_norm_N[data_Ctrl$Gene == GOI[iG]])
}

#Relative quantity of GOI of Treatment Groups compared to Control Group
for(iX in 1:nrow(data_corrected)){
  tmpGene = data_corrected$Gene[iX]
  data_corrected$Rel_quantity[iX] = data_corrected$HKG_norm_N[iX] / (tmp_Av_Ctrl$Ctrl_Av_N[tmp_Av_Ctrl$Gene == tmpGene])
}


################################################# Compare with 2^-DeltaDeltaCq method

data_corrected$DeltaDelta = deltadeltaCq(data_corrected_all, GOI, HKG, Control, Group)

#summary
data_summaryDeltaDelta = data.frame(matrix(ncol=4, nrow=3*length(GOI))) #Mean and SEM per GOI per Group
colnames(data_summaryDeltaDelta) = c('Group', 'Gene', 'Mean', 'SEM')

cnt=1
for(iT in 1:length(unique(data_corrected$Group))){
  
  for(iG in 1:length(GOI)){
    tmpMatrix = data_corrected$Group == unique(data_corrected$Group)[iT] & data_corrected$Gene == GOI[iG]
    data_summaryDeltaDelta$Group[cnt] = unique(data_corrected$Group)[iT]
    data_summaryDeltaDelta$Gene[cnt] = as.character(GOI[iG])
    data_summaryDeltaDelta$Mean[cnt] = mean(data_corrected$DeltaDelta[tmpMatrix])
    data_summaryDeltaDelta$SD[cnt] = sd(data_corrected$DeltaDelta[tmpMatrix])
    data_summaryDeltaDelta$SEM[cnt] = sd(data_corrected$DeltaDelta[tmpMatrix])/sqrt(sum(tmpMatrix))
    cnt=cnt+1
  }
}


################################################## Calculate Mean Difference & SEM
#summary
data_summary = data.frame(matrix(ncol=4, nrow=3*length(GOI))) #Mean and SEM per GOI per Group
colnames(data_summary) = c('Group', 'Gene', 'Mean', 'SEM')

cnt=1
for(iT in 1:length(unique(data_corrected$Group))){
  
  for(iG in 1:length(GOI)){
    tmpMatrix = data_corrected$Group == unique(data_corrected$Group)[iT] & data_corrected$Gene == GOI[iG]
    data_summary$Group[cnt] = unique(data_corrected$Group)[iT]
    data_summary$Gene[cnt] = as.character(GOI[iG])
    data_summary$Mean[cnt] = mean(data_corrected$Rel_quantity[tmpMatrix])
    data_summary$SD[cnt] = sd(data_corrected$Rel_quantity[tmpMatrix])
    data_summary$SEM[cnt] = sd(data_corrected$Rel_quantity[tmpMatrix])/sqrt(sum(tmpMatrix))
    cnt=cnt+1
  }
}



################################################### Output
#Save Corrected data
write.csv(data_corrected, paste(savepath,paste(part, 'data_corrected.csv', sep='_'),sep='/'))

write.csv(data_summary, paste(savepath,paste(part, 'data_summary.csv', sep='_'),sep='/'))

write.csv(data_summaryDeltaDelta, paste(savepath,paste(part, 'data_summaryDeltaDelta.csv', sep='_'),sep='/'))


# #final data, corrected to IR (only treatment groups)
# write.xlsx(data_corrected, paste(savepath,'data_corrected_Norm_to_IR_part.xlsx',sep='/'))
# 
# #exclusions
# write.xlsx(exclusions, paste(savepath,'excluded samples.xlsx',sep='/'))
# 
# #IC correction factors
# write.xlsx(IC_corrFactors, paste(savepath,'IC_correctino_factors_part.xlsx',sep='/'))
# 
# #Efficiencies
# write.xlsx(Efficiencies, paste(savepath,'Efficiencies.xlsx',sep='/'))
# 
# #final descriptive summary
# write.xlsx(data_summary, paste(savepath,'data_corrected_summary.xlsx',sep='/'))




################################################# Chain all data_corrected together and summarize
#...after all parts have been saved separately:
readpathData = 'C:/Users/lena_/Dropbox/studies/Osnabrück/Universität/Bachelorarbeit_DroBo/experiments/PCR/Mito_gene_expression/analysis/R_analysis'
x1=read.csv(paste(readpathData, 'part1_data_corrected.csv', sep='/'))
x2=read.csv(paste(readpathData, 'part2_data_corrected.csv', sep='/'))
x3=read.csv(paste(readpathData, 'part3_data_corrected.csv', sep='/'))
data_corrected_tot = rbind(x1, rbind(x2, x3))
data_corrected_tot = data_corrected_tot[,2:ncol(data_corrected_tot)]
write.xlsx(data_corrected_tot, paste(savepath,'data_corrected_ALL.xlsx',sep='/'))

allGOI = unique(data_corrected_tot$Gene)

#summary
data_summary_tot = data.frame(matrix(ncol=4, nrow=3*length(allGOI))) #Mean and SEM per GOI per Group
colnames(data_summary_tot) = c('Group', 'Gene', 'Mean', 'SEM')

cnt=1
for(iT in 1:length(unique(data_corrected_tot$Group))){
  
  for(iG in 1:length(allGOI)){
    tmpMatrix = data_corrected_tot$Group == unique(data_corrected_tot$Group)[iT] & data_corrected_tot$Gene == allGOI[iG]
    data_summary_tot$Group[cnt] = unique(data_corrected_tot$Group)[iT]
    data_summary_tot$Gene[cnt] = as.character(allGOI[iG])
    data_summary_tot$Mean[cnt] = mean(data_corrected_tot$Rel_quantity[tmpMatrix])
    data_summary_tot$SD[cnt] = sd(data_corrected_tot$Rel_quantity[tmpMatrix])
    data_summary_tot$SEM[cnt] = sd(data_corrected_tot$Rel_quantity[tmpMatrix])/sqrt(sum(tmpMatrix))
    cnt=cnt+1
  }
}


