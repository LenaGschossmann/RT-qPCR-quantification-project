### qPCR data quantification with the comparative delta-delta-2G method

library(ggplot2)
library(gridExtra)
library(cowplot)
library(grid)
library(ggsignif)
library(scales)
library(RColorBrewer)
# Sys.setenv(JAVA_HOME= 'C:\\Users\\lena_\\Downloads\\Matlab 2017a\\_temp_matlab_R2017a_win64\\sys\\java\\jre\\win64\\jre')
# library(rJava)
library(xlsx)

#### set some parameters and stuff
setwd('C:/Users/Gschossmann/Dropbox/studies/Osnabrück/Universität/Bachelorarbeit_DroBo/experiments/Analysis/PCR')
filepath = 'Y:/Lena/Bachelorthesis/Experiments/PCR/Mito_gene_expression/analysis'
# #testing:
# filepath = 'C:/Users/lena_/Dropbox/studies/Osnabrück/Universität/Bachelorarbeit_DroBo/experiments/Analysis/PCR/'


numRepl = 2 #number of technical replicates
HKG = c('HPRT', 'GAPDH') #number of housekeeping genes

acceptable_diff = 0.6 #between replicates

# import data from CFX manager exported excel table that has been collected in 1 excel
# the imported data table shall have following columns: plate, date, threshold, animal, group, well, Gene, Cq value
data_tot = read.csv(paste(filepath, 'data_tot1.csv', sep='/'), dec='.', sep=';')
colnames(data_tot) = c('Threshold', 'Plate', 'Well', 'Gene', 'Animal', 'Group', 'Cq')
data_tot = data_tot[which(data_tot$Animal != 'NTC'),]  #Kick out NTCs

# import the efficiencies
effs = read.csv(paste(filepath, 'efficiency calculation/efficiency_values.csv', sep='/'))
# #testing:
# effs = read.csv(paste(filepath, 'efficiency_values.csv', sep='/'))

##### get some fix variables extracted
plates= unique(data_tot$Plate)
animals = as.character(unique(data_tot$Animal))
# animals_perPlate
genes = as.character(unique(data_tot$Gene))
GOI = genes[!(genes %in% HKG)]
genes = append(HKG, as.character(GOI)) #so that the HKG come as first elements
Group = c('HR', 'LR')
control = 'IR'
#####

dataAveraged = data.frame(matrix(nrow=(nrow(data_tot)/numRepl), ncol=(ncol(data_tot)+1)))
colnames(dataAveraged) = c(colnames(data_tot), 'CV')

exclusions = data.frame(matrix(ncol=6, nrow=nrow(data_tot)))
colnames(exclusions) = c('Plate', 'Animal', 'Group', 'Gene', 'Cq', 'Included')

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
  
  for (iA in unique(data_tot$Animal[which(data_tot$Plate == iPl)])){
    
    for (iG in unique(data_tot$Gene[which(data_tot$Plate == iPl & data_tot$Animal == iA)])){
      
      if(iG %in% exclusions$Gene[exclusions$Animal == iA  & exclusions$Plate == iPl]){
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
            exclusions$Cq[iEx] = data_tot$Cq[tmpMatrix[iR]]
            exclusions$Included[iEx]= FALSE
            iEx = iEx+1
          }
          
          print(sprintf("Plate: %d | Animal: %s | Gene: %s", iPl, iA, iG))
          print(paste('Value:',data_tot$Cq[tmpMatrix]))
          print(paste('Difference:', round(max(tmp)-min(tmp), digits=4)))
          print(paste('Average:', mean(data_tot$Cq[data_tot$Gene == iG & data_tot$Animal != iA])))
          ui = readline(prompt='Type Number of Replicate (1. -> <1>) you want to keep, <0> if you dont want to keep any or <a> if you want to keep all:')
          input_ok = 0
          while(input_ok == 0){
            if(ui == '0'){
              input_ok = 1
              
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
outlier = data.frame(matrix(ncol=6, nrow=nrow(data_tot)))
colnames(outlier) = c('Plate', 'Animal', 'Group', 'Gene', 'Cq', 'Included')
iO=1
for(iT in unique(dataAveraged$Group[!(dataAveraged$Group == '') & is.na(dataAveraged$Cq) == FALSE])){
  for(iPl in plates){
    for(iA in unique(dataAveraged$Animal[dataAveraged$Plate == iPl & dataAveraged$Group == iT & is.na(dataAveraged$Cq) == FALSE])){
      excludeAnimal = 0
      for(iHKG in HKG){
        tmpMean = dataAv_prelimStats$Mean[dataAv_prelimStats$Gene == iHKG & dataAv_prelimStats$Treatment == iT]
        tmpSD = dataAv_prelimStats$SD[dataAv_prelimStats$Gene == iHKG & dataAv_prelimStats$Treatment == iT]
        if(abs(dataAveraged$Cq[dataAveraged$Animal == iA & dataAveraged$Gene == iHKG & is.na(dataAveraged$Cq) == FALSE]-tmpMean) >= 2*tmpSD){ #Apply Dixons outlier criterion on averaged replicates
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
        outlier$Included[iOrange]= FALSE
        dataAveraged[tmpOut,] = NA
        iO = 1+iOrange[length(iOrange)]
      }else{
        for(iG in GOI){
          tmpMean = dataAv_prelimStats$Mean[dataAv_prelimStats$Gene == iG & dataAv_prelimStats$Treatment == iT]
          tmpSD = dataAv_prelimStats$SD[dataAv_prelimStats$Gene == iG & dataAv_prelimStats$Treatment == iT]
          if(abs(dataAveraged$Cq[dataAveraged$Animal == iA & dataAveraged$Gene == iG & is.na(dataAveraged$Cq) == FALSE]-tmpMean) >= 2*tmpSD){ #Apply Dixons outlier criterion on averaged replicates
            tmpOut = dataAveraged$Animal == iA & dataAveraged$Plate == iPl & dataAveraged$Gene == iG & is.na(dataAveraged$Cq) == FALSE
            outlier$Plate[iO] = dataAveraged$Plate[tmpOut]
            outlier$Animal[iO] = as.character(dataAveraged$Animal[tmpOut])
            outlier$Group[iO] = as.character(dataAveraged$Group[tmpOut])
            outlier$Gene[iO] = as.character(dataAveraged$Gene[tmpOut])
            outlier$Cq[iO] = dataAveraged$Cq[tmpOut]
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

################################################### IC-Normalization
# 1) Average IC of different plates
# 2) Set  as 100%
# 3) Relate the in 4) calculated differences to the averaged IC differences --> Correction factor per plate


#split dataAveraged into data and tmp_IC
tmp_IC = dataAveraged[dataAveraged$Animal == 'IC',]
data_Genes = dataAveraged[dataAveraged$Animal != 'IC',]
IC_corrFactors = data.frame(matrix(ncol=(length(genes)+1), nrow=length(plates)))
colnames(IC_corrFactors) = append('Plate', as.character(genes))
IC_corrFactors$Plate[] = 1:length(plates)

av_corr_IC = data.frame(matrix(ncol=2, nrow=length(genes)))
colnames(av_corr_IC) = c('Gene', 'Av_Corr_Value')
av_corr_IC[,1] = genes

data_corrected = data_Genes

# 1)
for (iG in 1:length(genes)){ #average all IC values per gene
  tmpMatrix_genes = tmp_IC$Gene %in% genes[iG]
  av_corr_IC$Gene[iG] = genes[iG]
  av_corr_IC$Av_Corr_Value[iG] = mean(tmp_IC$Cq[tmpMatrix_genes])
}

# 3)
for (iPl in 1:length(plates)){
  tmpGenes = unique(tmp_IC$Gene[tmp_IC$Plate == plates[iPl]])
  for (iG in 1:length(tmpGenes)){
    IC_corrFactors[iPl,(iG+1)] = tmp_IC$Cq[tmp_IC$Plate == plates[iPl] & tmp_IC$Gene == tmpGenes[iG]]/av_corr_IC$Av_Corr_Value[iG]
    IC_corrFactors$Plate[iPl] = plates[iPl]
  }
}

#correct all animals by respective IC_correction-factor
for (iPl in 1:length(plates)){
  tmpGenes = unique(data_corrected$Gene[data_corrected$Plate == plates[iPl]])
  for (iA in unique(data_corrected$Animal[which(data_corrected$Plate == plates[iPl])])){
   for (iG in 1:length(tmpGenes)){
      tmpCq = data_corrected$Cq[data_corrected$Plate == plates[iPl] & data_corrected$Animal == iA & data_corrected$Gene == tmpGenes[iG]] * IC_corrFactors[iPl, (iG+1)] 
      data_corrected$IC_corr_Cq[data_corrected$Plate == plates[iPl] & data_corrected$Animal == iA & data_corrected$Gene == tmpGenes[iG]] = tmpCq
    }
  }
}


################################################## 1. Normalization (to HKGs)
#N = K*(1+Eff_ref)^Cq_ref/(1+Eff_sample)^Cq_sample  |  K: is going to chancel out in the second normalization step later, so wee dont regard it here
# N here is the ratio of initial amount of target gene over initial amount of reference gene
# tmpMatrix_HKG = data_corrected$Gene %in% HKG

for (iA in unique(data_corrected$Animal[])){
  #calculate mean of (1+E)^Cq_HKG for all HKG
  tmp_Cq_HKG = vector(length=length(HKG))
  for(iHKG in 1:length(HKG)){
    tmp_eff_HKG = effs$Efficiency[effs$Gene == HKG[iHKG]]
    tmp_Cq_HKG[iHKG] = (1+ tmp_eff_HKG)^(data_corrected$IC_corr_Cq[data_corrected$Gene == HKG[iHKG] & data_corrected$Animal == iA])
  }
  tmp_Av_Cq_HKG = mean(tmp_Cq_HKG)
  #calculate (1+E)^Cq_GOI for all GOI
  for (iG in 1:length(GOI)){
    tmp_eff_GOI = effs$Efficiency[effs$Gene == as.character(GOI[iG])]
    tmpX_GOI = (1+tmp_eff_GOI)^(data_corrected$IC_corr_Cq[data_corrected$Animal == iA & data_corrected$Gene == GOI[iG]])
    
    data_corrected$HKG_norm_N[data_corrected$Animal == iA & data_corrected$Gene == GOI[iG]] = tmp_Av_Cq_HKG/tmpX_GOI
  }
}

data_corrected = data_corrected[!(data_corrected$Gene %in% HKG),]

################################################## 2. Normalization (to Control Group)
#N = N_treatment * (N_ctrl)^-1 |  K: would going chancel out here
# for derivation see: Scheefe et al. (2006): Quantitative real-time RT-PCR data analysis: current concepts and the novel "gene expression's CT difference" formula.
# split data_corrected into Control (IR) and Treatment (LR, HR)
tmp_Ctrl = data_corrected[data_corrected$Group == 'IR',]
data_corrected = data_corrected[!(data_corrected$Group == 'IR'),]

tmp_Av_Ctrl = data.frame(matrix(ncol=2, nrow=length(GOI)))
colnames(tmp_Av_Ctrl) = c('Gene', 'Ctrl_Av_N')

#Average all HKG-normalized N values for Control Group for each GOI
for (iG in 1:length(GOI)){
  tmp_Av_Ctrl$Gene[iG] = as.character(GOI[iG])
  tmp_Av_Ctrl$Ctrl_Av_N[iG] = mean(tmp_Ctrl$HKG_norm_N[tmp_Ctrl$Gene == GOI[iG]])
}

#Relative quantity of GOI of Treatment Groups compared to Control Group
for(iX in 1:nrow(data_corrected)){
  tmpGene = data_corrected$Gene[iX]
  data_corrected$Rel_quantity[iX] = data_corrected$HKG_norm_N[iX] * (tmp_Av_Ctrl$Ctrl_Av_N[tmp_Av_Ctrl$Gene == tmpGene])^(-1)
}


################################################## Calculate Mean Difference & SEM

#summary
data_summary = data.frame(matrix(ncol=4, nrow=length(Group)*length(GOI))) #Mean and SEM per GOI per Group
colnames(data_summary) = c('Group', 'Gene', 'Mean', 'SEM')

cnt=1
for(iGr in 1:length(Group)){
  
  for(iG in 1:length(GOI)){
    tmpMatrix = data_corrected$Group == Group[iGr] & data_corrected$Gene == GOI[iG]
    data_summary$Group[cnt] = as.character(Group[iGr])
    data_summary$Gene[cnt] = as.character(GOI[iG])
    data_summary$Mean[cnt] = mean(data_corrected$Rel_quantity[tmpMatrix])
    data_summary$SD[cnt] = sd(data_corrected$Rel_quantity[tmpMatrix])
    data_summary$SEM[cnt] = sd(data_corrected$Rel_quantity[tmpMatrix])/sqrt(sum(tmpMatrix))
    cnt=cnt+1
  }
}

#Kick out outliers after Dixon's criterion (deviation from mean >2sd)
data_corrected_NoOut = data_corrected
for(idx in 1:nrow(data_corrected_NoOut)){
  tmpSD = data_corrected_NoOut$Rel_quantity[idx] - data_summary$Mean[data_summary$Gene == data_corrected_NoOut$Gene[idx] & data_summary$Group == data_corrected_NoOut$Group[idx]]
  SDx2 = 2*data_summary$SD[data_summary$Gene == data_corrected_NoOut$Gene[idx] & data_summary$Group == data_corrected_NoOut$Group[idx]]
  if(tmpSD >= SDx2){
    data_corrected_NoOut[idx,] = NA
  }
}
data_corrected_NoOut = na.omit(data_corrected_NoOut)

################################################## Plots

ggplot(data=data_corrected_NoOut, aes(x=data_corrected_NoOut$Gene, y=data_corrected_NoOut$Rel_quantity, fill=data_corrected_NoOut$Group))+
  stat_boxplot(stat='identity', position=position_dodge(.95))

ggplot(data=data_summary, aes(x=data_summary$Gene, y=data_summary$Mean, fill=data_summary$Group))+
  geom_bar(stat='identity', position=position_dodge(.95))+
  geom_errorbar(aes(ymin=data_summary$Mean-data_summary$SEM, ymax=data_summary$Mean+data_summary$SEM), width=.8, stat='identity', position=position_dodge(.95))




################################################## Stats

data_corrected$Rel_quantity[]=round(data_corrected$Rel_quantity[], digits=3)

#Kruskal Wallis
# sink('Kruskal-Wallis-Test.txt')
for (iG in 1:length(GOI)){
  varNames = GOI
  print(as.character(GOI[iG]))
  print(kruskal.test(data=data_corrected, data_corrected$Rel_quantity[data_corrected$Gene == as.character(GOI[iG])] ~ as.factor(data_corrected$Group[data_corrected$Gene == as.character(GOI[iG])])))
}
# sink()

#Man-Whitney-U Test
# sink('Man-Whitney-U-Test.txt')
for (iG in 1:length(GOI)){
  varNames = GOI
  print(as.character(GOI[iG]))
  print(wilcox.test(data=data_corrected_NoOut, data_corrected_NoOut$Rel_quantity[data_corrected_NoOut$Gene == as.character(GOI[iG])] ~ as.factor(data_corrected_NoOut$Group[data_corrected_NoOut$Gene == as.character(GOI[iG])])))
}
# sink()




################################ Output

print(paste('CV was smaller than it should be in following replicates:', tmp_warning_CV))








