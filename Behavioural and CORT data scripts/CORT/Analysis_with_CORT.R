##### Include CORT

# source('C:/Users/Gschossmann/Dropbox/studies/Osnabrück/Universität/Bachelorarbeit_DroBo/experiments/Analysis Code/PCR/qPCR_GeneExpressionAssay/getCORT.R')
source('C:/Users/lena_/Dropbox/studies/Osnabrück/Universität/Bachelorarbeit_DroBo/experiments/Analysis Code/CORT/getCORT.R')

savepath = 'C:/Users/lena_/Dropbox/studies/Osnabrück/Universität/Bachelorarbeit_DroBo/experiments/PCR/Mito_gene_expression/analysis/R_analysis'


############## Get CORT data
fileX = 'Cohort VR/Data Males Project 1-CORT.csv'
filepathCORT = paste('C:/Users/lena_/Dropbox/studies/Osnabrück/Universität/Bachelorarbeit_DroBo/experiments/Behavioural tests', fileX, sep='/')
cort = getCORT(filepathCORT)

############## Gene Expression data: get data_corrected_tot: data.frame with final data from all genes
#...after all parts have been saved separately:
x1=read.csv(paste(savepath, 'part1_data_corrected.csv', sep='/'))
x2=read.csv(paste(savepath, 'part2_data_corrected.csv', sep='/'))
x3=read.csv(paste(savepath, 'part3_data_corrected.csv', sep='/'))
data_corrected_tot = rbind(x1, rbind(x2, x3))

############## Behavioural Data
#xxx

##############write CORT data into other data.frame
#work with all data
data_with_CORT = data_corrected_tot
for(iR in 1:nrow(data_with_CORT)){
  data_with_CORT$CORT[iR] = cort$CORT[cort$ID == data_with_CORT$Animal[iR]]
}


#######################################  LinRegression & Pearson Correlation
tmpGOI = unique(data_with_CORT$Gene)
LinRegGOI = data.frame(matrix(ncol=4, nrow=length(tmpGOI)))
colnames(LinRegGOI) = c('Gene', 'Slope', 'Intercept', 'Observations')
corr_coeffsGOI = data.frame(matrix(ncol=4, nrow=length(tmpGOI)))
colnames(corr_coeffsGOI) = c('Gene', 'R', 'p', 'Observations')

for(iG in 1:length(tmpGOI)){
  tmpData = data_with_CORT[data_with_CORT$Gene == tmpGOI[iG],]
  #Pearson
  pearson_corr = rcorr(tmpData$Rel_quantity,tmpData$CORT, type='pearson')
  corr_coeffsGOI$Gene[iG] = as.character(tmpGOI[iG])
  corr_coeffsGOI$R[iG] = pearson_corr$r[1,2]
  corr_coeffsGOI$p[iG] = round(pearson_corr$P[1,2], digits=6)
  corr_coeffsGOI$Observations[iG] = nrow(tmpData)
  #LinReg
  mod = lm(Rel_quantity ~ CORT, data=tmpData)
  LinRegGOI$Gene[iG] = as.character(tmpGOI[iG])
  LinRegGOI$Slope[iG] = mod$coefficients[2]
  LinRegGOI$Intercept[iG] = mod$coefficients[1]
  LinRegGOI$Observations[iG] = nrow(tmpData)
}

#Output
savepath = 'C:/Users/lena_/Dropbox/studies/Osnabrück/Universität/Bachelorarbeit_DroBo/experiments/PCR/Mito_gene_expression/analysis/R_analysis'
write.xlsx(LinRegGOI, paste(savepath,'LinReg_CORT_GOI.xlsx',sep='/'))
write.xlsx(corr_coeffsGOI, paste(savepath,'Corr_coefficients_CORT_GOI.xlsx',sep='/'))
# 

############  SPLIT GROUPS - LinRegression & Pearson Correlation
tmpGOI = as.character(unique(data_with_CORT$Gene))
tmpGroups = as.character(unique(data_with_CORT$Group))
LinRegGOI_groups= data.frame(matrix(ncol=5, nrow=length(tmpGOI)*length(tmpGroups)))
colnames(LinRegGOI_groups) = c('Gene', 'Group', 'Slope', 'Intercept', 'Observations')
corr_coeffsGOI_groups = data.frame(matrix(ncol=5, nrow=length(tmpGOI)*length(tmpGroups)))
colnames(corr_coeffsGOI_groups) = c('Gene','Group', 'R', 'p', 'Observations')
idx=0
for(iG in tmpGOI){
  for(iGr in tmpGroups){
    idx=idx+1
    tmpData = data_with_CORT[data_with_CORT$Gene == iG & data_with_CORT$Group == iGr,]
    #Pearson
    if(nrow(tmpData) < 5){
      corr_coeffsGOI_groups$Gene[idx] = iG
      corr_coeffsGOI_groups$Group[idx] = iGr
      corr_coeffsGOI_groups$R[idx] = NA
      corr_coeffsGOI_groups$p[idx] = NA
      corr_coeffsGOI_groups$Observations[idx] = nrow(tmpData)
    }else{
      pearson_corr = rcorr(tmpData$Rel_quantity,tmpData$CORT, type='pearson')
      corr_coeffsGOI_groups$Gene[idx] = iG
      corr_coeffsGOI_groups$Group[idx] = iGr
      corr_coeffsGOI_groups$R[idx] = pearson_corr$r[1,2]
      corr_coeffsGOI_groups$p[idx] = round(pearson_corr$P[1,2], digits=6)
      corr_coeffsGOI_groups$Observations[idx] = nrow(tmpData)
    }
    #LinReg
    mod = lm(Rel_quantity ~ CORT, data=tmpData)
    LinRegGOI_groups$Gene[idx] = iG
    LinRegGOI_groups$Group[idx] = iGr
    LinRegGOI_groups$Slope[idx] = mod$coefficients[2]
    LinRegGOI_groups$Intercept[idx] = mod$coefficients[1]
    LinRegGOI_groups$Observations[idx] = nrow(tmpData)
  }
}

#Output
savepath = 'C:/Users/lena_/Dropbox/studies/Osnabrück/Universität/Bachelorarbeit_DroBo/experiments/PCR/Mito_gene_expression/analysis/R_analysis'
write.xlsx(LinRegGOI_groups, paste(savepath,'LinReg_CORT_GOI_groups.xlsx',sep='/'))
write.xlsx(corr_coeffsGOI_groups, paste(savepath,'Corr_coefficients_CORT_GOI_groups.xlsx',sep='/'))
# 