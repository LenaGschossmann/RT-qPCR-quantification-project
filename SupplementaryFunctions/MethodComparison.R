##################################################### Method comparison - Linear Regression
# source('C:/Users/Gschossmann/Dropbox/studies/Osnabrück/Universität/Bachelorarbeit_DroBo/experiments/Analysis Code/PCR/qPCR_GeneExpressionAssay/getCORT.R')
source('C:/Users/lena_/Dropbox/studies/Osnabrück/Universität/Bachelorarbeit_DroBo/experiments/Analysis Code/CORT/getCORT.R')

savepath = 'C:/Users/lena_/Dropbox/studies/Osnabrück/Universität/Bachelorarbeit_DroBo/experiments/PCR/Mito_gene_expression/analysis/R_analysis'



#...after all parts have been saved separately:
x1=read.csv(paste(savepath, 'part1_data_corrected.csv', sep='/'))
x2=read.csv(paste(savepath, 'part2_data_corrected.csv', sep='/'))
x3=read.csv(paste(savepath, 'part3_data_corrected.csv', sep='/'))
data_corrected_tot = rbind(x1, rbind(x2, x3))



# Pearson's R for each gene
statsData = data_corrected_tot
tmpGroups = unique(statsData$Group)
tmpGenes = unique(statsData$Gene)
corr_coeffs = data.frame(matrix(ncol=3, nrow=length(tmpGenes)))
colnames(corr_coeffs)[1:3] = c('Gene', 'R', 'p')

LinReg= data.frame(matrix(ncol=5, nrow=length(tmpGenes)*length(tmpGroups)))
colnames(LinReg) = c('Group', 'Gene', 'Slope', 'DevFrom1', 'Intercept')

idx=0
iRow = 1
for(iG in 1:length(tmpGenes)){
  idx = idx+1
  tmpStats = statsData[statsData$Gene == tmpGenes[iG],]
  pearson_corr = rcorr(tmpStats$Rel_quantity,tmpStats$DeltaDelta, type='pearson')
  corr_coeffs$Gene[idx] = as.character(tmpGenes[iG])
  corr_coeffs$R[idx] = pearson_corr$r[1,2]
  corr_coeffs$p[idx] = round(pearson_corr$P[1,2], digits=6)
  #Linear regression
  for(iGr in 1:length(tmpGroups)){
    tmpStatsGroup = tmpStats[tmpStats$Group == tmpGroups[iGr],]
    mod = lm(Rel_quantity ~ DeltaDelta, data=tmpStatsGroup)
    LinReg$Group[iRow] = as.character(tmpGroups[iGr])
    LinReg$Gene[iRow]= as.character(tmpGenes[iG])
    LinReg$Slope[iRow]= mod$coefficients[2]
    LinReg$Intercept[iRow]=mod$coefficients[1]
    LinReg$DevFrom1[iRow]=mod$coefficients[2]-1
    iCol = iRow = iRow+1
  }
}
print(corr_coeffs)
print(LinReg)

#Output
savepath = 'C:/Users/lena_/Dropbox/studies/Osnabrück/Universität/Bachelorarbeit_DroBo/experiments/PCR/Mito_gene_expression/analysis/R_analysis'
write.xlsx(LinReg, paste(savepath,'LinReg_MetComparison.xlsx',sep='/'))
write.xlsx(corr_coeffs, paste(savepath,'Corr_coefficients_MetComparison.xlsx',sep='/'))
# 
# #Kruskal-Wallis-Test
# print(kruskal.test(data=LinReg,LinReg$Intercept ~ as.factor(LinReg$Group)))


########################### Mean Square Deviation from Slope=1-Line




