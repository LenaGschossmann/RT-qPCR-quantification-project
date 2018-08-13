#### Check Efficiency of used Primers ####

# 1) fit linear regression line
# 2) compare slope of regression lines of GoI and HKG

library(ggplot2)
library(gridExtra)
library(cowplot)
library(grid)
library(ggsignif)
library(scales)
library(RColorBrewer)
library(xlsx)

#testing
#setwd('C:/Users/lena_/Dropbox/studies/Osnabrück/Universität/Bachelorarbeit_DroBo/experiments/Analysis')
#filepath = 'C:/Users/lena_/Dropbox/studies/Osnabrück/Universität/Bachelorarbeit_DroBo/experiments/PCR/Mito_gene_expression/Efficiency_mito_genes.csv'

#Set Color Palette
colors_Lena = c('firebrick', 'peachpuff4', 'plum3', 'red2', 'darkorchid2', 'coral1', 'bisque3', 'darkmagenta', 'maroon3','tan3', 'royalblue4')


setwd('C:/Users/Gschossmann/Dropbox/studies/Osnabrück/Universität/Bachelorarbeit_DroBo/experiments/Analysis/PCR')
filepath = 'Y:/Lena/Bachelorthesis/Experiments/PCR/Mito_gene_expression/analysis/efficiency calculation/PCR_mito_gene_efficiency.csv'

conc_cDNA = c(3,6,12)
conc_cDNA_log = log10(conc_cDNA)


raw_data = read.csv(filepath, dec='.', sep=';')
raw_data=na.omit(raw_data)
genes=unique(raw_data$Gene)
data = data.frame(matrix(nrow=nrow(raw_data)/2, ncol=3))
colnames(data)=c('Gene', 'cDNA', 'Mean_Cq')
linMod = data.frame(matrix(nrow= nrow(data)/3, ncol=4))
colnames(linMod)=c('Gene', 'Slope', 'Intercept', 'Efficiency')
residual_6ng = data.frame(matrix(ncol=2, nrow=nrow(data)/3))
colnames(residual_6ng)=c('Gene', 'Residual')

#build average of duplicates
cnt=1
for (iGene in genes){
  data$Gene[cnt:(cnt+2)]=iGene
  data$cDNA[cnt]= conc_cDNA_log[1]
  data$Mean_Cq[cnt]= mean(raw_data$Cq[which(raw_data$Gene == iGene & raw_data$cDNA == conc_cDNA[1])])
  data$cDNA[cnt+1]= conc_cDNA_log[2]
  data$Mean_Cq[cnt+1]= mean(raw_data$Cq[which(raw_data$Gene == iGene & raw_data$cDNA == conc_cDNA[2])])
  data$cDNA[cnt+2]= conc_cDNA_log[3]
  data$Mean_Cq[cnt+2]= mean(raw_data$Cq[which(raw_data$Gene == iGene & raw_data$cDNA == conc_cDNA[3])])
  cnt=cnt+3
  }

#linear regression
cnt = 1
for (idx in 1:nrow(linMod)){
  # data[cnt:(cnt+2),1]=as.character(raw_data[idx,1])
  
  tmp=lm(formula=Mean_Cq~cDNA, data=data[cnt:(cnt+2),])
  linMod$Gene[idx]= data$Gene[cnt]
  linMod$Slope[idx]= tmp$coefficients[2]
  linMod$Intercept[idx] = tmp$coefficients[1]
  
  linMod$Efficiency[idx] = 10^(-1/tmp$coefficients[2])-1
  
  residual_6ng$Gene[idx] = as.character(data$Gene[cnt])
  residual_6ng$Residual[idx] = tmp$residuals[2]
  
  cnt=cnt+3
}

ratio_GAPDH_to_HPRT = linMod$Slope[which(linMod=='GAPDH')]/linMod$Slope[which(linMod=='HPRT')]

ratio_genes_to_HKG = data.frame(matrix(nrow=nrow(linMod), ncol=2))
colnames(ratio_genes_to_HKG)=c('Gene', 'Ratio_Gene_to_HKG')
ratio_genes_to_HKG[,1]=linMod$Gene[]
ratio_genes_to_HKG[,2]= mean(linMod$Slope[which(linMod$Gene=='GAPDH')],linMod$Slope[which(linMod$Gene=='HPRT')]) /linMod$Slope[]

ratio_genes_to_GAPDH = data.frame(matrix(nrow=nrow(linMod), ncol=2))
colnames(ratio_genes_to_GAPDH)=c('Gene', 'Ratio_Gene_to_GAPDH')
ratio_genes_to_GAPDH[,1]=linMod$Gene[]
ratio_genes_to_GAPDH[,2]= linMod$Slope[which(linMod$Gene=='GAPDH')]/linMod$Slope[]

ratio_genes_to_HPRT = data.frame(matrix(nrow=nrow(linMod), ncol=2))
colnames(ratio_genes_to_HPRT)=c('Gene', 'Ratio_Gene_to_HPRT')
ratio_genes_to_HPRT[,1]=linMod$Gene[]
ratio_genes_to_HPRT[,2] = linMod$Slope[which(linMod$Gene=='HPRT')]/linMod$Slope[]


ggplot(data= data, aes(x=data[,2], y=data[,3], group=data[,1], color=data[,1]))+
  geom_point()+
  stat_smooth(method='lm', se=FALSE, formula=y~x)+
  scale_x_log10(breaks=trans_breaks('log10', function(x) 10^x), labels= trans_format('log10', math_format(10^.x)))+
  xlab('Amount of cDNA [ng] per 20ul')+ylab('Cq value')+scale_color_manual(name='Genes', values=colors_Lena)

#Output
print('Linear regression')
print(linMod)
print('Residuals: 6ng measurement and regression line')
print(residual_6ng)
print(sprintf('Ratio of GAPDH to HPRT: %f',ratio_GAPDH_to_HPRT))
print('Ratio of genes to HKP_mean:')
print(ratio_genes_to_HKG)
print('Ratio of genes to GAPDH:')
print(ratio_genes_to_GAPDH)
print('Ratio of genes to HPRT:')
print(ratio_genes_to_HPRT)

# Export as excel
write.csv(linMod, 'Y:/Lena/Bachelorthesis/Experiments/PCR/Mito_gene_expression/analysis/efficiency calculation/efficiency_values.csv')

