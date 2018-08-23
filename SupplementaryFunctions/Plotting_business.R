##### Plotting business

library(ggplot2)
library(gridExtra)
library(cowplot)
library(grid)
library(ggsignif)
library(scales)
library(RColorBrewer)

#specify colors
colHR = rgb(255, 0, 0, 255, names = 'HR', max=255)
colIR = rgb(5,190,120, 255, names= 'IR', max=255)
colLR = rgb(87,87,249, 255, names= 'LR', max=255)
cols = c(colHR, colIR, colLR)

#Square-ish
widthPlot = 12
heightPlot = 10

part='part2'
num=3

################################################## Mito gene expression
savepathPlot = 'C:/Users/lena_/Dropbox/studies/Osnabrück/Universität/Bachelorarbeit_DroBo/experiments/PCR/Mito_gene_expression/analysis/Plots/Gene_expression'

plotData = data_corrected

# tmpTtitle = 'Rel. expression according to the Delta-Delta-Cq Method'
tmpTitle = 'Rel. expression according to the GCD formula'

g=ggplot(data=plotData, aes(x=Group, y=Rel_quantity, fill=Group))+
  stat_boxplot(position=position_dodge(.95))+geom_point(aes(color=Group),alpha=.3, show.legend = FALSE)+
  xlab('') + ylab('Rel. gene expression [AU]') + theme(text=element_text(size = 15))+
  theme(axis.title.y=element_text(margin = margin(t=0, r=10, b=0, l=0))) +
  theme(axis.title.x=element_text(margin = margin(t=5, r=0, b=0, l=0)))+
  theme(legend.direction = 'vertical', legend.justification='center')+
  scale_fill_manual(values = cols, name='Line')+ ggtitle(tmpTitle)+
  facet_wrap(~Gene)
ggsave(paste(savepathPlot, paste(paste(part, 'BoxPlot', sep='_'), 'png', sep='.'), sep='/'), plot=g,units='cm', height=heightPlot,width=num*widthPlot, dpi=320)

for(iG in GOI){
  plotData = data_corrected[data_corrected$Gene == iG,]
  tmpTitle='Rel. expression'
  g=ggplot(data=plotData, aes(x=Group, y=Rel_quantity, fill=Group))+
    stat_boxplot(position=position_dodge(.95))+geom_point(aes(color=Group),alpha=.3, show.legend = FALSE)+
    xlab('') + ylab('Rel. gene expression [AU]') + theme(text=element_text(size = 15))+
    theme(axis.title.y=element_text(margin = margin(t=0, r=10, b=0, l=0))) +
    theme(axis.title.x=element_text(margin = margin(t=5, r=0, b=0, l=0)))+
    theme(legend.direction = 'vertical', legend.justification='center')+
    scale_fill_manual(values = cols, name='Line')+ ggtitle(paste(tmpTitle, iG, sep=': '))
  
  ggsave(paste(savepathPlot, paste(paste(part, paste(iG,'BoxPlot', sep='_'), sep='_'), 'png', sep='.'), sep='/'), plot=g,units='cm', height=heightPlot,width=widthPlot,dpi=320)
}



################################################ Method Comparison
savepathPlot = 'C:/Users/lena_/Dropbox/studies/Osnabrück/Universität/Bachelorarbeit_DroBo/experiments/PCR/Mito_gene_expression/analysis/Plots/MethodComparison'
plotData = data_corrected
g=ggplot(data=plotData, aes(x=Group, y=DeltaDelta-Rel_quantity, fill=Group))+
  geom_point()+geom_violin(alpha=0.2, linetype=0)+
  xlab('Genes') + ylab('Rel. expression difference [AU]') + theme(text=element_text(size = 15))+
  theme(axis.title.y=element_text(margin = margin(t=0, r=10, b=0, l=0))) +
  theme(axis.title.x=element_text(margin = margin(t=5, r=0, b=0, l=0)))+
  theme(legend.direction = 'vertical', legend.justification='center')+
  geom_hline(yintercept = 0)+
  scale_fill_manual(values = cols, name='Line')+ ggtitle('Difference in Rel. expression: 2^DeltaDelta Cq - GCD')+
  facet_wrap(~Gene)
ggsave(paste(savepathPlot, paste(paste(part, 'ViolinPlot', sep='_'), 'png', sep='.'), sep='/'), plot=g, units='cm',height=heightPlot,width=widthPlot*num,dpi=320)


for(iG in GOI){
  plotData = data_corrected[data_corrected$Gene == iG,]
  g1 = ggplot(data=plotData, aes(x=Rel_quantity, y=DeltaDelta, color=Group))+geom_point(size=1.3)+
    scale_x_continuous(limits=c(0.5,2))+scale_y_continuous(limits=c(0.5,2))+
    theme(text=element_text(size = 10))+coord_equal()+
    theme(axis.title.y=element_text(margin = margin(t=0, r=10, b=0, l=0))) +
    theme(axis.title.x=element_text(margin = margin(t=10, r=0, b=0, l=0)))+
    theme(axis.text = element_text(size=10))+
    theme(legend.direction = 'vertical', legend.justification='center', text=element_text(size=10), legend.key.height = unit(0.8,'line'))+
    scale_color_manual(values = cols, name='Line')+
    geom_line(aes(x=seq(from=0,to=2,by=2/(nrow(plotData)-1)), y=seq(from=0,to=2,by=2/(nrow(plotData)-1))),color='black')+
    ggtitle(iG)+theme(plot.title=element_text(size = 20))+
    xlab('GCD formula [AU]')+ylab('2^-DeltaDelta Cq formula [AU]')
  
  tmpPlotData = rbind(plotData, plotData)
  tmpPlotData$Method[1:nrow(plotData)] = 'GCD formula'
  tmpPlotData$Method[(nrow(plotData)+1):nrow(tmpPlotData)]= '2^-DeltaDelta Cq'
  tmpPlotData$Rel_quantity[(nrow(plotData)+1):nrow(tmpPlotData)]=plotData$DeltaDelta
  
  g2=ggplot(data=tmpPlotData, aes(x=Group, y=Rel_quantity, fill=Group))+
    stat_boxplot(position=position_dodge(.95))+geom_point(aes(color=Group),alpha=.3, show.legend = FALSE)+
    xlab('Genes') + ylab('Rel. gene expression [AU]') + theme(text=element_text(size = 15))+
    theme(axis.title.y=element_text(margin = margin(t=0, r=10, b=0, l=0))) +
    theme(axis.title.x=element_text(margin = margin(t=5, r=0, b=0, l=0)))+
    theme(legend.direction = 'vertical', legend.justification='center')+
    scale_fill_manual(values = cols, name='Line')+ ggtitle(tmpTitle)+
    facet_wrap(~Method)+ ggtitle(paste('Method Comparison -', iG, sep=' '))+theme(plot.title=element_text(size = 20))
  
  ggsave(paste(savepathPlot, paste(iG, 'png', sep='.'), sep='/'), plot=g1, units='cm', height =10,dpi=320)
  ggsave(paste(savepathPlot, paste(paste(iG, 'BoxPlot', sep='_'), 'png', sep='.'), sep='/'), plot=g2,units='cm', height=heightPlot,dpi=320)
}



# geom_signif(comparison=list(c('HR', 'LR')), annotations = '***', y_position=Y_sign1, tip_length = 0, vjust=0.1)+
# geom_signif(comparison=list(c('HR', 'IR')), annotations = 'NS.', y_position=Y_sign2, tip_length = 0, vjust=-0.1)+
# geom_signif(comparison=list(c('IR', 'LR')), annotations = '***', y_position=Y_sign3, tip_length = 0, vjust=0.1)
# 
# ggplot(data=plotData, aes(x=Group, y=Rel_quantity, fill=Group))+
#   geom_boxplot(aes(ymin = min(data_corrected$Rel_quantity),
#                    lower = quantile(data_corrected$Rel_quantity, 0.25),
#                    middle = median(data_corrected$Rel_quantity),
#                    upper = quantile(data_corrected$Rel_quantity, 0.75),
#                    ymax = max(data_corrected$Rel_quantity)),position=position_dodge(.95))+
#   scale_fill_manual(values = cols, name='Line')+
#   facet_wrap(~Gene)
#
# Barplot
ggplot(data=data_summary, aes(x=data_summary$Gene, y=data_summary$Mean, fill=data_summary$Group))+
  geom_bar(stat='identity', position=position_dodge(.95))+
  geom_errorbar(aes(ymin=data_summary$Mean-data_summary$SEM, ymax=data_summary$Mean+data_summary$SEM),
                width=.8, stat='identity', position=position_dodge(.95))+
  scale_fill_manual(values = cols, name='Line')

######## Further plots

# ##Stability of RefG
# ggplot(data=dataAveraged[dataAveraged$Gene %in% HKG,],
#        aes(x=(1:nrow(dataAveraged[dataAveraged$Gene %in% HKG,])),
#            y=dataAveraged$Cq[dataAveraged$Gene %in% HKG], group=dataAveraged$Gene[dataAveraged$Gene %in% HKG]))+
#   geom_line()

########################################## Efficiency histogram
savepathPlot = 'C:/Users/lena_/Dropbox/studies/Osnabrück/Universität/Bachelorarbeit_DroBo/experiments/PCR/Mito_gene_expression/analysis/Plots/Efficiencies_Dixon_corrected'

plotEffs= Efficiencies[Efficiencies$Outlier == 0,]
g1=ggplot(data=plotEffs, aes(x=Efficiency, fill=Group))+
  geom_density(color='black', alpha=.3)+ scale_fill_manual(values=cols, name='Line')+xlab('Efficiency')
g2=ggplot(data=plotEffs, aes(x=Efficiency, fill=Gene))+
  geom_density(color='black', alpha=.4)+ scale_fill_discrete(name='Gene')+xlab('Efficiency')


ggsave(paste(savepathPlot, paste(paste(part, 'Group', sep='_'), 'png', sep='.'), sep='/'), plot=g1,units='cm', height=10,dpi=320)
ggsave(paste(savepathPlot, paste(paste(part, 'Genes', sep='_'), 'png', sep='.'), sep='/'), plot=g2,units='cm', width=26, height=10,dpi=320)


########################################## Apoptosis Ratio
savepathPlot = 'C:/Users/lena_/Dropbox/studies/Osnabrück/Universität/Bachelorarbeit_DroBo/experiments/PCR/Mito_gene_expression/analysis/Plots/Gene_expression'

g=ggplot(data=data_Ratios, aes(x=data_Ratios$Group, y=data_Ratios$Bcl2_over_Bax, fill=data_Ratios$Group))+
  stat_boxplot(position=position_dodge(.95))+geom_point(aes(color=Group),alpha=.3, show.legend = FALSE)+
  xlab('') + ylab('Ratio Bcl2/Bax') + theme(text=element_text(size = 15))+
  theme(axis.title.y=element_text(margin = margin(t=0, r=10, b=0, l=0))) +
  theme(axis.title.x=element_text(margin = margin(t=5, r=0, b=0, l=0)))+
  theme(legend.direction = 'vertical', legend.justification='center')+
  scale_fill_manual(values = cols, name='Line')+
  ggtitle('Apoptosis Ratio')+theme(plot.title=element_text(size = 20))

ggsave(paste(savepathPlot, paste(paste('ApoptosisRatio', 'BoxPlot', sep='_'), 'png', sep='.'), sep='/'), plot=g, units='cm',height=heightPlot,width=widthPlot,dpi=320)




########################################## DNA damage
savepathPlot = 'C:/Users/lena_/Dropbox/studies/Osnabrück/Universität/Bachelorarbeit_DroBo/experiments/PCR/qPCR_DNA_damage/analysis/Plots/Gene_expression'

#mtDNA copy number
plotData = data_corrected[!is.na(data_corrected$mtCN),]
tmpTitle = 'Mitochondrial DNA Copy Number'
yname = 'mtDNA copy number [AU]'
namePlot = 'mtDNA_cn.png'
g=ggplot(data=plotData, aes(x=Group, y=mtCN, fill=Group))+
  stat_boxplot(position=position_dodge(.95))+geom_point(aes(color=Group),alpha=.3, show.legend = FALSE)+
  xlab('') + ylab(yname) + theme(text=element_text(size = 15))+
  theme(axis.title.y=element_text(margin = margin(t=0, r=10, b=0, l=0))) +
  theme(axis.title.x=element_text(margin = margin(t=5, r=0, b=0, l=0)))+
  theme(legend.direction = 'vertical', legend.justification='center')+
  scale_fill_manual(values = cols, name='Line')+ ggtitle(tmpTitle)

ggsave(paste(savepathPlot, namePlot, sep='/'), plot=g,units='cm', height=heightPlot, width=widthPlot, dpi=320)


#DNA damage
plotData = data_corrected[!is.na(data_corrected$LesionFreq),]
tmpTitle = 'Lesion Frequency per 10kb'
yname = 'lesion frequency'
namePlot = 'Les_freq.png'
g=ggplot(data=plotData, aes(x=Group, y=LesionFreq, fill=Group))+coord_equal()+
  stat_boxplot(position=position_dodge(.95))+geom_point(aes(color=Group),alpha=.3, show.legend = FALSE)+
  xlab('') + ylab(yname) + theme(text=element_text(size = 15))+
  theme(axis.title.y=element_text(margin = margin(t=0, r=10, b=0, l=0))) +
  theme(axis.title.x=element_text(margin = margin(t=5, r=0, b=0, l=0)))+
  theme(legend.direction = 'vertical', legend.justification='center')+
  scale_fill_manual(values = cols, name='Line')+ ggtitle(tmpTitle)+
  facet_wrap(~Gene)

ggsave(paste(savepathPlot, namePlot, sep='/'), plot=g,units='cm', height=heightPlot, width=num*widthPlot,dpi=320)




