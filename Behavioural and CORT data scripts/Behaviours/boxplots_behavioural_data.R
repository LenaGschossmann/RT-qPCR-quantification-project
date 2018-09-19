
# install.packages(rlang)
# install.packages(ggplot2)
library(ggplot2)
library(gridExtra)
library(cowplot)
library(grid)
library(ggsignif)


#specify colors
colHR = rgb(255, 0, 0, 255, names = 'HR', max=255)
colIR = rgb(5,190,120, 255, names= 'IR', max=255)
colLR = rgb(87,87,249, 255, names= 'LR', max=255)
cols = c(colHR, colIR, colLR)

widthPlot=16
heightPlot=12

###################################################### Boxplots cohorts I-IV

savepath = 'C:/Users/lena_/Dropbox/studies/Osnabrück/Universität/Bachelorarbeit_DroBo/experiments/Behavioural tests/CohortsI-IV'

#change Line code back to HR, IR, LR
raw_cleaned$Mouse_Line[raw_cleaned$Mouse_Line == 1] ='HR'
raw_cleaned$Mouse_Line[raw_cleaned$Mouse_Line == 2] ='IR'
raw_cleaned$Mouse_Line[raw_cleaned$Mouse_Line == 3] ='LR'
# colnames(raw_cleaned)[2] = 'Mouse_Line'

plotData = raw_cleaned

plotTitle='Average Relative Food Intake'
yname='Rel. Food Intake [g/g BW]'
saveVarName = c('food_intake_per_BW')
yIdx = which(variable.names(plotData) == saveVarName)
limY=c(0.05,0.2)

plt = ggplot(data = plotData, aes(x=Mouse_Line, y = plotData[,yIdx], group = Mouse_Line, fill = Mouse_Line))+
  geom_boxplot(outlier.color="black", outlier.stroke=1, outlier.shape=5, position=position_dodge(1))+
  scale_fill_manual(values = cols, name='Line', labels=c('HR', 'IR','LR'))+
  scale_y_continuous(limits = limY)+
  xlab('Line') + ylab(yname) + theme(text=element_text(size = 15))+
  theme(axis.title.y=element_text(margin = margin(t=0, r=10, b=0, l=0))) +
  theme(axis.title.x=element_text(margin = margin(t=5, r=0, b=0, l=0)))+
  theme(legend.direction = 'vertical', legend.justification='center')+
  geom_point(size=3, alpha=0.4)+ggtitle(plotTitle)+
  scale_color_manual(values = cols, name='Line', labels=c('HR', 'IR','LR'), guide=FALSE)+
  geom_hline(yintercept=65, linetype=4)
ggsave(paste(savepath,paste(saveVarName, 'png', sep='.'),sep='/'),plt, units='cm', height=heightPlot, width=widthPlot)

#Subline A
plotData = raw_cleaned[raw_cleaned$Mouse_SuBline == 1,]
plotTitleA = paste(plotTitle, 'Subline A', sep = ' | ')
plt = ggplot(data = plotData, aes(x=Mouse_Line, y = plotData[,yIdx], group = Mouse_Line, fill = Mouse_Line))+
  geom_boxplot(outlier.color="black", outlier.stroke=1, outlier.shape=5, position=position_dodge(1))+
  scale_fill_manual(values = cols, name='Line', labels=c('HR', 'IR','LR'))+
  xlab('Line') + ylab(yname) + theme(text=element_text(size = 15))+
  scale_y_continuous(limits = limY)+
  theme(axis.title.y=element_text(margin = margin(t=0, r=10, b=0, l=0))) +
  theme(axis.title.x=element_text(margin = margin(t=5, r=0, b=0, l=0)))+
  theme(legend.direction = 'vertical', legend.justification='center')+
  geom_point(size=3, alpha=0.4)+ggtitle(plotTitleA)+
  scale_color_manual(values = cols, name='Line', labels=c('HR', 'IR','LR'), guide=FALSE)
ggsave(paste(savepath,paste(paste(saveVarName, 'SublineA', sep='_'), 'png', sep='.'),sep='/'),plt, units='cm', height=heightPlot, width=widthPlot)

#Subline B
plotData = raw_cleaned[raw_cleaned$Mouse_SuBline == 2,]
plotTitleB = paste(plotTitle, 'Subline B', sep = ' | ')
plt = ggplot(data = plotData, aes(x=Mouse_Line, y = plotData[,yIdx], group = Mouse_Line, fill = Mouse_Line))+
  geom_boxplot(outlier.color="black", outlier.stroke=1, outlier.shape=5, position=position_dodge(1))+
  scale_fill_manual(values = cols, name='Line', labels=c('HR', 'IR','LR'))+
  scale_y_continuous(limits = limY)+
  xlab('Line') + ylab(yname) + theme(text=element_text(size = 15))+
  theme(axis.title.y=element_text(margin = margin(t=0, r=10, b=0, l=0))) +
  theme(axis.title.x=element_text(margin = margin(t=5, r=0, b=0, l=0)))+
  theme(legend.direction = 'vertical', legend.justification='center')+
  geom_point(size=3, alpha=0.4)+ggtitle(plotTitleB)+
  scale_color_manual(values = cols, name='Line', labels=c('HR', 'IR','LR'), guide=FALSE)
ggsave(paste(savepath,paste(paste(saveVarName, 'SublineB', sep='_'), 'png', sep='.'),sep='/'),plt, units='cm', height=heightPlot, width=widthPlot)


##################################################### SPT
# 
# plotData = raw_cleaned
# l=nrow(plotData)
# plotData= rbind(plotData, rbind(plotData, rbind(plotData, plotData)))
# 
# plotData$D[1:l]=1
# plotData$D[(l+1):(l+l)]=2
# plotData$D1_TOT[(l+1):(l+l)] = plotData$D2_TOT[1:l] #write 'per_D' values into 1 column (per_D1)
# plotData$D[(l+l+1):(l+l+l)]=3
# plotData$D1_TOT[(l+l+1):(l+l+l)] = plotData$D2_TOT[1:l]
# 
# yname='Average Sucrose Preference per Day'
# plt = ggplot(data = plotData, aes(x=D, y = per_AV, group=Mouse_Line,color = Mouse_Line))+
#   geom_point(aes(color=Mouse_Line),size=3, alpha=0.4, position=position_dodge(0.3))+
#   scale_x_continuous(breaks=c(1,2,3))+scale_y_continuous(limits=c(0,20))+
#   scale_fill_manual(values = cols, name='Line', labels=c('HR', 'IR','LR'))+
#   xlab('Line') + ylab(yname) + theme(text=element_text(size = 15))+
#   theme(axis.title.y=element_text(margin = margin(t=0, r=10, b=0, l=0))) +
#   theme(axis.title.x=element_text(margin = margin(t=5, r=0, b=0, l=0)))+
#   theme(legend.direction = 'vertical', legend.justification='center')+
#   scale_color_manual(values = cols, name='Line', labels=c('HR', 'IR','LR'), guide=FALSE)

#correlate Sucr Preference and Av. Water Intake
plotData = raw_cleaned[raw_cleaned$Mouse_Line == 'HR',]
yLimit=c(0,0.4)
ggplot(data = plotData, aes(x=Av_H_TOT, y = Av_rel_TOT, group = Mouse_Line, color = Mouse_Line))+
  scale_fill_manual(values = cols, name='Line', labels=c('HR', 'IR','LR'))+
  scale_y_continuous(limits=yLimit)+
  xlab('Av. Water') + ylab('Av. Sucrose Preference') + theme(text=element_text(size = 15))+
  theme(axis.title.y=element_text(margin = margin(t=0, r=10, b=0, l=0))) +
  theme(axis.title.x=element_text(margin = margin(t=5, r=0, b=0, l=0)))+
  theme(legend.direction = 'vertical', legend.justification='center')+
  geom_point(size=3, alpha=0.4)+
  scale_color_manual(values = cols, name='Line', labels=c('HR', 'IR','LR'), guide=FALSE)


###################################################### Boxplots cohort VR

savepath = 'C:/Users/lena_/Dropbox/studies/Osnabrück/Universität/Bachelorarbeit_DroBo/experiments/Behavioural tests/Cohort VR/plots'

plotData = raw_cleaned
plotTitle='Total distance'
yname = 'Total Distance [m]'
saveVarName = 'Tot_distance'
yIdx = which(variable.names(plotData) == saveVarName)
limY = c(0,75)

plt = ggplot(data = plotData, aes(x=Line, y = plotData[,yIdx], group = Line, fill = Line))+
  geom_boxplot(outlier.color="black", outlier.stroke=1, outlier.shape=5, position=position_dodge(1))+
  scale_fill_manual(values = cols, name='Line', labels=c('HR', 'IR','LR'))+
  scale_y_continuous(limits = limY)+
  xlab('Line') + ylab(yname) + theme(text=element_text(size = 15))+
  theme(axis.title.y=element_text(margin = margin(t=0, r=10, b=0, l=0))) +
  theme(axis.title.x=element_text(margin = margin(t=5, r=0, b=0, l=0)))+
  theme(legend.direction = 'vertical', legend.justification='center')+
  geom_point(size=3, alpha=0.4)+ggtitle(plotTitle)+
  scale_color_manual(values = cols, name='Line', labels=c('HR', 'IR','LR'), guide=FALSE)
ggsave(paste(savepath,paste(saveVarName, 'png', sep='.'),sep='/'),plt, units='cm', height=heightPlot, width=widthPlot)


############################################################## CORT
plotData = raw_cleaned
plotCORT = cort

#write CORT data into other data.frame
for(iR in 1:nrow(plotData)){
  plotData$CORT[iR] = plotCORT$CORT[plotCORT$ID == plotData$ID[iR]]
}

plotData = plotData[order(plotData$CORT, decreasing = FALSE),]

ggplot(plotData, aes(x=CORT, y=AV_BW, color=Line))+geom_point()


# limits Rel Food: 0.1,0.15
# BW: 23, 55
# tot food: 4,10
#
#
#

