library(ggplot2)
library(gridExtra)
library(cowplot)
library(grid)
library(ggsignif)

# library(Rmisc)

setwd('C:/Users/lena_/Documents/Studium_2/VB/Bachelorarbeit/experiments/Analysis')

filepath = 'C:/Users/lena_/Documents/Studium_2/VB/Bachelorarbeit/experiments/Animals_Information/Data Males VR Project 1-Metabolism.csv'

raw_data = read.csv(filepath)

mouseLines = c('HR', 'IR', 'LR')

descriptStats = data.frame(matrix(ncol = length(raw_data), nrow= 3*5))
names = variable.names(raw_data[1:length(raw_data)])
names[1] = 'Line'
names[2] = 'DescriptiveParam'
colnames(descriptStats) = names

cleanedStats = data.frame(matrix(ncol = length(raw_data), nrow= 3*5))
names = variable.names(raw_data[1:length(raw_data)])
names[1] = 'Line'
names[2] = 'DescriptiveParam'
colnames(cleanedStats) = names

#simles descriptive stats table
count = c(1:5)
for (idx in 1:3){
  tmp_data = raw_data[which(raw_data$Line == mouseLines[idx]),]
  #mean
  descriptStats[count[1],3:length(descriptStats)] = apply(tmp_data[,3:length(tmp_data)], 2, function(x) round(mean(x, na.rm = TRUE), digits=3))
  #standard deviation
  descriptStats[count[2],3:length(descriptStats)] = apply(tmp_data[,3:length(tmp_data)], 2, function(x) round(sd(x, na.rm = TRUE), digits=3))
  #SEM
  descriptStats[count[3],3:length(descriptStats)] = apply(tmp_data[,3:length(tmp_data)], 2, function(x) round(sd(x, na.rm = TRUE)/sqrt(sum(!is.na(x))), digits=3))
  #AV+sd
  descriptStats[count[4],3:length(descriptStats)] = apply(tmp_data[,3:length(tmp_data)], 2, function(x) round((mean(x, na.rm = TRUE)+ (2*sd(x, na.rm = TRUE))), digits=3))
  #AV-sd
  descriptStats[count[5],3:length(descriptStats)] = apply(tmp_data[,3:length(tmp_data)], 2, function(x) round((mean(x, na.rm = TRUE)-(2*sd(x, na.rm = TRUE))), digits=3))
  
  
  descriptStats[count[1:5], 2] = c('Average', 'Standard_Deviation', 'SEM', 'AV+sd', 'AV-sd')
  descriptStats[count[1:5],1] = rep(mouseLines[idx], 5)
  count = c((count[5]+1):(count[5]+5))
}


#Outlier detection after Dixon
#specify variables that are essential -> when outlier herein, animal should be removed (completely/from test)
exclusion_criteria = c('AV_BW')
FST_exclusion = c('Dur_floating')
OF_exclusion = c('thigmotaxis')

FST = 18:22 #Variable Indices corresponding to FST
OF = 23:31 #Variable Indices corresponding to OF

outlierTable = raw_data
outlierTable[,3:length(outlierTable)]=0

raw_cleaned = raw_data

elim = 0
elimFST = 0
elimOF = 0

#outer loop -> animals
for (idxA in 1:nrow(raw_data)){
  
  LineID = raw_data$Line[idxA]
  tmp_stats_AVplus = descriptStats[descriptStats$Line == LineID & descriptStats$DescriptiveParam == 'AV+sd',]
  tmp_stats_AVminus = descriptStats[descriptStats$Line == LineID & descriptStats$DescriptiveParam == 'AV-sd',]

  
  #inner loop -> variables
  for (idxV in 3:(length(raw_data))){
    isOutlier = 0
    #if Value is NA -> ignore
    if(is.na(raw_data[idxA, idxV]) == 1){
      break
    }
    #check for each variable if animal shows values outside of mean+/- 2*sd
    if(raw_data[idxA, idxV] > tmp_stats_AVplus[idxV] | raw_data[idxA, idxV] < tmp_stats_AVminus[idxV]){
        outlierTable[idxA, idxV] = 1
        isOutlier = 1
    }
    #check if animal needs to be removed completely
    if(variable.names(raw_data)[idxV] %in% exclusion_criteria & isOutlier ==1){
      raw_cleaned[idxA:(nrow(raw_cleaned)-1),] = raw_cleaned[(idxA+1):nrow(raw_cleaned), ]
      raw_cleaned = raw_cleaned[1:(nrow(raw_cleaned)-1),]
      # raw_cleaned[idxA,] = NA
      elim = elim+1
      break
    }
    #check if animal needs to be removed from FST
    if(variable.names(raw_data)[idxV] %in% FST_exclusion  & isOutlier ==1){
      raw_cleaned[idxA, FST] = NA
      elimFST = elimFST+1
    }
    #check if animal needs to be removed from OF
    if(variable.names(raw_data)[idxV] %in% OF_exclusion  & isOutlier ==1){
      raw_cleaned[idxA, OF] = NA
      elimOF = elimOF+1
    }
  }
}

# nHR = sum(raw_cleaned$Line == 'HR')
# nIR = sum(raw_cleaned$Line == 'IR')
# nLR = sum(raw_cleaned$Line == 'LR')


#### Run descriptive Stats again

count = c(1:5)
for (idx in 1:3){
  tmp_data = raw_cleaned[which(raw_cleaned$Line == mouseLines[idx]),]
  #mean
  cleanedStats[count[1],3:length(descriptStats)] = apply(tmp_data[,3:length(tmp_data)], 2, function(x) round(mean(x, na.rm = TRUE), digits=3))
  #standard deviation
  cleanedStats[count[2],3:length(descriptStats)] = apply(tmp_data[,3:length(tmp_data)], 2, function(x) round(sd(x, na.rm = TRUE), digits=3))
  #SEM
  cleanedStats[count[3],3:length(descriptStats)] = apply(tmp_data[,3:length(tmp_data)], 2, function(x) round(sd(x, na.rm = TRUE)/sqrt(sum(!is.na(x))), digits=3))
  #AV+sd
  cleanedStats[count[4],3:length(descriptStats)] = apply(tmp_data[,3:length(tmp_data)], 2, function(x) round((mean(x, na.rm = TRUE)+ 2*sd(x, na.rm = TRUE)), digits=3))
  #AV-sd
  cleanedStats[count[5],3:length(descriptStats)] = apply(tmp_data[,3:length(tmp_data)], 2, function(x) round((mean(x, na.rm = TRUE)-2*sd(x, na.rm = TRUE)), digits=3))
  
  
  cleanedStats[count[1:5], 2] = c('Average', 'Standard_Deviation', 'SEM', 'AV+sd', 'AV-sd')
  cleanedStats[count[1:5],1] = rep(mouseLines[idx], 5)
  count = c((count[5]+1):(count[5]+5))
}




#### Plots ####

#specify colors
colHR = rgb(255, 0, 0, 255, names = 'HR', max=255)
colIR = rgb(5,190,120, 255, names= 'IR', max=255)
colLR = rgb(87,87,249, 255, names= 'LR', max=255)
cols = c(colHR, colIR, colLR)

##

subset_AV = cleanedStats$DescriptiveParam == 'Average'
subset_AVminsd = cleanedStats$DescriptiveParam == 'AV-sd'
subset_AVplussd = cleanedStats$DescriptiveParam == 'AV+sd'
subset_SEM = cleanedStats$DescriptiveParam == 'SEM'
factorLine = factor(cleanedStats$Line[subset_AV])

# factorLine = factor(raw_cleaned$Line)
# colorLines = 1:length(raw_cleaned$Line)
# colorLines[raw_cleaned$Line == 'HR'] = colHR
# colorLines[raw_cleaned$Line == 'IR'] = colIR
# colorLines[raw_cleaned$Line == 'LR'] = colLR


##Morphometrics

tmp_Plot1 = cleanedStats$AV_food_intake
#significance bars y-value
Y_sign2= max(tmp_Plot1[subset_AV][factorLine!= 'LR']+tmp_Plot1[subset_SEM][factorLine!= 'LR'])+ (min(tmp_Plot1[subset_AV]+tmp_Plot1[subset_SEM])/20)
Y_sign3= max(tmp_Plot1[subset_AV][factorLine!= 'HR']+tmp_Plot1[subset_SEM][factorLine!= 'HR']) + (min(tmp_Plot1[subset_AV]+tmp_Plot1[subset_SEM])/20)
Y_sign1= max(Y_sign2, Y_sign3) + (min(tmp_Plot1[subset_AV]+tmp_Plot1[subset_SEM])/15)
#
  plot1 = ggplot(data = descriptStats[subset_AV,], aes(x=factorLine, y = tmp_Plot1[subset_AV], fill = factorLine))+
    geom_bar(stat="identity", position="dodge", width=0.7 )+
    geom_errorbar(aes(ymin=tmp_Plot1[subset_AV]- tmp_Plot1[subset_SEM],
                      ymax=tmp_Plot1[subset_AV]+tmp_Plot1[subset_SEM]), width=0.2, position=position_dodge(0.9))+
    scale_fill_manual(values = cols, name='Line')+
    xlab('Line') + ylab('Average food intake [g]') + theme(text=element_text(size = 10))+
    theme(axis.title.y=element_text(margin = margin(t=0, r=10, b=0, l=0))) +
    theme(axis.title.x=element_text(margin = margin(t=5, r=0, b=0, l=0)))+
    theme(legend.direction = 'horizontal', legend.justification='center')+
    geom_signif(comparison=list(c('HR', 'LR')), annotations = '***', y_position=Y_sign1, tip_length = 0, vjust=0.1)+
    geom_signif(comparison=list(c('HR', 'IR')), annotations = 'NS.', y_position=Y_sign2, tip_length = 0, vjust=-0.1)+
    geom_signif(comparison=list(c('IR', 'LR')), annotations = '***', y_position=Y_sign3, tip_length = 0, vjust=0.1)
    
   

tmp_Plot2 = cleanedStats$AV_BW
Y_sign2= max(tmp_Plot2[subset_AV][factorLine!= 'LR']+tmp_Plot2[subset_SEM][factorLine!= 'LR'])+ (min(tmp_Plot2[subset_AV]+tmp_Plot2[subset_SEM])/20)
Y_sign3= max(tmp_Plot2[subset_AV][factorLine!= 'HR']+tmp_Plot2[subset_SEM][factorLine!= 'HR']) + (min(tmp_Plot2[subset_AV]+tmp_Plot2[subset_SEM])/20)
Y_sign1= max(Y_sign2, Y_sign3) + (min(tmp_Plot2[subset_AV]+tmp_Plot2[subset_SEM])/10)
#
  plot2 = ggplot(data = descriptStats[subset_AV,], aes(x=factorLine, y = tmp_Plot2[subset_AV], fill = factorLine))+
    geom_bar(stat="identity", position="dodge", width=0.7 )+
    geom_errorbar(aes(ymin=tmp_Plot2[subset_AV]- tmp_Plot2[subset_SEM],
                      ymax=tmp_Plot2[subset_AV]+tmp_Plot2[subset_SEM]), width=0.2, position=position_dodge(0.9))+
    scale_fill_manual(values = cols, guide=FALSE)+
    xlab('Line') + ylab('Average bodyweight [g]') + theme(text=element_text(size = 10))+
    theme(axis.title.y=element_text(margin = margin(t=0, r=10, b=0, l=0))) +
    theme(axis.title.x=element_text(margin = margin(t=5, r=0, b=0, l=0)))+
    geom_signif(comparison=list(c('HR', 'LR')), annotations = 'T', y_position=Y_sign1, tip_length = 0, vjust=-0.1)+
    geom_signif(comparison=list(c('HR', 'IR')), annotations = '***', y_position=Y_sign2, tip_length = 0, vjust=0.1)+
    geom_signif(comparison=list(c('IR', 'LR')), annotations = 'NS.', y_position=Y_sign3, tip_length = 0, vjust=-0.1)
  
  
tmp_Plot3 = cleanedStats$food_intake_per_BW
Y_sign2= max(tmp_Plot3[subset_AV][factorLine!= 'LR']+tmp_Plot3[subset_SEM][factorLine!= 'LR'])+ (min(tmp_Plot3[subset_AV]+tmp_Plot3[subset_SEM])/20)
Y_sign3= max(tmp_Plot3[subset_AV][factorLine!= 'HR']+tmp_Plot3[subset_SEM][factorLine!= 'HR']) + (min(tmp_Plot3[subset_AV]+tmp_Plot3[subset_SEM])/20)
Y_sign1= max(Y_sign2, Y_sign3) + (min(tmp_Plot3[subset_AV]+tmp_Plot3[subset_SEM])/15)
#
  plot3 = ggplot(data = descriptStats[subset_AV,], aes(x=factorLine, y = tmp_Plot3[subset_AV], fill = factorLine))+
    geom_bar(stat="identity", position="dodge", width=0.7 )+
    geom_errorbar(aes(ymin=tmp_Plot3[subset_AV]- tmp_Plot3[subset_SEM],
                      ymax=tmp_Plot3[subset_AV]+tmp_Plot3[subset_SEM]), width=0.2, position=position_dodge(0.9))+
    scale_fill_manual(values = cols, guide=FALSE)+
    xlab('Line') + ylab('Average food intake per bodyweight [g/g]') + theme(text=element_text(size = 10))+
    theme(axis.title.y=element_text(margin = margin(t=0, r=10, b=0, l=0))) +
    theme(axis.title.x=element_text(margin = margin(t=5, r=0, b=0, l=0)))+
    geom_signif(comparison=list(c('HR', 'LR')), annotations = 'NS.', y_position=Y_sign1, tip_length = 0, vjust=-0.1)+
    geom_signif(comparison=list(c('HR', 'IR')), annotations = 'T', y_position=Y_sign2, tip_length = 0, vjust=-0.1)+
    geom_signif(comparison=list(c('IR', 'LR')), annotations = '***', y_position=Y_sign3, tip_length = 0, vjust=0.1)


tmp_Plot4 = cleanedStats$cal_intake_per_24h
Y_sign2= max(tmp_Plot4[subset_AV][factorLine!= 'LR']+tmp_Plot4[subset_SEM][factorLine!= 'LR'])+ (min(tmp_Plot4[subset_AV]+tmp_Plot4[subset_SEM])/20)
Y_sign3= max(tmp_Plot4[subset_AV][factorLine!= 'HR']+tmp_Plot4[subset_SEM][factorLine!= 'HR']) + (min(tmp_Plot4[subset_AV]+tmp_Plot4[subset_SEM])/20)
Y_sign1= max(Y_sign2, Y_sign3) + (min(tmp_Plot4[subset_AV]+tmp_Plot4[subset_SEM])/15)
#
  plot4 = ggplot(data = descriptStats[subset_AV,], aes(x=factorLine, y = tmp_Plot4[subset_AV], fill = factorLine))+
    geom_bar(stat="identity", position="dodge", width=0.7 )+
    geom_errorbar(aes(ymin=tmp_Plot4[subset_AV]- tmp_Plot4[subset_SEM],
                      ymax=tmp_Plot4[subset_AV]+tmp_Plot4[subset_SEM]), width=0.2, position=position_dodge(0.9))+
    scale_fill_manual(values = cols, guide=FALSE)+
    xlab('Line') + ylab('Calorie intake in 24h [kcal]') + theme(text=element_text(size = 10))+
    theme(axis.title.y=element_text(margin = margin(t=0, r=10, b=0, l=0))) +
    theme(axis.title.x=element_text(margin = margin(t=5, r=0, b=0, l=0)))+
    geom_signif(comparison=list(c('HR', 'LR')), annotations = '***', y_position=Y_sign1, tip_length = 0, vjust=0.1)+
    geom_signif(comparison=list(c('HR', 'IR')), annotations = 'NS.', y_position=Y_sign2, tip_length = 0, vjust=-0.1)+
    geom_signif(comparison=list(c('IR', 'LR')), annotations = '***', y_position=Y_sign3, tip_length = 0, vjust=0.1)
  

tmp_Plot5 = cleanedStats$cal_intake_per_BW
Y_sign2= max(tmp_Plot5[subset_AV][factorLine!= 'LR']+tmp_Plot5[subset_SEM][factorLine!= 'LR'])+ (min(tmp_Plot5[subset_AV]+tmp_Plot5[subset_SEM])/20)
Y_sign3= max(tmp_Plot5[subset_AV][factorLine!= 'HR']+tmp_Plot5[subset_SEM][factorLine!= 'HR']) + (min(tmp_Plot5[subset_AV]+tmp_Plot5[subset_SEM])/20)
Y_sign1= max(Y_sign2, Y_sign3) + (min(tmp_Plot5[subset_AV]+tmp_Plot5[subset_SEM])/15)
#
  plot5 = ggplot(data = descriptStats[subset_AV,], aes(x=factorLine, y = tmp_Plot5[subset_AV], fill = factorLine))+
    geom_bar(stat="identity", position="dodge", width=0.7 )+
    geom_errorbar(aes(ymin=tmp_Plot5[subset_AV]- tmp_Plot5[subset_SEM],
                      ymax=tmp_Plot5[subset_AV]+tmp_Plot5[subset_SEM]), width=0.2, position=position_dodge(0.9))+
    scale_fill_manual(values = cols, guide=FALSE)+
    xlab('Line') + ylab('Calorie intake per BW [kcal/g]') + theme(text=element_text(size = 10))+
    theme(axis.title.y=element_text(margin = margin(t=0, r=10, b=0, l=0))) +
    theme(axis.title.x=element_text(margin = margin(t=5, r=0, b=0, l=0)))+
    geom_signif(comparison=list(c('HR', 'LR')), annotations = 'NS.', y_position=Y_sign1, tip_length = 0, vjust=-0.1)+
    geom_signif(comparison=list(c('HR', 'IR')), annotations = '*', y_position=Y_sign2, tip_length = 0, vjust=0.1)+
    geom_signif(comparison=list(c('IR', 'LR')), annotations = '***', y_position=Y_sign3, tip_length = 0, vjust=0.1)
  
  
tmp_Plot6 = cleanedStats$Total.WAT_mg_per_gBW
Y_sign2= max(tmp_Plot6[subset_AV][factorLine!= 'LR']+tmp_Plot6[subset_SEM][factorLine!= 'LR'])+ (min(tmp_Plot6[subset_AV]+tmp_Plot6[subset_SEM])/20)
Y_sign3= max(tmp_Plot6[subset_AV][factorLine!= 'HR']+tmp_Plot6[subset_SEM][factorLine!= 'HR']) + (min(tmp_Plot6[subset_AV]+tmp_Plot6[subset_SEM])/20)
Y_sign1= max(Y_sign2, Y_sign3) + (min(tmp_Plot6[subset_AV]+tmp_Plot6[subset_SEM])/10)
#
  plot6 = ggplot(data = descriptStats[subset_AV,], aes(x=factorLine, y = tmp_Plot6[subset_AV], fill = factorLine))+
    geom_bar(stat="identity", position="dodge", width=0.7 )+
    geom_errorbar(aes(ymin=tmp_Plot6[subset_AV]- tmp_Plot6[subset_SEM],
                      ymax=tmp_Plot6[subset_AV]+tmp_Plot6[subset_SEM]), width=0.2, position=position_dodge(0.9))+
    scale_fill_manual(values = cols, guide=FALSE)+
    xlab('Line') + ylab('Total WAT/BW [mg/g]') + theme(text=element_text(size = 10))+
    theme(axis.title.y=element_text(margin = margin(t=0, r=10, b=0, l=0))) +
    theme(axis.title.x=element_text(margin = margin(t=5, r=0, b=0, l=0)))+
    geom_signif(comparison=list(c('HR', 'LR')), annotations = 'NS.', y_position=Y_sign1, tip_length = 0, vjust=-0.1)+
    geom_signif(comparison=list(c('HR', 'IR')), annotations = 'NS.', y_position=Y_sign2, tip_length = 0, vjust=-0.1)+
    geom_signif(comparison=list(c('IR', 'LR')), annotations = '*', y_position=Y_sign3, tip_length = 0, vjust=0.1)
  
#add legend
  g_legend<-function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)}

  mylegend = g_legend(plot1)

Morphometrics_Plot= grid.arrange(arrangeGrob(plot1+theme(legend.position='none'),plot2,plot3,plot4,plot5,plot6, nrow=2, ncol=3),
                           mylegend, heights=c(10,1),top='Morphometrics')

# grid_arrange_shared_legend(plot1, plot2,plot3, plot4, plot5, plot6, ncol=3, nrow=2, position='bottom', plotTitle='Morphometrics')




##FST
rm(tmp_Plot1, tmp_Plot2, tmp_Plot3, tmp_Plot4, tmp_Plot5, tmp_Plot6)

tmp_Plot1 = cleanedStats$Dur_struggling
Y_sign2= max(tmp_Plot1[subset_AV][factorLine!= 'LR']+tmp_Plot1[subset_SEM][factorLine!= 'LR'])+ (min(tmp_Plot1[subset_AV]+tmp_Plot1[subset_SEM])/20)
Y_sign3= max(tmp_Plot1[subset_AV][factorLine!= 'HR']+tmp_Plot1[subset_SEM][factorLine!= 'HR']) + (min(tmp_Plot1[subset_AV]+tmp_Plot1[subset_SEM])/20)
Y_sign1= max(Y_sign2, Y_sign3) + (min(tmp_Plot1[subset_AV]+tmp_Plot1[subset_SEM])/5)
#
plot1 = ggplot(data = descriptStats[subset_AV,], aes(x=factorLine, y = tmp_Plot1[subset_AV], fill = factorLine))+
  geom_bar(stat="identity", position="dodge", width=0.7 )+
  geom_errorbar(aes(ymin=tmp_Plot1[subset_AV]- tmp_Plot1[subset_SEM],
                    ymax=tmp_Plot1[subset_AV]+tmp_Plot1[subset_SEM]), width=0.2, position=position_dodge(0.9))+
  scale_fill_manual(values = cols, name='Line')+
  xlab('Line') + ylab('Duration struggling [s]') + theme(text=element_text(size = 10))+
  theme(axis.title.y=element_text(margin = margin(t=0, r=10, b=0, l=0))) +
  theme(axis.title.x=element_text(margin = margin(t=5, r=0, b=0, l=0)))+
  theme(legend.direction = 'horizontal', legend.justification='center')+
  geom_signif(comparison=list(c('HR', 'LR')), annotations = '***', y_position=Y_sign1, tip_length = 0, vjust=0.1)+
  geom_signif(comparison=list(c('HR', 'IR')), annotations = '***', y_position=Y_sign2, tip_length = 0, vjust=0.1)+
  geom_signif(comparison=list(c('IR', 'LR')), annotations = 'NS.', y_position=Y_sign3, tip_length = 0, vjust=-0.1)


tmp_Plot2 = cleanedStats$Lat_firstBout
Y_sign2= max(tmp_Plot2[subset_AV][factorLine!= 'LR']+tmp_Plot2[subset_SEM][factorLine!= 'LR'])+ (min(tmp_Plot2[subset_AV]+tmp_Plot2[subset_SEM])/20)
Y_sign3= max(tmp_Plot2[subset_AV][factorLine!= 'HR']+tmp_Plot2[subset_SEM][factorLine!= 'HR']) + (min(tmp_Plot2[subset_AV]+tmp_Plot2[subset_SEM])/20)
Y_sign1= max(Y_sign2, Y_sign3) + (min(tmp_Plot2[subset_AV]+tmp_Plot2[subset_SEM])/5)
#
plot2 = ggplot(data = descriptStats[subset_AV,], aes(x=factorLine, y = tmp_Plot2[subset_AV], fill = factorLine))+
  geom_bar(stat="identity", position="dodge", width=0.7 )+
  geom_errorbar(aes(ymin=tmp_Plot2[subset_AV]- tmp_Plot2[subset_SEM],
                    ymax=tmp_Plot2[subset_AV]+tmp_Plot2[subset_SEM]), width=0.2, position=position_dodge(0.9))+
  scale_fill_manual(values = cols, guide=FALSE)+
  xlab('Line') + ylab('Latency to first bout [s]') + theme(text=element_text(size = 10))+
  theme(axis.title.y=element_text(margin = margin(t=0, r=10, b=0, l=0))) +
  theme(axis.title.x=element_text(margin = margin(t=5, r=0, b=0, l=0)))+
  geom_signif(comparison=list(c('HR', 'LR')), annotations = '***', y_position=Y_sign1, tip_length = 0, vjust=0.1)+
  geom_signif(comparison=list(c('HR', 'IR')), annotations = '***', y_position=Y_sign2, tip_length = 0, vjust=0.1)+
  geom_signif(comparison=list(c('IR', 'LR')), annotations = 'NS.', y_position=Y_sign3, tip_length = 0, vjust=-0.1)


tmp_Plot3 = cleanedStats$Dur_swimming
Y_sign2= max(tmp_Plot3[subset_AV][factorLine!= 'LR']+tmp_Plot3[subset_SEM][factorLine!= 'LR'])+ (min(tmp_Plot3[subset_AV]+tmp_Plot3[subset_SEM])/20)
Y_sign3= max(tmp_Plot3[subset_AV][factorLine!= 'HR']+tmp_Plot3[subset_SEM][factorLine!= 'HR']) + (min(tmp_Plot3[subset_AV]+tmp_Plot3[subset_SEM])/20)
Y_sign1= max(Y_sign2, Y_sign3) + (min(tmp_Plot3[subset_AV]+tmp_Plot3[subset_SEM])/10)
#
plot3 = ggplot(data = descriptStats[subset_AV,], aes(x=factorLine, y = tmp_Plot3[subset_AV], fill = factorLine))+
  geom_bar(stat="identity", position="dodge", width=0.7 )+
  geom_errorbar(aes(ymin=tmp_Plot3[subset_AV]- tmp_Plot3[subset_SEM],
                    ymax=tmp_Plot3[subset_AV]+tmp_Plot3[subset_SEM]), width=0.2, position=position_dodge(0.9))+
  scale_fill_manual(values = cols, guide=FALSE)+
  xlab('Line') + ylab('Duration swimming [s]') + theme(text=element_text(size = 10))+
  theme(axis.title.y=element_text(margin = margin(t=0, r=10, b=0, l=0))) +
  theme(axis.title.x=element_text(margin = margin(t=5, r=0, b=0, l=0)))+
  geom_signif(comparison=list(c('HR', 'LR')), annotations = 'NS.', y_position=Y_sign1, tip_length = 0, vjust=-0.1)+
  geom_signif(comparison=list(c('HR', 'IR')), annotations = 'NS.', y_position=Y_sign2, tip_length = 0, vjust=-0.1)+
  geom_signif(comparison=list(c('IR', 'LR')), annotations = 'NS.', y_position=Y_sign3, tip_length = 0, vjust=-0.1)


tmp_Plot4 = cleanedStats$Dur_floating
Y_sign2= max(tmp_Plot4[subset_AV][factorLine!= 'LR']+tmp_Plot4[subset_SEM][factorLine!= 'LR'])+ (max(tmp_Plot4[subset_AV]+tmp_Plot4[subset_SEM])/20)
Y_sign3= max(tmp_Plot4[subset_AV][factorLine!= 'HR']+tmp_Plot4[subset_SEM][factorLine!= 'HR']) + (max(tmp_Plot4[subset_AV]+tmp_Plot4[subset_SEM])/20)
Y_sign1= max(Y_sign2, Y_sign3) + (max(tmp_Plot4[subset_AV]+tmp_Plot4[subset_SEM])/15)
#
plot4 = ggplot(data = descriptStats[subset_AV,], aes(x=factorLine, y = tmp_Plot4[subset_AV], fill = factorLine))+
  geom_bar(stat="identity", position="dodge", width=0.7 )+
  geom_errorbar(aes(ymin=tmp_Plot4[subset_AV]- tmp_Plot4[subset_SEM],
                    ymax=tmp_Plot4[subset_AV]+tmp_Plot4[subset_SEM]), width=0.2, position=position_dodge(0.9))+
  scale_fill_manual(values = cols, guide=FALSE)+
  xlab('Line') + ylab('Duration floating [s]') + theme(text=element_text(size = 10))+
  theme(axis.title.y=element_text(margin = margin(t=0, r=10, b=0, l=0))) +
  theme(axis.title.x=element_text(margin = margin(t=5, r=0, b=0, l=0)))+
  geom_signif(comparison=list(c('HR', 'LR')), annotations = '***', y_position=Y_sign1, tip_length = 0, vjust=0.1)+
  geom_signif(comparison=list(c('HR', 'IR')), annotations = '***', y_position=Y_sign2, tip_length = 0, vjust=0.1)+
  geom_signif(comparison=list(c('IR', 'LR')), annotations = 'NS.', y_position=Y_sign3, tip_length = 0, vjust=-0.1)


tmp_Plot5 = cleanedStats$Dur_actCoping
Y_sign2= max(tmp_Plot5[subset_AV][factorLine!= 'LR']+tmp_Plot5[subset_SEM][factorLine!= 'LR'])+ (min(tmp_Plot5[subset_AV]+tmp_Plot5[subset_SEM])/20)
Y_sign3= max(tmp_Plot5[subset_AV][factorLine!= 'HR']+tmp_Plot5[subset_SEM][factorLine!= 'HR']) + (min(tmp_Plot5[subset_AV]+tmp_Plot5[subset_SEM])/20)
Y_sign1= max(Y_sign2, Y_sign3) + (min(tmp_Plot5[subset_AV]+tmp_Plot5[subset_SEM])/10)
#
plot5 = ggplot(data = descriptStats[subset_AV,], aes(x=factorLine, y = tmp_Plot5[subset_AV], fill = factorLine))+
  geom_bar(stat="identity", position="dodge", width=0.7 )+
  geom_errorbar(aes(ymin=tmp_Plot5[subset_AV]- tmp_Plot5[subset_SEM],
                    ymax=tmp_Plot5[subset_AV]+tmp_Plot5[subset_SEM]), width=0.2, position=position_dodge(0.9))+
  scale_fill_manual(values = cols, guide=FALSE)+
  xlab('Line') + ylab('Duration of active coping [s]') + theme(text=element_text(size = 10))+
  theme(axis.title.y=element_text(margin = margin(t=0, r=10, b=0, l=0))) +
  theme(axis.title.x=element_text(margin = margin(t=5, r=0, b=0, l=0)))+
  geom_signif(comparison=list(c('HR', 'LR')), annotations = '***', y_position=Y_sign1, tip_length = 0, vjust=0.1)+
  geom_signif(comparison=list(c('HR', 'IR')), annotations = '***', y_position=Y_sign2, tip_length = 0, vjust=0.1)+
  geom_signif(comparison=list(c('IR', 'LR')), annotations = 'NS.', y_position=Y_sign3, tip_length = 0, vjust=-0.1)


#add legend
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend = g_legend(plot1)

FST_Plot= grid.arrange(arrangeGrob(plot1+theme(legend.position='none'),plot3,plot5,plot2,plot4, nrow=2, ncol=3),
                                 mylegend, heights=c(10,1),top='Forced Swimming Test')




##OF
rm(tmp_Plot1, tmp_Plot2, tmp_Plot3, tmp_Plot4, tmp_Plot5, tmp_Plot6)

tmp_Plot1 = cleanedStats$thigmotaxis
Y_sign2= max(tmp_Plot1[subset_AV][factorLine!= 'LR']+tmp_Plot1[subset_SEM][factorLine!= 'LR'])+ (min(tmp_Plot1[subset_AV]+tmp_Plot1[subset_SEM])/20)
Y_sign3= max(tmp_Plot1[subset_AV][factorLine!= 'HR']+tmp_Plot1[subset_SEM][factorLine!= 'HR']) + (min(tmp_Plot1[subset_AV]+tmp_Plot1[subset_SEM])/20)
Y_sign1= max(Y_sign2, Y_sign3) + (min(tmp_Plot1[subset_AV]+tmp_Plot1[subset_SEM])/15)
#
plot1 = ggplot(data = descriptStats[subset_AV,], aes(x=factorLine, y = tmp_Plot1[subset_AV], fill = factorLine))+
  geom_bar(stat="identity", position="dodge", width=0.7 )+
  geom_errorbar(aes(ymin=tmp_Plot1[subset_AV]- tmp_Plot1[subset_SEM],
                    ymax=tmp_Plot1[subset_AV]+tmp_Plot1[subset_SEM]), width=0.2, position=position_dodge(0.9))+
  scale_fill_manual(values = cols, name='Line')+
  xlab('Line') + ylab('Distance Periphery/Total') + theme(text=element_text(size = 10))+
  theme(axis.title.y=element_text(margin = margin(t=0, r=10, b=0, l=0))) +
  theme(axis.title.x=element_text(margin = margin(t=5, r=0, b=0, l=0)))+
  theme(legend.direction = 'horizontal', legend.justification='center')+
  geom_signif(comparison=list(c('HR', 'LR')), annotations = 'NS.', y_position=Y_sign1, tip_length = 0, vjust=-0.1)+
  geom_signif(comparison=list(c('HR', 'IR')), annotations = 'NS.', y_position=Y_sign2, tip_length = 0, vjust=-0.1)+
  geom_signif(comparison=list(c('IR', 'LR')), annotations = 'NS.', y_position=Y_sign3, tip_length = 0, vjust=-0.1)


tmp_Plot2 = cleanedStats$Tot_distance
Y_sign2= max(tmp_Plot2[subset_AV][factorLine!= 'LR']+tmp_Plot2[subset_SEM][factorLine!= 'LR'])+ (min(tmp_Plot2[subset_AV]+tmp_Plot2[subset_SEM])/20)
Y_sign3= max(tmp_Plot2[subset_AV][factorLine!= 'HR']+tmp_Plot2[subset_SEM][factorLine!= 'HR']) + (min(tmp_Plot2[subset_AV]+tmp_Plot2[subset_SEM])/20)
Y_sign1= max(Y_sign2, Y_sign3) + (min(tmp_Plot2[subset_AV]+tmp_Plot2[subset_SEM])/10)
#
plot2 = ggplot(data = descriptStats[subset_AV,], aes(x=factorLine, y = tmp_Plot2[subset_AV], fill = factorLine))+
  geom_bar(stat="identity", position="dodge", width=0.7 )+
  geom_errorbar(aes(ymin=tmp_Plot2[subset_AV]- tmp_Plot2[subset_SEM],
                    ymax=tmp_Plot2[subset_AV]+tmp_Plot2[subset_SEM]), width=0.2, position=position_dodge(0.9))+
  scale_fill_manual(values = cols, guide=FALSE)+
  xlab('Line') + ylab('Total distance [m]') + theme(text=element_text(size = 10))+
  theme(axis.title.y=element_text(margin = margin(t=0, r=10, b=0, l=0))) +
  theme(axis.title.x=element_text(margin = margin(t=5, r=0, b=0, l=0)))+
  geom_signif(comparison=list(c('HR', 'LR')), annotations = '***', y_position=Y_sign1, tip_length = 0, vjust=0.1)+
  geom_signif(comparison=list(c('HR', 'IR')), annotations = 'NS.', y_position=Y_sign2, tip_length = 0, vjust=-0.1)+
  geom_signif(comparison=list(c('IR', 'LR')), annotations = '***', y_position=Y_sign3, tip_length = 0, vjust=0.1)


tmp_Plot3 = cleanedStats$Rearing_freq
Y_sign2= max(tmp_Plot3[subset_AV][factorLine!= 'LR']+tmp_Plot3[subset_SEM][factorLine!= 'LR'])+ (min(tmp_Plot3[subset_AV]+tmp_Plot3[subset_SEM])/20)
Y_sign3= max(tmp_Plot3[subset_AV][factorLine!= 'HR']+tmp_Plot3[subset_SEM][factorLine!= 'HR']) + (min(tmp_Plot3[subset_AV]+tmp_Plot3[subset_SEM])/20)
Y_sign1= max(Y_sign2, Y_sign3) + (min(tmp_Plot3[subset_AV]+tmp_Plot3[subset_SEM])/10)
#
plot3 = ggplot(data = descriptStats[subset_AV,], aes(x=factorLine, y = tmp_Plot3[subset_AV], fill = factorLine))+
  geom_bar(stat="identity", position="dodge", width=0.7 )+
  geom_errorbar(aes(ymin=tmp_Plot3[subset_AV]- tmp_Plot3[subset_SEM],
                    ymax=tmp_Plot3[subset_AV]+tmp_Plot3[subset_SEM]), width=0.2, position=position_dodge(0.9))+
  scale_fill_manual(values = cols, guide=FALSE)+
  xlab('Line') + ylab('Rearing frequenc') + theme(text=element_text(size = 10))+
  theme(axis.title.y=element_text(margin = margin(t=0, r=10, b=0, l=0))) +
  theme(axis.title.x=element_text(margin = margin(t=5, r=0, b=0, l=0)))+
  geom_signif(comparison=list(c('HR', 'LR')), annotations = 'NS.', y_position=Y_sign1, tip_length = 0, vjust=-0.1)+
  geom_signif(comparison=list(c('HR', 'IR')), annotations = 'NS.', y_position=Y_sign2, tip_length = 0, vjust=-0.1)+
  geom_signif(comparison=list(c('IR', 'LR')), annotations = 'NS.', y_position=Y_sign3, tip_length = 0, vjust=-0.1)


tmp_Plot4 = cleanedStats$Dur_selfGrooming
Y_sign2= max(tmp_Plot4[subset_AV][factorLine!= 'LR']+tmp_Plot4[subset_SEM][factorLine!= 'LR'])+ (min(tmp_Plot4[subset_AV]+tmp_Plot4[subset_SEM])/20)
Y_sign3= max(tmp_Plot4[subset_AV][factorLine!= 'HR']+tmp_Plot4[subset_SEM][factorLine!= 'HR']) + (min(tmp_Plot4[subset_AV]+tmp_Plot4[subset_SEM])/20)
Y_sign1= max(Y_sign2, Y_sign3) + (min(tmp_Plot4[subset_AV]+tmp_Plot4[subset_SEM])/5)
#
plot4 = ggplot(data = descriptStats[subset_AV,], aes(x=factorLine, y = tmp_Plot4[subset_AV], fill = factorLine))+
  geom_bar(stat="identity", position="dodge", width=0.7 )+
  geom_errorbar(aes(ymin=tmp_Plot4[subset_AV]- tmp_Plot4[subset_SEM],
                    ymax=tmp_Plot4[subset_AV]+tmp_Plot4[subset_SEM]), width=0.2, position=position_dodge(0.9))+
  scale_fill_manual(values = cols, guide=FALSE)+
  xlab('Line') + ylab('Duration of Self-Grooming [s]') + theme(text=element_text(size = 10))+
  theme(axis.title.y=element_text(margin = margin(t=0, r=10, b=0, l=0))) +
  theme(axis.title.x=element_text(margin = margin(t=5, r=0, b=0, l=0)))+
  geom_signif(comparison=list(c('HR', 'LR')), annotations = 'NS.', y_position=Y_sign1, tip_length = 0, vjust=-0.1)+
  geom_signif(comparison=list(c('HR', 'IR')), annotations = 'NS.', y_position=Y_sign2, tip_length = 0, vjust=-0.1)+
  geom_signif(comparison=list(c('IR', 'LR')), annotations = 'NS.', y_position=Y_sign3, tip_length = 0, vjust=-0.1)


tmp_Plot5 = cleanedStats$Ctr_time
Y_sign2= max(tmp_Plot5[subset_AV][factorLine!= 'LR']+tmp_Plot5[subset_SEM][factorLine!= 'LR'])+ (min(tmp_Plot5[subset_AV]+tmp_Plot5[subset_SEM])/20)
Y_sign3= max(tmp_Plot5[subset_AV][factorLine!= 'HR']+tmp_Plot5[subset_SEM][factorLine!= 'HR']) + (min(tmp_Plot5[subset_AV]+tmp_Plot5[subset_SEM])/20)
Y_sign1= max(Y_sign2, Y_sign3) + (min(tmp_Plot5[subset_AV]+tmp_Plot5[subset_SEM])/10)
#
plot5 = ggplot(data = descriptStats[subset_AV,], aes(x=factorLine, y = tmp_Plot5[subset_AV], fill = factorLine))+
  geom_bar(stat="identity", position="dodge", width=0.7 )+
  geom_errorbar(aes(ymin=tmp_Plot5[subset_AV]- tmp_Plot5[subset_SEM],
                    ymax=tmp_Plot5[subset_AV]+tmp_Plot5[subset_SEM]), width=0.2, position=position_dodge(0.9))+
  scale_fill_manual(values = cols, guide=FALSE)+
  xlab('Line') + ylab('Time in center [s]') + theme(text=element_text(size = 10))+
  theme(axis.title.y=element_text(margin = margin(t=0, r=10, b=0, l=0))) +
  theme(axis.title.x=element_text(margin = margin(t=5, r=0, b=0, l=0)))+
  geom_signif(comparison=list(c('HR', 'LR')), annotations = 'NS.', y_position=Y_sign1, tip_length = 0, vjust=-0.1)+
  geom_signif(comparison=list(c('HR', 'IR')), annotations = 'NS.', y_position=Y_sign2, tip_length = 0, vjust=-0.1)+
  geom_signif(comparison=list(c('IR', 'LR')), annotations = 'NS.', y_position=Y_sign3, tip_length = 0, vjust=-0.1)


tmp_Plot6 = cleanedStats$Periphery_time
Y_sign2= max(tmp_Plot6[subset_AV][factorLine!= 'LR']+tmp_Plot6[subset_SEM][factorLine!= 'LR'])+ (min(tmp_Plot6[subset_AV]+tmp_Plot6[subset_SEM])/20)
Y_sign3= max(tmp_Plot6[subset_AV][factorLine!= 'HR']+tmp_Plot6[subset_SEM][factorLine!= 'HR']) + (min(tmp_Plot6[subset_AV]+tmp_Plot6[subset_SEM])/20)
Y_sign1= max(Y_sign2, Y_sign3) + (min(tmp_Plot6[subset_AV]+tmp_Plot6[subset_SEM])/15)
#
plot6 = ggplot(data = descriptStats[subset_AV,], aes(x=factorLine, y = tmp_Plot6[subset_AV], fill = factorLine))+
  geom_bar(stat="identity", position="dodge", width=0.7 )+
  geom_errorbar(aes(ymin=tmp_Plot6[subset_AV]- tmp_Plot6[subset_SEM],
                    ymax=tmp_Plot6[subset_AV]+tmp_Plot6[subset_SEM]), width=0.2, position=position_dodge(0.9))+
  scale_fill_manual(values = cols, guide=FALSE)+
  xlab('Line') + ylab('Time in periphery [s]') + theme(text=element_text(size = 10))+
  theme(axis.title.y=element_text(margin = margin(t=0, r=10, b=0, l=0))) +
  theme(axis.title.x=element_text(margin = margin(t=5, r=0, b=0, l=0)))
 
  
#add legend
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend = g_legend(plot1)

OF_Plot= grid.arrange(arrangeGrob(plot1+theme(legend.position='none'),plot2,plot3,plot4,plot5,plot6, nrow=2, ncol=3),
                                 mylegend, heights=c(10,1),top='Open Field Test')





#########################################################################################################

##### Stat Tests

#Kruskal Wallis
sink('Kruskal-Wallis-Test.txt')
for (idx in 3:length(raw_cleaned)){
  varNames = variable.names(raw_cleaned)
  print(varNames[idx])
  print(kruskal.test(data=raw_cleaned, raw_cleaned[,idx] ~ raw_cleaned$Line))
}
sink()

#Man-Whitney-U Test
sink('Man-Whitney-U-Test.txt')
for (idx in 3:length(raw_cleaned)){
  varNames = variable.names(raw_cleaned)
  
  for(groups in list(c('HR', 'IR'), c('IR', 'LR'), c('HR', 'LR'))){
    print(paste(varNames[idx], ':', groups[1], 'vs', groups[2], sep=' '))
    print(wilcox.test(data=raw_cleaned[raw_cleaned$Line %in% groups], raw_cleaned[raw_cleaned$Line %in% groups,idx] ~ raw_cleaned$Line[raw_cleaned$Line %in% groups]))
  }
}
sink()


  