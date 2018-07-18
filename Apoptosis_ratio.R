
#need to run qPCR_quantification.R first!!

ini_data = data_corrected_Treat


data_Ratios = data.frame(matrix(ncol = 4, nrow=nrow(ini_data)))
colnames(data_Ratios) = c('Plate', 'Animal', 'Group', 'Bcl2_over_Bax')
tmpGenes = c('Bax', 'Bcl2')
cntR = 1
for(iA in unique(ini_data$Animal)){
  if(sum(tmpGenes %in% ini_data$Gene[ini_data$Animal == iA]) == 2 ){
    tmpBax = ini_data[ini_data$Animal == iA & ini_data$Gene == 'Bax', ncol(ini_data)]
    tmpBcl2 = ini_data[ini_data$Animal == iA & ini_data$Gene == 'Bcl2', ncol(ini_data)]
    data_Ratios$Animal[cntR] = iA
    data_Ratios$Plate[cntR] = ini_data$Plate[ini_data$Animal == iA & ini_data$Gene == 'Bax']
    data_Ratios$Group[cntR] = ini_data$Group[ini_data$Animal == iA & ini_data$Gene == 'Bax']
    data_Ratios$Bcl2_over_Bax[cntR] = tmpBcl2 / tmpBax
    cntR = cntR+1
  }
}

data_Ratios = na.omit(data_Ratios)

########################### Plot

#specify colors
colHR = rgb(255, 0, 0, 255, names = 'HR', max=255)
colIR = rgb(5,190,120, 255, names= 'IR', max=255)
colLR = rgb(87,87,249, 255, names= 'LR', max=255)
cols = c(colHR, colIR, colLR)

ggplot(data=data_Ratios, aes(x=data_Ratios$Group, y=data_Ratios$Bcl2_over_Bax, fill=data_Ratios$Group))+
  stat_boxplot(position=position_dodge(.95))+
  xlab('Groups') + ylab('Apoptosis Ratio: Bcl2/Bax') + theme(text=element_text(size = 15))+
  theme(axis.title.y=element_text(margin = margin(t=0, r=10, b=0, l=0))) +
  theme(axis.title.x=element_text(margin = margin(t=5, r=0, b=0, l=0)))+
  theme(legend.direction = 'vertical', legend.justification='center')+
  scale_fill_manual(values = cols, name='Line')

print(kruskal.test(data=data_Ratios,data_Ratios$Bcl2_over_Bax ~ as.factor(data_Ratios$Group)))

for(iGr in Group){
  tmp_data_Ratios = data_Ratios[data_Ratios$Group == Control | data_Ratios$Group == iGr,]
  print(paste(iGr, 'vs.', Control))
  print(wilcox.test(data=tmp_data_Ratios, tmp_data_Ratios$Bcl2_over_Bax ~ as.factor(tmp_data_Ratios$Group)))
}

print(paste(Group[1], 'vs.', Group[2]))
print(wilcox.test(data=data_Ratios[data_Ratios$Group == Group[1] | data_Ratios$Group == Group[2],],
                  data_Ratios$Bcl2_over_Bax[data_Ratios$Group == Group[1] | data_Ratios$Group == Group[2]] ~
                    as.factor(data_Ratios$Group[data_Ratios$Group == Group[1] | data_Ratios$Group == Group[2]])))




