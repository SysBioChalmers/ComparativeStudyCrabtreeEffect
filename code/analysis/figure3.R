# Load required libraries
library(tidyverse)
library(RColorBrewer)
library(ggpubr)

# Load proteomics data
# set working directory
setwd('../../data/experimentalData')

# Read proteomics data data in g/gDW
protData_abs_Sce <- read.delim('proteomicsData_g_per_gDW_Scer.txt',stringsAsFactors = FALSE)
protData_abs_Spo <- read.delim('proteomicsData_g_per_gDW_Spombe.txt',stringsAsFactors = FALSE)
protData_abs_Sstip <- read.delim('proteomicsData_g_per_gDW_Pstip.txt',stringsAsFactors = FALSE)
protData_abs_Kma <- read.delim('proteomicsData_g_per_gDW_Kmarx.txt',stringsAsFactors = FALSE)

# Extract data for the complexes of the electron transport chain
# read file with GO term annotation
GOannotation_Sce <- read.delim('../GOterms/Sce_GOterms.txt')
GOannotation_Spo <- read.delim('../GOterms/Spo_GOterms.txt')
GOannotation_Sstip <- read.delim('../GOterms/Pst_GOterms.txt')
GOannotation_Kma <- read.delim('../GOterms/Kma_GOterms.txt')

# Define required functions
# map proteins to GO terms (since multiple GO terms per gene, map GO terms to proteins in proteomics data)
map2GO <- function(proteomicsData,GOannotation,GOterms){
  if(length(GOterms)>1){
    #GOdata <- data.frame(matrix(ncol = ncol(proteomicsData), nrow = 0))
    proteinsInGOterms <- data.frame(matrix(ncol = 3, nrow = 0))
    for(GOterm in GOterms){
      proteinsMapped <- GOannotation[grep(GOterm,GOannotation$GOterm),1:3]
      proteinsInGOterms <- rbind(proteinsInGOterms,proteinsMapped)
    }
  }else{
    proteinsInGOterms <- GOannotation[grep(GOterms,GOannotation$GOterm),1:3]
  }
  #proteinsInGOterms <- data.frame(proteinsInGOterms)
  colnames(proteinsInGOterms) <- c('Accession','GOterm','GOtermName')
  GOdata <- merge(proteomicsData,proteinsInGOterms)
  # Remove duplicated entries
  GOdata <- GOdata[!duplicated(GOdata$Accession),]
  return(GOdata)
}

mapGenesToProtData <- function(geneList,protData){
  data <- data.frame(matrix(ncol = 7,nrow = 0))
  for(gene in geneList){
    geneData = protData[grep(gene,protData$GeneName),]
    if(!nrow(geneData)==0){
      data <- rbind(data,geneData)
    }
  }
  return(data)
}

# fig 3A

# Read usage data
usages_Kma <- read.delim('../../results/enzymeUsage/enzymeUsage_Kmarx.txt',stringsAsFactors = F)
usages_Sce <- read.delim('../../results/enzymeUsage/enzymeUsage_Scer.txt',stringsAsFactors = F)

# Extract usages (exclude proteins with 0 usage)
capUse_Kma <- usages_Kma[usages_Kma$capUsage != 0,1:4]
capUse_Kma[,4] <- capUse_Kma[,4]*100
capUse_Sce <- usages_Sce[usages_Sce$capUsage != 0,1:4]
capUse_Sce[,4] <- capUse_Sce[,4]*100

# Read GO term annotation for selected terms
GO_Kma <- read.delim('../GOterms/GOannotation_Kma.txt',stringsAsFactors = F)
GO_Sce <- read.delim('../GOterms/GOannotation_Sce.txt',stringsAsFactors = F)

# Keep data from relevant GO terms for Kma
capUse_Kma <- capUse_Kma[capUse_Kma$protID %in% GO_Kma$Entry,]
idx <- match(capUse_Kma$protID,GO_Kma$Entry)
capUse_Kma$GOterm <- GO_Kma$System[idx]
capUse_Kma$protName <- GO_Kma$Gene.names...primary..[idx]
capUse_Kma <- capUse_Kma %>% mutate_if(is.numeric,round,digits = 3)
# Save to file
write.table(capUse_Kma,'capUsage_Kma.txt',sep = '\t',row.names = F)

# Keep relevant GO terms for Sce
capUse_Sce <- capUse_Sce[capUse_Sce$protID %in% GO_Sce$Entry,]
idx <- match(capUse_Sce$protID,GO_Sce$Entry)
capUse_Sce$GOterm <- GO_Sce$System[idx]
capUse_Sce <- capUse_Sce %>% mutate_if(is.numeric,round,digits = 3)
# Save to file
write.table(capUse_Sce,'capUsage_Sce.txt',sep = '\t',row.names = F)

# Combine usages into one data frame for plotting
Kma <- capUse_Kma[,4:5]
Kma$GOterm <- factor(Kma$GOterm,levels = c('Glycolysis','TCA cycle','ETC',
                                           'ATP synthase','PP shunt',
                                           'Amino acid metabolism','Fatty acid metabolism',
                                           'Nucleotide biosynthesis','Sterol biosynthesis'))
Kma$Organism <- c(rep('Kma',nrow(Kma)))

Sce <- capUse_Sce[,4:5]
Sce$GOterm <- factor(Sce$GOterm,levels = c('Glycolysis','TCA cycle','ETC',
                                           'ATP synthase','PP shunt',
                                           'Amino acid metabolism','Fatty acid metabolism',
                                           'Nucleotide biosynthesis','Sterol biosynthesis'))
Sce$Organism <- c(rep('Sce',nrow(Sce)))

# Combine data for both species
combined <- rbind(Sce,Kma)
combined$Organism <- factor(combined$Organism,levels = c('Kma','Sce'))
dataToPlot <- combined[combined$GOterm %in% c('Glycolysis','PP shunt','TCA cycle','ETC','ATP synthase'),]
fig3A <- ggplot(dataToPlot,aes(x = Organism, y = capUsage, fill = Organism)) +
  geom_boxplot(alpha=0.9) +
  facet_grid(. ~ GOterm) +
  labs(x = '', y = 'Capacity usage (%)') +
  theme_classic() +
  theme(text = element_text(size = 8, color = 'black', face = 'bold'),
        axis.text = element_text(size = 8, color = 'black', face = 'bold'),
        axis.title.y.right = element_blank(),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        panel.spacing = unit(3, "mm"),
        #line = element_line(size = 0.15, color = 'black'),
        #axis.line = element_line(size = 0.15, color = 'black'),
        panel.border = element_rect(size = 1, color = 'black',fill = NA),
        legend.position = 'none') +
  scale_y_continuous(sec.axis = dup_axis(breaks = 0))
fig3A
#ggsave('../../results/fig3A.pdf',
 #      fig3A, device = 'pdf',height = 7, width = 7, units = 'cm',dpi = 300)

# fig3B

# Define functions used
map2GO <- function(proteomicsData,GOannotation,GOterms){
  if(length(GOterms)>1){
    #GOdata <- data.frame(matrix(ncol = ncol(proteomicsData), nrow = 0))
    proteinsInGOterms <- data.frame(matrix(ncol = 3, nrow = 0))
    for(GOterm in GOterms){
      proteinsMapped <- GOannotation[grep(GOterm,GOannotation$GOterm),1:3]
      proteinsInGOterms <- rbind(proteinsInGOterms,proteinsMapped)
    }
  }else{
    proteinsInGOterms <- GOannotation[grep(GOterms,GOannotation$GOterm),1:3]
  }
  #proteinsInGOterms <- data.frame(proteinsInGOterms)
  colnames(proteinsInGOterms) <- c('Accession','GOterm','GOtermName')
  GOdata <- merge(proteomicsData,proteinsInGOterms)
  # Remove duplicated entries
  GOdata <- GOdata[!duplicated(GOdata$Accession),]
  return(GOdata)
}

# read proteomics data (allocation [g/g protein])
protData_Sce <- read.delim('proteomicsData_massPercent_Scer.txt',stringsAsFactors = FALSE)
protData_Spo <- read.delim('proteomicsData_massPercent_Spombe.txt',stringsAsFactors = FALSE)
protData_Sstip <- read.delim('proteomicsData_massPercent_Pstip.txt',stringsAsFactors = FALSE)
protData_Kma <- read.delim('proteomicsData_massPercent_Kmarx.txt',stringsAsFactors = FALSE)

# read file with GO term annotation
GOannotation_Sce <- read.delim('../GOterms/Sce_GOterms.txt')
GOannotation_Spo <- read.delim('../GOterms/Spo_GOterms.txt')
GOannotation_Sstip <- read.delim('../GOterms/Pst_GOterms.txt')
GOannotation_Kma <- read.delim('../GOterms/Kma_GOterms.txt')

# S. cerevisiae
TCAcycle_Sce <- map2GO(protData_Sce,GOannotation_Sce,'GO:0006099')
TCAcycle_Sce$GOtermName <- 'TCA cycle'
# Remove ICL1, which is part of the glyoxylate cycle
TCAcycle_Sce <- TCAcycle_Sce[-grep('ICL1',TCAcycle_Sce$PrimaryGeneName),]
TCAcycle_Sce <- data.frame(cbind(t(colSums(TCAcycle_Sce[,5:7])),TCAcycle_Sce[1,9]))
colnames(TCAcycle_Sce) <- c('repl1','repl2','repl3','Process')
ETC_Sce <- map2GO(protData_Sce,GOannotation_Sce,c('GO:0005750','GO:0005751',
                                                  'GO:0004129','GO:0008137',
                                                  'GO:0003954','GO:0006122',
                                                  'GO:0006123'))
ETC_Sce <- ETC_Sce[!duplicated(ETC_Sce$Accession),]
# Remove NAR1 and COQ11, which are not actually part of the ETC
ETC_Sce <- ETC_Sce[-grep('NAR1',ETC_Sce$PrimaryGeneName),]
ETC_Sce <- ETC_Sce[-grep('COQ11',ETC_Sce$PrimaryGeneName),]
ETC_Sce$GOtermName <- 'ETC'
ETC_Sce <- data.frame(cbind(t(colSums(ETC_Sce[,5:7])),ETC_Sce[1,9]))
colnames(ETC_Sce) <- c('repl1','repl2','repl3','Process')
ATPsynth_Sce <- map2GO(protData_Sce,GOannotation_Sce,'GO:0015986')
ATPsynth_Sce$GOtermName <- 'ATP synthase'
ATPsynth_Sce <- data.frame(cbind(t(colSums(ATPsynth_Sce[,5:7])),ATPsynth_Sce[1,9]))
colnames(ATPsynth_Sce) <- c('repl1','repl2','repl3','Process')
oxphos_Sce <- rbind(TCAcycle_Sce,ETC_Sce,ATPsynth_Sce)
oxphos_Sce$Organism <- c(rep('Sce',3))
oxphos_Sce[,1:3] <- sapply(oxphos_Sce[,1:3],as.numeric)
oxphos_Sce$meanAllocation <- rowMeans(oxphos_Sce[,1:3])
oxphos_Sce$sd <- apply(oxphos_Sce[,1:3], 1, sd)

# S. pombe
TCAcycle_Spo <- map2GO(protData_Spo,GOannotation_Spo,'GO:0006099')
TCAcycle_Spo$GOtermName <- 'TCA cycle'
TCAcycle_Spo <- data.frame(cbind(t(colSums(TCAcycle_Spo[,5:7])),TCAcycle_Spo[1,9]))
colnames(TCAcycle_Spo) <- c('repl1','repl2','repl3','Process')
ETC_Spo <- map2GO(protData_Spo,GOannotation_Spo,c('GO:0005750','GO:0005751',
                                                  'GO:0004129','GO:0008137',
                                                  'GO:0003954','GO:0006122',
                                                  'GO:0006123'))
ETC_Spo <- ETC_Spo[!duplicated(ETC_Spo$Accession),]
ETC_Spo <- data.frame(cbind(t(colSums(ETC_Spo[,5:7])),ETC_Spo[1,9]))
colnames(ETC_Spo) <- c('repl1','repl2','repl3','Process')
ETC_Spo$Process <- 'ETC'
ATPsynth_Spo <- map2GO(protData_Spo,GOannotation_Spo,'GO:0015986')
ATPsynth_Spo$GOtermName <- 'ATP synthase'
ATPsynth_Spo <- data.frame(cbind(t(colSums(ATPsynth_Spo[,5:7])),ATPsynth_Spo[1,9]))
colnames(ATPsynth_Spo) <- c('repl1','repl2','repl3','Process')
oxphos_Spo <- rbind(TCAcycle_Spo,ETC_Spo,ATPsynth_Spo)
oxphos_Spo$Organism <- c(rep('Spo',3))
oxphos_Spo[,1:3] <- sapply(oxphos_Spo[,1:3],as.numeric)
oxphos_Spo$meanAllocation <- rowMeans(oxphos_Spo[,1:3])
oxphos_Spo$sd <- apply(oxphos_Spo[,1:3], 1, sd)

# S. stipitis
TCAcycle_Sstip <- map2GO(protData_Sstip,GOannotation_Sstip,'GO:0006099')
TCAcycle_Sstip$GOtermName <- 'TCA cycle'
TCAcycle_Sstip <- data.frame(cbind(t(colSums(TCAcycle_Sstip[,5:7])),TCAcycle_Sstip[1,9]))
colnames(TCAcycle_Sstip) <- c('repl1','repl2','repl3','Process')
ETC_Sstip <- map2GO(protData_Sstip,GOannotation_Sstip,c('GO:0005750','GO:0005751',
                                                        'GO:0004129','GO:0008137',
                                                        'GO:0003954','GO:0006122',
                                                        'GO:0006123','GO:0042773'))
ETC_Sstip <- ETC_Sstip[!duplicated(ETC_Sstip$Accession),]
ETC_Sstip$GOtermName <- 'ETC'
ETC_Sstip <- data.frame(cbind(t(colSums(ETC_Sstip[,5:7])),ETC_Sstip[1,9]))
colnames(ETC_Sstip) <- c('repl1','repl2','repl3','Process')
ATPsynth_Sstip <- map2GO(protData_Sstip,GOannotation_Sstip,'GO:0015986')
ATPsynth_Sstip$GOtermName <- 'ATP synthase'
ATPsynth_Sstip <- data.frame(cbind(t(colSums(ATPsynth_Sstip[,5:7])),ATPsynth_Sstip[1,9]))
colnames(ATPsynth_Sstip) <- c('repl1','repl2','repl3','Process')
oxphos_Sstip <- rbind(TCAcycle_Sstip,ETC_Sstip,ATPsynth_Sstip)
oxphos_Sstip$Organism <- c(rep('Sstip',3))
oxphos_Sstip[,1:3] <- sapply(oxphos_Sstip[,1:3],as.numeric)
oxphos_Sstip$meanAllocation <- rowMeans(oxphos_Sstip[,1:3])
oxphos_Sstip$sd <- apply(oxphos_Sstip[,1:3], 1, sd)

# K. marxianus
TCAcycle_Kma <- map2GO(protData_Kma,GOannotation_Kma,'GO:0006099')
TCAcycle_Kma$GOtermName <- 'TCA cycle'
TCAcycle_Kma <- data.frame(cbind(t(colSums(TCAcycle_Kma[,5:7])),TCAcycle_Kma[1,9]))
colnames(TCAcycle_Kma) <- c('repl1','repl2','repl3','Process')
ETC_Kma <- map2GO(protData_Kma,GOannotation_Kma,c('GO:0005750','GO:0005751',
                                                  'GO:0004129','GO:0008137',
                                                  'GO:0003954','GO:0006122',
                                                  'GO:0006123','GO:0042773'))
ETC_Kma <- ETC_Kma[!duplicated(ETC_Kma$Accession),]
ETC_Kma$GOtermName <- 'ETC'
ETC_Kma <- data.frame(cbind(t(colSums(ETC_Kma[,5:7])),ETC_Kma[1,9]))
colnames(ETC_Kma) <- c('repl1','repl2','repl3','Process')
ATPsynth_Kma <- map2GO(protData_Kma,GOannotation_Kma,'GO:0015986')
ATPsynth_Kma$GOtermName <- 'ATP synthase'
ATPsynth_Kma <- data.frame(cbind(t(colSums(ATPsynth_Kma[,5:7])),ATPsynth_Kma[1,9]))
colnames(ATPsynth_Kma) <- c('repl1','repl2','repl3','Process')
oxphos_Kma <- rbind(TCAcycle_Kma,ETC_Kma,ATPsynth_Kma)
oxphos_Kma$Organism <- c(rep('Kma',3))
oxphos_Kma[,1:3] <- sapply(oxphos_Kma[,1:3],as.numeric)
oxphos_Kma$meanAllocation <- rowMeans(oxphos_Kma[,1:3])
oxphos_Kma$sd <- apply(oxphos_Kma[,1:3], 1, sd)

# Combine data
oxphos <- rbind(oxphos_Sce,oxphos_Spo,oxphos_Sstip,oxphos_Kma)
#oxphos$Organism <- factor(oxphos$Organism,levels = c('Sce','Spo','Kma','Sstip'))
oxphos$rGlu <- c(rep(18.8,3),rep(8.9,3),rep(4.5,3),rep(5.6,3))
ATPS <- oxphos[grepl('ATP synthase',oxphos$Process),]
df1 <- melt(ATPS[,c(1:4,8)],id.vars = c('Process','rGlu'),measure.vars = c('repl1','repl2','repl3'))
ETC <- oxphos[grepl('ETC',oxphos$Process),]
df2 <- melt(ETC[,c(1:4,8)],id.vars = c('Process','rGlu'),measure.vars = c('repl1','repl2','repl3'))
TCA <- oxphos[grepl('TCA cycle',oxphos$Process),]
df3 <- melt(TCA[,c(1:4,8)],id.vars = c('Process','rGlu'),measure.vars = c('repl1','repl2','repl3'))

# Plot figure
fig3B <- ggplot() +
  geom_point(data=ATPS,aes(x=rGlu,y=meanAllocation*100),color='#9970AB',shape=22,fill='#9970AB') +
  geom_errorbar(data=ATPS,aes(x=rGlu,ymin=(meanAllocation*100)-sd,ymax=(meanAllocation*100)+sd),
                width=.2,position=position_dodge(.9)) +
  geom_jitter(data=df1,aes(x=rGlu,y=value*100),size=2,
              position = position_jitter(0.2),color='#9970AB') +
  geom_point(data=ETC,aes(x=rGlu,y=meanAllocation*100),color='#6BAED6',shape=22,fill='#6BAED6') +
  geom_errorbar(data=ETC,aes(x=rGlu,ymin=(meanAllocation*100)-sd,ymax=(meanAllocation*100)+sd),
                width=.2,position=position_dodge(.9)) +
  geom_jitter(data=df2,aes(x=rGlu,y=value*100),size=2,
              position = position_jitter(0.2),color='#6BAED6') +
  geom_point(data=TCA,aes(x=rGlu,y=meanAllocation*100),color='#FD8D3C',shape=22,fill='#FD8D3C') +
  geom_errorbar(data=TCA,aes(x=rGlu,ymin=(meanAllocation*100)-sd,ymax=(meanAllocation*100)+sd),
                width=.2,position=position_dodge(.9)) +
  geom_jitter(data=df3,aes(x=rGlu,y=value*100),size=2,
              position = position_jitter(0.2),color='#FD8D3C') +
  theme_classic() +
  theme(axis.title = element_text(size = 8,face = 'bold',color = 'black'),
        axis.text = element_text(size = 8,face = 'bold',color = 'black'),
        axis.text.x = element_text(size = 8,face = 'bold',color = 'black'),
        legend.title = element_blank(),
        legend.text = element_text(size = 8,face = 'bold',color = 'black'),
        legend.position = c(0.7,0.9)) +
  labs(x='rGlu (mmol/gDW/h)',y='Proteome allocation (%)') +
  #scale_color_manual(values = c('#9970AB','#5AAE61','#6BAED6','#A6761D')) +
  #scale_fill_manual(values = c('#9970AB','#5AAE61','#6BAED6','#A6761D')) +
  scale_y_continuous(limits = c(0,3.25),breaks = seq(0,3,0.5),expand = c(0,0)) +
  scale_x_continuous(breaks=c(4.5,5.6,8.9,18.8))
fig3B
# Save plot
#ggsave('../../results/fig3B.pdf',
 #      fig3B,width=7,height = 7,units = 'cm',device = 'pdf',dpi=300)

# fig 3C
# NADH dehydrogenases
# S. cerevisiae
NADH_dehydrogenase_abs_Sce <- map2GO(protData_abs_Sce,GOannotation_Sce,c('GO:0008137','GO:0003954'))
NADH_dehydrogenase_abs_Sce <- NADH_dehydrogenase_abs_Sce[which(NADH_dehydrogenase_abs_Sce$PrimaryGeneName == 'NDE1' | NADH_dehydrogenase_abs_Sce$PrimaryGeneName == 'NDI1'),]
NADH_dehydrogenase_abs_Sce$complex <- c(rep('NADH',nrow(NADH_dehydrogenase_abs_Sce)))
NADH_dehydrogenase_abs_Sce$PrimaryGeneName <- paste(NADH_dehydrogenase_abs_Sce$PrimaryGeneName,'_Sce',sep = '')
# S. pombe
NADH_dehydrogenase_abs_Spo <- map2GO(protData_abs_Spo,GOannotation_Spo,c('GO:0008137','GO:0003954'))
NADH_dehydrogenase_abs_Spo <- NADH_dehydrogenase_abs_Spo[1:2,]
NADH_dehydrogenase_abs_Spo$complex <- c(rep('NADH',nrow(NADH_dehydrogenase_abs_Spo)))
NADH_dehydrogenase_abs_Spo$PrimaryGeneName <- c('nde1','ndi1')
# K. marxianus
NADH_dehydrogenase_abs_Kma <- map2GO(protData_abs_Kma,GOannotation_Kma,c('GO:0008137','GO:0003954'))
NADH_dehydrogenase_abs_Kma <- NADH_dehydrogenase_abs_Kma[1:2,]
NADH_dehydrogenase_abs_Kma$complex <- c(rep('NADH',nrow(NADH_dehydrogenase_abs_Kma)))
NADH_dehydrogenase_abs_Kma$PrimaryGeneName <- paste(NADH_dehydrogenase_abs_Kma$PrimaryGeneName,'_Kma',sep = '')
# S. stipitis
NADH_dehydrogenase_abs_Sstip <- map2GO(protData_abs_Sstip,GOannotation_Sstip,c('GO:0008137','GO:0003954'))
NADH_dehydrogenase_abs_Sstip <- NADH_dehydrogenase_abs_Sstip[grep('NDE1',NADH_dehydrogenase_abs_Sstip$PrimaryGeneName),]
complexIgenes <- c('PICST_85822','PICST_68160','PICST_74163','PICST_62206','PICST_81272',
                   'PICST_58506','PICST_75518','PICST_88630','PICST_31688','PICST_76559',
                   'PICST_65836','PICST_46630','PICST_84570','PICST_82538','PICST_76018',
                   'PICST_44890','PICST_76226','PICST_63166','PICST_63376','PICST_80638',
                   'PICST_90638','PICST_79236')
complexI_abs_Sstip <- mapGenesToProtData(complexIgenes,protData_abs_Sstip)
NADH_dehydrogenase_abs_Sstip <- rbind(NADH_dehydrogenase_abs_Sstip[,1:7],complexI_abs_Sstip)
NADH_dehydrogenase_abs_Sstip$complex <- c(rep('NADH',nrow(NADH_dehydrogenase_abs_Sstip)))

# Complex II
CII_abs_Sce <- map2GO(protData_abs_Sce,GOannotation_Sce,c('GO:0005749'))
CII_abs_Sce$complex <- c(rep('CII',nrow(CII_abs_Sce)))

CII_abs_Spo <- map2GO(protData_abs_Spo,GOannotation_Spo,c('GO:0005749'))

CII_abs_Kma <- map2GO(protData_abs_Kma,GOannotation_Kma,c('GO:0005749'))
CII_abs_Kma <- CII_abs_Kma[which(CII_abs_Kma$PrimaryGeneName != 'FMP21' & CII_abs_Kma$PrimaryGeneName != 'TIM18'),]

CII_abs_Sstip <- map2GO(protData_abs_Sstip,GOannotation_Sstip,c('GO:0005749'))
CII_abs_Sstip <- CII_abs_Sstip[-grep('A3LUM4',CII_abs_Sstip$Accession),]

# Complex III
CIII_abs_Sce <- map2GO(protData_abs_Sce,GOannotation_Sce,c('GO:0005750','GO:0006122'))
CIII_abs_Sce <- CIII_abs_Sce[-grep('MAS1',CIII_abs_Sce$PrimaryGeneName),]

CIII_abs_Spo <- map2GO(protData_abs_Spo,GOannotation_Spo,c('GO:0005750','GO:0006122'))

CIII_abs_Kma <- map2GO(protData_abs_Kma,GOannotation_Kma,c('GO:0005750','GO:0006122'))
CIII_abs_Kma <- CIII_abs_Kma[-grep('MAS1',CIII_abs_Kma$PrimaryGeneName),]

CIII_abs_Sstip <- map2GO(protData_abs_Sstip,GOannotation_Sstip,c('GO:0005750','GO:0006122'))
CIII_abs_Sstip <- CIII_abs_Sstip[-grep('MAS1',CIII_abs_Sstip$PrimaryGeneName),]

# Complex IV
CIV_abs_Sce <- map2GO(protData_abs_Sce,GOannotation_Sce,c('GO:0005751','GO:0006123','GO:0004129'))

CIV_abs_Spo <- map2GO(protData_abs_Spo,GOannotation_Spo,c('GO:0005751','GO:0006123','GO:0004129'))
CIV_abs_Spo <- CIV_abs_Spo[-grep('P05511',CIV_abs_Spo$Accession),]

CIV_abs_Kma <- map2GO(protData_abs_Kma,GOannotation_Kma,c('GO:0005751','GO:0006123','GO:0004129'))
CIV_abs_Kma <- CIV_abs_Kma[-grep('W0T4B8',CIV_abs_Kma$Accession),]

CIV_abs_Sstip <- map2GO(protData_abs_Sstip,GOannotation_Sstip,c('GO:0005751','GO:0006123','GO:0004129'))

# ATP synthase
ATPsynth_abs_Sce <- map2GO(protData_abs_Sce,GOannotation_Sce,c('GO:0015986'))
ATPsynth_abs_Sce <- ATPsynth_abs_Sce[-grep('VPH1',ATPsynth_abs_Sce$PrimaryGeneName),]
ATPsynth_abs_Sce <- ATPsynth_abs_Sce[-grep('STV1',ATPsynth_abs_Sce$PrimaryGeneName),]
ATPsynth_abs_Sce <- ATPsynth_abs_Sce[-grep('Q3E824',ATPsynth_abs_Sce$Accession),]

ATPsynth_abs_Spo <- map2GO(protData_abs_Spo,GOannotation_Spo,c('GO:0015986'))
ATPsynth_abs_Spo <- ATPsynth_abs_Spo[-grep('vph1',ATPsynth_abs_Spo$PrimaryGeneName),]

ATPsynth_abs_Kma <- map2GO(protData_abs_Kma,GOannotation_Kma,c('GO:0015986'))
ATPsynth_abs_Kma <- ATPsynth_abs_Kma[-grep('VPH1',ATPsynth_abs_Kma$PrimaryGeneName),]
ATPsynth_abs_Kma <- ATPsynth_abs_Kma[-grep('STV1',ATPsynth_abs_Kma$PrimaryGeneName),]

ATPsynth_abs_Sstip <- map2GO(protData_abs_Sstip,GOannotation_Sstip,c('GO:0015986'))
ATPsynth_abs_Sstip <- ATPsynth_abs_Sstip[-grep('VPH1',ATPsynth_abs_Sstip$PrimaryGeneName),]

# Combine all data
combined_abs_Sce <- rbind(data.frame(t(colSums(NADH_dehydrogenase_abs_Sce[,5:7]))),
                                         data.frame(t(colSums(CII_abs_Sce[,5:7]))),
                                         data.frame(t(colSums(CIII_abs_Sce[,5:7]))),
                                         data.frame(t(colSums(CIV_abs_Sce[,5:7]))),
                                         data.frame(t(colSums(ATPsynth_abs_Sce[,5:7]))))
combined_abs_Sce$complex <- c('NADH','CII','CIII','CIV','ATP synthase')
combined_abs_Sce$organism <- c(rep('Sce',5))
combined_abs_Sce$meanAbundance <- rowMeans(combined_abs_Sce[,1:3])
combined_abs_Sce$sd <- apply(combined_abs_Sce[,1:3], 1, sd)
# S. pombe
combined_abs_Spo <- rbind(data.frame(t(colSums(NADH_dehydrogenase_abs_Spo[,5:7]))),
                          data.frame(t(colSums(CII_abs_Spo[,5:7]))),
                          data.frame(t(colSums(CIII_abs_Spo[,5:7]))),
                          data.frame(t(colSums(CIV_abs_Spo[,5:7]))),
                          data.frame(t(colSums(ATPsynth_abs_Spo[,5:7]))))
combined_abs_Spo$complex <- c('NADH','CII','CIII','CIV','ATP synthase')
combined_abs_Spo$organism <- c(rep('Spo',5))
combined_abs_Spo$meanAbundance <- rowMeans(combined_abs_Spo[,1:3])
combined_abs_Spo$sd <- apply(combined_abs_Spo[,1:3], 1, sd)
# K. marxianus
combined_abs_Kma <- rbind(data.frame(t(colSums(NADH_dehydrogenase_abs_Kma[,5:7]))),
                          data.frame(t(colSums(CII_abs_Kma[,5:7]))),
                          data.frame(t(colSums(CIII_abs_Kma[,5:7]))),
                          data.frame(t(colSums(CIV_abs_Kma[,5:7]))),
                          data.frame(t(colSums(ATPsynth_abs_Kma[,5:7]))))
combined_abs_Kma$complex <- c('NADH','CII','CIII','CIV','ATP synthase')
combined_abs_Kma$organism <- c(rep('Kma',5))
combined_abs_Kma$meanAbundance <- rowMeans(combined_abs_Kma[,1:3])
combined_abs_Kma$sd <- apply(combined_abs_Kma[,1:3], 1, sd)
# S. stipitis
combined_abs_Sstip <- rbind(data.frame(t(colSums(NADH_dehydrogenase_abs_Sstip[,5:7]))),
                            data.frame(t(colSums(CII_abs_Sstip[,5:7]))),
                            data.frame(t(colSums(CIII_abs_Sstip[,5:7]))),
                            data.frame(t(colSums(CIV_abs_Sstip[,5:7]))),
                            data.frame(t(colSums(ATPsynth_abs_Sstip[,5:7]))))
combined_abs_Sstip$complex <- c('NADH','CII','CIII','CIV','ATP synthase')
combined_abs_Sstip$organism <- c(rep('Sstip',5))
combined_abs_Sstip$meanAbundance <- rowMeans(combined_abs_Sstip[,1:3])
combined_abs_Sstip$sd <- apply(combined_abs_Sstip[,1:3], 1, sd)

combined_abs_All <- rbind(combined_abs_Sce,combined_abs_Spo,combined_abs_Kma,combined_abs_Sstip)
combined_abs_All$organism <- factor(combined_abs_All$organism,levels = c('Spo','Sce','Kma','Sstip'))
combined_abs_All$complex <- factor(combined_abs_All$complex,levels = c('NADH','CII','CIII','CIV','ATP synthase'))
colnames(combined_abs_All)[1:3] <- c('repl1','repl2','repl3')
df1 <- melt(combined_abs_All[grepl('NADH'),1:5],id.vars = c('complex','organism'),measure.vars = c('repl1','repl2','repl3'))

# Plot results
fig3C <- ggplot(combined_abs_All,aes(x=complex,y=meanAbundance*1000,fill=organism)) +
  geom_bar(position=position_dodge2(), stat="identity",color='black',alpha=0.85) +
  geom_errorbar(aes(ymin=(meanAbundance-sd)*1000,ymax=(meanAbundance+sd)*1000),width=0.4,
                position=position_dodge(0.9)) +
  geom_jitter(data=df1,aes(x=complex,y=value*1000,color=organism),
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),size=1) +
  scale_y_continuous(breaks = seq(0,12.5,1),limits = c(0,13),expand = c(0,0)) +
  scale_fill_manual(values=c('#CBBBA0','#000000','#1D71B8','#878787')) +
  scale_color_manual(values=c('#CBBBA0','#000000','#1D71B8','#878787')) +
  theme_classic() +
  theme(text = element_text(size = 8, color = 'black',face = 'bold'),
        legend.position = c(0.2,0.8),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10,face = 'bold',color = 'black'),
        axis.text.x = element_text(size = 7,face = 'bold',color = 'black'),
        axis.text.y = element_text(size = 8,face = 'bold',color = 'black'),
        axis.ticks = element_line(color = 'black')) +
  ylab('Abundance (mg/gDW)')
fig3C
# Save figure
#ggsave('../../results/fig3C.pdf',
 #      fig3C,width=7,height = 7,units = 'cm',device = 'pdf',dpi=300)




