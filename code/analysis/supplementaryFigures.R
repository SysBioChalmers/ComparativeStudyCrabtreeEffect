# Load required packages
library(tidyverse)
library(RColorBrewer)
library(reshape2)
library(ggpubr)
library(readr)
library(readxl)

# Define required functions
mapGenesToProtData <- function(geneList,protData){
  data <- data.frame(matrix(ncol = 7,nrow = 0))
  for(gene in geneList){
    geneData <- protData[grep(gene,protData$PrimaryGeneName),]
    if(!nrow(geneData)==0){
      data <- rbind(data,geneData)
    }
  }
  return(data)
}
map2GO <- function(proteomicsData,GOannotation,GOterms){
  if(length(GOterms)>1){
    proteinsInGOterms <- data.frame(matrix(ncol = 3, nrow = 0))
    for(GOterm in GOterms){
      proteinsMapped <- GOannotation[grep(GOterm,GOannotation$GOterm),1:3]
      proteinsInGOterms <- rbind(proteinsInGOterms,proteinsMapped)
    }
  }else{
    proteinsInGOterms <- GOannotation[grep(GOterms,GOannotation$GOterm),1:3]
  }
  colnames(proteinsInGOterms) <- c('Accession','GOterm','GOtermName')
  GOdata <- merge(proteomicsData,proteinsInGOterms)
  # Remove duplicated entries
  GOdata <- GOdata[!duplicated(GOdata$Accession),]
  return(GOdata)
}
setwd('../../data/experimentalData/')
# read proteomics data (allocation [g/g protein])
protData_Sce <- read.delim('proteomicsData_massPercent_Scer.txt',stringsAsFactors = F)
protData_Spo <- read.delim('proteomicsData_massPercent_Spombe.txt',stringsAsFactors = F)
protData_Sstip <- read.delim('proteomicsData_massPercent_Pstip.txt',stringsAsFactors = F)
protData_Kma <- read.delim('proteomicsData_massPercent_Kmarx.txt',stringsAsFactors = F)

# Read absolute proteome data [g/gDW]
protData_Sce_absolute <- read.delim('proteomicsData_g_per_gDW_Scer.txt',stringsAsFactors = F)
protData_Spo_absolute <- read.delim('proteomicsData_g_per_gDW_Spombe.txt',stringsAsFactors = F)
protData_Sstip_absolute <- read.delim('proteomicsData_g_per_gDW_Pstip.txt',stringsAsFactors = F)
protData_Kma_absolute <- read.delim('proteomicsData_g_per_gDW_Kmarx.txt',stringsAsFactors = F)

# read file with GO term annotation
GOannotation_Sce <- read.delim('../GOterms/Sce_GOterms.txt')
GOannotation_Spo <- read.delim('../GOterms/Spo_GOterms.txt')
GOannotation_Sstip <- read.delim('../GOterms/Pst_GOterms.txt')
GOannotation_Kma <- read.delim('../GOterms/Kma_GOterms.txt')

# fig S3
# Manually add annotation in Spo
protData_Spo_absolute[grep('SPCC191.02c',protData_Spo_absolute$GeneName),3] <- 'acs1'
protData_Spo_absolute[grep('O42873',protData_Spo_absolute$Accession),3] <- 'pdc201'
protData_Spo_absolute[grep('Q92345',protData_Spo_absolute$Accession),3] <- 'pdc101'
protData_Spo_absolute[grep('O14293',protData_Spo_absolute$Accession),3] <- 'atd1'
protData_Spo_absolute[grep('Q9URW9',protData_Spo_absolute$Accession),3] <- 'atd2'
protData_Spo_absolute[grep('P78733',protData_Spo_absolute$Accession),3] <- 'atd2'

# PDC
pdc_Spo <- c('pdc101','pdc201')
pdc_data_Spo <- mapGenesToProtData(pdc_Spo,protData_Spo_absolute)
pdc_data_Spo$organism <- c(rep('Spo',2))
pdc_Sce <- c('PDC1','PDC2','PDC5','PDC6')
pdc_data_Sce <- mapGenesToProtData(pdc_Sce,protData_Sce_absolute)
pdc_data_Sce$organism <- c(rep('Sce',3))
pdc_Kma <- c('PDC1')
pdc_data_Kma <- mapGenesToProtData(pdc_Kma,protData_Kma_absolute)
pdc_data_Kma$organism <- c('Kma')
pdc_Sstip <- c('PDC1','PDC2')
pdc_data_Sstip <- mapGenesToProtData(pdc_Sstip,protData_Sstip_absolute)
pdc_data_Sstip$organism <- c(rep('Sstip',2))
pdc_combined <- data.frame(rbind(pdc_data_Spo,pdc_data_Sce,pdc_data_Sstip,pdc_data_Kma))
pdc_combined$organism <- factor(pdc_combined$organism,levels=c('Spo','Sce','Kma','Sstip'))
# Calculate mean and sd
pdc_combined$mean <- rowMeans(pdc_combined[,5:7])
pdc_combined$sd <- apply(pdc_combined[,5:7],1,sd)
pdc_combined$PrimaryGeneName <- c('pdc101','pdc201','PDC1_Sce','PDC2_Sce','PDC6',
                             'PDC1_Kma','PDC1_Sstip','PDC2_Sstip')
pdc_combined$PrimaryGeneName <- factor(pdc_combined$PrimaryGeneName,levels = c('pdc101','pdc201',
                                                                     'PDC1_Sce','PDC2_Sce','PDC6',
                                                                     'PDC1_Kma',
                                                                     'PDC1_Sstip','PDC2_Sstip'))
df1 <- melt(pdc_combined[,c(3,5:8)],id.vars = c('PrimaryGeneName','organism'))

# plot data
figS3A <- ggplot(pdc_combined, aes(x = PrimaryGeneName, y = mean*1000, fill=organism)) +
  geom_bar(stat="identity", color="black", 
           position=position_dodge(),alpha=0.85) +
  geom_errorbar(aes(ymin=1000*(mean-sd), ymax=1000*(mean+sd)), width=.2,
                position=position_dodge(.9)) +
  geom_jitter(data=df1,aes(x=,y=value*1000,color=organism),
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),size=1) +
  scale_fill_manual(values=c('#CBBBA0','#000000','#1D71B8','#878787')) +
  scale_color_manual(values=c('#CBBBA0','#000000','#1D71B8','#878787')) +
  theme_classic() +
  theme(plot.title = element_text(size = 8,face = 'bold',color = 'black',hjust = 0.5),
        axis.text.x = element_text(angle = 90),
        axis.text = element_text(size = 8,face = 'bold',color='black'),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 8,face = 'bold',color='black'),
        legend.title = element_blank(),
        panel.border = element_rect(size = 1.1, color = 'black',fill = NA)) +
  scale_y_continuous(expand = c(0,0),limits = c(0,13),breaks = seq(0,13,1)) +
  scale_x_discrete(labels = c('pdc101','pdc201','PDC1','PDC2','PDC6','PDC1','PDC1','PDC2')) + 
  ylab('Abundance (mg/gDW)') +
  ggtitle('PDC')
figS3A

# ADH
adh_Spo <- c('adh1','adh4','adh8','fmd1','fmd3')
adh_data_Spo <- mapGenesToProtData(adh_Spo,protData_Spo_absolute)
adh_data_Spo$organism <- c(rep('Spo',2))
adh_Sce <- c('ADH1','ADH2','ADH3','ADH4','ADH5','ADH6')
adh_data_Sce <- mapGenesToProtData(adh_Sce,protData_Sce_absolute)
adh_data_Sce$organism <- c(rep('Sce',5))
adh_Kma <- c('adh','ADH1','ADH2','ADH3','ADH4a','ADH4b','ADH6')
adh_data_Kma <- mapGenesToProtData(adh_Kma,protData_Kma_absolute)
adh_data_Kma$organism <- c(rep('Kma',6))
adh_Sstip <- c('ADH1','ADH2','ADH3','ADH4','ADH5')
adh_data_Sstip <- mapGenesToProtData(adh_Sstip,protData_Sstip_absolute)
adh_data_Sstip$organism <- c(rep('Sstip',2))

adh_combined <- rbind(adh_data_Spo,adh_data_Sce,
                 adh_data_Kma,adh_data_Sstip)
adh_combined$organism <- factor(adh_combined$organism,levels = c('Spo','Sce','Kma','Sstip'))
adh_combined$PrimaryGeneName <- c('adh1','adh4','ADH1_Sce','ADH3_Sce','ADH4_Sce','ADH5_Sce','ADH6_Sce',
                             'adh','ADH1_Kma','ADH2_Kma','ADH3_Kma','ADH4b','ADH6_Kma',
                             'ADH11_Sstip','ADH5_Sstip')
adh_combined$PrimaryGeneName <- factor(adh_combined$PrimaryGeneName,levels = c('adh1','adh4',
                                                                     'ADH1_Sce','ADH3_Sce','ADH4_Sce','ADH5_Sce','ADH6_Sce',
                                                                     'adh','ADH1_Kma','ADH2_Kma','ADH3_Kma','ADH4b','ADH6_Kma',
                                                                     'ADH1_Sstip','ADH5_Sstip'))
# Calculate mean and standard deviation
adh_combined$mean = rowMeans(adh_combined[,5:7])
adh_combined$sd <- apply(adh_combined[,5:7],1,sd)
df1 <- melt(adh_combined[,c(3,5:8)],id.vars = c('PrimaryGeneName','organism'))

# plot data
figS3B <- ggplot(adh_combined, aes(x = PrimaryGeneName, y = mean*1000, fill=organism)) +
  geom_bar(stat="identity", color="black", 
           position=position_dodge(),alpha=0.85) +
  geom_errorbar(aes(ymin=1000*(mean-sd), ymax=1000*(mean+sd)), width=.2,
                position=position_dodge(.9)) +
  geom_jitter(data=df1,aes(x=,y=value*1000,color=organism),
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),size=1) +
  scale_fill_manual(values=c('#CBBBA0','#000000','#1D71B8','#878787')) +
  scale_color_manual(values=c('#CBBBA0','#000000','#1D71B8','#878787')) +
  theme_classic() +
  theme(plot.title = element_text(size = 8,face = 'bold',color = 'black',hjust = 0.5),
        axis.text.x = element_text(angle = 90),
        axis.text = element_text(size = 8,face = 'bold',color='black'),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 8,face = 'bold',color='black'),
        legend.title = element_blank(),
        panel.border = element_rect(size = 1.1, color = 'black',fill = NA)) +
  scale_y_continuous(expand = c(0,0),limits = c(0,13),breaks = seq(0,13,1)) +
  scale_x_discrete(labels = c('adh1','adh4','ADH1','ADH3','ADH4','ADH5','ADH6',
                              'adh','ADH1','ADH2','ADH3','ADH4b','ADH6','ADH1','ADH5')) +
  ylab('Abundance (mg/gDW)') +
  ggtitle('ADH')
figS3B

# ALD
ald_Spo <- c('atd1','atd2')
ald_data_Spo <- mapGenesToProtData(ald_Spo,protData_Spo_absolute)
ald_data_Spo$organism <- c(rep('Spo',2))
ald_Sce <- c('ALD2','ALD3','ALD4','ALD5','ALD6')
ald_data_Sce <- mapGenesToProtData(ald_Sce,protData_Sce_absolute)
ald_data_Sce$organism <- c(rep('Sce',5))
ald_Kma <- c('ALD2','ALD4','ALD5','ALD6')
ald_data_Kma <- mapGenesToProtData(ald_Kma,protData_Kma_absolute)
ald_data_Kma$organism <- c(rep('Kma',4))
ald_Sstip <- c('ALD2','ALD3','ALD4','ALD5','ALD6','ALD7')
ald_data_Sstip <- mapGenesToProtData(ald_Sstip,protData_Sstip_absolute)
ald_data_Sstip$organism <- c(rep('Sstip',5))

ald_combined <- rbind(ald_data_Spo,ald_data_Sce,
                 ald_data_Kma,ald_data_Sstip)
ald_combined$organism <- factor(ald_combined$organism,levels = c('Spo','Sce','Kma','Sstip'))
ald_combined$PrimaryGeneName <- c('atd1','atd2','ALD2_Sce','ALD3_Sce','ALD4_Sce','ALD5_Sce','ALD6_Sce',
                             'ALD2_Kma','ALD4_Kma','ALD5','ALD6_Kma','ALD2_Sstip','ALD3_Sstip',
                             'ALD5_Sstip','ALD6_Sstip','ALD7')
ald_combined$PrimaryGeneName <- factor(ald_combined$PrimaryGeneName,levels = c('atd1','atd2','ALD2_Sce','ALD3_Sce','ALD4_Sce','ALD5_Sce','ALD6_Sce',
                                                                     'ALD2_Kma','ALD4_Kma','ALD5','ALD6_Kma','ALD2_Sstip','ALD3_Sstip',
                                                                     'ALD5_Sstip','ALD6_Sstip','ALD7'))

# Calculate mean and standard deviation
ald_combined$mean = rowMeans(ald_combined[,5:7])
ald_combined$sd <- apply(ald_combined[,5:7],1,sd)
df1 <- melt(ald_combined[,c(3,5:8)],id.vars = c('PrimaryGeneName','organism'))

# plot data
figS3C <- ggplot(ald_combined, aes(x = PrimaryGeneName, y = mean*1000, fill=organism)) +
  geom_bar(stat="identity", color="black", 
           position=position_dodge(),alpha=0.85) +
  geom_errorbar(aes(ymin=1000*(mean-sd), ymax=1000*(mean+sd)), width=.2,
                position=position_dodge(.9)) +
  geom_jitter(data=df1,aes(x=,y=value*1000,color=organism),
             position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),size=1) +
  scale_fill_manual(values=c('#CBBBA0','#000000','#1D71B8','#878787')) +
  scale_color_manual(values=c('#CBBBA0','#000000','#1D71B8','#878787')) +
  theme_classic() +
  theme(plot.title = element_text(size = 8,face = 'bold',color = 'black',hjust = 0.5),
        axis.text.x = element_text(angle = 90),
        axis.text = element_text(size = 8,face = 'bold',color='black'),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 8,face = 'bold',color='black'),
        legend.title = element_blank(),
        panel.border = element_rect(size = 1.1, color = 'black',fill = NA)) +
  scale_y_continuous(expand = c(0,0),limits = c(0,1.75),breaks = seq(0,1.75,0.25)) +
  scale_x_discrete(labels = c('atd1','atd2','ALD2','ALD3','ALD4','ALD5','ALD6','ALD2','ALD4',
                              'ALD5','ALD6','ALD2','ALD3','ALD5','ALD6','ALD7')) +
  ylab('Abundance (mg/gDW)') +
  ggtitle('ALD')
figS3C

# ACS
acs_Spo <- c('acs1')
acs_data_Spo <- mapGenesToProtData(acs_Spo,protData_Spo_absolute)
acs_data_Spo$organism <- c('Spo')
acs_Sce <- c('ACS1','ACS2')
acs_data_Sce <- mapGenesToProtData(acs_Sce,protData_Sce_absolute)
acs_data_Sce$organism <- c('Sce')
acs_Kma <- c('ACS1','ACS2')
acs_data_Kma <- mapGenesToProtData(acs_Kma,protData_Kma_absolute)
acs_data_Kma$organism <- c(rep('Kma',2))
acs_Sstip <- c('ACS1','ACS2')
acs_data_Sstip <- mapGenesToProtData(acs_Sstip,protData_Sstip_absolute)
acs_data_Sstip$organism <- c(rep('Sstip',2))

acs_combined <- rbind(acs_data_Spo,acs_data_Sce,
                 acs_data_Kma,acs_data_Sstip)
acs_combined$organism <- factor(acs_combined$organism,levels = c('Spo','Sce','Kma','Sstip'))
acs_combined$PrimaryGeneName <- c('acs11','ACS2_Sce','ACS1_Kma','ACS2_Kma','ACS1_Sstip','ACS2_Sstip')
acs_combined$PrimaryGeneName <- factor(acs_combined$PrimaryGeneName,levels = c('acs11','ACS2_Sce','ACS1_Kma','ACS2_Kma','ACS1_Sstip','ACS2_Sstip'))

# Calculate mean and standard deviation
acs_combined$mean <- rowMeans(acs_combined[,5:7])
acs_combined$sd <- apply(acs_combined[,5:7],1,sd)
df1 <- melt(acs_combined[,c(3,5:8)],id.vars = c('PrimaryGeneName','organism'))

# plot data
figS3D <- ggplot(acs_combined, aes(x = PrimaryGeneName, y = mean*1000, fill=organism)) +
  geom_bar(stat="identity", color="black", 
           position=position_dodge(),alpha=0.85) +
  geom_errorbar(aes(ymin=1000*(mean-sd), ymax=1000*(mean+sd)), width=.2,
                position=position_dodge(.9)) +
  geom_jitter(data=df1,aes(x=,y=value*1000,color=organism),
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),size=1) +
  scale_fill_manual(values=c('#CBBBA0','#000000','#1D71B8','#878787')) +
  scale_color_manual(values=c('#CBBBA0','#000000','#1D71B8','#878787')) +
  theme_classic() +
  theme(plot.title = element_text(size = 8,face = 'bold',color = 'black',hjust = 0.5),
        axis.text.x = element_text(angle = 90),
        axis.text = element_text(size = 8,face = 'bold',color='black'),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 8,face = 'bold',color='black'),
        legend.title = element_blank(),
        panel.border = element_rect(size = 1.1, color = 'black',fill = NA)) +
  scale_y_continuous(expand = c(0,0),limits = c(0,0.9),breaks = seq(0,0.9,0.1)) +
  scale_x_discrete(labels = c('acs1','ACS2','ACS1','ACS2','ACS1','ACS2')) +
  ylab('Abundance (mg/gDW)') +
  ggtitle('ACS')
figS3D

# Combine into one figure
figS3 <- ggarrange(figS3A,figS3B,figS3C,figS3D,nrow=2,ncol=2,
                   labels=c('A','B','C','D'),
                   common.legend = TRUE,legend = 'bottom')
figS3

# Save figure
#ggsave('../../results/figS3.pdf',
 #      figS3,width=16,height = 16,units = 'cm',device = 'pdf',dpi=300)

# Fig S4
# Load list of GO terms related to amino acid biosynthesis
AAchildTerms <- read.delim('../GOterms/childTerms_AAbiosynthesis.tsv',
                           stringsAsFactors = FALSE,header = FALSE)
AA_GOterms <-  append(AAchildTerms$V1,'GO:0008652')

# fig S4A
# S. cerevisiae
lipidBiosynthesis_Sce <- map2GO(protData_Sce,GOannotation_Sce,c('GO:0008610','GO:0006633','GO:0009922',
                                                                'GO:0006636','GO:0008654','GO:0008299',
                                                                'GO:0006696'))
lipidBiosynthesis_Sce$GOtermName <- 'lipid biosynthesis'
nucleotideBioSynthesis_Sce <- map2GO(protData_Sce,GOannotation_Sce,c('GO:0009165','GO:0006164','GO:0006221'))
nucleotideBioSynthesis_Sce$GOtermName <- 'nucleotide biosynthesis'
AA_cerevisiae <- map2GO(protData_Sce,GOannotation_Sce,AAchildTerms$V1)
biosynthesis_Sce <- rbind(data.frame(t(colSums(lipidBiosynthesis_Sce[,5:7]))),
                          data.frame(t(colSums(nucleotideBioSynthesis_Sce[,5:7]))),
                          data.frame(t(colSums(AA_cerevisiae[,5:7]))))
biosynthesis_Sce$Process <- c('lipid biosynthesis','nucleotide biosynthesis','amino acid biosynthesis')
biosynthesis_Sce$Organism <- c(rep('Sce',3))
# Calculate mean and sd
biosynthesis_Sce$meanAllocation = rowMeans(biosynthesis_Sce[,1:3])
biosynthesis_Sce$sd <- apply(biosynthesis_Sce[,1:3],1,sd)

# S. pombe
lipidBiosynthesis_Spo <- map2GO(protData_Spo,GOannotation_Spo,c('GO:0008610','GO:0006633','GO:0009922',
                                                                'GO:0006636','GO:0008654','GO:0008299',
                                                                'GO:0006696'))
lipidBiosynthesis_Spo$GOtermName <- 'lipid biosynthesis'
nucleotideBioSynthesis_Spo <- map2GO(protData_Spo,GOannotation_Spo,c('GO:0009165','GO:0006164','GO:0006221'))
nucleotideBioSynthesis_Spo$GOtermName <- 'nucleotide biosynthesis'
AA_pombe <- map2GO(protData_Spo,GOannotation_Spo,AAchildTerms$V1)
biosynthesis_Spo <- rbind(data.frame(t(colSums(lipidBiosynthesis_Spo[,5:7]))),
                          data.frame(t(colSums(nucleotideBioSynthesis_Spo[,5:7]))),
                          data.frame(t(colSums(AA_pombe[,5:7]))))
biosynthesis_Spo$Process <- c('lipid biosynthesis','nucleotide biosynthesis','amino acid biosynthesis')
biosynthesis_Spo$Organism <- c(rep('Spo',3))
# Calculate mean and sd
biosynthesis_Spo$meanAllocation = rowMeans(biosynthesis_Spo[,1:3])
biosynthesis_Spo$sd <- apply(biosynthesis_Spo[,1:3],1,sd)

# K. marxianus
lipidBiosynthesis_Kma <- map2GO(protData_Kma,GOannotation_Kma,c('GO:0008610','GO:0006633','GO:0009922',
                                                                'GO:0006636','GO:0008654','GO:0008299',
                                                                'GO:0006696'))
lipidBiosynthesis_Kma$GOtermName <- 'lipid biosynthesis'
nucleotideBioSynthesis_Kma <- map2GO(protData_Kma,GOannotation_Kma,c('GO:0009165','GO:0006164','GO:0006221'))
nucleotideBioSynthesis_Kma$GOtermName <- 'nucleotide biosynthesis'
AA_marxianus <- map2GO(protData_Kma,GOannotation_Kma,AAchildTerms$V1)
biosynthesis_Kma <- rbind(data.frame(t(colSums(lipidBiosynthesis_Kma[,5:7]))),
                          data.frame(t(colSums(nucleotideBioSynthesis_Kma[,5:7]))),
                          data.frame(t(colSums(AA_marxianus[,5:7]))))
biosynthesis_Kma$Process <- c('lipid biosynthesis','nucleotide biosynthesis','amino acid biosynthesis')
biosynthesis_Kma$Organism <- c(rep('Kma',3))
# Calculate mean and sd
biosynthesis_Kma$meanAllocation = rowMeans(biosynthesis_Kma[,1:3])
biosynthesis_Kma$sd <- apply(biosynthesis_Kma[,1:3],1,sd)

# S. stipitis
lipidBiosynthesis_Sstip <- map2GO(protData_Sstip,GOannotation_Sstip,c('GO:0008610','GO:0006633','GO:0009922',
                                                                      'GO:0006636','GO:0008654','GO:0008299',
                                                                      'GO:0006696'))
lipidBiosynthesis_Sstip$GOtermName <- 'lipid biosynthesis'
nucleotideBioSynthesis_Sstip <- map2GO(protData_Sstip,GOannotation_Sstip,c('GO:0009165','GO:0006164','GO:0006221'))
nucleotideBioSynthesis_Sstip$GOtermName <- 'nucleotide biosynthesis'
AA_stipitis <- map2GO(protData_Sstip,GOannotation_Sstip,AAchildTerms$V1)
biosynthesis_Sstip <- rbind(data.frame(t(colSums(lipidBiosynthesis_Sstip[,5:7]))),
                            data.frame(t(colSums(nucleotideBioSynthesis_Sstip[,5:7]))),
                            data.frame(t(colSums(AA_stipitis[,5:7]))))
biosynthesis_Sstip$Process <- c('lipid biosynthesis','nucleotide biosynthesis','amino acid biosynthesis')
biosynthesis_Sstip$Organism <- c(rep('Sstip',3))
# Calculate mean and sd
biosynthesis_Sstip$meanAllocation = rowMeans(biosynthesis_Sstip[,1:3])
biosynthesis_Sstip$sd <- apply(biosynthesis_Sstip[,1:3],1,sd)

# Combine data
biosynthesis_combined <- rbind(biosynthesis_Sce,biosynthesis_Spo,biosynthesis_Sstip,biosynthesis_Kma)
biosynthesis_combined$Organism <- factor(biosynthesis_combined$Organism,levels = c('Spo','Sce','Kma','Sstip'))
df1 <- melt(biosynthesis_combined[,1:5],id.vars = c('Process','Organism'))

# Plot data
figS4A <- ggplot(biosynthesis_combined,aes(x=reorder(Process,meanAllocation),y=meanAllocation*100,
                                      fill=Organism)) +
  geom_bar(position=position_dodge2(), stat="identity",color='black',alpha=0.85) +
  geom_errorbar(aes(ymin=100*(meanAllocation-sd),ymax=100*(meanAllocation+sd)),width=0.4,
                position=position_dodge(0.9)) +
  geom_jitter(data=df1,aes(x=Process,y=value*100,color=Organism),
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),size=1) +
  scale_y_continuous(breaks = seq(0,10,1),limits = c(0,10),expand = c(0,0)) +
  scale_fill_manual(values=c('#CBBBA0','#1D1D1B','#1D71B8','#878787')) +
  scale_color_manual(values=c('#CBBBA0','#1D1D1B','#1D71B8','#878787')) +
  theme_classic() +
  theme(text = element_text(size = 8, color = 'black',face = 'bold'),
        legend.position = c(0.2,0.8),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10,face = 'bold',color = 'black'),
        axis.text.x = element_text(size = 7,face = 'bold',color = 'black'),
        axis.text.y = element_text(size = 8,face = 'bold',color = 'black')) +
  ylab('Proteome allocation (%)')
figS4A

# Save figure
#ggsave('../../results/figS4A.pdf',
 #      figS4A, device = 'pdf',height = 6, width = 6, units = 'cm',dpi = 300)

# figS4B
# S. cerevisiae
lipidBiosynthesis_Sce_abs <- map2GO(protData_Sce_absolute,GOannotation_Sce,c('GO:0008610','GO:0006633','GO:0009922',
                                                                             'GO:0006636','GO:0008654','GO:0008299',
                                                                             'GO:0006696'))
lipidBiosynthesis_Sce_abs$GOtermName <- 'lipid biosynthesis'
nucleotideBioSynthesis_Sce_abs <- map2GO(protData_Sce_absolute,GOannotation_Sce,c('GO:0009165','GO:0006164','GO:0006221'))
nucleotideBioSynthesis_Sce_abs$GOtermName <- 'nucleotide biosynthesis'
AA_cerevisiae_absolute <- map2GO(protData_Sce_absolute,GOannotation_Sce,AAchildTerms$V1)
biosynthesis_Sce_abs <- rbind(data.frame(t(colSums(lipidBiosynthesis_Sce_abs[,5:7]))),
                              data.frame(t(colSums(nucleotideBioSynthesis_Sce_abs[,5:7]))),
                              data.frame(t(colSums(AA_cerevisiae_absolute[,5:7]))))
biosynthesis_Sce_abs$Process <- c('lipid biosynthesis','nucleotide biosynthesis','amino acid biosynthesis')
biosynthesis_Sce_abs$Organism <- c(rep('Sce',3))
# Calculate mean and sd
biosynthesis_Sce_abs$meanAllocation = rowMeans(biosynthesis_Sce_abs[,1:3])
biosynthesis_Sce_abs$sd <- apply(biosynthesis_Sce_abs[,1:3],1,sd)

# S. pombe
lipidBiosynthesis_Spo_abs <- map2GO(protData_Spo_absolute,GOannotation_Spo,c('GO:0008610','GO:0006633','GO:0009922',
                                                                             'GO:0006636','GO:0008654','GO:0008299',
                                                                             'GO:0006696'))
lipidBiosynthesis_Spo_abs$GOtermName <- 'lipid biosynthesis'
nucleotideBioSynthesis_Spo_abs <- map2GO(protData_Spo_absolute,GOannotation_Spo,c('GO:0009165','GO:0006164','GO:0006221'))
nucleotideBioSynthesis_Spo_abs$GOtermName <- 'nucleotide biosynthesis'
AA_pombe_absolute <- map2GO(protData_Spo_absolute,GOannotation_Spo,AAchildTerms$V1)
biosynthesis_Spo_abs <- rbind(data.frame(t(colSums(lipidBiosynthesis_Spo_abs[,5:7]))),
                              data.frame(t(colSums(nucleotideBioSynthesis_Spo_abs[,5:7]))),
                              data.frame(t(colSums(AA_pombe_absolute[,5:7]))))
biosynthesis_Spo_abs$Process <- c('lipid biosynthesis','nucleotide biosynthesis','amino acid biosynthesis')
biosynthesis_Spo_abs$Organism <- c(rep('Spo',3))
# Calculate mean and sd
biosynthesis_Spo_abs$meanAllocation = rowMeans(biosynthesis_Spo_abs[,1:3])
biosynthesis_Spo_abs$sd <- apply(biosynthesis_Spo_abs[,1:3],1,sd)

# K. marxianus
lipidBiosynthesis_Kma_abs <- map2GO(protData_Kma_absolute,GOannotation_Kma,c('GO:0008610','GO:0006633','GO:0009922',
                                                                             'GO:0006636','GO:0008654','GO:0008299',
                                                                             'GO:0006696'))
lipidBiosynthesis_Kma_abs$GOtermName <- 'lipid biosynthesis'
nucleotideBioSynthesis_Kma_abs <- map2GO(protData_Kma_absolute,GOannotation_Kma,c('GO:0009165','GO:0006164','GO:0006221'))
nucleotideBioSynthesis_Kma_abs$GOtermName <- 'nucleotide biosynthesis'
AA_marxianus_absolute <- map2GO(protData_Kma_absolute,GOannotation_Kma,AAchildTerms$V1)
biosynthesis_Kma_abs <- rbind(data.frame(t(colSums(lipidBiosynthesis_Kma_abs[,5:7]))),
                              data.frame(t(colSums(nucleotideBioSynthesis_Kma_abs[,5:7]))),
                              data.frame(t(colSums(AA_marxianus_absolute[,5:7]))))
biosynthesis_Kma_abs$Process <- c('lipid biosynthesis','nucleotide biosynthesis','amino acid biosynthesis')
biosynthesis_Kma_abs$Organism <- c(rep('Kma',3))
# Calculate mean and sd
biosynthesis_Kma_abs$meanAllocation = rowMeans(biosynthesis_Kma_abs[,1:3])
biosynthesis_Kma_abs$sd <- apply(biosynthesis_Kma_abs[,1:3],1,sd)

# S. stipitis
lipidBiosynthesis_Sstip_abs <- map2GO(protData_Sstip_absolute,GOannotation_Sstip,c('GO:0008610','GO:0006633','GO:0009922',
                                                                                   'GO:0006636','GO:0008654','GO:0008299',
                                                                                   'GO:0006696'))
lipidBiosynthesis_Sstip_abs$GOtermName <- 'lipid biosynthesis'
nucleotideBioSynthesis_Sstip_abs <- map2GO(protData_Sstip_absolute,GOannotation_Sstip,c('GO:0009165','GO:0006164','GO:0006221'))
nucleotideBioSynthesis_Sstip_abs$GOtermName <- 'nucleotide biosynthesis'
AA_stipitis_absolute <- map2GO(protData_Sstip_absolute,GOannotation_Sstip,AAchildTerms$V1)
biosynthesis_Sstip_abs <- rbind(data.frame(t(colSums(lipidBiosynthesis_Sstip_abs[,5:7]))),
                                data.frame(t(colSums(nucleotideBioSynthesis_Sstip_abs[,5:7]))),
                                data.frame(t(colSums(AA_stipitis_absolute[,5:7]))))
biosynthesis_Sstip_abs$Process <- c('lipid biosynthesis','nucleotide biosynthesis','amino acid biosynthesis')
biosynthesis_Sstip_abs$Organism <- c(rep('Sstip',3))
# Calculate mean and sd
biosynthesis_Sstip_abs$meanAllocation = rowMeans(biosynthesis_Sstip_abs[,1:3])
biosynthesis_Sstip_abs$sd <- apply(biosynthesis_Sstip_abs[,1:3],1,sd)

# Combine data
biosynthesis_combined_abs <- rbind(biosynthesis_Sce_abs,biosynthesis_Spo_abs,biosynthesis_Sstip_abs,biosynthesis_Kma_abs)
biosynthesis_combined_abs$Organism <- factor(biosynthesis_combined_abs$Organism,levels = c('Spo','Sce','Kma','Sstip'))
df1 <- melt(biosynthesis_combined_abs[,1:5],id.vars = c('Process','Organism'))

# Plot data
figS4B <- ggplot(biosynthesis_combined_abs,aes(x=reorder(Process,meanAllocation),y=meanAllocation*1000,
                                          fill=Organism)) +
  geom_bar(position=position_dodge2(), stat="identity",color='black',alpha=0.85) +
  geom_errorbar(aes(ymin=1000*(meanAllocation-sd),ymax=1000*(meanAllocation+sd)),width=0.4,
                position=position_dodge(0.9)) +
  geom_jitter(data=df1,aes(x=Process,y=value*1000,color=Organism),
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),size=1) +
  scale_y_continuous(breaks = seq(0,45,5),limits = c(0,46),expand = c(0,0)) +
  scale_fill_manual(values=c('#CBBBA0','#1D1D1B','#1D71B8','#878787')) +
  scale_color_manual(values=c('#CBBBA0','#1D1D1B','#1D71B8','#878787')) +
  theme_classic() +
  theme(text = element_text(size = 8, color = 'black',face = 'bold'),
        legend.position = c(0.2,0.8),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10,face = 'bold',color = 'black'),
        axis.text.x = element_text(size = 7,face = 'bold',color = 'black'),
        axis.text.y = element_text(size = 8,face = 'bold',color = 'black')) +
  ylab('Abundance (mg/gDW)')
figS4B

# Save figure
#ggsave('../../results/figS4B.pdf',
 #      figS4B, device = 'pdf',height = 6, width = 6, units = 'cm',dpi = 300)

# fig S4C

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

# Keep relevant GO terms for Sce
capUse_Sce <- capUse_Sce[capUse_Sce$protID %in% GO_Sce$Entry,]
idx <- match(capUse_Sce$protID,GO_Sce$Entry)
capUse_Sce$GOterm <- GO_Sce$System[idx]
capUse_Sce <- capUse_Sce %>% mutate_if(is.numeric,round,digits = 3)

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
dataToPlot <- combined[combined$GOterm %in% c('Amino acid metabolism','Fatty acid metabolism',
                                                           'Nucleotide biosynthesis','Sterol biosynthesis'),]

figS4C <- ggplot(dataToPlot,aes(x = Organism, y = capUsage, fill = Organism)) +
  geom_boxplot(alpha=0.85) +
  facet_grid(. ~ GOterm) +
  labs(x = '', y = 'Capacity usage (%)') +
  theme_classic() +
  theme(text = element_text(size = 8, color = 'black', face = 'bold'),
        axis.text = element_text(size = 8, color = 'black', face = 'bold'),
        axis.title.y.right = element_blank(),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        panel.spacing = unit(3, "mm"),
        panel.border = element_rect(size = 1, color = 'black',fill = NA),
        legend.position = 'none') +
  scale_fill_manual(values=c('#1D71B8','#000000')) +
  scale_y_continuous(sec.axis = dup_axis(breaks = 0))
figS4C

# Save figure
#ggsave('../../results/figS4C.pdf',
 #      figS4C, device = 'pdf',height = 6, width = 6, units = 'cm',dpi = 300)

# Fig S5
# Read data on PPP genes
PPP_proteins <- read.delim('PPP_all.txt')
PPP_proteins <- PPP_proteins[,c(1,3,8)]

glycolysis_Sce_absolute <- map2GO(protData_Sce_absolute,GOannotation_Sce,c('GO:0006096'))
# remove pyruvate dehydrogenase genes (PDA1 and PDB1) from glycolysis GO term
glycolysis_Sce_absolute <- glycolysis_Sce_absolute[-grep('PD',glycolysis_Sce_absolute$PrimaryGeneName),]
mitochondrial_Sce_absolute <- mapMitoProteins(protData_Sce_absolute,GOannotation_Sce)
mitochondrial_Sce_absolute <- mitochondrial_Sce_absolute[!mitochondrial_Sce_absolute$Accession %in% glycolysis_Sce_absolute$Accession[grep('glycolytic process',glycolysis_Sce_absolute$GOtermName)],]
# load child terms of amino acid biosynthesis
AAchildTerms <- read.delim('../GOterms/childTerms_AAbiosynthesis.tsv',
                           stringsAsFactors = FALSE,header = FALSE)
AA_cerevisiae_absolute <- map2GO(protData_Sce_absolute,GOannotation_Sce,AAchildTerms$V1)
PPP_proteins_Sce_absolute <- PPP_proteins[grep('Sce',PPP_proteins$organism),]
PPP_data_Sce_absolute <- mapUniprotToProtData(PPP_proteins_Sce_absolute$Accession,protData_Sce_absolute)
combined_Sce_absolute <- rbind(data.frame(t(colSums(glycolysis_Sce_absolute[,5:7]))),
                               data.frame(t(colSums(PPP_data_Sce_absolute[,5:7]))),
                               data.frame(t(colSums(AA_cerevisiae_absolute[,5:7]))),
                               data.frame(t(colSums(mitochondrial_Sce_absolute[,5:7]))))
combined_Sce_absolute$Process <- c('Glycolysis','PPP','Amino acid biosynthesis','Mitochondrion')
combined_Sce_absolute$Organism <- c(rep('Sce',4))
combined_Sce_absolute$meanAbundance <- rowMeans(combined_Sce_absolute[,1:3])
combined_Sce_absolute$sd <- apply(combined_Sce_absolute[,1:3], 1, sd)
#combined_Sce_absolute <- combined_Sce_absolute[,4:7]

# S. pombe
glycolysis_Spo_absolute <- map2GO(protData_Spo_absolute,GOannotation_Spo,c('GO:0006096'))
glycolysis_Spo_absolute <- glycolysis_Spo_absolute[-grep('pda1',glycolysis_Spo_absolute$PrimaryGeneName),]
glycolysis_Spo_absolute <- glycolysis_Spo_absolute[-grep('pdb1',glycolysis_Spo_absolute$PrimaryGeneName),]
mitochondrial_Spo_absolute <- mapMitoProteins(protData_Spo_absolute,GOannotation_Spo)
mitochondrial_Spo_absolute <- mitochondrial_Spo_absolute[!mitochondrial_Spo_absolute$Accession %in% glycolysis_Spo_absolute$Accession[grep('glycolytic process',glycolysis_Spo_absolute$GOtermName)],]
AA_pombe_absolute <- map2GO(protData_Spo_absolute,GOannotation_Spo,AAchildTerms$V1)
PPP_proteins_Spo_absolute <- PPP_proteins[grep('Spo',PPP_proteins$organism),]
PPP_data_Spo_absolute <- mapUniprotToProtData(PPP_proteins_Spo_absolute$Accession,protData_Spo_absolute)
combined_Spo_absolute <- rbind(data.frame(t(colSums(glycolysis_Spo_absolute[,5:7]))),
                               data.frame(t(colSums(PPP_data_Spo_absolute[,5:7]))),
                               data.frame(t(colSums(AA_pombe_absolute[,5:7]))),
                               data.frame(t(colSums(mitochondrial_Spo_absolute[,5:7]))))
combined_Spo_absolute$Process <- c('Glycolysis','PPP','Amino acid biosynthesis','Mitochondrion')
combined_Spo_absolute$Organism <- c(rep('Spo',4))
combined_Spo_absolute$meanAbundance <- rowMeans(combined_Spo_absolute[,1:3])
combined_Spo_absolute$sd <- apply(combined_Spo_absolute[,1:3], 1, sd)
#combined_Spo_absolute <- combined_Spo_absolute[,4:7]

# S.stipitis
glycolysis_Sstip_absolute <- map2GO(protData_Sstip_absolute,GOannotation_Sstip,c('GO:0006096'))
glycolysis_Sstip_absolute <- glycolysis_Sstip_absolute[-grep('PD',glycolysis_Sstip_absolute$PrimaryGeneName),]
mitochondrial_Sstip_absolute <- mapMitoProteins(protData_Sstip_absolute,GOannotation_Sstip)
mitochondrial_Sstip_absolute <- mitochondrial_Sstip_absolute[!mitochondrial_Sstip_absolute$Accession %in% glycolysis_Sstip_absolute$Accession[grep('glycolytic process',glycolysis_Sstip_absolute$GOtermName)],]
AA_stipitis_absolute <- map2GO(protData_Sstip_absolute,GOannotation_Sstip,AAchildTerms$V1)
PPP_proteins_Sstip_absolute <- PPP_proteins[grep('Sstip',PPP_proteins$organism),]
PPP_data_Sstip_absolute <- mapUniprotToProtData(PPP_proteins_Sstip_absolute$Accession,protData_Sstip_absolute)
combined_Sstip_absolute <- rbind(data.frame(t(colSums(glycolysis_Sstip_absolute[,5:7]))),
                                 data.frame(t(colSums(PPP_data_Sstip_absolute[,5:7]))),
                                 data.frame(t(colSums(AA_stipitis_absolute[,5:7]))),
                                 data.frame(t(colSums(mitochondrial_Sstip_absolute[,5:7]))))
combined_Sstip_absolute$Process <- c('Glycolysis','PPP','Amino acid biosynthesis','Mitochondrion')
combined_Sstip_absolute$Organism <- c(rep('Sstip',4))
combined_Sstip_absolute$meanAbundance <- rowMeans(combined_Sstip_absolute[,1:3])
combined_Sstip_absolute$sd <- apply(combined_Sstip_absolute[,1:3], 1, sd)
#combined_Sstip_absolute <- combined_Sstip_absolute[,4:7]

# K. marxianus
glycolysis_Kma_absolute <- map2GO(protData_Kma_absolute,GOannotation_Kma,c('GO:0006096'))
glycolysis_Kma_absolute <- glycolysis_Kma_absolute[-grep('PDA1',glycolysis_Kma_absolute$PrimaryGeneName),]
glycolysis_Kma_absolute <- glycolysis_Kma_absolute[-grep('PDB1',glycolysis_Kma_absolute$PrimaryGeneName),]
mitochondrial_Kma_absolute <- mapMitoProteins(protData_Kma_absolute,GOannotation_Kma)
mitochondrial_Kma_absolute <- mitochondrial_Kma_absolute[!mitochondrial_Kma_absolute$Accession %in% glycolysis_Kma_absolute$GOtermName,]
AA_marxianus_absolute <- map2GO(protData_Kma_absolute,GOannotation_Kma,AAchildTerms$V1)
PPP_proteins_Kma_absolute <- PPP_proteins[grep('Kma',PPP_proteins$organism),]
PPP_data_Kma_absolute <- mapUniprotToProtData(PPP_proteins_Kma_absolute$Accession,protData_Kma_absolute)
combined_Kma_absolute <- rbind(data.frame(t(colSums(glycolysis_Kma_absolute[,5:7]))),
                               data.frame(t(colSums(PPP_data_Kma_absolute[,5:7]))),
                               data.frame(t(colSums(AA_marxianus_absolute[,5:7]))),
                               data.frame(t(colSums(mitochondrial_Kma_absolute[,5:7]))))
combined_Kma_absolute$Process <- c('Glycolysis','PPP','Amino acid biosynthesis','Mitochondrion')
combined_Kma_absolute$Organism <- c(rep('Kma',4))
combined_Kma_absolute$meanAbundance <- rowMeans(combined_Kma_absolute[,1:3])
combined_Kma_absolute$sd <- apply(combined_Kma_absolute[,1:3], 1, sd)
#combined_Kma_absolute <- combined_Kma_absolute[,4:7]

# summarize the processes and plot results
combined_All_absolute <- rbind(combined_Sce_absolute,combined_Spo_absolute,combined_Sstip_absolute,combined_Kma_absolute)
combined_All_absolute$Organism <- factor(combined_All_absolute$Organism,levels = c('Spo','Sce','Kma','Sstip'))
combined_All_absolute$rGlu <- c(rep(18.8,4),rep(8.9,4),rep(4.5,4),rep(5.6,4))
colnames(combined_All_absolute)[1:3] <- c('repl1','repl2','repl3')
AA <- combined_All_absolute[grepl('Amino acid biosynthesis',combined_All_absolute$Process),]
df1 <- melt(AA[,c(1:4,8)],id.vars = c('Process','rGlu'),measure.vars = c('repl1','repl2','repl3'))
glycolysis <- combined_All_absolute[grepl('Glycolysis',combined_All_absolute$Process),]
df2 <- melt(glycolysis[,c(1:4,8)],id.vars = c('Process','rGlu'),measure.vars = c('repl1','repl2','repl3'))
mito <- combined_All_absolute[grepl('Mitochondrion',combined_All_absolute$Process),]
df3 <- melt(mito[,c(1:4,8)],id.vars = c('Process','rGlu'),measure.vars = c('repl1','repl2','repl3'))
PPP <- combined_All_absolute[grepl('PPP',combined_All_absolute$Process),]
df4 <- melt(PPP[,c(1:4,8)],id.vars = c('Process','rGlu'),measure.vars = c('repl1','repl2','repl3'))

# Plot figure
figS5 <- ggplot() +
  geom_point(data=AA,aes(x=rGlu,y=meanAbundance*1000),color='#9970AB',shape=22,fill='#9970AB') +
  geom_errorbar(data=AA,aes(x=rGlu,ymin=(meanAbundance-sd)*1000,ymax=(meanAbundance+sd)*1000),
                width=.2,position=position_dodge(.9)) +
  geom_jitter(data=df1,aes(x=rGlu,y=value*1000),
              position = position_jitter(0.2),color='#9970AB',size=2) +
  geom_point(data=glycolysis,aes(x=rGlu,y=meanAbundance*1000),color='#5AAE61',shape=22,fill='#5AAE61') +
  geom_errorbar(data=glycolysis,aes(x=rGlu,ymin=(meanAbundance-sd)*1000,ymax=(meanAbundance+sd)*1000),
                width=.2,position=position_dodge(.9)) +
  geom_jitter(data=df2,aes(x=rGlu,y=value*1000),
              position = position_jitter(0.2),color='#5AAE61',size=2) +
  geom_point(data=mito,aes(x=rGlu,y=meanAbundance*1000),color='#6BAED6',shape=22,fill='#6BAED6') +
  geom_errorbar(data=mito,aes(x=rGlu,ymin=(meanAbundance-sd)*1000,ymax=(meanAbundance+sd)*1000),
                width=.2,position=position_dodge(.9)) +
  geom_jitter(data=df3,aes(x=rGlu,y=value*1000),
              position = position_jitter(0.2),color='#6BAED6',size=2) +
  geom_point(data=PPP,aes(x=rGlu,y=meanAbundance*1000),color='#A6761D',shape=22,fill='#A6761D') +
  geom_errorbar(data=PPP,aes(x=rGlu,ymin=(meanAbundance-sd)*1000,ymax=(meanAbundance+sd)*1000),
                width=.2,position=position_dodge(.9)) +
  geom_jitter(data=df4,aes(x=rGlu,y=value*1000),
              position = position_jitter(0.2),color='#A6761D',size=2) +
  theme_classic() +
  theme(axis.title = element_text(size = 8,face = 'bold',color = 'black'),
        axis.text = element_text(size = 8,face = 'bold',color = 'black'),
        axis.text.x = element_text(size = 8,face = 'bold',color = 'black'),
        legend.title = element_blank(),
        legend.text = element_text(size = 8,face = 'bold',color = 'black'),
        legend.position = c(0.7,0.9)) +
  labs(x='rGlu (mmol/gDW/h)',y='Proteome allocation (%)') +
  scale_color_manual(values = c('#9970AB','#5AAE61','#6BAED6','#A6761D')) +
  scale_fill_manual(values = c('#9970AB','#5AAE61','#6BAED6','#A6761D')) +
  scale_y_continuous(limits = c(0,125),breaks = seq(0,120,10),expand = c(0,0)) +
  scale_x_continuous(breaks=c(4.5,5.6,8.9,18.8))
figS5
# Save plot
#ggsave('../../results/figS5.pdf',
 #      figS5,width=6,height = 6,units = 'cm',device = 'pdf',dpi=300)


# Fig S6
# Load data on glycolytic proteins
glycolysisData <- read.delim('glycolysis_combined.txt',sep = '\t',stringsAsFactors = F)
# Convert rxn and organism columns to factor
glycolysisData$Rxn <- as.factor(glycolysisData$Rxn)
glycolysisData$organism <- factor(glycolysisData$organism,levels = c('Spo','Sce','Kma','Sstip'))

# Calculate mean and standard deviation
glycolysisData$meanAbundance = rowMeans(glycolysisData[,4:6])
glycolysisData$sd <- apply(glycolysisData[,4:6],1,sd)

# fig S6A Glucose phosphorylation
hexokinase <- glycolysisData[grep('Glucose phosphorylation',glycolysisData$Rxn),]
hexokinase$gene = factor(hexokinase$gene, levels = c('hxk1','hxk2','HXK1_Sce','HXK2_Sce','GLK1_Sce',
                                                     'GLK1_Kma','RAG5','NAG5','GLK1_Sstip','HXK1_Sstip'))
df1 <- melt(hexokinase[,c(1,4:7)],id.vars = c('gene','organism'))

figS6A <- ggplot(hexokinase, aes(x = gene, y = meanAbundance*1000, fill=organism)) +
  geom_bar(position=position_dodge2(), stat="identity",color='black',alpha=0.85) +
  geom_errorbar(aes(ymin=(meanAbundance-sd)*1000,ymax=(meanAbundance+sd)*1000),width=0.4,
                position=position_dodge(0.9)) +
  geom_jitter(data=df1,aes(x=gene,y=value*1000,color=organism),
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),size=0.5) +
  scale_fill_manual(values=c('#CBBBA0','#000000','#1D71B8','#878787')) +
  scale_color_manual(values=c('#CBBBA0','#000000','#1D71B8','#878787')) +
  theme_classic() +
  theme(plot.title = element_text(size = 10,face = 'bold',color = 'black',hjust = 0.5),
        axis.text.x = element_text(angle = 90),
        axis.text = element_text(size = 7,face = 'bold',color='black'),
        axis.title = element_blank(),
        panel.border = element_rect(size = 0.5, color = 'black',fill = NA),
        legend.title = element_blank()) +
  scale_y_continuous(expand = c(0,0),limits = c(0,1.12),breaks = seq(0,1,0.25)) +
  scale_x_discrete(labels = c('hxk1','hxk2','HXK1','HXK2','GLK1','GLK1',
                              'RAG5','NAG5','GLK1','HXK1')) +
  ggtitle('GLK')
figS6A

# fig S6B Glucose 6-phosphate isomerase
g6p_isomerase <- glycolysisData[grep('Glucose-6-phosphate isomerase',glycolysisData$Rxn),]
g6p_isomerase$gene <- factor(g6p_isomerase$gene, levels = c('pgi1','PGI1_Sce','RAG2','PGI1_Sstip'))
# Calculate mean and standard deviation
g6p_isomerase$meanAbundance = rowMeans(g6p_isomerase[,4:6])
g6p_isomerase$sd <- apply(g6p_isomerase[,4:6],1,sd)
df1 <- melt(g6p_isomerase[,c(1,4:7)],id.vars = c('gene','organism'))

figS6B <- ggplot(g6p_isomerase, aes(x = gene, y = meanAbundance*1000, fill=organism)) +
  geom_bar(position=position_dodge2(), stat="identity",color='black',alpha=0.85) +
  geom_errorbar(aes(ymin=(meanAbundance-sd)*1000,ymax=(meanAbundance+sd)*1000),width=0.4,
                position=position_dodge(0.9)) +
  geom_jitter(data=df1,aes(x=gene,y=value*1000,color=organism),
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),size=0.5) +
  scale_fill_manual(values=c('#CBBBA0','#000000','#1D71B8','#878787')) +
  scale_color_manual(values=c('#CBBBA0','#000000','#1D71B8','#878787')) +
  theme_classic() +
  theme(plot.title = element_text(size = 10,face = 'bold',color = 'black',hjust = 0.5),
        axis.text.x = element_text(angle = 90),
        axis.text = element_text(size = 7,face = 'bold',color='black'),
        axis.title = element_blank(),
        panel.border = element_rect(size = 0.5, color = 'black',fill = NA),
        legend.title = element_blank()) +
  scale_y_continuous(expand = c(0,0),limits = c(0,1.6),breaks = seq(0,1.6,0.2)) +
  scale_x_discrete(labels=c('pgi1','PGI1','RAG2','PGI1')) +
  ggtitle('PGI')
figS6B

# fig S6C phosphofructokinase
pfk <- glycolysisData[grep('Phosphofructokinase',glycolysisData$Rxn),]
pfk$gene <- factor(pfk$gene,levels = c('pfk1','PFK1_Sce','PFK2_Sce','PFK1_Kma','PFK2_Kma',
                                       'PFK1_Sstip','PFK2_Sstip'))
# calculate mean abundance and sd
pfk$meanAbundance <- rowMeans(pfk[,4:6])
pfk$sd <- apply(pfk[,4:6],1,sd)
df1 <- melt(pfk[,c(1,4:7)],id.vars = c('gene','organism'))

figS6C <- ggplot(pfk, aes(x = gene, y = meanAbundance*1000, fill=organism)) +
  geom_bar(position=position_dodge2(), stat="identity",color='black',alpha=0.85) +
  geom_errorbar(aes(ymin=(meanAbundance-sd)*1000,ymax=(meanAbundance+sd)*1000),width=0.4,
                position=position_dodge(0.9)) +
  geom_jitter(data=df1,aes(x=gene,y=value*1000,color=organism),
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),size=0.5) +
  scale_fill_manual(values=c('#CBBBA0','#000000','#1D71B8','#878787')) +
  scale_color_manual(values=c('#CBBBA0','#000000','#1D71B8','#878787')) +
  theme_classic() +
  theme(plot.title = element_text(size = 10,face = 'bold',color = 'black',hjust = 0.5),
        axis.text.x = element_text(angle = 90),
        axis.text = element_text(size = 7,face = 'bold',color='black'),
        axis.title = element_blank(),
        panel.border = element_rect(size = 0.5, color = 'black',fill = NA),
        legend.title = element_blank()) +
  scale_y_continuous(expand = c(0,0),limits = c(0,2.4),breaks = seq(0,2.4,0.25)) +
  scale_x_discrete(labels = c('pfk1','PFK1','PFK2','PFK1','PFK2','PFK1','PFK2')) +
  ggtitle('PFK')
figS6C

# fig S6D fructose 1,6-bisphosphate aldolase
fba <- glycolysisData[grep('Fructose-bisphosphate aldolase',glycolysisData$Rxn),]
fba$gene <- factor(fba$gene,levels = c('fba1','FBA1_Sce','FBA1_Kma','FBA1_Sstip'))
# calculate mean and sd
fba$meanAllocation <- rowMeans(fba[,4:6])
fba$sd <- apply(fba[,4:6],1,sd)
df1 <- melt(fba[,c(1,4:7)],id.vars = c('gene','organism'))

figS6D <- ggplot(fba, aes(x = gene, y = meanAbundance*1000, fill=organism)) +
  geom_bar(position=position_dodge2(), stat="identity",color='black',alpha=0.85) +
  geom_errorbar(aes(ymin=(meanAbundance-sd)*1000,ymax=(meanAbundance+sd)*1000),width=0.4,
                position=position_dodge(0.9)) +
  geom_jitter(data=df1,aes(x=gene,y=value*1000,color=organism),
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),size=0.5) +
  scale_fill_manual(values=c('#CBBBA0','#000000','#1D71B8','#878787')) +
  scale_color_manual(values=c('#CBBBA0','#000000','#1D71B8','#878787')) +
  theme_classic() +
  theme(plot.title = element_text(size = 10,face = 'bold',color = 'black',hjust = 0.5),
        axis.text.x = element_text(angle = 90),
        axis.text = element_text(size = 7,face = 'bold',color='black'),
        axis.title = element_blank(),
        panel.border = element_rect(size = 0.5, color = 'black',fill = NA),
        legend.title = element_blank()) +
  scale_y_continuous(expand = c(0,0),limits = c(0,6.25),breaks = seq(0,6,0.5)) +
  scale_x_discrete(labels=c('fba1','FBA1','FBA1','FBA1')) +
  ggtitle('FBA')
figS6D

# fig S6E triose phosphate isomerase
tpi <- glycolysisData[grep('Triose phosphate isomerase',glycolysisData$Rxn),]
tpi$gene <- factor(tpi$gene,levels = c('tpi1','TPI1_Sce','TPI1_Kma','TPI1_Sstip'))
# calculate mean and sd
tpi$meanAbundance <- rowMeans(tpi[,4:6])
tpi$sd <- apply(tpi[,4:6],1,sd)
df1 <- melt(tpi[,c(1,4:7)],id.vars = c('gene','organism'))

figS6E <- ggplot(tpi, aes(x = gene, y = meanAbundance*1000, fill=organism)) +
  geom_bar(position=position_dodge2(), stat="identity",color='black',alpha=0.85) +
  geom_errorbar(aes(ymin=(meanAbundance-sd)*1000,ymax=(meanAbundance+sd)*1000),width=0.4,
                position=position_dodge(0.9)) +
  geom_jitter(data=df1,aes(x=gene,y=value*1000,color=organism),
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),size=0.5) +
  scale_fill_manual(values=c('#CBBBA0','#000000','#1D71B8','#878787')) +
  scale_color_manual(values=c('#CBBBA0','#000000','#1D71B8','#878787')) +
  theme_classic() +
  theme(plot.title = element_text(size = 10,face = 'bold',color = 'black',hjust = 0.5),
        axis.text.x = element_text(angle = 90),
        axis.text = element_text(size = 7,face = 'bold',color='black'),
        axis.title = element_blank(),
        panel.border = element_rect(size = 0.5, color = 'black',fill = NA),
        legend.title = element_blank()) +
  scale_y_continuous(expand = c(0,0),limits = c(0,4),breaks = seq(0,4,0.5)) +
  scale_x_discrete(labels=c('tpi1','TPI1','TPI1','TPI1')) +
  ggtitle('TPI')
figS6E

# fig S6F Glyceraldehyde 3-phosphate dehudrogenase
GA3PDH <- glycolysisData[grep('Glyceraldehyde 3-phosphate hydrogenase',glycolysisData$Rxn),]
GA3PDH$gene <- factor(GA3PDH$gene,levels = c('tdh1','gpd3','TDH1_Sce','TDH2_Sce','TDH3_Sce',
                                             'GAP1','GAP3','GPD2','TDH1_Sstip','TDH2_Sstip','TDH3_Sstip'))
# calculate mean and sd
GA3PDH$meanAbundance <- rowMeans(GA3PDH[,4:6])
GA3PDH$sd <- apply(GA3PDH[,4:6],1,sd)
df1 <- melt(GA3PDH[,c(1,4:7)],id.vars = c('gene','organism'))

figS6F <- ggplot(GA3PDH, aes(x = gene, y = meanAbundance*1000, fill=organism)) +
  geom_bar(position=position_dodge2(), stat="identity",color='black',alpha=0.85) +
  geom_errorbar(aes(ymin=(meanAbundance-sd)*1000,ymax=(meanAbundance+sd)*1000),width=0.4,
                position=position_dodge(0.9)) +
  geom_jitter(data=df1,aes(x=gene,y=value*1000,color=organism),
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),size=0.5) +
  scale_fill_manual(values=c('#CBBBA0','#000000','#1D71B8','#878787')) +
  scale_color_manual(values=c('#CBBBA0','#000000','#1D71B8','#878787')) +
  theme_classic() +
  theme(plot.title = element_text(size = 10,face = 'bold',color = 'black',hjust = 0.5),
        axis.text.x = element_text(angle = 90),
        axis.text = element_text(size = 7,face = 'bold',color='black'),
        axis.title = element_blank(),
        panel.border = element_rect(size = 0.5, color = 'black',fill = NA),
        legend.title = element_blank()) +
  scale_y_continuous(expand = c(0,0),limits = c(0,17.5),breaks = seq(0,17,1.5)) +
  scale_x_discrete(labels=c('tdh1','gpd3','TDH1','TDH2','TDH3',
                            'GAP1','GAP3','GPD2','TDH1','TDH2','TDH3')) +
  ggtitle('GA3PDH')
figS6F

# fig S6G phosphoglycerate kinase
pgk <- glycolysisData[grep('Phosphoglycerate kinase',glycolysisData$Rxn),]
pgk$gene <- factor(pgk$gene,levels = c('pgk1','PGK1_Sce','PGK','PGK1_Sstip'))
# calculate mean and sd
pgk$meanAbundance <- rowMeans(pgk[,4:6])
pgk$sd <- apply(pgk[,4:6],1,sd)
df1 <- melt(pgk[,c(1,4:7)],id.vars = c('gene','organism'))

figS6G <- ggplot(pgk, aes(x = gene, y = meanAbundance*1000, fill=organism)) +
  geom_bar(position=position_dodge2(), stat="identity",color='black',alpha=0.85) +
  geom_errorbar(aes(ymin=(meanAbundance-sd)*1000,ymax=(meanAbundance+sd)*1000),width=0.4,
                position=position_dodge(0.9)) +
  geom_jitter(data=df1,aes(x=gene,y=value*1000,color=organism),
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),size=0.5) +
  scale_fill_manual(values=c('#CBBBA0','#000000','#1D71B8','#878787')) +
  scale_color_manual(values=c('#CBBBA0','#000000','#1D71B8','#878787')) +
  theme_classic() +
  theme(plot.title = element_text(size = 10,face = 'bold',color = 'black',hjust = 0.5),
        axis.text.x = element_text(angle = 90),
        axis.text = element_text(size = 7,face = 'bold',color='black'),
        axis.title = element_blank(),
        panel.border = element_rect(size = 0.5, color = 'black',fill = NA),
        legend.title = element_blank()) +
  scale_y_continuous(expand = c(0,0),limits = c(0,11.5),breaks = seq(0,11.5,1.5)) +
  scale_x_discrete(labels=c('pgk1','PGK1','PGK','PGK1')) +
  ggtitle('PGK')
figS6G

# fig S6H phosphoglycerate mutase
gpm <-  glycolysisData[grep('Phosphoglycerate mutase',glycolysisData$Rxn),]
gpm$gene <- factor(gpm$gene,levels = c('gpm1','GPM1_Sce','GPM2_Sce','GPM3_Sce',
                                       'GPM1_Kma','GPM3_Kma','GPM1.1','GPM1.2'))
# calculate mean and sd
gpm$meanAbundance <- rowMeans(gpm[,4:6])
gpm$sd <- apply(gpm[,4:6],1,sd)
df1 <- melt(gpm[,c(1,4:7)],id.vars = c('gene','organism'))

figS6H <- ggplot(gpm, aes(x = gene, y = meanAbundance*1000, fill=organism)) +
  geom_bar(position=position_dodge2(), stat="identity",color='black',alpha=0.85) +
  geom_errorbar(aes(ymin=(meanAbundance-sd)*1000,ymax=(meanAbundance+sd)*1000),width=0.4,
                position=position_dodge(0.9)) +
  geom_jitter(data=df1,aes(x=gene,y=value*1000,color=organism),
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),size=0.5) +
  scale_fill_manual(values=c('#CBBBA0','#000000','#1D71B8','#878787')) +
  scale_color_manual(values=c('#CBBBA0','#000000','#1D71B8','#878787')) +
  theme_classic() +
  theme(plot.title = element_text(size = 10,face = 'bold',color = 'black',hjust = 0.5),
        axis.text.x = element_text(angle = 90),
        axis.text = element_text(size = 7,face = 'bold',color='black'),
        axis.title = element_blank(),
        panel.border = element_rect(size = 0.5, color = 'black',fill = NA),
        legend.title = element_blank()) +
  scale_y_continuous(expand = c(0,0),limits = c(0,3.8),breaks = seq(0,3.5,0.5)) +
  scale_x_discrete(labels=c('gpm1','GPM1','GPM2','GPM3',
                            'GPM1','GPM3','GPM1.1','GPM1.2')) +
  ggtitle('GPM')
figS6H

# fig S6I enolase
enolase <- glycolysisData[grep('Enolase',glycolysisData$Rxn),]
enolase$gene <- factor(enolase$gene, levels = c('eno101','eno102','ENO1','ENO2','ENO','ENO1_Sstip'))
# calculate mean and sd
enolase$meanAbundance <- rowMeans(enolase[4:6])
enolase$sd <- apply(enolase[,4:6],1,sd)
df1 <- melt(enolase[,c(1,4:7)],id.vars = c('gene','organism'))

figS6I <- ggplot(enolase, aes(x = gene, y = meanAbundance*1000, fill=organism)) +
  geom_bar(position=position_dodge2(), stat="identity",color='black',alpha=0.85) +
  geom_errorbar(aes(ymin=(meanAbundance-sd)*1000,ymax=(meanAbundance+sd)*1000),width=0.4,
                position=position_dodge(0.9)) +
  geom_jitter(data=df1,aes(x=gene,y=value*1000,color=organism),
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),size=0.5) +
  scale_fill_manual(values=c('#CBBBA0','#000000','#1D71B8','#878787')) +
  scale_color_manual(values=c('#CBBBA0','#000000','#1D71B8','#878787')) +
  theme_classic() +
  theme(plot.title = element_text(size = 10,face = 'bold',color = 'black',hjust = 0.5),
        axis.text.x = element_text(angle = 90),
        axis.text = element_text(size = 7,face = 'bold',color='black'),
        axis.title = element_blank(),
        panel.border = element_rect(size = 0.5, color = 'black',fill = NA),
        legend.title = element_blank()) +
  scale_y_continuous(expand = c(0,0),limits = c(0,17),breaks = seq(0,17,1)) +
  scale_x_discrete(labels=c('eno101','eno102','ENO1','ENO2','ENO','ENO1')) +
  ggtitle('ENO')
figS6I

# fig S6J pyruvate kinase
pyk <- glycolysisData[grep('Pyruvate kinase',glycolysisData$Rxn),]
pyk$gene <- factor(pyk$gene, levels = c('pyk1','CDC19','PYK2','PYK1_Kma','PYK1_Sstip'))
# calculate mean and sd
pyk$meanAbundance <- rowMeans(pyk[,4:6])
pyk$sd <- apply(pyk[,4:6],1,sd)
df1 <- melt(pyk[,c(1,4:7)],id.vars = c('gene','organism'))

figS6J <- ggplot(pyk, aes(x = gene, y = meanAbundance*1000, fill=organism)) +
  geom_bar(position=position_dodge2(), stat="identity",color='black',alpha=0.85) +
  geom_errorbar(aes(ymin=(meanAbundance-sd)*1000,ymax=(meanAbundance+sd)*1000),width=0.4,
                position=position_dodge(0.9)) +
  geom_jitter(data=df1,aes(x=gene,y=value*1000,color=organism),
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),size=0.5) +
  scale_fill_manual(values=c('#CBBBA0','#000000','#1D71B8','#878787')) +
  scale_color_manual(values=c('#CBBBA0','#000000','#1D71B8','#878787')) +
  theme_classic() +
  theme(plot.title = element_text(size = 10,face = 'bold',color = 'black',hjust = 0.5),
        axis.text.x = element_text(angle = 90),
        axis.text = element_text(size = 7,face = 'bold',color='black'),
        axis.title = element_blank(),
        panel.border = element_rect(size = 0.5, color = 'black',fill = NA),
        legend.title = element_blank()) +
  scale_y_continuous(expand = c(0,0),limits = c(0,10),breaks = seq(0,9.5,1.5)) +
  scale_x_discrete(labels=c('pyk1','CDC19','PYK2','PYK1','PYK1')) +
  ggtitle('PYK')
figS6J

# combine in one plot
figS6 <- ggarrange(figS6A,figS6B,figS6C,figS6D,figS6E,figS6F,figS6G,figS6H,figS6I,figS6J,
                   ncol = 5,nrow = 2, common.legend = TRUE,legend = 'bottom')
figS6

# Save figure
#ggsave('../../results/figS6.pdf',
 #      figS6,width=16,height = 16,units = 'cm',device = 'pdf',dpi=300)

# fig S7
PPP_data <- read.delim('PPP_all.txt',sep = '\t',stringsAsFactors = FALSE)
# Convert rxn and organism columns to factor
PPP_data$Rxn <- as.factor(PPP_data$Rxn)
PPP_data$organism <- factor(PPP_data$organism,levels = c('Spo','Sce','Kma','Sstip'))
# Calculate mean and standard deviation
PPP_data$meanAbundance = rowMeans(PPP_data[,5:7])
PPP_data$sd <- apply(PPP_data[,5:7],1,sd)

# fig S7A
# zwf
zwf <- PPP_data[grep('ZWF',PPP_data$Rxn),]
zwf$PrimaryGeneName <- factor(zwf$PrimaryGeneName,levels = zwf$PrimaryGeneName)
# calculate mean and sd
zwf$meanAbundance <- rowMeans(zwf[,5:7])
zwf$sd <- apply(zwf[,5:7],1,sd)
df1 <- melt(zwf[,c(3,5:8)],id.vars = c('PrimaryGeneName','organism'))

figS7A <- ggplot(zwf, aes(x = PrimaryGeneName, y = meanAbundance*1000, fill=organism)) +
  geom_bar(position=position_dodge2(), stat="identity",color='black',alpha=0.85) +
  geom_errorbar(aes(ymin=1000*(meanAbundance-sd), ymax=1000*(meanAbundance+sd)), width=.2,
                position=position_dodge(.9)) +
  geom_jitter(data=df1,aes(x=PrimaryGeneName,y=value*1000,color=organism),
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),size=0.5) +
  scale_fill_manual(values=c('#CBBBA0','#000000','#1D71B8','#878787')) +
  scale_color_manual(values=c('#CBBBA0','#000000','#1D71B8','#878787')) +
  theme_classic() +
  theme(plot.title = element_text(size = 10,face = 'bold',color = 'black',hjust = 0.5),
        axis.text.x = element_text(angle = 30,vjust=0.6),
        axis.text = element_text(size = 7,face = 'bold',color='black'),
        axis.title = element_blank(),
        panel.border = element_rect(size = 1.1, color = 'black',fill = NA),
        legend.title = element_blank()) +
  scale_y_continuous(expand = c(0,0),limits = c(0,0.25),breaks = seq(0,0.25,0.05)) +
  scale_x_discrete(labels=c('zwf1','SPAC3C7.13c','SPCC794.01c','ZWF1','ZWF','ZWF1')) +
  ggtitle('ZWF')
figS7A

# fig S7B
# pgl
pgl <- PPP_data[grep('PGL',PPP_data$Rxn),]
pgl$PrimaryGeneName <- factor(pgl$PrimaryGeneName,levels = pgl$PrimaryGeneName)
# calculate mean and sd
pgl$meanAbundance <- rowMeans(pgl[,5:7])
pgl$sd <- apply(pgl[,5:7],1,sd)
df1 <- melt(pgl[,c(3,5:8)],id.vars=c('PrimaryGeneName','organism'))

figS7B <- ggplot(pgl, aes(x = PrimaryGeneName, y = meanAbundance*1000, fill=organism)) +
  geom_bar(position=position_dodge2(), stat="identity",color='black',alpha=0.85) +
  geom_errorbar(aes(ymin=1000*(meanAbundance-sd), ymax=1000*(meanAbundance+sd)), width=.2,
                position=position_dodge(.9)) +
  geom_jitter(data=df1,aes(x=PrimaryGeneName,y=value*1000,color=organism),
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),size=0.5) +
  scale_fill_manual(values=c('#CBBBA0','#000000','#1D71B8','#878787')) +
  scale_color_manual(values=c('#CBBBA0','#000000','#1D71B8','#878787')) +
  theme_classic() +
  theme(plot.title = element_text(size = 10,face = 'bold',color = 'black',hjust = 0.5),
        axis.text.x = element_text(angle = 30,vjust=0.6),
        axis.text = element_text(size = 7,face = 'bold',color='black'),
        axis.title = element_blank(),
        panel.border = element_rect(size = 1.1, color = 'black',fill = NA),
        legend.title = element_blank()) +
  scale_y_continuous(expand = c(0,0),limits = c(0,0.1),breaks = seq(0,0.1,0.025)) +
  scale_x_discrete(labels=c('SPCC16C4.10','SOL2','SOL3','SOL1','SOL4',
                            'SOL3','SOL1','SOL2','SOL1')) +
  ggtitle('PGL')
figS7B

# fig S7C
# gnd
gnd <- PPP_data[grep('GND',PPP_data$Rxn),]
gnd$PrimaryGeneName <- factor(gnd$PrimaryGeneName,levels = gnd$PrimaryGeneName)
# calculate mean and sd
gnd$meanAbundance <- rowMeans(gnd[,5:7])
gnd$sd <- apply(gnd[,5:7],1,sd)
df1 <- melt(gnd[,c(3,5:8)],id.vars = c('PrimaryGeneName','organism'))

figS7C <- ggplot(gnd, aes(x = PrimaryGeneName, y = meanAbundance*1000, fill=organism)) +
  geom_bar(position=position_dodge2(), stat="identity",color='black',alpha=0.85) +
  geom_errorbar(aes(ymin=1000*(meanAbundance-sd), ymax=1000*(meanAbundance+sd)), width=.2,
                position=position_dodge(.9)) +
  geom_jitter(data=df1,aes(x=PrimaryGeneName,y=value*1000,color=organism),
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),size=0.5) +
  scale_fill_manual(values=c('#CBBBA0','#000000','#1D71B8','#878787')) +
  scale_color_manual(values=c('#CBBBA0','#000000','#1D71B8','#878787')) +
  theme_classic() +
  theme(plot.title = element_text(size = 10,face = 'bold',color = 'black',hjust = 0.5),
        axis.text.x = element_text(angle = 30,vjust=0.6),
        axis.text = element_text(size = 7,face = 'bold',color='black'),
        axis.title = element_blank(),
        panel.border = element_rect(size = 1.1, color = 'black',fill = NA),
        legend.title = element_blank()) +
  scale_y_continuous(expand = c(0,0),limits = c(0,1.95),breaks = seq(0,1.8,0.2)) +
  scale_x_discrete(labels=c('SPBC660.16','GND1','GND1','GND1')) +
  ggtitle('GND')
figS7C

# fig S7D
# rki
rki <- PPP_data[grep('RKI',PPP_data$Rxn),]
rki$PrimaryGeneName <- factor(rki$PrimaryGeneName,levels = rki$PrimaryGeneName)
# calculate mean and sd
rki$meanAbundance <- rowMeans(rki[,5:7])
rki$sd <- apply(rki[,5:7],1,sd)
df1 <- melt(rki[,c(3,5:8)],id.vars = c('PrimaryGeneName','organism'))

figS7D <- ggplot(rki, aes(x = PrimaryGeneName, y = meanAbundance*1000, fill=organism)) +
  geom_bar(position=position_dodge2(), stat="identity",color='black',alpha=0.85) +
  geom_errorbar(aes(ymin=1000*(meanAbundance-sd), ymax=1000*(meanAbundance+sd)), width=.2,
                position=position_dodge(.9)) +
  geom_jitter(data=df1,aes(x=PrimaryGeneName,y=value*1000,color=organism),
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),size=0.5) +
  scale_fill_manual(values=c('#CBBBA0','#000000','#1D71B8','#878787')) +
  scale_color_manual(values=c('#CBBBA0','#000000','#1D71B8','#878787')) +
  theme_classic() +
  theme(plot.title = element_text(size = 10,face = 'bold',color = 'black',hjust = 0.5),
        axis.text.x = element_text(angle = 0),
        axis.text = element_text(size = 7,face = 'bold',color='black'),
        axis.title = element_blank(),
        panel.border = element_rect(size = 1.1, color = 'black',fill = NA),
        legend.title = element_blank()) +
  scale_y_continuous(expand = c(0,0),limits = c(0,0.06),breaks = seq(0,0.06,0.01)) +
  scale_x_discrete(labels=c('rki1','RKI1','RKI1','RKI1')) +
  ggtitle('RKI')
figS7D

# fig S7E
# rpe
rpe <- PPP_data[grep('RPE',PPP_data$Rxn),]
rpe$PrimaryGeneName <- factor(rpe$PrimaryGeneName,levels = rpe$PrimaryGeneName)
# calculate mean and sd
rpe$meanAbundance <- rowMeans(rpe[,5:7])
rpe$sd <- apply(rpe[,5:7],1,sd)
df1 <- melt(rpe[,c(3,5:8)],id.vars = c('PrimaryGeneName','organism'))

figS7E <- ggplot(rpe, aes(x = PrimaryGeneName, y = meanAbundance*1000, fill=organism)) +
  geom_bar(position=position_dodge2(), stat="identity",color='black',alpha=0.85) +
  geom_errorbar(aes(ymin=1000*(meanAbundance-sd), ymax=1000*(meanAbundance+sd)), width=.2,
                position=position_dodge(.9)) +
  geom_jitter(data=df1,aes(x=PrimaryGeneName,y=value*1000,color=organism),
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),size=0.5) +
  scale_fill_manual(values=c('#CBBBA0','#000000','#1D71B8','#878787')) +
  scale_color_manual(values=c('#CBBBA0','#000000','#1D71B8','#878787')) +
  theme_classic() +
  theme(plot.title = element_text(size = 10,face = 'bold',color = 'black',hjust = 0.5),
        axis.text.x = element_text(angle = 30,vjust=0.6),
        axis.text = element_text(size = 7,face = 'bold',color='black'),
        axis.title = element_blank(),
        panel.border = element_rect(size = 1.1, color = 'black',fill = NA),
        legend.title = element_blank()) +
  scale_y_continuous(expand = c(0,0),limits = c(0,0.045),breaks = seq(0,0.045,0.005)) +
  scale_x_discrete(labels=c('SPAC31G5.05c','RPE1','RPE1','RPE1')) +
  ggtitle('RPE')
figS7E

# fig S7F
# tkl
tkl <- PPP_data[grep('TKL',PPP_data$Rxn),]
tkl$PrimaryGeneName <- factor(tkl$PrimaryGeneName,levels = tkl$PrimaryGeneName)
# calculate mean and sd
tkl$meanAbundance <- rowMeans(tkl[,5:7])
tkl$sd <- apply(tkl[5:7],1,sd)
df1 <- melt(tkl[,c(3,5:8)],id.vars = c('PrimaryGeneName','organism'))

figS7F <- ggplot(tkl, aes(x = PrimaryGeneName, y = meanAbundance*1000, fill=organism)) +
  geom_bar(position=position_dodge2(), stat="identity",color='black',alpha=0.85) +
  geom_errorbar(aes(ymin=1000*(meanAbundance-sd), ymax=1000*(meanAbundance+sd)), width=.2,
                position=position_dodge(.9)) +
  geom_jitter(data=df1,aes(x=PrimaryGeneName,y=value*1000,color=organism),
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),size=0.5) +
  scale_fill_manual(values=c('#CBBBA0','#000000','#1D71B8','#878787')) +
  scale_color_manual(values=c('#CBBBA0','#000000','#1D71B8','#878787')) +
  theme_classic() +
  theme(plot.title = element_text(size = 10,face = 'bold',color = 'black',hjust = 0.5),
        axis.text.x = element_text(angle = 30,vjust=0.6),
        axis.text = element_text(size = 7,face = 'bold',color='black'),
        axis.title = element_blank(),
        panel.border = element_rect(size = 1.1, color = 'black',fill = NA),
        legend.title = element_blank()) +
  scale_y_continuous(expand = c(0,0),limits = c(0,3.82),breaks = seq(0,3.5,0.5)) +
  scale_x_discrete(labels=c('SPBC2G5.05','TKL1','TKL2','TKL1','TKT')) +
  ggtitle('TKL')
figS7F

# fig S7G
# tal
tal <- PPP_data[grep('TAL',PPP_data$Rxn),]
tal$PrimaryGeneName <- factor(tal$PrimaryGeneName,levels = tal$PrimaryGeneName)
# calculate mean and sd
tal$meanAbundance <- rowMeans(tal[,5:7])
tal$sd <- apply(tal[,5:7],1,sd)
df1 <- melt(tal[,c(3,5:8)],id.vars = c('PrimaryGeneName','organism'))

figS7G <- ggplot(tal, aes(x = PrimaryGeneName, y = meanAbundance*1000, fill=organism)) +
  geom_bar(position=position_dodge2(), stat="identity",color='black',alpha=0.85) +
  geom_errorbar(aes(ymin=1000*(meanAbundance-sd), ymax=1000*(meanAbundance+sd)), width=.2,
                position=position_dodge(.9)) +
  geom_jitter(data=df1,aes(x=PrimaryGeneName,y=value*1000,color=organism),
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),size=0.5) +
  scale_fill_manual(values=c('#CBBBA0','#000000','#1D71B8','#878787')) +
  scale_color_manual(values=c('#CBBBA0','#000000','#1D71B8','#878787')) +
  theme_classic() +
  theme(plot.title = element_text(size = 10,face = 'bold',color = 'black',hjust = 0.5),
        axis.text.x = element_text(angle = 0),
        axis.text = element_text(size = 7,face = 'bold',color='black'),
        axis.title = element_blank(),
        panel.border = element_rect(size = 1.1, color = 'black',fill = NA),
        legend.title = element_blank()) +
  scale_y_continuous(expand = c(0,0),limits = c(0,1.6),breaks = seq(0,1.6,0.25)) +
  scale_x_discrete(labels=c('tal1','TAL1','NQM1','TAL1','TAL1')) +
  ggtitle('TAL')
figS7G

# Combine into one plot
figS7 <- ggarrange(figS7A,figS7B,figS7C,figS7D,figS7E,figS7F,figS7G,
                   nrow = 3,ncol = 3,common.legend = TRUE,legend = 'bottom')
figS7

# Save figure
#ggsave('../../results/figS7.pdf',
 #      figS7,width=16,height = 14,units = 'cm',device = 'pdf',dpi=300)

# Fig S8
# S. cerevisiae
TCAcycle_Sce_absolute <- map2GO(protData_Sce_absolute,GOannotation_Sce,'GO:0006099')
TCAcycle_Sce_absolute$GOtermName <- 'TCA cycle'
# Remove ICL1, which is part of the glyoxylate cycle
TCAcycle_Sce_absolute <- TCAcycle_Sce_absolute[-grep('ICL1',TCAcycle_Sce_absolute$PrimaryGeneName),]
TCAcycle_Sce_absolute <- data.frame(cbind(t(colSums(TCAcycle_Sce_absolute[,5:7])),TCAcycle_Sce_absolute[1,9]))
colnames(TCAcycle_Sce_absolute) <- c('repl1','repl2','repl3','Process')
ETC_Sce_absolute <- map2GO(protData_Sce_absolute,GOannotation_Sce,c('GO:0005750','GO:0005751',
                                                                    'GO:0004129','GO:0008137',
                                                                    'GO:0003954','GO:0006122',
                                                                    'GO:0006123'))
ETC_Sce_absolute <- ETC_Sce_absolute[!duplicated(ETC_Sce_absolute$Accession),]
# Remove NAR1 and COQ11, which are not actually part of the ETC
ETC_Sce_absolute <- ETC_Sce_absolute[-grep('NAR1',ETC_Sce_absolute$PrimaryGeneName),]
ETC_Sce_absolute <- ETC_Sce_absolute[-grep('COQ11',ETC_Sce_absolute$PrimaryGeneName),]
ETC_Sce_absolute$GOtermName <- 'ETC'
ETC_Sce_absolute <- data.frame(cbind(t(colSums(ETC_Sce_absolute[,5:7])),ETC_Sce_absolute[1,9]))
colnames(ETC_Sce_absolute) <- c('repl1','repl2','repl3','Process')
ATPsynth_Sce_absolute <- map2GO(protData_Sce_absolute,GOannotation_Sce,'GO:0015986')
ATPsynth_Sce_absolute$GOtermName <- 'ATP synthase'
ATPsynth_Sce_absolute <- data.frame(cbind(t(colSums(ATPsynth_Sce_absolute[,5:7])),ATPsynth_Sce_absolute[1,9]))
colnames(ATPsynth_Sce_absolute) <- c('repl1','repl2','repl3','Process')
oxphos_Sce_absolute <- rbind(TCAcycle_Sce_absolute,ETC_Sce_absolute,ATPsynth_Sce_absolute)
oxphos_Sce_absolute$Organism <- c(rep('Sce',3))
oxphos_Sce_absolute[,1:3] <- sapply(oxphos_Sce_absolute[,1:3],as.numeric)
oxphos_Sce_absolute$meanAllocation <- rowMeans(oxphos_Sce_absolute[,1:3])
oxphos_Sce_absolute$sd <- apply(oxphos_Sce_absolute[,1:3], 1, sd)

# S. pombe
TCAcycle_Spo_absolute <- map2GO(protData_Spo_absolute,GOannotation_Spo,'GO:0006099')
TCAcycle_Spo_absolute$GOtermName <- 'TCA cycle'
TCAcycle_Spo_absolute <- data.frame(cbind(t(colSums(TCAcycle_Spo_absolute[,5:7])),TCAcycle_Spo_absolute[1,9]))
colnames(TCAcycle_Spo_absolute) <- c('repl1','repl2','repl3','Process')
ETC_Spo_absolute <- map2GO(protData_Spo_absolute,GOannotation_Spo,c('GO:0005750','GO:0005751',
                                                                    'GO:0004129','GO:0008137',
                                                                    'GO:0003954','GO:0006122',
                                                                    'GO:0006123'))
ETC_Spo_absolute <- ETC_Spo_absolute[!duplicated(ETC_Spo_absolute$Accession),]

ETC_Spo_absolute$GOtermName <- 'ETC'
ETC_Spo_absolute <- data.frame(cbind(t(colSums(ETC_Spo_absolute[,5:7])),ETC_Spo_absolute[1,9]))
colnames(ETC_Spo_absolute) <- c('repl1','repl2','repl3','Process')
ATPsynth_Spo_absolute <- map2GO(protData_Spo_absolute,GOannotation_Spo,'GO:0015986')
ATPsynth_Spo_absolute$GOtermName <- 'ATP synthase'
ATPsynth_Spo_absolute <- data.frame(cbind(t(colSums(ATPsynth_Spo_absolute[,5:7])),ATPsynth_Spo_absolute[1,9]))
colnames(ATPsynth_Spo_absolute) <- c('repl1','repl2','repl3','Process')
oxphos_Spo_absolute <- rbind(TCAcycle_Spo_absolute,ETC_Spo_absolute,ATPsynth_Spo_absolute)
oxphos_Spo_absolute$Organism <- c(rep('Spo',3))
oxphos_Spo_absolute[,1:3] <- sapply(oxphos_Spo_absolute[,1:3],as.numeric)
oxphos_Spo_absolute$meanAllocation <- rowMeans(oxphos_Spo_absolute[,1:3])
oxphos_Spo_absolute$sd <- apply(oxphos_Spo_absolute[,1:3], 1, sd)

# S. stipitis
TCAcycle_Sstip_absolute <- map2GO(protData_Sstip_absolute,GOannotation_Sstip,'GO:0006099')
TCAcycle_Sstip_absolute$GOtermName <- 'TCA cycle'
TCAcycle_Sstip_absolute <- data.frame(cbind(t(colSums(TCAcycle_Sstip_absolute[,5:7])),TCAcycle_Sstip_absolute[1,9]))
colnames(TCAcycle_Sstip_absolute) <- c('repl1','repl2','repl3','Process')
ETC_Sstip_absolute <- map2GO(protData_Sstip_absolute,GOannotation_Sstip,c('GO:0005750','GO:0005751',
                                                                          'GO:0004129','GO:0008137',
                                                                          'GO:0003954','GO:0006122',
                                                                          'GO:0006123','GO:0042773'))
ETC_Sstip_absolute <- ETC_Sstip_absolute[!duplicated(ETC_Sstip_absolute$Accession),]

ETC_Sstip_absolute$GOtermName <- 'ETC'
ETC_Sstip_absolute <- data.frame(cbind(t(colSums(ETC_Sstip_absolute[,5:7])),ETC_Sstip_absolute[1,9]))
colnames(ETC_Sstip_absolute) <- c('repl1','repl2','repl3','Process')
ATPsynth_Sstip_absolute <- map2GO(protData_Sstip_absolute,GOannotation_Sstip,'GO:0015986')
ATPsynth_Sstip_absolute$GOtermName <- 'ATP synthase'
ATPsynth_Sstip_absolute <- data.frame(cbind(t(colSums(ATPsynth_Sstip_absolute[,5:7])),ATPsynth_Sstip_absolute[1,9]))
colnames(ATPsynth_Sstip_absolute) <- c('repl1','repl2','repl3','Process')
oxphos_Sstip_absolute <- rbind(TCAcycle_Sstip_absolute,ETC_Sstip_absolute,ATPsynth_Sstip_absolute)
oxphos_Sstip_absolute$Organism <- c(rep('Sstip',3))
oxphos_Sstip_absolute[,1:3] <- sapply(oxphos_Sstip_absolute[,1:3],as.numeric)
oxphos_Sstip_absolute$meanAllocation <- rowMeans(oxphos_Sstip_absolute[,1:3])
oxphos_Sstip_absolute$sd <- apply(oxphos_Sstip_absolute[,1:3], 1, sd)

# K. marxianus
TCAcycle_Kma_absolute <- map2GO(protData_Kma_absolute,GOannotation_Kma,'GO:0006099')
TCAcycle_Kma_absolute$GOtermName <- 'TCA cycle'
TCAcycle_Kma_absolute <- data.frame(cbind(t(colSums(TCAcycle_Kma_absolute[,5:7])),TCAcycle_Kma_absolute[1,9]))
colnames(TCAcycle_Kma_absolute) <- c('repl1','repl2','repl3','Process')
ETC_Kma_absolute <- map2GO(protData_Kma_absolute,GOannotation_Kma,c('GO:0005750','GO:0005751',
                                                                    'GO:0004129','GO:0008137',
                                                                    'GO:0003954','GO:0006122',
                                                                    'GO:0006123','GO:0042773'))
ETC_Kma_absolute <- ETC_Kma_absolute[!duplicated(ETC_Kma_absolute$Accession),]
ETC_Kma_absolute$GOtermName <- 'ETC'
ETC_Kma_absolute <- data.frame(cbind(t(colSums(ETC_Kma_absolute[,5:7])),ETC_Kma_absolute[1,9]))
colnames(ETC_Kma_absolute) <- c('repl1','repl2','repl3','Process')
ATPsynth_Kma_absolute <- map2GO(protData_Kma_absolute,GOannotation_Kma,'GO:0015986')
ATPsynth_Kma_absolute$GOtermName <- 'ATP synthase'
ATPsynth_Kma_absolute <- data.frame(cbind(t(colSums(ATPsynth_Kma_absolute[,5:7])),ATPsynth_Kma_absolute[1,9]))
colnames(ATPsynth_Kma_absolute) <- c('repl1','repl2','repl3','Process')
oxphos_Kma_absolute <- rbind(TCAcycle_Kma_absolute,ETC_Kma_absolute,ATPsynth_Kma_absolute)
oxphos_Kma_absolute$Organism <- c(rep('Kma',3))
oxphos_Kma_absolute[,1:3] <- sapply(oxphos_Kma_absolute[,1:3],as.numeric)
oxphos_Kma_absolute$meanAllocation <- rowMeans(oxphos_Kma_absolute[,1:3])
oxphos_Kma_absolute$sd <- apply(oxphos_Kma_absolute[,1:3], 1, sd)

# Combine data
oxphos_abs <- rbind(oxphos_Sce_absolute,oxphos_Spo_absolute,oxphos_Sstip_absolute,oxphos_Kma_absolute)
#oxphos$Organism <- factor(oxphos$Organism,levels = c('Sce','Spo','Kma','Sstip'))
oxphos_abs$rGlu <- c(rep(18.8,3),rep(8.9,3),rep(4.5,3),rep(5.6,3))
ATPS_abs <- oxphos_abs[grepl('ATP synthase',oxphos_abs$Process),]
df1 <- melt(ATPS_abs[,c(1:4,8)],id.vars = c('Process','rGlu'),measure.vars = c('repl1','repl2','repl3'))
ETC_abs <- oxphos_abs[grepl('ETC',oxphos_abs$Process),]
df2 <- melt(ETC_abs[,c(1:4,8)],id.vars = c('Process','rGlu'),measure.vars = c('repl1','repl2','repl3'))
TCA_abs <- oxphos_abs[grepl('TCA cycle',oxphos_abs$Process),]
df3 <- melt(TCA_abs[,c(1:4,8)],id.vars = c('Process','rGlu'),measure.vars = c('repl1','repl2','repl3'))

# Plot figure
figS8A <- ggplot() +
  geom_point(data=ATPS_abs,aes(x=rGlu,y=meanAllocation*1000),color='#9970AB',shape=22,fill='#9970AB') +
  #geom_errorbar(data=ATPS_abs,aes(x=rGlu,ymin=(meanAllocation-sd)*1000,ymax=(meanAllocation+sd)*1000),
  #              width=.2,position=position_dodge(.9)) +
  geom_jitter(data=df1,aes(x=rGlu,y=value*1000),size=2,
              position = position_jitter(0.2),color='#9970AB') +
  geom_point(data=ETC_abs,aes(x=rGlu,y=meanAllocation*1000),color='#6BAED6',shape=22,fill='#6BAED6') +
  #geom_errorbar(data=ETC_abs,aes(x=rGlu,ymin=(meanAllocation-sd)*1000,ymax=(meanAllocation+sd)*1000),
  #              width=.2,position=position_dodge(.9)) +
  geom_jitter(data=df2,aes(x=rGlu,y=value*1000),size=2,
              position = position_jitter(0.2),color='#6BAED6') +
  geom_point(data=TCA_abs,aes(x=rGlu,y=meanAllocation*1000),color='#FD8D3C',shape=22,fill='#FD8D3C') +
  #geom_errorbar(data=TCA_abs,aes(x=rGlu,ymin=(meanAllocation-sd)*1000,ymax=(meanAllocation+sd)*1000),
  #              width=.2,position=position_dodge(.9)) +
  geom_jitter(data=df3,aes(x=rGlu,y=value*1000),size=2,
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
  scale_y_continuous(limits = c(0,13),breaks = seq(0,13,1),expand = c(0,0)) +
  scale_x_continuous(breaks=c(4.5,5.6,8.9,18.8))
figS8A
# Save plot
#ggsave('../../results/figS8A.pdf',
 #      figS8A,width=8,height = 8,units = 'cm',device = 'pdf',dpi=300)

# Fig S8B
# S. cerevisiae
NADH_dehydrogenase_Sce <- map2GO(protData_Sce,GOannotation_Sce,c('GO:0008137','GO:0003954'))
NADH_dehydrogenase_Sce <- NADH_dehydrogenase_Sce[which(NADH_dehydrogenase_Sce$PrimaryGeneName == 'NDE1' | NADH_dehydrogenase_Sce$PrimaryGeneName == 'NDI1'),]
NADH_dehydrogenase_Sce$complex <- c(rep('NADH',nrow(NADH_dehydrogenase_Sce)))
NADH_dehydrogenase_Sce$PrimaryGeneName <- paste(NADH_dehydrogenase_Sce$PrimaryGeneName,'_Sce',sep = '')
# S. pombe
NADH_dehydrogenase_Spo <- map2GO(protData_Spo,GOannotation_Spo,c('GO:0008137','GO:0003954'))
NADH_dehydrogenase_Spo <- NADH_dehydrogenase_Spo[1:2,]
NADH_dehydrogenase_Spo$complex <- c(rep('NADH',nrow(NADH_dehydrogenase_Spo)))
NADH_dehydrogenase_Spo$PrimaryGeneName <- c('nde1','ndi1')

# K. marxianus
NADH_dehydrogenase_Kma <- map2GO(protData_Kma,GOannotation_Kma,c('GO:0008137','GO:0003954'))
NADH_dehydrogenase_Kma <- NADH_dehydrogenase_Kma[1:2,]
NADH_dehydrogenase_Kma$complex <- c(rep('NADH',nrow(NADH_dehydrogenase_Kma)))
NADH_dehydrogenase_Kma$PrimaryGeneName <- paste(NADH_dehydrogenase_Kma$PrimaryGeneName,'_Kma',sep = '')
# S. stipitis
NADH_dehydrogenase_Sstip <- map2GO(protData_Sstip,GOannotation_Sstip,c('GO:0008137','GO:0003954'))
NADH_dehydrogenase_Sstip <- NADH_dehydrogenase_Sstip[grep('NDE1',NADH_dehydrogenase_Sstip$PrimaryGeneName),]
complexIgenes <- c('PICST_85822','PICST_68160','PICST_74163','PICST_62206','PICST_81272',
                   'PICST_58506','PICST_75518','PICST_88630','PICST_31688','PICST_76559',
                   'PICST_65836','PICST_46630','PICST_84570','PICST_82538','PICST_76018',
                   'PICST_44890','PICST_76226','PICST_63166','PICST_63376','PICST_80638',
                   'PICST_90638','PICST_79236')
complexI_Sstip <- mapGenesToProtData(complexIgenes,protData_Sstip)
NADH_dehydrogenase_Sstip <- rbind(NADH_dehydrogenase_Sstip[,1:7],complexI_Sstip)
NADH_dehydrogenase_Sstip$complex <- c(rep('NADH',nrow(NADH_dehydrogenase_Sstip)))

# Complex II
CII_Sce <- map2GO(protData_Sce,GOannotation_Sce,c('GO:0005749'))
CII_Sce$complex <- c(rep('CII',nrow(CII_Sce)))

CII_Spo <- map2GO(protData_Spo,GOannotation_Spo,c('GO:0005749'))

CII_Kma <- map2GO(protData_Kma,GOannotation_Kma,c('GO:0005749'))
CII_Kma <- CII_Kma[which(CII_Kma$PrimaryGeneName != 'FMP21' & CII_Kma$PrimaryGeneName != 'TIM18'),]

CII_Sstip <- map2GO(protData_Sstip,GOannotation_Sstip,c('GO:0005749'))
CII_Sstip <- CII_Sstip[-grep('A3LUM4',CII_Sstip$Accession),]

# Complex III
CIII_Sce <- map2GO(protData_Sce,GOannotation_Sce,c('GO:0005750','GO:0006122'))
CIII_Sce <- CIII_Sce[-grep('MAS1',CIII_Sce$PrimaryGeneName),]

CIII_Spo <- map2GO(protData_Spo,GOannotation_Spo,c('GO:0005750','GO:0006122'))

CIII_Kma <- map2GO(protData_Kma,GOannotation_Kma,c('GO:0005750','GO:0006122'))
CIII_Kma <- CIII_Kma[-grep('MAS1',CIII_Kma$PrimaryGeneName),]

CIII_Sstip <- map2GO(protData_Sstip,GOannotation_Sstip,c('GO:0005750','GO:0006122'))
CIII_Sstip <- CIII_Sstip[-grep('MAS1',CIII_Sstip$PrimaryGeneName),]

# Complex IV
CIV_Sce <- map2GO(protData_Sce,GOannotation_Sce,c('GO:0005751','GO:0006123','GO:0004129'))

CIV_Spo <- map2GO(protData_Spo,GOannotation_Spo,c('GO:0005751','GO:0006123','GO:0004129'))
CIV_Spo <- CIV_Spo[-grep('P05511',CIV_Spo$Accession),]

CIV_Kma <- map2GO(protData_Kma,GOannotation_Kma,c('GO:0005751','GO:0006123','GO:0004129'))
CIV_Kma <- CIV_Kma[-grep('W0T4B8',CIV_Kma$Accession),]

CIV_Sstip <- map2GO(protData_Sstip,GOannotation_Sstip,c('GO:0005751','GO:0006123','GO:0004129'))

# ATP synthase
ATPsynth_Sce <- map2GO(protData_Sce,GOannotation_Sce,c('GO:0015986'))
ATPsynth_Sce <- ATPsynth_Sce[-grep('VPH1',ATPsynth_Sce$PrimaryGeneName),]
ATPsynth_Sce <- ATPsynth_Sce[-grep('STV1',ATPsynth_Sce$PrimaryGeneName),]
ATPsynth_Sce <- ATPsynth_Sce[-grep('Q3E824',ATPsynth_Sce$Accession),]

ATPsynth_Spo <- map2GO(protData_Spo,GOannotation_Spo,c('GO:0015986'))
ATPsynth_Spo <- ATPsynth_Spo[-grep('vph1',ATPsynth_Spo$PrimaryGeneName),]

ATPsynth_Kma <- map2GO(protData_Kma,GOannotation_Kma,c('GO:0015986'))
ATPsynth_Kma <- ATPsynth_Kma[-grep('VPH1',ATPsynth_Kma$PrimaryGeneName),]
ATPsynth_Kma <- ATPsynth_Kma[-grep('STV1',ATPsynth_Kma$PrimaryGeneName),]

ATPsynth_Sstip <- map2GO(protData_Sstip,GOannotation_Sstip,c('GO:0015986'))
ATPsynth_Sstip <- ATPsynth_Sstip[-grep('VPH1',ATPsynth_Sstip$PrimaryGeneName),]

# Collect all data
# S. cerevisiae
combined_Sce <- rbind(data.frame(t(colSums(NADH_dehydrogenase_Sce[,5:7]))),
                      data.frame(t(colSums(CII_Sce[,5:7]))),
                      data.frame(t(colSums(CIII_Sce[,5:7]))),
                      data.frame(t(colSums(CIV_Sce[,5:7]))),
                      data.frame(t(colSums(ATPsynth_Sce[,5:7]))))
combined_Sce$complex <- c('NADH','CII','CIII','CIV','ATP synthase')
combined_Sce$organism <- c(rep('Sce',5))
combined_Sce$meanAllocation <- rowMeans(combined_Sce[,1:3])
combined_Sce$sd <- apply(combined_Sce[,1:3], 1, sd)

# S. pombe
combined_Spo <- rbind(data.frame(t(colSums(NADH_dehydrogenase_Spo[,5:7]))),
                      data.frame(t(colSums(CII_Spo[,5:7]))),
                      data.frame(t(colSums(CIII_Spo[,5:7]))),
                      data.frame(t(colSums(CIV_Spo[,5:7]))),
                      data.frame(t(colSums(ATPsynth_Spo[,5:7]))))
combined_Spo$complex <- c('NADH','CII','CIII','CIV','ATP synthase')
combined_Spo$organism <- c(rep('Spo',5))
combined_Spo$meanAllocation <- rowMeans(combined_Spo[,1:3])
combined_Spo$sd <- apply(combined_Spo[,1:3], 1, sd)

# K. marxianus
combined_Kma <- rbind(data.frame(t(colSums(NADH_dehydrogenase_Kma[,5:7]))),
                      data.frame(t(colSums(CII_Kma[,5:7]))),
                      data.frame(t(colSums(CIII_Kma[,5:7]))),
                      data.frame(t(colSums(CIV_Kma[,5:7]))),
                      data.frame(t(colSums(ATPsynth_Kma[,5:7]))))
combined_Kma$complex <- c('NADH','CII','CIII','CIV','ATP synthase')
combined_Kma$organism <- c(rep('Kma',5))
combined_Kma$meanAllocation <- rowMeans(combined_Kma[,1:3])
combined_Kma$sd <- apply(combined_Kma[,1:3], 1, sd)

# S. stipitis
combined_Sstip <- rbind(data.frame(t(colSums(NADH_dehydrogenase_Sstip[,5:7]))),
                        data.frame(t(colSums(CII_Sstip[,5:7]))),
                        data.frame(t(colSums(CIII_Sstip[,5:7]))),
                        data.frame(t(colSums(CIV_Sstip[,5:7]))),
                        data.frame(t(colSums(ATPsynth_Sstip[,5:7]))))
combined_Sstip$complex <- c('NADH','CII','CIII','CIV','ATP synthase')
combined_Sstip$organism <- c(rep('Sstip',5))
combined_Sstip$meanAllocation <- rowMeans(combined_Sstip[,1:3])
combined_Sstip$sd <- apply(combined_Sstip[,1:3], 1, sd)

combined_All <- rbind(combined_Sce,combined_Spo,combined_Kma,combined_Sstip)
combined_All$organism <- factor(combined_All$organism,levels = c('Spo','Sce','Kma','Sstip'))
combined_All$complex <- factor(combined_All$complex,levels = c('NADH','CII','CIII','CIV','ATP synthase'))
df1 <- melt(combined_All[,1:5],id.vars = c('complex','organism'))

# Plot results
figS8B <- ggplot(combined_All,aes(x=complex,y=meanAllocation*100,fill=organism)) +
  geom_bar(position=position_dodge2(), stat="identity",color='black',alpha=0.85) +
  geom_errorbar(aes(ymin=(meanAllocation-sd)*100,ymax=(meanAllocation+sd)*100),width=0.4,
                position=position_dodge(0.9)) +
  geom_jitter(data=df1,aes(x=complex,y=value*100,color=organism),
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),size=1) +
  scale_y_continuous(breaks = seq(0,3,0.5),limits = c(0,3.1),expand = c(0,0)) +
  scale_fill_manual(values=c('#CBBBA0','#1D1D1B','#1D71B8','#878787')) +
  scale_color_manual(values=c('#CBBBA0','#1D1D1B','#1D71B8','#878787')) +
  theme_classic() +
  theme(text = element_text(size = 8, color = 'black',face = 'bold'),
        legend.position = c(0.2,0.8),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10,face = 'bold',color = 'black'),
        axis.text.x = element_text(size = 7,face = 'bold',color = 'black'),
        axis.text.y = element_text(size = 8,face = 'bold',color = 'black')) +
  ylab('Proteome allocation (%)')
figS8B

# Save figure
#ggsave('../../results/figS8B.pdf',
 #      figS8B,width=8,height = 8,units = 'cm',device = 'pdf',dpi=300)

# Fig S9
TCA_data <- read.delim('TCA_cycle_combined.txt',sep = '\t',stringsAsFactors = F)
# Convert rxn and organism columns to factor
TCA_data$Rxn <- as.factor(TCA_data$Rxn)
TCA_data$organism <- factor(TCA_data$organism,levels = c('Spo','Sce','Kma','Sstip'))
# Calculate mean and standard deviation
TCA_data$meanAbundance = rowMeans(TCA_data[,5:7])
TCA_data$sd <- apply(TCA_data[,5:7],1,sd)

# fig S9A
mapGenesToProtData <- function(geneList,protData){
  data <- data.frame(matrix(ncol = 7,nrow = 0))
  for(gene in geneList){
    geneData <- protData[grep(gene,protData$PrimaryGeneName),]
    if(!nrow(geneData)==0){
      data <- rbind(data,geneData)
    }
  }
  return(data)
}
# pdh
protData_Spo_absolute[grep('O94709',protData_Spo_absolute$Accession),3] <- 'pdx1'
pdh_Spo <- c('pda1','pdb1','pdx1','lat1')
pdh_data_Spo <- mapGenesToProtData(pdh_Spo,protData_Spo_absolute)
pdh_data_Spo$organism <- c(rep('Spo',4))
pdh_Sce <- c('LAT1','LPD1','PDA1','PDB1','PDX1')
pdh_data_Sce <- mapGenesToProtData(pdh_Sce,protData_Sce_absolute)
pdh_data_Sce$organism <- c(rep('Sce',5))
pdh_Kma <- c('LAT1','LPD1','PDA1','PDB1','PDX1')
pdh_data_Kma <- mapGenesToProtData(pdh_Kma,protData_Kma_absolute)
pdh_data_Kma$organism <- c(rep('Kma',5))
pdh_Sstip <- c('LAT1','LPD1','PDA1','PDB1','PDX1')
pdh_data_Sstip <- mapGenesToProtData(pdh_Sstip,protData_Sstip_absolute)
pdh_data_Sstip$organism <- c(rep('Sstip',5))
pdh_all <- rbind(pdh_data_Spo,pdh_data_Sce,pdh_data_Kma,pdh_data_Sstip)
pdh_all$meanAbundance <- rowMeans(pdh_all[,5:7])
pdh_all$sd <- apply(pdh_all[,5:7],1,sd)
#names(pdh_all)[names(pdh_all)=='GeneName'] <- 'Rxn'
pdh_all$Rxn <- 'PDH'
pdh_all$organism <- factor(pdh_all$organism,levels = c('Spo','Sce','Kma','Sstip'))
pdh_all$PrimaryGeneName <- c('pda1','pdb1','pdx1','lat1',
                             'LAT1_Sce','LPD1_Sce','PDA1_Sce','PDB1_Sce','PDX1_Sce',
                             'LAT1_Kma','LPD1_Kma','PDA1_Kma','PDB1_Kma','PDX1_Kma',
                             'LAT1_Sstip','LPD1_Sstip','PDA1_Sstip','PDB1_Sstip','PDX1_Sstip')
pdh_all$PrimaryGeneName <- factor(pdh_all$PrimaryGeneName,levels=c('pda1','pdb1','pdx1','lat1',
                                                                   'LAT1_Sce','LPD1_Sce','PDA1_Sce','PDB1_Sce','PDX1_Sce',
                                                                   'LAT1_Kma','LPD1_Kma','PDA1_Kma','PDB1_Kma','PDX1_Kma',
                                                                   'LAT1_Sstip','LPD1_Sstip','PDA1_Sstip','PDB1_Sstip','PDX1_Sstip'))
df1 <- melt(pdh_all[,c(3,5:8)],id.vars = c('PrimaryGeneName','organism'))

figS9A <- ggplot(pdh_all, aes(x = PrimaryGeneName, y = meanAbundance*1000, fill=organism)) +
  geom_bar(position=position_dodge2(), stat="identity",color='black',alpha=0.85) +
  geom_errorbar(aes(ymin=1000*(meanAbundance-sd), ymax=1000*(meanAbundance+sd)), width=.2,
                position=position_dodge(.9)) +
  geom_jitter(data=df1,aes(x=PrimaryGeneName,y=value*1000,color=organism),
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),size=0.5) +
  scale_fill_manual(values=c('#CBBBA0','#000000','#1D71B8','#878787')) +
  scale_color_manual(values=c('#CBBBA0','#000000','#1D71B8','#878787')) +
  theme_classic() +
  theme(plot.title = element_text(size = 10,face = 'bold',color = 'black',hjust = 0.5),
        axis.text.x = element_text(angle = 90),
        axis.text = element_text(size = 7,face = 'bold',color='black'),
        axis.title = element_blank(),
        panel.border = element_rect(size = 1.1, color = 'black',fill = NA),
        legend.title = element_blank(),
        axis.ticks = element_line(color='black')) +
  scale_y_continuous(expand = c(0,0),limits = c(0,1),breaks = seq(0,1,0.25)) +
  scale_x_discrete(labels=c('pda1','pdb1','pdx1','lat1',
                            'LAT1','LPD1','PDA1','PDB1','PDX1',
                            'LAT1','LPD1','PDA1','PDB1','PDX1',
                            'LAT1','LPD1','PDA1','PDB1','PDX1')) +
  ggtitle('PDH')
figS9A
# save figure
#ggsave('../../results/figS9A.pdf',
 #      figS9A,width=7,height = 7,units = 'cm',device = 'pdf',dpi=300)

# fig S9B
# cit
cit <- TCA_data[grep('CIT',TCA_data$Rxn),]
cit$PrimaryGeneName <- factor(cit$PrimaryGeneName, levels = cit$PrimaryGeneName)
df1 <- melt(cit[,c(3,5:8)],id.vars = c('PrimaryGeneName','organism'))

figS9B <- ggplot(cit, aes(x = PrimaryGeneName, y = meanAbundance*1000, fill=organism)) +
  geom_bar(position=position_dodge2(), stat="identity",color='black',alpha=0.85) +
  geom_errorbar(aes(ymin=1000*(meanAbundance-sd), ymax=1000*(meanAbundance+sd)), width=.2,
                position=position_dodge(.9)) +
  geom_jitter(data=df1,aes(x=PrimaryGeneName,y=value*1000,color=organism),
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),size=0.5) +
  scale_fill_manual(values=c('#CBBBA0','#000000','#1D71B8','#878787')) +
  scale_color_manual(values=c('#CBBBA0','#000000','#1D71B8','#878787')) +
  theme_classic() +
  theme(plot.title = element_text(size = 10,face = 'bold',color = 'black',hjust = 0.5),
        axis.text.x = element_text(angle = 90),
        axis.text = element_text(size = 7,face = 'bold',color='black'),
        axis.title = element_blank(),
        panel.border = element_rect(size = 1.1, color = 'black',fill = NA),
        axis.ticks = element_line(color='black'),
        legend.title = element_blank()) +
  scale_y_continuous(expand = c(0,0),limits = c(0,1.85),breaks = seq(0,1.75,0.25)) +
  scale_x_discrete(labels = c('cit1','CIT1','CIT2','CIT1','CIT3','CIT1')) +
  ggtitle('CIT')
figS9B

# fig S9C
aco <- TCA_data[grep('ACO',TCA_data$Rxn),]
aco$PrimaryGeneName <- factor(aco$PrimaryGeneName,levels = aco$PrimaryGeneName)
df1 <- melt(aco[,c(3,5:8)],id.vars = c('PrimaryGeneName','organism'))

figS9C <- ggplot(aco, aes(x = PrimaryGeneName, y = meanAbundance*1000, fill=organism)) +
  geom_bar(position=position_dodge2(), stat="identity",color='black',alpha=0.85) +
  geom_errorbar(aes(ymin=1000*(meanAbundance-sd), ymax=1000*(meanAbundance+sd)), width=.2,
                position=position_dodge(.9)) +
  geom_jitter(data=df1,aes(x=PrimaryGeneName,y=value*1000,color=organism),
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),size=0.5) +
  scale_fill_manual(values=c('#CBBBA0','#000000','#1D71B8','#878787')) +
  scale_color_manual(values=c('#CBBBA0','#000000','#1D71B8','#878787')) +
  theme_classic() +
  theme(plot.title = element_text(size = 10,face = 'bold',color = 'black',hjust = 0.5),
        axis.text.x = element_text(angle = 90),
        axis.text = element_text(size = 7,face = 'bold',color='black'),
        axis.title = element_blank(),
        panel.border = element_rect(size = 1.1, color = 'black',fill = NA),
        axis.ticks = element_line(color='black'),
        legend.title = element_blank()) +
  scale_y_continuous(expand = c(0,0),limits = c(0,3),breaks = seq(0,3,0.5)) +
  scale_x_discrete(labels = c('aco1','aco2','ACO1','ACO2','ACO2a','ACO2b',
                              'ACO2','ACO1')) +
  ggtitle('ACO')
figS9C

# fig S9D
# idh
idh <- TCA_data[grep('IDH',TCA_data$Rxn),]
idh$PrimaryGeneName <- factor(idh$PrimaryGeneName,levels = idh$PrimaryGeneName)
df1 <- melt(idh[,c(3,5:8)],id.vars = c('PrimaryGeneName','organism'))

figS9D <- ggplot(idh, aes(x = PrimaryGeneName, y = meanAbundance*1000, fill=organism)) +
  geom_bar(position=position_dodge2(), stat="identity",color='black',alpha=0.85) +
  geom_errorbar(aes(ymin=1000*(meanAbundance-sd), ymax=1000*(meanAbundance+sd)), width=.2,
                position=position_dodge(.9)) +
  geom_jitter(data=df1,aes(x=PrimaryGeneName,y=value*1000,color=organism),
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),size=0.5) +
  scale_fill_manual(values=c('#CBBBA0','#000000','#1D71B8','#878787')) +
  scale_color_manual(values=c('#CBBBA0','#000000','#1D71B8','#878787')) +
  theme_classic() +
  theme(plot.title = element_text(size = 10,face = 'bold',color = 'black',hjust = 0.5),
        axis.text.x = element_text(angle = 90),
        axis.text = element_text(size = 7,face = 'bold',color='black'),
        axis.title = element_blank(),
        panel.border = element_rect(size = 1.1, color = 'black',fill = NA),
        axis.ticks = element_line(color='black'),
        legend.title = element_blank()) +
  scale_y_continuous(expand = c(0,0),limits = c(0,1.1),breaks = seq(0,1,0.25)) +
  scale_x_discrete(labels = c('idh1','idp1','idh2','IDP1','IDH2','IDH1','IDP2',
                              'IDH2','IDP2','IDP1','IDH1','IDH2','IDP2','IDP1','IDH1')) +
  ggtitle('IDH')
figS9D

# fig S9E
# kgd
kgd <- TCA_data[grep('KGD',TCA_data$Rxn),]
kgd$PrimaryGeneName <- factor(kgd$PrimaryGeneName,levels = kgd$PrimaryGeneName)
df1 <- melt(kgd[,c(3,5:8)],id.vars = c('PrimaryGeneName','organism'))

figS9E <- ggplot(kgd, aes(x = PrimaryGeneName, y = meanAbundance*1000, fill=organism)) +
  geom_bar(position=position_dodge2(), stat="identity",color='black',alpha=0.85) +
  geom_errorbar(aes(ymin=1000*(meanAbundance-sd), ymax=1000*(meanAbundance+sd)), width=.2,
                position=position_dodge(.9)) +
  geom_jitter(data=df1,aes(x=PrimaryGeneName,y=value*1000,color=organism),
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),size=0.5) +
  scale_fill_manual(values=c('#CBBBA0','#000000','#1D71B8','#878787')) +
  scale_color_manual(values=c('#CBBBA0','#000000','#1D71B8','#878787')) +
  theme_classic() +
  theme(plot.title = element_text(size = 10,face = 'bold',color = 'black',hjust = 0.5),
        axis.text.x = element_text(angle = 90),
        axis.text = element_text(size = 7,face = 'bold',color='black'),
        axis.title = element_blank(),
        panel.border = element_rect(size = 1.1, color = 'black',fill = NA),
        axis.ticks = element_line(color='black'),
        legend.title = element_blank()) +
  scale_y_continuous(expand = c(0,0),limits = c(0,0.75),breaks = seq(0,0.75,0.25)) +
  scale_x_discrete(labels = c('kgd1','kgd2','KGD2','KGD1','KGD2','KGD1','KGD2','KGD1')) +
  ggtitle('KGD')
figS9E

# fig S9F
# lsc
lsc <- TCA_data[grep('succinate-CoA ligase',TCA_data$Rxn),]
lsc$PrimaryGeneName <- factor(lsc$PrimaryGeneName,levels = lsc$PrimaryGeneName)
df1 <- melt(lsc[,c(3,5:8)],id.vars = c('PrimaryGeneName','organism'))

figS9F <- ggplot(lsc, aes(x = PrimaryGeneName, y = meanAbundance*1000, fill=organism)) +
  geom_bar(position=position_dodge2(), stat="identity",color='black',alpha=0.85) +
  geom_errorbar(aes(ymin=1000*(meanAbundance-sd), ymax=1000*(meanAbundance+sd)), width=.2,
                position=position_dodge(.9)) +
  geom_jitter(data=df1,aes(x=PrimaryGeneName,y=value*1000,color=organism),
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),size=0.5) +
  scale_fill_manual(values=c('#CBBBA0','#000000','#1D71B8','#878787')) +
  scale_color_manual(values=c('#CBBBA0','#000000','#1D71B8','#878787')) +
  theme_classic() +
  theme(plot.title = element_text(size = 10,face = 'bold',color = 'black',hjust = 0.5),
        axis.text.x = element_text(angle = 90),
        axis.text = element_text(size = 7,face = 'bold',color='black'),
        axis.title = element_blank(),
        panel.border = element_rect(size = 1.1, color = 'black',fill = NA),
        axis.ticks = element_line(color='black'),
        legend.title = element_blank()) +
  scale_y_continuous(expand = c(0,0),limits = c(0,0.8),breaks = seq(0,0.75,0.25)) +
  scale_x_discrete(labels = c('17c','lsc2','LSC2','LSC1','LSC2','LSC1',
                              'LSC1','LSC2')) +
  ggtitle('succinate-CoA ligase')
figS9F

# fig S9G
# sdh
sdh <- TCA_data[grep('SDH',TCA_data$Rxn),]
sdh$PrimaryGeneName <- factor(sdh$PrimaryGeneName,levels = sdh$PrimaryGeneName)
df1 <- melt(sdh[,c(3,5:8)],id.vars = c('PrimaryGeneName','organism'))

figS9G <- ggplot(sdh, aes(x = PrimaryGeneName, y = meanAbundance*1000, fill=organism)) +
  geom_bar(position=position_dodge2(), stat="identity",color='black',alpha=0.85) +
  geom_errorbar(aes(ymin=1000*(meanAbundance-sd), ymax=1000*(meanAbundance+sd)), width=.2,
                position=position_dodge(.9)) +
  geom_jitter(data=df1,aes(x=PrimaryGeneName,y=value*1000,color=organism),
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),size=0.5) +
  scale_fill_manual(values=c('#CBBBA0','#000000','#1D71B8','#878787')) +
  scale_color_manual(values=c('#CBBBA0','#000000','#1D71B8','#878787')) +
  theme_classic() +
  theme(plot.title = element_text(size = 10,face = 'bold',color = 'black',hjust = 0.5),
        axis.text.x = element_text(angle = 90),
        axis.text = element_text(size = 7,face = 'bold',color='black'),
        axis.title = element_blank(),
        panel.border = element_rect(size = 1.1, color = 'black',fill = NA),
        axis.ticks = element_line(color='black'),
        legend.title = element_blank()) +
  scale_y_continuous(expand = c(0,0),limits = c(0,0.6),breaks = seq(0,0.6,0.1)) +
  scale_x_discrete(labels=c('sdh2','sdh4','sdh1','SDH2','SDH3','SDH4','SDH1',
                            'SDH1','3012','SDH4','SDH2','SDH1','SDH6',
                            'SDH4','SDH2')) +
  ggtitle('SDH')
figS9G

# fig S9H
# fum
fum <- TCA_data[grep('FUM',TCA_data$Rxn),]
fum$PrimaryGeneName <- factor(fum$PrimaryGeneName,levels = fum$PrimaryGeneName)
df1 <- melt(fum[,c(3,5:8)],id.vars = c('PrimaryGeneName','organism'))

figS9H <- ggplot(fum, aes(x = PrimaryGeneName, y = meanAbundance*1000, fill=organism)) +
  geom_bar(position=position_dodge2(), stat="identity",color='black',alpha=0.85) +
  geom_errorbar(aes(ymin=1000*(meanAbundance-sd), ymax=1000*(meanAbundance+sd)), width=.2,
                position=position_dodge(.9)) +
  geom_jitter(data=df1,aes(x=PrimaryGeneName,y=value*1000,color=organism),
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),size=0.5) +
  scale_fill_manual(values=c('#CBBBA0','#000000','#1D71B8','#878787')) +
  scale_color_manual(values=c('#CBBBA0','#000000','#1D71B8','#878787')) +
  theme_classic() +
  theme(plot.title = element_text(size = 10,face = 'bold',color = 'black',hjust = 0.5),
        axis.text.x = element_text(angle = 90),
        axis.text = element_text(size = 7,face = 'bold',color='black'),
        axis.title = element_blank(),
        panel.border = element_rect(size = 1.1, color = 'black',fill = NA),
        axis.ticks = element_line(color='black'),
        legend.title = element_blank()) +
  scale_y_continuous(expand = c(0,0),limits = c(0,0.7),breaks = seq(0,0.7,0.1)) +
  scale_x_discrete(labels=c('fum1','FUM1','FUM1','FUM2','FUM1')) +
ggtitle('FUM')
figS9H

# fig S9I
# mdh
mdh <- TCA_data[grep('MDH',TCA_data$Rxn),]
mdh$PrimaryGeneName <- factor(mdh$PrimaryGeneName,levels = mdh$PrimaryGeneName)
df1 <- melt(mdh[,c(3,5:8)],id.vars = c('PrimaryGeneName','organism'))

figS9I <- ggplot(mdh, aes(x = PrimaryGeneName, y = meanAbundance*1000, fill=organism)) +
  geom_bar(position=position_dodge2(), stat="identity",color='black',alpha=0.85) +
  geom_errorbar(aes(ymin=1000*(meanAbundance-sd), ymax=1000*(meanAbundance+sd)), width=.2,
                position=position_dodge(.9)) +
  geom_jitter(data=df1,aes(x=PrimaryGeneName,y=value*1000,color=organism),
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),size=0.5) +
  scale_fill_manual(values=c('#CBBBA0','#000000','#1D71B8','#878787')) +
  scale_color_manual(values=c('#CBBBA0','#000000','#1D71B8','#878787')) +
  theme_classic() +
  theme(plot.title = element_text(size = 10,face = 'bold',color = 'black',hjust = 0.5),
        axis.text.x = element_text(angle = 90),
        axis.text = element_text(size = 7,face = 'bold',color='black'),
        axis.title = element_blank(),
        panel.border = element_rect(size = 1.1, color = 'black',fill = NA),
        axis.ticks = element_line(color='black'),
        legend.title = element_blank()) +
  scale_y_continuous(expand = c(0,0),limits = c(0,2.35),breaks = seq(0,2,0.5)) +
  scale_x_discrete(labels=c('MDH1','MDH1','MDH3','MDH2','MDH3','MDH1',
                            'MDH1','MDHM','MDH2')) +
ggtitle('MDH')
figS9I

# Arrange into a single figure
figS9 <- ggarrange(figS9A,figS9B,figS9C,figS9D,figS9E,
                   figS9F,figS9G,figS9H,figS9I,
                   ncol = 3,nrow = 3,
                   common.legend = TRUE,legend = 'bottom')
figS9

# Save figure
#ggsave('../../results/figS9.pdf',
 #      figS9,width=16,height = 16,units = 'cm',device = 'pdf',dpi=300)

# Fig S10
mapGenesToProtData <- function(geneList,protData){
  data <- data.frame(matrix(ncol = 7,nrow = 0))
  for(gene in geneList){
    geneData <- protData[grep(gene,protData$GeneName),]
    if(!nrow(geneData)==0){
      data <- rbind(data,geneData)
    }
  }
  return(data)
}
# Read data in g/gDW
protData_abs_Sce <- read.delim('proteomicsData_g_per_gDW_Scer.txt',stringsAsFactors = FALSE)
protData_abs_Spo <- read.delim('proteomicsData_g_per_gDW_Spombe.txt',stringsAsFactors = FALSE)
protData_abs_Sstip <- read.delim('proteomicsData_g_per_gDW_Pstip.txt',stringsAsFactors = FALSE)
protData_abs_Kma <- read.delim('proteomicsData_g_per_gDW_Kmarx.txt',stringsAsFactors = FALSE)

# NADH dehydrogenases
# S. stipitis
NADH_dehydrogenase_abs_Sstip <- map2GO(protData_abs_Sstip,GOannotation_Sstip,c('GO:0008137','GO:0003954'))
NADH_dehydrogenase_abs_Sstip <- NADH_dehydrogenase_abs_Sstip[grep('NDE1',NADH_dehydrogenase_abs_Sstip$PrimaryGeneName),]
NADH_dehydrogenase_abs_Sstip <- NADH_dehydrogenase_abs_Sstip[,1:7]
NADH_dehydrogenase_abs_Sstip$component <- 'NDE1'
complexIgenes <- c('PICST_85822','PICST_68160','PICST_74163','PICST_62206','PICST_81272',
                   'PICST_58506','PICST_75518','PICST_88630','PICST_31688','PICST_76559',
                   'PICST_65836','PICST_46630','PICST_84570','PICST_82538','PICST_76018',
                   'PICST_44890','PICST_76226','PICST_63166','PICST_63376','PICST_80638',
                   'PICST_90638','PICST_79236')
complexI_abs_Sstip <- mapGenesToProtData(complexIgenes,protData_abs_Sstip)
NADH_dehydrogenase_abs_Sstip <- rbind(NADH_dehydrogenase_abs_Sstip[,1:7],complexI_abs_Sstip)
NADH_dehydrogenase_abs_Sstip$complex <- c(rep('NADH',nrow(NADH_dehydrogenase_abs_Sstip)))

# Complex IV
CIV_abs_Sstip <- map2GO(protData_abs_Sstip,GOannotation_Sstip,c('GO:0005751','GO:0006123','GO:0004129'))

NDE1 <- c('PICST_66598')
NDI1 <- c('PICST_58800')
AOX <- c('PICST_67332')
AOXdata <- mapGenesToProtData(AOX,protData_abs_Sstip)
AOXdata$component <- 'AOX'
CIV_abs_Sstip$component <- 'Complex IV'
CIV_abs_Sstip <- CIV_abs_Sstip[,c(1,2,3,4,5,6,7,10)]
complexI_abs_Sstip$component <- 'Complex I'
#NADH_dehydrogenase_abs_Sstip$component <- 'NDE1'
#NADH_dehydrogenase_abs_Sstip <- NADH_dehydrogenase_abs_Sstip[,c(1,2,3,4,5,6,7,10)]

toPlot <- rbind(NADH_dehydrogenase_abs_Sstip,complexI_abs_Sstip,AOXdata,CIV_abs_Sstip)
# Calculate mean abundance and standard deviation
toPlot$meanAbundance <- rowMeans(toPlot[,5:7])
toPlot$sd <- apply(toPlot[,5:7], 1, sd)
toPlot$component <- factor(toPlot$component,levels=c('NDE1','Complex I','AOX','Complex IV'))
# Prepare data for individual process jitter plots
df1 <- melt(toPlot[grepl('NDE1',toPlot$component),c(3,5:7)],id.vars = 'PrimaryGeneName')
df2 <- melt(toPlot[which(toPlot$component=='Complex I'),4:7],id.vars = c('GeneName'))
df3 <- melt(toPlot[grepl('AOX',toPlot$component),c(3,5:7)],id.vars = c('PrimaryGeneName'))
df4 <- melt(toPlot[which(toPlot$component=='Complex IV'),4:7],id.vars = c('GeneName'))

# Create bar plot for individual components
figS10A <- ggplot(toPlot[grepl('NDE1',toPlot$component),],
                    aes(x=PrimaryGeneName,y=meanAbundance*1000,fill='#DEEBF7')) +
  geom_bar(position=position_dodge2(), stat="identity",color='black',alpha=0.85) +
  geom_errorbar(aes(ymin=(meanAbundance-sd)*1000,ymax=(meanAbundance+sd)*1000),width=0.2,
                position=position_dodge(0.9)) +
  geom_jitter(data=df1,aes(x=PrimaryGeneName,y=value*1000,color='#878787'),
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),size=1) +
  scale_y_continuous(breaks = seq(0,0.6,0.1),limits = c(0,0.65),expand = c(0,0)) +
  scale_fill_manual(values=c('#DEEBF7')) +
  scale_color_manual(values=c('#878787')) +
  theme_classic() +
  theme(text = element_text(size = 8, color = 'black',face = 'bold'),
        legend.position = 'None',
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10,face = 'bold',color = 'black'),
        axis.text.x = element_text(angle = 0,size = 7,face = 'bold',color = 'black'),
        axis.text.y = element_text(size = 8,face = 'bold',color = 'black'),
        axis.ticks = element_line(color = 'black'),
        plot.title = element_text(hjust=0.5)) +
  ggtitle('External NADH dehydrogenase') +
  ylab('Abundance (mg/gDW)')
figS10A

# Save figure
#ggsave('../../results/figS10A.pdf',
 #      figS10A,width=3.5,height = 6,units = 'cm',device = 'pdf',dpi=300)

figS10B <- ggplot(data=toPlot[which(toPlot$component=='Complex I'),],aes(x=GeneName,y=meanAbundance*1000,fill='#9ECAE1')) +
  geom_bar(position=position_dodge2(), stat="identity",color='black',alpha=0.85) +
  geom_errorbar(aes(x=GeneName,ymin=(meanAbundance-sd)*1000,ymax=(meanAbundance+sd)*1000),width=0.4,
                position=position_dodge(0.9)) +
  geom_jitter(data=df2,aes(x=GeneName,y=value*1000,color='#878787'),
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),size=0.5) +
  scale_y_continuous(breaks = seq(0,0.35,0.05),limits = c(0,0.375),expand = c(0,0)) +
  scale_fill_manual(values=c('#9ECAE1')) +
  scale_color_manual(values=c('#878787')) +
  theme_classic() +
  theme(text = element_text(size = 8, color = 'black',face = 'bold'),
        legend.position = 'None',
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10,face = 'bold',color = 'black'),
        axis.text.x = element_text(angle = 45,size = 7,face = 'bold',color = 'black',vjust=0.6),
        axis.text.y = element_text(size = 8,face = 'bold',color = 'black'),
        axis.ticks = element_line(color = 'black'),
        plot.title = element_text(hjust=0.5)) +
  ggtitle('Complex I') +
  ylab('Abundance (mg/gDW)')
figS10B

# Save figure
#ggsave('../../results/figS10B.pdf',
 #      figS10B,width=9.4,height = 7.5,units = 'cm',device = 'pdf',dpi=300)


figS10C <- ggplot(toPlot[grepl('AOX',toPlot$component),],
                   aes(x=PrimaryGeneName,y=meanAbundance*1000,fill='#4292C6')) +
  geom_bar(position=position_dodge2(), stat="identity",color='black',alpha=0.85) +
  geom_errorbar(aes(ymin=(meanAbundance-sd)*1000,ymax=(meanAbundance+sd)*1000),width=0.2,
                position=position_dodge(0.9)) +
  geom_jitter(data=df3,aes(x=PrimaryGeneName,y=value*1000,color='#878787')) +
  scale_y_continuous(breaks = seq(0,0.035,0.005),limits = c(0,0.038),expand = c(0,0)) +
  scale_fill_manual(values=c('#4292C6')) +
  scale_color_manual(values=c('#878787')) +
  theme_classic() +
  theme(text = element_text(size = 8, color = 'black',face = 'bold'),
        legend.position = 'None',
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10,face = 'bold',color = 'black'),
        axis.text.x = element_text(size = 7,face = 'bold',color = 'black'),
        axis.text.y = element_text(size = 8,face = 'bold',color = 'black'),
        axis.ticks = element_line(color = 'black'),
        plot.title = element_text(hjust = 0.5)) +
  ggtitle('Alternative oxidase (STO1)') +
  ylab('Abundance (mg/gDW)')
figS10C

# Save figure
#ggsave('../../results/figS10C.pdf',
 #      figS10C,width=3.3,height = 7,units = 'cm',device = 'pdf',dpi=300)

figS10D <- ggplot(data=toPlot[which(toPlot$component == 'Complex IV'),],
                   aes(x=GeneName,y=meanAbundance*1000,fill='#08519C')) +
  geom_bar(position=position_dodge2(), stat="identity",color='black',alpha=0.85) +
  geom_errorbar(aes(x=GeneName,ymin=(meanAbundance-sd)*1000,ymax=(meanAbundance+sd)*1000),width=0.4,
                position=position_dodge(0.9)) +
  geom_jitter(data=df4,aes(x=GeneName,y=value*1000,color='#878787')) +
  scale_y_continuous(breaks = seq(0,0.35,0.05),limits = c(0,0.375),expand = c(0,0)) +
  scale_fill_manual(values=c('#08519C')) +
  scale_color_manual(values=c('#878787')) +
  theme_classic() +
  theme(text = element_text(size = 8, color = 'black',face = 'bold'),
        legend.position = 'None',
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10,face = 'bold',color = 'black'),
        axis.text.x = element_text(angle = 90,size = 7,face = 'bold',color = 'black'),
        axis.text.y = element_text(size = 8,face = 'bold',color = 'black'),
        axis.ticks = element_line(color = 'black'),
        plot.title = element_text(hjust=0.5)) +
  ggtitle('Complex IV') +
  ylab('Abundance (mg/gDW)')
figS10D

# Save figure
ggsave('../../results/figS10D.pdf',
       figS10D,width=9.4,height = 6.2,units = 'cm',device = 'pdf',dpi=300)

# fig S11C
yeast <- rep(c("Scer", "Kmarx", "Spombe", "Pstip"), each=3) #all lists are with this ordered 
crabtree <- rep(c("pos", "neg", "pos", "neg"), each=3)
gProt.gDW <- c(0.526511468, 0.446126383, 0.476327579, #Scer
               0.364441271, 0.386171605, 0.367812889, #Kmarx
               0.323023467, 0.370122217, 0.296523312, #Spombe
               0.320718791, 0.438759556, 0.32447428) #Pstip
growth <- c(0.4131,0.43657,0.4182, #Scer
            0.44954, 0.42759, 0.44624, #Kmarx
            0.22454, 0.21923, 0.22382, #Spombe
            0.46499, 0.4661, 0.46809) #Pstip

Scer <- read_excel("Scer.xlsx")
Kmarx <- read_excel("Kmarx.xlsx")
Spombe <- read_excel("Spombe.xlsx")
Pstip <- read_excel("Pstip.xlsx")

#pull out RPs from Uniprot:
uniprot_Scer <- read_excel("../databases/20200731 uniprot_scer.xlsx")
uniprot_Kmarx <- read_excel("../databases/20200731 uniprot_kmarx.xlsx")
uniprot_Spombe <- read_excel("../databases/20200731 uniprot_spombe.xlsx")
uniprot_Pstip <- read_excel("../databases/20200731 uniprot_pstip.xlsx")


RP_Scer <- uniprot_Scer[grepl("ribosom", uniprot_Scer$`Protein names`) | grepl("Ribosom", uniprot_Scer$`Protein names`),c(1,3)]
RP_Scer <- RP_Scer[!grepl("mitochondr", RP_Scer$`Protein names`) & !grepl("Mitochondr", RP_Scer$`Protein names`) & 
                     !grepl("biogenesis", RP_Scer$`Protein names`) & !grepl("export", RP_Scer$`Protein names`) & 
                     !grepl("RNA", RP_Scer$`Protein names`),]
RP_Scer <- RP_Scer[order(RP_Scer$`Protein names`),]
RP_Scer <- RP_Scer[-c(134:135,137:160,164:165),]

RP_Spombe <- uniprot_Spombe[grepl("ribosom", uniprot_Spombe$`Protein names`)| grepl("Ribosom", uniprot_Spombe$`Protein names`),c(1,3)]
RP_Spombe <- RP_Spombe[!grepl("mitochondr", RP_Spombe$`Protein names`) &  !grepl("Mitochondr", RP_Spombe$`Protein names`) & 
                         !grepl("biogenesis", RP_Spombe$`Protein names`) & 
                         !grepl("export", RP_Spombe$`Protein names`) & 
                         !grepl("RNA", RP_Spombe$`Protein names`),]
RP_Spombe <- RP_Spombe[order(RP_Spombe$`Protein names`),]
RP_Spombe <- RP_Spombe[-c(139,141:152,157),]

RP_Kmarx <- uniprot_Kmarx[grepl("ribosom", uniprot_Kmarx$`Protein names`) | grepl("Ribosom", uniprot_Kmarx$`Protein names`),c(1,3)]
RP_Kmarx <- RP_Kmarx[!grepl("mitochondr", RP_Kmarx$`Protein names`) & 
                       !grepl("Mitochondr", RP_Kmarx$`Protein names`) & 
                       !grepl("biogenesis", RP_Kmarx$`Protein names`) & 
                       !grepl("export", RP_Kmarx$`Protein names`) & 
                       !grepl("RNA", RP_Kmarx$`Protein names`) &
                       !grepl("37S", RP_Kmarx$`Protein names`) & 
                       !grepl("54S", RP_Kmarx$`Protein names`) & 
                       !grepl("30S", RP_Kmarx$`Protein names`),]
RP_Kmarx <- RP_Kmarx[order(RP_Kmarx$`Protein names`),]
RP_Kmarx <- RP_Kmarx[-c(66:70, 82:91),]


RP_Pstip <- uniprot_Pstip[grepl("ribosom", uniprot_Pstip$`Protein names`) |
                            grepl("Ribosom", uniprot_Pstip$`Protein names`) ,c(1,3)]
RP_Pstip <- RP_Pstip[!grepl("mitochondr", RP_Pstip$`Protein names`) & 
                       !grepl("Mitochondr", RP_Pstip$`Protein names`) & 
                       !grepl("biogenesis", RP_Pstip$`Protein names`) & 
                       !grepl("export", RP_Pstip$`Protein names`) & 
                       !grepl("RNA", RP_Pstip$`Protein names`) &
                       !grepl("37S", RP_Pstip$`Protein names`) & 
                       !grepl("54S", RP_Pstip$`Protein names`),]
RP_Pstip <- RP_Pstip[order(RP_Pstip$`Protein names`),]
RP_Pstip <- RP_Pstip[-c(49:50,52, 78,79,89:92),]






# pull out RP abundances:


RP_Scer <- merge(Scer, RP_Scer, by.x="Accession", by.y="Entry", all.y=T)
RP_Scer <- RP_Scer[,c(1,15,3:14)]
RP_Scer <- RP_Scer[order(RP_Scer$`Protein names`),]
rownames(RP_Scer) <- 1:nrow(RP_Scer)
nrow(na.omit(RP_Scer))

RP_Spombe <- merge(Spombe, RP_Spombe, by.x="Accession", by.y="Entry", all.y=T)
RP_Spombe <- RP_Spombe[,c(1,15,3:14)]
RP_Spombe <- RP_Spombe[order(RP_Spombe$`Protein names`),]
rownames(RP_Spombe) <- 1:nrow(RP_Spombe)
nrow(na.omit(RP_Spombe))

RP_Kmarx <- merge(Kmarx, RP_Kmarx, by.x="Accession", by.y="Entry", all.y=T)
RP_Kmarx <- RP_Kmarx[,c(1,15,3:14)]
RP_Kmarx <- RP_Kmarx[order(RP_Kmarx$`Protein names`),]
rownames(RP_Kmarx) <- 1:nrow(RP_Kmarx)
nrow(na.omit(RP_Kmarx))

RP_Pstip <- merge(Pstip, RP_Pstip, by.x="Accession", by.y="Entry", all.y=T)
RP_Pstip <- RP_Pstip[,c(1,15,3:14)]
RP_Pstip <- RP_Pstip[order(RP_Pstip$`Protein names`),]
rownames(RP_Pstip) <- 1:nrow(RP_Pstip)
nrow(na.omit(RP_Pstip))


gRibo.gProt <- c(colSums(RP_Scer[,9:11], na.rm=T), colSums(RP_Kmarx[,9:11], na.rm=T), 
                 colSums(RP_Spombe[,9:11], na.rm=T), colSums(RP_Pstip[,9:11], na.rm=T))



RP_Scer[,15:17] <- NA
colnames(RP_Scer)[15:17] <- paste("nmol/mgDW", 1:3, sep="_")
RP_Scer[,15] <- RP_Scer[,6] * gProt.gDW[1]/1000
RP_Scer[,16] <- RP_Scer[,7] * gProt.gDW[2]/1000
RP_Scer[,17] <- RP_Scer[,8] * gProt.gDW[3]/1000


RP_Kmarx[,15:17] <- NA
colnames(RP_Kmarx)[15:17] <- paste("nmol/mgDW", 1:3, sep="_")
RP_Kmarx[,15] <- RP_Kmarx[,6] * gProt.gDW[4]/1000
RP_Kmarx[,16] <- RP_Kmarx[,7] * gProt.gDW[5]/1000
RP_Kmarx[,17] <- RP_Kmarx[,8] * gProt.gDW[6]/1000

RP_Spombe[,15:17] <- NA
colnames(RP_Spombe)[15:17] <- paste("nmol/mgDW", 1:3, sep="_")
RP_Spombe[,15] <- RP_Spombe[,6] * gProt.gDW[7]/1000
RP_Spombe[,16] <- RP_Spombe[,7] * gProt.gDW[8]/1000
RP_Spombe[,17] <- RP_Spombe[,8] * gProt.gDW[9]/1000

RP_Pstip[,15:17] <- NA
colnames(RP_Pstip)[15:17] <- paste("nmol/mgDW", 1:3, sep="_")
RP_Pstip[,15] <- RP_Pstip[,6] * gProt.gDW[10]/1000
RP_Pstip[,16] <- RP_Pstip[,7] * gProt.gDW[11]/1000
RP_Pstip[,17] <- RP_Pstip[,8] * gProt.gDW[12]/1000

nmolRibo.mgDW <- c(colMeans(RP_Scer[,15:17], na.rm=T), colMeans(RP_Kmarx[,15:17], na.rm=T), 
                   colMeans(RP_Spombe[,15:17], na.rm=T), colMeans(RP_Pstip[,15:17], na.rm=T))

# now can plot Supplemental C

# suppl figure C
# Prepare data for plotting
dataToPlot <- cbind(data.frame(growth),data.frame(gRibo.gProt))
colnames(dataToPlot)[2] <- 'gRibo_gProt'
dataToPlot$organism <- c(rep('C',3),rep('M',3),rep('P',3),rep('S',3))
dataToPlot$group <- c(rep('Crab-pos',3),rep('Crab-neg',3),rep('Crab-pos',3),rep('Crab-neg',3))

supp_fig4C <- ggplot(data = dataToPlot,aes(x=growth,y=gRibo_gProt,color=group)) +
  geom_text(aes(label=organism),size=5,) +
  theme_classic() +
  theme(axis.text.x = element_text(size=8,color='black',face='bold'),
        axis.title.x = element_text(size=8,color='black',face='bold'),
        axis.text.y = element_text(size=8,color='black',face='bold'),
        axis.title.y = element_text(size=8,color='black',face='bold'),
        legend.position = 'None') + 
  ylab('Ribosomes (g/g protein)') +
  xlab('Growth rate (h-1)') +
  scale_y_continuous(limits = c(0,0.3),breaks=seq(0,0.3,0.05),expand = c(0,0)) +
  scale_x_continuous(limits=c(0,0.55),breaks=seq(0,0.5,0.1),expand = c(0,0)) +
  scale_color_manual(values=c('blue','orange')) 
supp_fig4C

# Construct tables
data_supp_RY_A <- data.frame(matrix(nrow=0,ncol=5))
data_supp_RY_A <- rbind(data.frame(t(c('# RP (Uniprot)',137,143,78,86))),
                        data.frame(t(c('# RP detected',101,110,76,84))),
                        data.frame(t(c('% RP detected',74,77,97,98))))
colnames(data_supp_RY_A) <- c('','Sce','Spo','Kma','Sstip')
data_supp_RY_A.t <- ggtexttable(data_supp_RY_A,rows = NULL)
data_supp_RY_A.t

data_supp_RY_B <- data.frame(matrix(nrow=0,ncol=5))
data_supp_RY_B <- rbind(data.frame(t(c('# RP subunits (paralogs binned)',78,81,78,86))),
                        data.frame(t(c('# RP subunits detected',77,78,76,84))),
                        data.frame(t(c('% RP subunits detected',99,96,97,98))))
colnames(data_supp_RY_B) <- c('','Sce','Spo','Kma','Sstip')
data_supp_RY_B.t <- ggtexttable(data_supp_RY_B,rows = NULL)
data_supp_RY_B.t

data_supp_RY_D <- data.frame(matrix(nrow=0,ncol=5))
data_supp_RY_D <- rbind(data.frame(t(c('# MRPs (Uniprot)',77,71,66,73))),
                        data.frame(t(c('# MRPs detected',75,64,65,72))),
                        data.frame(t(c('% MRPs detected',97,90,98,99))))
colnames(data_supp_RY_D) <- c('','Sce','Spo','Kma','Sstip')
data_supp_RY_D.t <- ggtexttable(data_supp_RY_D,rows = NULL)
data_supp_RY_D.t

# Construct and save figure
figS11 <- ggarrange(data_supp_RY_A.t,data_supp_RY_B.t,supp_fig4C,data_supp_RY_D.t,
                       nrow = 2,ncol = 2)
figS11
#ggsave('../../results/figS11.pdf',
 #      figS11,width=16,height = 16,units = 'cm',device = 'pdf',dpi=300)

