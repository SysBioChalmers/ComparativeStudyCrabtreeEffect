### Figure 2

# load required packages
library(tidyverse)
library(RColorBrewer)
library(dplyr)
library(reshape2)

# set working directory
setwd('../../data/experimentalData')

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

# Define function for extracting abundance data given a list of proteins
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

# Define function for mapping mitochondrial proteins based on child terms to GO:0005739 mitochondrion
mapMitoProteins <- function(proteomicsData,GOannotation){
  mitoProteins <- GOannotation[grep('mitochondri',GOannotation$GOname),]
  # Remove terms related to cellular component, keeping only biological process terms
  #mitoProteins <- mitoProteins[-grep('GO:0005739',mitoProteins$GOterm),]
  mitoProteins <- mitoProteins[!duplicated(mitoProteins$UniprotID),]
  mitoProteins <- data.frame(mitoProteins[,1:3])
  colnames(mitoProteins) <- c('Accession','GOterm','GOtermName')
  mitochondrialProteins <- merge(proteomicsData,mitoProteins)
  return(mitochondrialProteins)
}

mapUniprotToProtData <- function(UniprotIDs,protData){
  data <- data.frame(matrix(ncol = 7,nrow = 0))
  for(ID in UniprotIDs){
    proteinData <- protData[grep(ID,protData$Accession),]
    if(!nrow(proteinData)==0){
      data <- rbind(data,proteinData)
    }
  }
  return(data)
}

# fig 2B
# Read data on PPP genes
PPP_proteins <- read.delim('PPP_all.txt')
PPP_proteins <- PPP_proteins[,c(1,3,8)]

# S. cerevisiae
glycolysis_Sce <- map2GO(protData_Sce,GOannotation_Sce,c('GO:0006096'))
# remove pyruvate dehydrogenase genes (PDA1 and PDB1) from glycolysis GO term
glycolysis_Sce <- glycolysis_Sce[-grep('PD',glycolysis_Sce$PrimaryGeneName),]
mitochondrial_Sce <- mapMitoProteins(protData_Sce,GOannotation_Sce)
mitochondrial_Sce <- mitochondrial_Sce[!mitochondrial_Sce$Accession %in% glycolysis_Sce$Accession[grep('glycolytic process',glycolysis_Sce$GOtermName)],]
# load child terms of amino acid biosynthesis
AAchildTerms <- read.delim('../GOterms/childTerms_AAbiosynthesis.tsv',
                           stringsAsFactors = FALSE,header = FALSE)
AA_cerevisiae <- map2GO(protData_Sce,GOannotation_Sce,AAchildTerms$V1)
PPP_proteins_Sce <- PPP_proteins[grep('Sce',PPP_proteins$organism),]
PPP_data_Sce <- mapUniprotToProtData(PPP_proteins_Sce$Accession, protData_Sce)
combined_Sce <- rbind(data.frame(t(colSums(glycolysis_Sce[,5:7]))),
                      data.frame(t(colSums(PPP_data_Sce[,5:7]))),
                      data.frame(t(colSums(AA_cerevisiae[,5:7]))),
                      data.frame(t(colSums(mitochondrial_Sce[,5:7]))))
combined_Sce$Process <- c('Glycolysis','PPP','Amino acid biosynthesis','Mitochondrion')
combined_Sce$Organism <- c(rep('Sce',4))
# Add process representing other processes
row.names(combined_Sce) <- NULL
combined_Sce$meanAllocation <- rowMeans(combined_Sce[,1:3])
combined_Sce$sd <- apply(combined_Sce[,1:3], 1, sd)
#combined_Sce <- combined_Sce[,4:7]

# S. pombe
glycolysis_Spo <- map2GO(protData_Spo,GOannotation_Spo,c('GO:0006096'))
glycolysis_Spo <- glycolysis_Spo[-grep('pda1',glycolysis_Spo$PrimaryGeneName),]
glycolysis_Spo <- glycolysis_Spo[-grep('pdb1',glycolysis_Spo$PrimaryGeneName),]
mitochondrial_Spo <- mapMitoProteins(protData_Spo,GOannotation_Spo)
mitochondrial_Spo <- mitochondrial_Spo[!mitochondrial_Spo$Accession %in% glycolysis_Spo$Accession[grep('glycolytic process',glycolysis_Spo$GOtermName)],]
AA_pombe <- map2GO(protData_Spo,GOannotation_Spo,AAchildTerms$V1)
PPP_proteins_Spo <- PPP_proteins[grep('Spo',PPP_proteins$organism),]
PPP_data_Spo <- mapUniprotToProtData(PPP_proteins_Spo$Accession, protData_Spo)
combined_Spo <- rbind(data.frame(t(colSums(glycolysis_Spo[,5:7]))),
                      data.frame(t(colSums(PPP_data_Spo[,5:7]))),
                      data.frame(t(colSums(AA_pombe[,5:7]))),
                      data.frame(t(colSums(mitochondrial_Spo[,5:7]))))
combined_Spo$Process <- c('Glycolysis','PPP','Amino acid biosynthesis','Mitochondrion')
combined_Spo$Organism <- c(rep('Spo',4))
# Add process representing other processes
combined_Spo$meanAllocation <- rowMeans(combined_Spo[,1:3])
combined_Spo$sd <- apply(combined_Spo[,1:3], 1, sd)
#combined_Spo <- combined_Spo[,4:7]

# S.stipitis
glycolysis_Sstip <- map2GO(protData_Sstip,GOannotation_Sstip,c('GO:0006096'))
glycolysis_Sstip <- glycolysis_Sstip[-grep('PD',glycolysis_Sstip$PrimaryGeneName),]
mitochondrial_Sstip <- mapMitoProteins(protData_Sstip,GOannotation_Sstip)
mitochondrial_Sstip <- mitochondrial_Sstip[!mitochondrial_Sstip$Accession %in% glycolysis_Sstip$Accession[grep('glycolytic process',glycolysis_Sstip$GOtermName)],]
AA_stipitis <- map2GO(protData_Sstip,GOannotation_Sstip,AAchildTerms$V1)
PPP_proteins_Sstip <- PPP_proteins[grep('Sstip',PPP_proteins$organism),]
PPP_data_Sstip <- mapUniprotToProtData(PPP_proteins_Sstip$Accession, protData_Sstip)
combined_Sstip <- rbind(data.frame(t(colSums(glycolysis_Sstip[,5:7]))),
                        data.frame(t(colSums(PPP_data_Sstip[,5:7]))),
                        data.frame(t(colSums(AA_stipitis[,5:7]))),
                        data.frame(t(colSums(mitochondrial_Sstip[,5:7]))))
combined_Sstip$Process <- c('Glycolysis','PPP','Amino acid biosynthesis','Mitochondrion')
combined_Sstip$Organism <- c(rep('Sstip',4))
# Add process representing other processes
combined_Sstip$meanAllocation <- rowMeans(combined_Sstip[,1:3])
combined_Sstip$sd <- apply(combined_Sstip[,1:3], 1, sd)
#combined_Sstip <- combined_Sstip[,4:7]

# K. marxianus
glycolysis_Kma <- map2GO(protData_Kma,GOannotation_Kma,c('GO:0006096'))
glycolysis_Kma <- glycolysis_Kma[-grep('PDA1',glycolysis_Kma$PrimaryGeneName),]
glycolysis_Kma <- glycolysis_Kma[-grep('PDB1',glycolysis_Kma$PrimaryGeneName),]
mitochondrial_Kma <- mapMitoProteins(protData_Kma,GOannotation_Kma)
mitochondrial_Kma <- mitochondrial_Kma[!mitochondrial_Kma$Accession %in% glycolysis_Kma$Accession[grep('glycolytic process',glycolysis_Kma$GOtermName)],]
AA_marxianus <- map2GO(protData_Kma,GOannotation_Kma,AAchildTerms$V1)
PPP_proteins_Kma <- PPP_proteins[grep('Kma',PPP_proteins$organism),]
PPP_data_Kma <- mapUniprotToProtData(PPP_proteins_Kma$Accession, protData_Kma)
#combined_Kma <- rbind(glycolysis_ribosomes_Kma,mitochondrial_Kma,AA_marxianus)
combined_Kma <- rbind(data.frame(t(colSums(glycolysis_Kma[,5:7]))),
                      data.frame(t(colSums(PPP_data_Kma[,5:7]))),
                      data.frame(t(colSums(AA_marxianus[,5:7]))),
                      data.frame(t(colSums(mitochondrial_Kma[,5:7]))))
combined_Kma$Process <- c('Glycolysis','PPP','Amino acid biosynthesis','Mitochondrion')
combined_Kma$Organism <- c(rep('Kma',4))
# Add process representing other processes
combined_Kma$meanAllocation <- rowMeans(combined_Kma[,1:3])
combined_Kma$sd <- apply(combined_Kma[,1:3], 1, sd)
#combined_Kma <- combined_Kma[,4:7]
# summarize the processes and plot results
combined_All <- rbind(combined_Sce,combined_Spo,combined_Sstip,combined_Kma)
combined_All$Organism <- factor(combined_All$Organism,levels = c('Spo','Sce','Kma','Sstip'))
combined_All$rGlu <- c(rep(18.8,4),rep(8.9,4),rep(4.5,4),rep(5.6,4))
colnames(combined_All)[1:3] <- c('repl1','repl2','repl3')
AA <- combined_All[grepl('Amino acid biosynthesis',combined_All$Process),]
df1 <- melt(AA[,c(1:4,8)],id.vars = c('Process','rGlu'),measure.vars = c('repl1','repl2','repl3'))
glycolysis <- combined_All[grepl('Glycolysis',combined_All$Process),]
df2 <- melt(glycolysis[,c(1:4,8)],id.vars = c('Process','rGlu'),measure.vars = c('repl1','repl2','repl3'))
mito <- combined_All[grepl('Mitochondrion',combined_All$Process),]
df3 <- melt(mito[,c(1:4,8)],id.vars = c('Process','rGlu'),measure.vars = c('repl1','repl2','repl3'))
PPP <- combined_All[grepl('PPP',combined_All$Process),]
df4 <- melt(PPP[,c(1:4,8)],id.vars = c('Process','rGlu'),measure.vars = c('repl1','repl2','repl3'))

# Plot figure
fig2B <- ggplot() +
  geom_point(data=AA,aes(x=rGlu,y=meanAllocation*100),color='#9970AB',shape=22,fill='#9970AB') +
  geom_errorbar(data=AA,aes(x=rGlu,ymin=(meanAllocation-sd)*100,ymax=(meanAllocation+sd)*100),
                width=.2,position=position_dodge(.9)) +
  geom_jitter(data=df1,aes(x=rGlu,y=value*100),
              position = position_jitter(0.2),color='#9970AB',size=2) +
  geom_point(data=glycolysis,aes(x=rGlu,y=meanAllocation*100),color='#5AAE61',shape=22,fill='#5AAE61') +
  geom_errorbar(data=glycolysis,aes(x=rGlu,ymin=(meanAllocation-sd)*100,ymax=(meanAllocation+sd)*100),
                width=.2,position=position_dodge(.9)) +
  geom_jitter(data=df2,aes(x=rGlu,y=value*100),
              position = position_jitter(0.2),color='#5AAE61',size=2) +
  geom_point(data=mito,aes(x=rGlu,y=meanAllocation*100),color='#6BAED6',shape=22,fill='#6BAED6') +
  geom_errorbar(data=mito,aes(x=rGlu,ymin=(meanAllocation-sd)*100,ymax=(meanAllocation+sd)*100),
                width=.2,position=position_dodge(.9)) +
  geom_jitter(data=df3,aes(x=rGlu,y=value*100),
              position = position_jitter(0.2),color='#6BAED6',size=2) +
  geom_point(data=PPP,aes(x=rGlu,y=meanAllocation*100),color='#A6761D',shape=22,fill='#A6761D') +
  geom_errorbar(data=PPP,aes(x=rGlu,ymin=(meanAllocation-sd)*100,ymax=(meanAllocation+sd)*100),
                width=.2,position=position_dodge(.9)) +
  geom_jitter(data=df4,aes(x=rGlu,y=value*100),
              position = position_jitter(0.2),color='#A6761D',size=2) +
  #geom_point(size = 2) +
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
  scale_y_continuous(limits = c(0,30),breaks = seq(0,30,5),expand = c(0,0)) +
  scale_x_continuous(breaks=c(4.5,5.6,8.9,18.8))
fig2B
# Save plot
#ggsave('../../results/fig2B.pdf',
 #      fig2B,width=6,height = 6,units = 'cm',device = 'pdf',dpi=300)

# fig 2C
protData_fmol_per_ug_Sce <- read.delim('proteomicsData_fmol_per_ug_protein_Scer.txt',
                                       sep = '\t',stringsAsFactors = FALSE)
protData_fmol_per_ug_Spo <- read.delim('proteomicsData_fmol_per_ug_protein_Spombe.txt',
                                       sep = '\t',stringsAsFactors = FALSE)
protData_fmol_per_ug_Kma <- read.delim('/proteomicsData_fmol_per_ug_protein_Kmarx.txt',
                                       sep = '\t',stringsAsFactors = FALSE)
protData_fmol_per_ug_Sstip <- read.delim('proteomicsData_fmol_per_ug_protein_Pstip.txt',
                                         sep = '\t',stringsAsFactors = FALSE)

GlcTransporters_Spo <- c('ght1','ght2','ght3','ght4','ght5','ght6','ght7','ght8')
GlcTransporters_data_Spo <- mapGenesToProtData(GlcTransporters_Spo,protData_fmol_per_ug_Spo)
GlcTransporters_data_Spo$organism <- c(rep('Spo',4))
GlcTransporters_Sce <- c('HXT1','HXT2','HXT3','HXT4','HXT5','HXT6','HXT7','HXT8','HXT9',
                         'HXT10','HXT11','HXT13','HXT14','HXT15','HXT16','HXT17')
GlcTransporters_data_Sce <- mapGenesToProtData(GlcTransporters_Sce,protData_fmol_per_ug_Sce)
GlcTransporters_data_Sce$organism <- c(rep('Sce',7))
GlcTransporters_Kma <- c('RAG1','HGT1','KHT2')
GlcTransporters_data_Kma <- mapGenesToProtData(GlcTransporters_Kma,protData_fmol_per_ug_Kma)
# Two entries identified for KHT2, take average of those
GlcTransporters_data_Kma[2,5:7] <- colMeans(GlcTransporters_data_Kma[2:3,5:7])
GlcTransporters_data_Kma <- GlcTransporters_data_Kma[1:2,]
GlcTransporters_data_Kma$organism <- c(rep('Kma',2))
GlcTransporters_Sstip <- c('HGT1','HGT2','SUT1','SUT2','SUT3')
GlcTransporters_data_Sstip <- mapGenesToProtData(GlcTransporters_Sstip,protData_fmol_per_ug_Sstip)
GlcTransporters_data_Sstip$organism <- c(rep('Sstip',2))
GlcTransporters_All <- rbind(GlcTransporters_data_Spo,GlcTransporters_data_Sce,
                             GlcTransporters_data_Kma,GlcTransporters_data_Sstip)
GlcTransporters_All$organism <- factor(GlcTransporters_All$organism,levels = c('Spo','Sce','Kma','Sstip'))
GlcTransporters_All$PrimaryGeneName <- factor(GlcTransporters_All$PrimaryGeneName,levels = c('ght2','ght4','ght5','ght8',
                                                                                             'HXT1','HXT2',
                                                                                             'HXT3','HXT4',
                                                                                             'HXT5','HXT6',
                                                                                             'HXT7','KHT2',
                                                                                             'RAG1','HGT2',
                                                                                             'SUT1'))

# Calculate mean and standard deviation
GlcTransporters_All$meanAbundance = rowMeans(GlcTransporters_All[,5:7])
GlcTransporters_All$sd <- apply(GlcTransporters_All[,5:7],1,sd)
df1 <- GlcTransporters_All[,c(3,5:8)]
colnames(df1)[2:4] <- c('repl1','repl2','repl3')
df1 <- melt(df1,id.vars = c('PrimaryGeneName','organism'),measure.vars = c('repl1','repl2','repl3'))

# Plot abundances
fig2C <- ggplot() +
  geom_bar(data=GlcTransporters_All,aes(x=PrimaryGeneName,y=meanAbundance,fill=organism),stat="identity", color="black", 
           position=position_dodge(),alpha=0.8) +
  geom_errorbar(data=GlcTransporters_All,aes(x=PrimaryGeneName,ymin=meanAbundance-sd, ymax=meanAbundance+sd), width=.2,
                position=position_dodge(.9)) +
  geom_jitter(data=df1,aes(x=PrimaryGeneName,y=value,
                           color=organism),position = position_jitter(0.2)) +
  scale_fill_manual(values=c('#CBBBA0','#000000','#1D71B8','#878787')) +
  scale_color_manual(values=c('#CBBBA0','#000000','#1D71B8','#878787')) +
  theme_classic() +
  theme(plot.title = element_text(size = 10,face = 'bold',color = 'black',hjust = 0.5),
        axis.text.x = element_text(angle = 90),
        axis.text = element_text(size = 7,face = 'bold',color='black'),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 8,face = 'bold',color = 'black'),
        panel.border = element_rect(size = 1.1, color = 'black',fill = NA),
        legend.title = element_blank(),
        legend.text = element_text(size = 8,face = 'bold',color = 'black'),
        legend.position = c(0.8,0.8),
        axis.ticks = element_line(color = 'black')) +
  scale_y_continuous(expand = c(0,0),limits = c(0,30),breaks = seq(0,30,5)) +
  #ggtitle('Glucose transporters') +
  ylab('Abundance (fmol/ug protein)')
fig2C

# Save figure
# Save plot
#ggsave('../../results/fig2C.pdf',
 #      fig2C,width=6,height = 6,units = 'cm',device = 'pdf',dpi=300)

# fig 2D
# Load data on ATP production
ATPdata <- read.delim('../../results/modelSimulation/ATPproduction.txt',stringsAsFactors = F)

# Plot net ATP production
fluxData <- ATPdata[,c(1,4,5)]
fluxData <- melt(fluxData,id='Organism')
fluxData$Organism <- factor(fluxData$Organism,levels=c('Sce','Kma'))

fig2D <- ggplot(data = fluxData[grep('Net.ATP.production.mmol.gDW.h',fluxData$variable),],
                      aes(x=reorder(Organism,value),y=value,fill=Organism)) +
  geom_bar(position = 'dodge',stat = 'identity',color = 'black') +
  theme_classic() +
  theme(axis.title.y = element_text(size = 8,face = 'bold',color = 'black'),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 8,face = 'bold',color = 'black'),
        axis.text.x = element_text(size = 8,face = 'bold',color = 'black'),
        legend.title = element_blank(),
        legend.text = element_text(size = 8,face = 'bold',color = 'black'),
        legend.position = 'None') +
  scale_fill_manual(values = c('#1D1D1B','#1D71B8')) +
  ylab('ATP production (mmol ATP/gDW/h)') +
  scale_y_continuous(expand = c(0,0),breaks = seq(0,35,5),limits = c(0,37)) #+
#scale_x_discrete(labels=c('','HY'))
fig2D
# Save figure
#ggsave('../../results/fig2D.pdf',
 #      fig2D, device = 'pdf',height = 6, width = 5, units = 'cm',dpi = 300)

# fig 2E
yieldData <- ATPdata[,c(1,5)]
yieldData <- melt(yieldData,id=c('Organism'))
yieldData$Organism <- factor(yieldData$Organism,levels = c('Kma','Sce'))

fig2E <- ggplot(data = yieldData,aes(x=Organism,y=value,fill=Organism)) +
  geom_bar(position = 'dodge',stat = 'identity',color = 'black') +
  theme_classic() +
  theme(axis.title.y = element_text(size = 8,face = 'bold',color = 'black'),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 8,face = 'bold',color = 'black'),
        axis.text.x = element_text(size = 8,face = 'bold',color = 'black'),
        legend.title = element_blank(),
        legend.text = element_text(size = 8,face = 'bold',color = 'black'),
        legend.position = 'None') +
  scale_fill_manual(values = c('#1D71B8','#1D1D1B')) +
  ylab('ATP yield (mmol ATP/mmol glucose)') +
  scale_y_continuous(expand = c(0,0),breaks = seq(0,6,0.5),limits = c(0,6))# +
#scale_x_discrete(labels=c('LY','HY'))
fig2E
# Save figure
#ggsave('../../results/fig2E.pdf',
 #      fig2E, device = 'pdf',height = 6, width = 5, units = 'cm',dpi = 300)





