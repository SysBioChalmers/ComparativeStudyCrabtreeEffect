# Figure 1
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

# Define functions used for analysis
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

# fig 1A-D
# Load HPLC datasets

# Scheffersomyces stipitis
filename <- 'HPLC_stipitis.csv'
data_stipitis <- read.csv(filename,header = TRUE,sep = ';',dec = ',')

data_stipitis <- data_stipitis %>%
  rowwise() %>%
  mutate(meanOD = mean(c(OD_repl1,OD_repl2,OD_repl3))) %>%
  mutate(sdOD = sd(c(OD_repl1,OD_repl2,OD_repl3))) %>%
  mutate(meanGlucose = mean(c(Glucose_repl1,Glucose_repl2,Glucose_repl3))) %>%
  mutate(sdGlucose = sd(c(Glucose_repl1,Glucose_repl2,Glucose_repl3))) %>%
  mutate(meanGlycerol = mean(c(Glycerol_repl1,Glycerol_repl2,Glycerol_repl3))) %>%
  mutate(sdGlycerol = sd(c(Glycerol_repl1,Glycerol_repl2,Glycerol_repl3))) %>%
  mutate(meanAcetate = mean(c(Acetate_repl1,Acetate_repl2,Acetate_repl3))) %>%
  mutate(sdAcetate = sd(c(Acetate_repl1,Acetate_repl2,Acetate_repl3))) %>%
  mutate(meanEthanol = mean(c(Ethanol_repl1,Ethanol_repl2,Ethanol_repl3))) %>%
  mutate(sdEthanol = sd(c(Ethanol_repl1,Ethanol_repl2,Ethanol_repl3))) %>%
  mutate(meanPyruvate = mean(c(Pyruvate_repl1,Pyruvate_repl2,Pyruvate_repl3))) %>%
  mutate(sdPyruvate = sd(c(Pyruvate_repl1,Pyruvate_repl2,Pyruvate_repl3))) %>%
  mutate(meanSuccinate = mean(c(Succinate_repl1,Succinate_repl2,Succinate_repl3))) %>%
  mutate(sdSuccinate = sd(c(Succinate_repl1,Succinate_repl2,Succinate_repl3)))

# Plot metabolite concentration over time
param_enlarge <- 20

ggplot(data_stipitis,aes(x = Time)) +
  theme_bw() +
  geom_point(aes(y=meanOD,color = 'OD')) +
  geom_line(aes(y=meanOD,color ='OD')) +
  geom_errorbar(aes(ymin=meanOD-sdOD,ymax=meanOD+sdOD),size=0.5,width=0.3) +
  geom_point(aes(y=meanGlucose,color='Glucose')) +
  geom_line(aes(y=meanGlucose,color='Glucose')) +
  geom_errorbar(aes(ymin=meanGlucose-sdGlucose,ymax=meanGlucose+sdGlucose),size=0.5,width=0.3) +
  geom_point(aes(y=param_enlarge*meanEthanol,color='Ethanol')) +
  geom_line(aes(y=param_enlarge*meanEthanol,color='Ethanol')) +
  geom_errorbar(aes(ymin=param_enlarge*(meanEthanol-sdEthanol),ymax=param_enlarge*(meanEthanol+sdEthanol)),size=0.5,width=0.3) +
  geom_point(aes(y=param_enlarge*meanGlycerol,color='Glycerol')) +
  geom_line(aes(y=param_enlarge*meanGlycerol,color='Glycerol')) +
  geom_errorbar(aes(ymin=param_enlarge*(meanGlycerol-sdGlycerol),ymax=param_enlarge*(meanGlycerol+sdGlycerol)),size=0.5,width=0.3) +
  geom_point(aes(y=param_enlarge*meanAcetate,color='Acetate')) +
  geom_line(aes(y=param_enlarge*meanAcetate,color='Acetate')) +
  geom_errorbar(aes(ymin=param_enlarge*(meanAcetate-sdAcetate),ymax=param_enlarge*(meanAcetate+sdAcetate)),size=0.5,width=0.3) +
  geom_point(aes(y=param_enlarge*meanPyruvate,color='Pyruvate')) +
  geom_line(aes(y=param_enlarge*meanPyruvate,color='Pyruvate')) +
  geom_errorbar(aes(ymin=param_enlarge*(meanPyruvate-sdPyruvate),ymax=param_enlarge*(meanPyruvate+sdPyruvate)),size=0.5,width=0.3) +
  geom_point(aes(y=param_enlarge*meanSuccinate,color='Succinate')) +
  geom_line(aes(y=param_enlarge*meanSuccinate,color='Succinate')) +
  geom_errorbar(aes(ymin=param_enlarge*(meanSuccinate-sdSuccinate),ymax=param_enlarge*(meanSuccinate+sdSuccinate)),size=0.5,width=0.3) +
  scale_y_continuous(name='OD, Glucose (g/l)',
                     sec.axis = sec_axis(~./param_enlarge,breaks = seq(0,2,0.25),name='Acetate, Ethanol, Glycerol, Pyruvate, Succinate (g/l)'),
                     limits = c(-0.4,23),breaks = seq(0,22,2)) +
  scale_x_continuous(limits = c(0,16),breaks = seq(0,16,2)) +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.border = element_rect(color = 'black',size = 0.5,fill = NA),
        axis.text = element_text(size=8,face = 'bold',color = 'black'),
        axis.title = element_text(size=8,face = 'bold',color = 'black'),
        legend.position = 'top',
        legend.title = element_blank(),
        axis.ticks = element_line(color = 'black')) +
  xlab('Time (h)')

# Kluyveromyces marxianus
filename <- 'HPLC_marxianus.csv'
data_marxianus <- read.csv(filename,header = TRUE,sep = ';',dec = ',')

data_marxianus <- data_marxianus %>%
  rowwise() %>%
  mutate(meanOD = mean(c(OD_repl1,OD_repl2,OD_repl3))) %>%
  mutate(sdOD = sd(c(OD_repl1,OD_repl2,OD_repl3))) %>%
  mutate(meanGlucose = mean(c(Glucose_repl1,Glucose_repl2,Glucose_repl3))) %>%
  mutate(sdGlucose = sd(c(Glucose_repl1,Glucose_repl2,Glucose_repl3))) %>%
  mutate(meanGlycerol = mean(c(Glycerol_repl1,Glycerol_repl2,Glycerol_repl3))) %>%
  mutate(sdGlycerol = sd(c(Glycerol_repl1,Glycerol_repl2,Glycerol_repl3))) %>%
  mutate(meanAcetate = mean(c(Acetate_repl1,Acetate_repl2,Acetate_repl3))) %>%
  mutate(sdAcetate = sd(c(Acetate_repl1,Acetate_repl2,Acetate_repl3))) %>%
  mutate(meanEthanol = mean(c(Ethanol_repl1,Ethanol_repl2,Ethanol_repl3))) %>%
  mutate(sdEthanol = sd(c(Ethanol_repl1,Ethanol_repl2,Ethanol_repl3))) %>%
  mutate(meanPyruvate = mean(c(Pyruvate_repl1,Pyruvate_repl2,Pyruvate_repl3))) %>%
  mutate(sdPyruvate = sd(c(Pyruvate_repl1,Pyruvate_repl2,Pyruvate_repl3))) %>%
  mutate(meanSuccinate = mean(c(Succinate_repl1,Succinate_repl2,Succinate_repl3))) %>%
  mutate(sdSuccinate = sd(c(Succinate_repl1,Succinate_repl2,Succinate_repl3)))

# Plot metabolite concentration over time
param_enlarge <- 10

ggplot(data_marxianus,aes(x = Time)) +
  theme_bw() +
  geom_point(aes(y=meanOD,color = 'OD')) +
  geom_line(aes(y=meanOD,color ='OD')) +
  geom_errorbar(aes(ymin=meanOD-sdOD,ymax=meanOD+sdOD),size=0.5,width=0.3) +
  geom_point(aes(y=meanGlucose,color='Glucose')) +
  geom_line(aes(y=meanGlucose,color='Glucose')) +
  geom_errorbar(aes(ymin=meanGlucose-sdGlucose,ymax=meanGlucose+sdGlucose),size=0.5,width=0.3) +
  geom_point(aes(y=meanEthanol,color='Ethanol')) +
  geom_line(aes(y=meanEthanol,color='Ethanol')) +
  geom_errorbar(aes(ymin=meanEthanol-sdEthanol,ymax=meanEthanol+sdEthanol),size=0.5,width=0.3) +
  geom_point(aes(y=param_enlarge*meanGlycerol,color='Glycerol')) +
  geom_line(aes(y=param_enlarge*meanGlycerol,color='Glycerol')) +
  geom_errorbar(aes(ymin=param_enlarge*(meanGlycerol-sdGlycerol),ymax=param_enlarge*(meanGlycerol+sdGlycerol)),size=0.5,width=0.3) +
  geom_point(aes(y=param_enlarge*meanAcetate,color='Acetate')) +
  geom_line(aes(y=param_enlarge*meanAcetate,color='Acetate')) +
  geom_errorbar(aes(ymin=param_enlarge*(meanAcetate-sdAcetate),ymax=param_enlarge*(meanAcetate+sdAcetate)),size=0.5,width=0.3) +
  geom_point(aes(y=param_enlarge*meanPyruvate,color='Pyruvate')) +
  geom_line(aes(y=param_enlarge*meanPyruvate,color='Pyruvate')) +
  geom_errorbar(aes(ymin=param_enlarge*(meanPyruvate-sdPyruvate),ymax=param_enlarge*(meanPyruvate+sdPyruvate)),size=0.5,width=0.3) +
  geom_point(aes(y=param_enlarge*meanSuccinate,color='Succinate')) +
  geom_line(aes(y=param_enlarge*meanSuccinate,color='Succinate')) +
  geom_errorbar(aes(ymin=param_enlarge*(meanSuccinate-sdSuccinate),ymax=param_enlarge*(meanSuccinate+sdSuccinate)),size=0.5,width=0.3) +
  scale_y_continuous(name='OD, Glucose (g/l)',
                     sec.axis = sec_axis(~./param_enlarge,breaks = seq(0,2,0.25),name='Acetate, Ethanol, Glycerol, Pyruvate, Succinate (g/l)'),
                     limits = c(0,12),breaks = seq(0,22,2),expand = c(0,0.1)) +
  scale_x_continuous(limits = c(0,16),breaks = seq(0,16,2),expand = c(0,0)) +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.border = element_rect(color = 'black',size = 0.5,fill = NA),
        axis.text = element_text(size=8,face = 'bold',color = 'black'),
        axis.title = element_text(size=8,face = 'bold',color = 'black'),
        legend.position = 'top',
        legend.title = element_blank(),
        axis.ticks = element_line(color = 'black')) +
  xlab('Time (h)')

# Schizosaccharomyces pombe
filename <- 'HPLC_pombe.csv'
data_pombe <- read.csv(filename,header = TRUE,sep = ';',dec = ',')

data_pombe <- data_pombe %>%
  rowwise() %>%
  mutate(meanOD = mean(c(OD_repl1,OD_repl2,OD_repl3))) %>%
  mutate(sdOD = sd(c(OD_repl1,OD_repl2,OD_repl3))) %>%
  mutate(meanGlucose = mean(c(Glucose_repl1,Glucose_repl2,Glucose_repl3))) %>%
  mutate(sdGlucose = sd(c(Glucose_repl1,Glucose_repl2,Glucose_repl3))) %>%
  mutate(meanGlycerol = mean(c(Glycerol_repl1,Glycerol_repl2,Glycerol_repl3))) %>%
  mutate(sdGlycerol = sd(c(Glycerol_repl1,Glycerol_repl2,Glycerol_repl3))) %>%
  mutate(meanAcetate = mean(c(Acetate_repl1,Acetate_repl2,Acetate_repl3))) %>%
  mutate(sdAcetate = sd(c(Acetate_repl1,Acetate_repl2,Acetate_repl3))) %>%
  mutate(meanEthanol = mean(c(Ethanol_repl1,Ethanol_repl2,Ethanol_repl3))) %>%
  mutate(sdEthanol = sd(c(Ethanol_repl1,Ethanol_repl2,Ethanol_repl3))) %>%
  mutate(meanPyruvate = mean(c(Pyruvate_repl1,Pyruvate_repl2,Pyruvate_repl3))) %>%
  mutate(sdPyruvate = sd(c(Pyruvate_repl1,Pyruvate_repl2,Pyruvate_repl3))) %>%
  mutate(meanSuccinate = mean(c(Succinate_repl1,Succinate_repl2,Succinate_repl3))) %>%
  mutate(sdSuccinate = sd(c(Succinate_repl1,Succinate_repl2,Succinate_repl3)))

# Plot metabolite concentration over time
param_enlarge <- 15

ggplot(data_pombe,aes(x = Time)) +
  theme_bw() +
  geom_point(aes(y=meanOD,color = 'OD')) +
  geom_line(aes(y=meanOD,color ='OD')) +
  geom_errorbar(aes(ymin=meanOD-sdOD,ymax=meanOD+sdOD),size=0.5,width=0.3) +
  geom_point(aes(y=meanGlucose,color='Glucose')) +
  geom_line(aes(y=meanGlucose,color='Glucose')) +
  geom_errorbar(aes(ymin=meanGlucose-sdGlucose,ymax=meanGlucose+sdGlucose),size=0.5,width=0.3) +
  geom_point(aes(y=meanEthanol,color='Ethanol')) +
  geom_line(aes(y=meanEthanol,color='Ethanol')) +
  geom_errorbar(aes(ymin=meanEthanol-sdEthanol,ymax=meanEthanol+sdEthanol),size=0.5,width=0.3) +
  geom_point(aes(y=param_enlarge*meanGlycerol,color='Glycerol')) +
  geom_line(aes(y=param_enlarge*meanGlycerol,color='Glycerol')) +
  geom_errorbar(aes(ymin=param_enlarge*(meanGlycerol-sdGlycerol),ymax=param_enlarge*(meanGlycerol+sdGlycerol)),size=0.5,width=0.3) +
  geom_point(aes(y=param_enlarge*meanAcetate,color='Acetate')) +
  geom_line(aes(y=param_enlarge*meanAcetate,color='Acetate')) +
  geom_errorbar(aes(ymin=param_enlarge*(meanAcetate-sdAcetate),ymax=param_enlarge*(meanAcetate+sdAcetate)),size=0.5,width=0.3) +
  geom_point(aes(y=param_enlarge*meanPyruvate,color='Pyruvate')) +
  geom_line(aes(y=param_enlarge*meanPyruvate,color='Pyruvate')) +
  geom_errorbar(aes(ymin=param_enlarge*(meanPyruvate-sdPyruvate),ymax=param_enlarge*(meanPyruvate+sdPyruvate)),size=0.5,width=0.3) +
  geom_point(aes(y=param_enlarge*meanSuccinate,color='Succinate')) +
  geom_line(aes(y=param_enlarge*meanSuccinate,color='Succinate')) +
  geom_errorbar(aes(ymin=param_enlarge*(meanSuccinate-sdSuccinate),ymax=param_enlarge*(meanSuccinate+sdSuccinate)),size=0.5,width=0.3) +
  scale_y_continuous(name = 'OD, Glucose, Ethanol (g/l)',
                     sec.axis = sec_axis(~./param_enlarge,breaks = seq(0,3,0.25),name='Acetate, Glycerol, Pyruvate, Succinate (g/l)',
                     ),limits = c(-0.1,22),breaks = seq(0,22,2),expand = c(0,0)) +
  scale_x_continuous(limits = c(-0.1,20),breaks = seq(0,20,2),expand=c(0,0)) +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.border = element_rect(color = 'black',size = 0.5,fill = NA),
        axis.text = element_text(size=8,face = 'bold',color = 'black'),
        axis.title = element_text(size=8,face = 'bold',color = 'black'),
        legend.position = 'top',
        legend.title = element_blank(),
        axis.ticks = element_line(color = 'black')) +
  xlab('Time (h)')

# Saccharomyces cerevisiae
filename <- 'HPLC_cerevisiae.csv'
data_cerevisiae <- read.csv(filename,header = TRUE,sep = ';',dec = ',')

data_cerevisiae <- data_cerevisiae %>%
  rowwise() %>%
  mutate(meanOD = mean(c(OD_repl1,OD_repl2,OD_repl3))) %>%
  mutate(sdOD = sd(c(OD_repl1,OD_repl2,OD_repl3))) %>%
  mutate(meanGlucose = mean(c(Glucose_repl1,Glucose_repl2,Glucose_repl3))) %>%
  mutate(sdGlucose = sd(c(Glucose_repl1,Glucose_repl2,Glucose_repl3))) %>%
  mutate(meanGlycerol = mean(c(Glycerol_repl1,Glycerol_repl2,Glycerol_repl3))) %>%
  mutate(sdGlycerol = sd(c(Glycerol_repl1,Glycerol_repl2,Glycerol_repl3))) %>%
  mutate(meanAcetate = mean(c(Acetate_repl1,Acetate_repl2,Acetate_repl3))) %>%
  mutate(sdAcetate = sd(c(Acetate_repl1,Acetate_repl2,Acetate_repl3))) %>%
  mutate(meanEthanol = mean(c(Ethanol_repl1,Ethanol_repl2,Ethanol_repl3))) %>%
  mutate(sdEthanol = sd(c(Ethanol_repl1,Ethanol_repl2,Ethanol_repl3))) %>%
  mutate(meanPyruvate = mean(c(Pyruvate_repl1,Pyruvate_repl2,Pyruvate_repl3))) %>%
  mutate(sdPyruvate = sd(c(Pyruvate_repl1,Pyruvate_repl2,Pyruvate_repl3))) %>%
  mutate(meanSuccinate = mean(c(Succinate_repl1,Succinate_repl2,Succinate_repl3))) %>%
  mutate(sdSuccinate = sd(c(Succinate_repl1,Succinate_repl2,Succinate_repl3)))

# Plot metabolite concentration over time
param_enlarge <- 20

ggplot(data_cerevisiae,aes(x = Time)) +
  theme_bw() +
  geom_point(aes(y=meanOD,color = 'OD')) +
  geom_line(aes(y=meanOD,color ='OD')) +
  geom_errorbar(aes(ymin=meanOD-sdOD,ymax=meanOD+sdOD),size=0.5,width=0.3) +
  geom_point(aes(y=meanGlucose,color='Glucose')) +
  geom_line(aes(y=meanGlucose,color='Glucose')) +
  geom_errorbar(aes(ymin=meanGlucose-sdGlucose,ymax=meanGlucose+sdGlucose),size=0.5,width=0.3) +
  geom_point(aes(y=meanEthanol,color='Ethanol')) +
  geom_line(aes(y=meanEthanol,color='Ethanol')) +
  geom_errorbar(aes(ymin=meanEthanol-sdEthanol,ymax=meanEthanol+sdEthanol),size=0.5,width=0.3) +
  geom_point(aes(y=param_enlarge*meanGlycerol,color='Glycerol')) +
  geom_line(aes(y=param_enlarge*meanGlycerol,color='Glycerol')) +
  geom_errorbar(aes(ymin=param_enlarge*(meanGlycerol-sdGlycerol),ymax=param_enlarge*(meanGlycerol+sdGlycerol)),size=0.5,width=0.3) +
  geom_point(aes(y=param_enlarge*meanAcetate,color='Acetate')) +
  geom_line(aes(y=param_enlarge*meanAcetate,color='Acetate')) +
  geom_errorbar(aes(ymin=param_enlarge*(meanAcetate-sdAcetate),ymax=param_enlarge*(meanAcetate+sdAcetate)),size=0.5,width=0.3) +
  geom_point(aes(y=param_enlarge*meanPyruvate,color='Pyruvate')) +
  geom_line(aes(y=param_enlarge*meanPyruvate,color='Pyruvate')) +
  geom_errorbar(aes(ymin=param_enlarge*(meanPyruvate-sdPyruvate),ymax=param_enlarge*(meanPyruvate+sdPyruvate)),size=0.5,width=0.3) +
  geom_point(aes(y=param_enlarge*meanSuccinate,color='Succinate')) +
  geom_line(aes(y=param_enlarge*meanSuccinate,color='Succinate')) +
  geom_errorbar(aes(ymin=param_enlarge*(meanSuccinate-sdSuccinate),ymax=param_enlarge*(meanSuccinate+sdSuccinate)),size=0.5,width=0.3) +
  scale_y_continuous(name='OD, Glucose, Ethanol (g/l)',
                     sec.axis = sec_axis(~./param_enlarge,breaks = seq(0,2,0.25),name='Acetate, Glycerol, Pyruvate, Succinate (g/l)'),
                     limits = c(-0.1,22),breaks = seq(0,22,2),expand = c(0,0)) +
  scale_x_continuous(limits = c(0,16.2),breaks = seq(0,16,2),expand = c(0,0)) +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.border = element_rect(color = 'black',size = 0.5,fill = NA),
        axis.text = element_text(size=8,face = 'bold',color = 'black'),
        axis.title = element_text(size=8,face = 'bold',color = 'black'),
        legend.position = 'top',
        legend.title = element_blank(),
        axis.ticks = element_line(color = 'black')) +
  xlab('Time (h)')

# fig 1E
# load data on exchange fluxes
filename <- 'fluxProfile.csv'
fluxData <- read.csv(filename,sep = ';')

# Calculate the carbon recovery
# Initiate a data frame for storing results
C_recovery <-  data.frame(matrix(nrow = 4,ncol = 10))
colnames(C_recovery) <- c('Organism','CO2','Biomass','Ethanol','Acetate','Glycerol',
                          'Pyruvate','Succinate','Glucose','Total recovered')
C_recovery$Organism <- fluxData$Organism
# (q(EtOH)*2+q(Acetate)*2+q(Glycerol)*3+q(CO2)*1+q(Biomass)*1)/(-q(Glucose)*6)
MW_biomass <- 1000*1/23.58 # g/mol -> mmol/g, assumed to be similar for all 4 species

species <- fluxData$Organism

# Go through data for each species and calculate the carbon revovery (C-mol/C-mol glucose)
for(organism in species){
  rowInd <- which(fluxData$Organism == organism)
  qGlucose <- 6*fluxData$qGlucose[rowInd]
  C_recovery$CO2[rowInd] <- fluxData$CER[rowInd]/qGlucose
  C_recovery$Biomass[rowInd] <- MW_biomass*fluxData$gRate[rowInd]/qGlucose
  C_recovery$Ethanol[rowInd] <- 2*fluxData$qEthanol[rowInd]/qGlucose
  C_recovery$Acetate[rowInd] <- 2*fluxData$qAcetate[rowInd]/qGlucose
  C_recovery$Glycerol[rowInd] <- 3*fluxData$qGlycerol[rowInd]/qGlucose
  C_recovery$Pyruvate[rowInd] <- 3*fluxData$qPyruvate[rowInd]/qGlucose
  C_recovery$Succinate[rowInd] <- 4*fluxData$qSuccinate[rowInd]/qGlucose
  C_recovery$Glucose[rowInd] <- qGlucose
}

C_recovery$`Total recovered` <- rowSums(C_recovery[,2:8])
C_recovery_normalized <- C_recovery %>%
  rowwise %>%
  mutate(CO2 = CO2/`Total recovered`,
         Biomass = Biomass/`Total recovered`,
         Ethanol = Ethanol/`Total recovered`,
         Acetate = Acetate/`Total recovered`,
         Glycerol = Glycerol/`Total recovered`,
         Pyruvate = Pyruvate/`Total recovered`,
         Succinate = Succinate/`Total recovered`)

# Create dataset for plotting
dataToPlot <- C_recovery_normalized[,1:8]
dataToPlot <- melt(dataToPlot,id.vars = 'Organism')
dataToPlot$Organism <- factor(dataToPlot$Organism,levels = c('Sstip','Kma','Spo','Sce'))
#dataToPlot$Organism <- factor(dataToPlot$Organism,levels = c('Spo','Sce','Kma','Sstip'))
dataToPlot$variable <- factor(dataToPlot$variable,levels = c('Biomass','CO2','Ethanol',
                                                             'Glycerol','Acetate','Succinate',
                                                             'Pyruvate'))

# Create a stacked bar plot
fig1E <- ggplot(dataToPlot,aes(x=Organism,y=value*100,fill=forcats::fct_rev(variable))) +
  geom_bar(position = 'stack',stat = 'identity') +
  theme_classic() +
  theme(axis.text = element_text(size=8,color = 'black',face = 'bold'),
        axis.title.y = element_text(size = 8,color = 'black',face = 'bold'),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position = 'top',
        legend.text = element_text(size = 8,color = 'black',face = 'bold')) +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_brewer(palette = 'Paired',direction = -1) +
  ylab('Fraction of glucose consumed (%)')
fig1E
# Save figure
#ggsave('../../results/fig1E.pdf',
 #      fig1E, device = 'pdf',height = 7, width = 7, units = 'cm',dpi = 300)

# fig 1F
glc_uptake_replicates <- read.delim('glcUptakeRate.csv',sep = ';')
# Calculate mean and sd
glc_uptake_replicates$mean <- rowMeans(glc_uptake_replicates[,1:3])
glc_uptake_replicates$sd <- apply(glc_uptake_replicates[,1:3], 1, sd)
df1 <- glc_uptake_replicates[,1:4]
df1 <- melt(glc_uptake_replicates,id.vars = 'organism',measure.vars = c('repl1','repl2','repl3'))
df1$organism <- factor(df1$organism,levels = c('Sstip','Kma','Spo','Sce'))

fig1F <- ggplot() +
  geom_bar(data=glc_uptake_replicates,aes(x=reorder(organism,mean),y=mean,fill=organism),
           stat='identity',color='black',position=position_dodge(),alpha=0.875) +
  geom_errorbar(data=glc_uptake_replicates,aes(x=organism,ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(.9)) +
  geom_jitter(data = df1,aes(x=organism,y=value,color=organism),
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8)) +
  theme_classic() +
  theme(plot.title = element_text(size = 10,face = 'bold',color = 'black',hjust = 0.5),
        axis.text = element_text(size = 8,face = 'bold',color = 'black'),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 8,face = 'bold',color = 'black'),
        panel.border = element_rect(size = 1.1, color = 'black',fill = NA),
        legend.title = element_blank(),
        legend.position = c(0.15,0.85),
        legend.text = element_text(size = 8,face = 'bold',color = 'black')) +
  scale_fill_manual(values=c('#1D71B8','#000000','#CBBBA0','#878787')) +
  scale_color_manual(values=c('#1D71B8','#000000','#CBBBA0','#878787')) +
  scale_y_continuous(limits = c(0,21.2),breaks = seq(0,20,2), expand = c(0,0)) +
  ylab('Glucose uptake rate (mmol/gDW/h)')
fig1F
# Save figure
#ggsave('../../results/fig1F.pdf',
 #      fig1F, device = 'pdf',height = 7, width = 7, units = 'cm',dpi = 300)

# fig 1G
data <- read.delim('Yxs_gRate_Ptot.csv',sep = ';',dec = '.')
Yxs <- data[grepl('Yxs',data$Parameter),]
Yxs$GUR <- c(18.8,8.9,5.6,4.5)
Yxs$organism <- factor(Yxs$organism,levels = c('Sstip','Kma','Spo','Sce'))
df1 <- melt(Yxs[,c(2:5,8)],id.vars = c('organism','GUR'),measure.vars = c('repl1','repl2','repl3'))
Ptot <- data[grepl('Ptot',data$Parameter),]
Ptot$GUR <- c(18.8,8.9,5.6,4.5)
Ptot$organism <- factor(Ptot$organism,levels = c('Sstip','Kma','Spo','Sce'))
df2 <- melt(Ptot[,c(2:5,8)],id.vars = c('organism','GUR'),measure.vars = c('repl1','repl2','repl3'))
gRate <- data[grepl('gRate',data$Parameter),]
gRate$GUR <- c(18.8,8.9,5.6,4.5)
gRate$organism <- factor(gRate$organism,levels = c('Sstip','Kma','Spo','Sce'))
df3 <- melt(gRate[,c(2:5,8)],id.vars = c('organism','GUR'),measure.vars = c('repl1','repl2','repl3'))

coeff <- 1
fig1G <- ggplot() +
  geom_errorbar(data=Yxs,aes(x=GUR,ymin=mean-sd,ymax=mean+sd),
                width=.2,position=position_dodge(.9)) +
  geom_point(data=Yxs,aes(x=GUR,y=mean),color='#08519C',shape=22,
             fill = '#08519C') +
  geom_jitter(data=df1,aes(x=GUR,y=value),color='#08519C',
              position = position_jitter(0.2),alpha=0.85) +
  geom_errorbar(data=Ptot,aes(x=GUR,ymin=mean-sd,ymax=mean+sd),
                width=.2,position=position_dodge(.9)) +
  geom_point(data=Ptot,aes(x=GUR,y=mean),color='#6BAED6',shape=22,
             fill = '#6BAED6') +
  geom_jitter(data=df2,aes(x=GUR,y=value),color='#6BAED6',
              position = position_jitter(0.2),alpha=0.85) +
  geom_errorbar(data=gRate,aes(x=GUR,ymin=mean-sd,ymax=mean+sd),
                width=.2,position=position_dodge(.9)) +
  geom_point(data=gRate,aes(x=GUR,y=mean),color='black',shape=22,
             fill = 'black',alpha=0.85) +
  geom_jitter(data=df3,aes(x=GUR,y=value),color='black',
              position = position_jitter(0.2),alpha=0.85) +
  scale_y_continuous(name = 'Total protein content (g/gDW), growth rate (h-1)',
                     sec.axis = sec_axis(~.*coeff,name = 'Biomass yield (gDW/g glucose)')) +
  scale_x_continuous(name = 'rGlu (mmol/gDW/h)',breaks=c(4.5,5.6,8.9,18.8)) +
  theme_classic() +
  theme(axis.title = element_text(size = 8,face = 'bold',color = 'black'),
        axis.text = element_text(size = 8,face = 'bold',color = 'black'))
fig1G
#ggsave('../../results/fig1G.pdf',
 #      fig1G,width=7,height = 7,units = 'cm',device = 'pdf',dpi=300)

# fig 1H
data <- read.delim('OUR_CER_RQ.csv',sep = ';',dec = ',')
OUR <- data[grepl('OUR',data$Parameter),]
OUR$GUR <- c(18.8,8.9,5.6,4.5)
OUR$organism <- factor(OUR$organism,levels = c('Sstip','Kma','Spo','Sce'))
df1 <- melt(OUR[,c(2:5,8)],id.vars = c('organism','GUR'),measure.vars = c('repl1','repl2','repl3'))
df1$organism <- factor(df1$organism,levels=c('Sstip','Kma','Spo','Sce'))
CER <- data[grepl('CER',data$Parameter),]
CER$GUR <- c(18.8,8.9,5.6,4.5)
CER$organism <- factor(CER$organism,levels = c('Sstip','Kma','Spo','Sce'))
df2 <- melt(CER[,c(2:5,8)],id.vars = c('organism','GUR'),measure.vars = c('repl1','repl2','repl3'))
RQ <- data[grepl('RQ',data$Parameter),]
RQ$GUR <- c(18.8,8.9,5.6,4.5)
RQ$organism <- factor(RQ$organism,levels = c('Sstip','Kma','Spo','Sce'))
df3 <- melt(RQ[,c(2:5,8)],id.vars = c('organism','GUR'),measure.vars = c('repl1','repl2','repl3'))

# Plot OUR, CER vs glucose uptake
coeff <- 0.1
fig1H <- ggplot() +
  geom_errorbar(data=OUR,aes(x=GUR,ymin=mean-sd,ymax=mean+sd),
                width=.2,position=position_dodge(.9)) +
  geom_point(data=OUR,aes(x=GUR,y=mean),color='#6BAED6',shape=22,
             fill = '#6BAED6') +
  geom_jitter(data=df1,aes(x=GUR,y=value),color='#6BAED6',
              position = position_jitter(0.2),alpha=0.85) +
  geom_errorbar(data=CER,aes(x=GUR,ymin=mean-sd,ymax=mean+sd),
                width=.2,position=position_dodge(.9)) +
  geom_point(data=CER,aes(x=GUR,y=mean),color='#08519C',shape=22,
             fill = '#08519C') +
  geom_jitter(data=df2,aes(x=GUR,y=value),color='#08519C',
              position = position_jitter(0.2),alpha=0.85) +
  geom_errorbar(data=RQ,aes(x=GUR,ymin=mean-sd,ymax=mean+sd),
                width=.2,position=position_dodge(.9)) +
  geom_point(data=RQ,aes(x=GUR,y=mean),color='black',shape=22,
             fill = 'black') +
  geom_jitter(data=df3,aes(x=GUR,y=value),color='black',
              position = position_jitter(0.2),alpha=0.85) +
  scale_y_continuous(name = 'OUR, CER (mmol/gDW/h), RQ',limits = c(0,29),
                     breaks = c(0,1,2,3,4,5,6,7,8,9,10,15,20,25), expand = c(0,0)) +
  scale_x_continuous(name = 'rGlu (mmol/gDW/h)',breaks=c(4.5,5.6,8.9,18.8)) +
  theme_classic() +
  theme(axis.title = element_text(size = 8,face = 'bold',color = 'black'),
        axis.text = element_text(size = 8,face = 'bold',color = 'black'),
        legend.position = 'top')
fig1H
# Save figure
#ggsave('../../results/fig1G.pdf',
 #      fig1H, device = 'pdf',height = 7, width = 7, units = 'cm',dpi = 300)

