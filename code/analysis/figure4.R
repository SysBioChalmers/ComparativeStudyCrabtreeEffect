# Load required packages
library(readr)
library(readxl)
library(ggplot2)
library(ggpubr)

# load data
setwd('../../data/experimentalData')

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

# cytoplasmic ribosomes

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

# for Scer and Spombe, bin the paralogs:
test <- RP_Scer
for(i in 1:nrow(test)){
  if(grepl("-A", test[i,"Protein names"])){
    test[i,6:17] <- colSums(test[c(i,i+1),6:17], na.rm = T)
    test[i+1,] <- NA
  }
}
test[120,] <- NA #this row is across-the-board zero
test[,3:5] <- "x" #fill these so the below na.omit doesn't kill everything
RP_Scer_bin <- na.omit(test)
RP_Scer_bin <- RP_Scer_bin[,c(2,6:17)]


test <- RP_Spombe
for(i in 1:nrow(test)){
  if(grepl("-A", test[i,"Protein names"])){
    test[i,6:17] <- colSums(test[c(i,i+1),6:17], na.rm = T)
    test[i+1,] <- NA
  }
}
for(i in 1:nrow(test)){
  if(grepl("-C", test[i,"Protein names"])){
    test[i-2,6:17] <- colSums(test[c(i,i-1, i-2),6:17], na.rm = T)
    test[i,] <- NA
    test[i-1,] <- NA
  }
}

test[140,] <- NA #paralog with row 141
test[143,] <- NA #paralog with row 142
#test[58:59,1] <- NA
test[125,] <- NA #across-the-board zero
test[,3:5] <- "x" #fill these so the below na.omit doesn't kill everything
RP_Spombe_bin <- na.omit(test)
RP_Spombe_bin <- RP_Spombe_bin[,c(2,6:17)]


nmolRibo.mgDW.binned <- c(colMeans(RP_Scer_bin[,11:13], na.rm=T), colMeans(RP_Kmarx[,15:17], na.rm=T), 
                          colMeans(RP_Spombe_bin[,11:13], na.rm=T), colMeans(RP_Pstip[,15:17], na.rm=T))

ribo.efficiency <- gProt.gDW*growth/nmolRibo.mgDW.binned
#now can plot Fig 4A

# pull out mitoribosomes
# for Scer, Spombe, Kmarx this is from Uniprot
# for Pstip, if pull from Uniprot, we only find 21 MRPs. so for Pstip this is done based on orthologs with Scer

MRP_Scer <- uniprot_Scer[grep("ribosom",uniprot_Scer$`Protein names`),]
MRP_Scer <- MRP_Scer[grep("mito", MRP_Scer$`Protein names`),]
MRP_Scer <- MRP_Scer[order(MRP_Scer$`Protein names`),]
MRP_Scer[76,] #not a mitoribo subunit
MRP_Scer[79,] #not a mitoribo subunit
MRP_Scer <- MRP_Scer[-c(76,79),]
MRP_Scer <- merge(MRP_Scer[,c(1,3)], Scer[,c(1,3:14)], by.x="Entry", by.y="Accession", all.x=T)


MRP_Spombe <- uniprot_Spombe[grep("ribosom",uniprot_Spombe$`Protein names`),]
MRP_Spombe <- MRP_Spombe[grep("mito", MRP_Spombe$`Protein names`),]
MRP_Spombe <- MRP_Spombe[order(MRP_Spombe$`Protein names`),]
MRP_Spombe[72:75,] #not a mitoribo subunit
MRP_Spombe <- MRP_Spombe[1:71,]
MRP_Spombe <- merge(MRP_Spombe[,c(1,3)], Spombe[,c(1,3:14)], by.x="Entry", by.y="Accession", all.x=T)


MRP_Kmarx <- uniprot_Kmarx[grep("ribosom",uniprot_Kmarx$`Protein names`),]
MRP_Kmarx <- MRP_Kmarx[grep("37S|54S|mito|Mito|30S", MRP_Kmarx$`Protein names`),]
MRP_Kmarx <- MRP_Kmarx[order(MRP_Kmarx$`Protein names`),]
MRP_Kmarx <- merge(MRP_Kmarx[,c(1,3)], Kmarx[,c(1,3:14)], by.x="Entry", by.y="Accession", all.x=T)


#MRP_Pstip <- uniprot_Pstip[grep("ribosom|Ribosom",uniprot_Pstip$`Protein names`),]
#MRP_Pstip <- MRP_Pstip[grep("mito|Mito|37S|54S", MRP_Pstip$`Protein names`),]
#MRP_Pstip <- MRP_Pstip[order(MRP_Pstip$`Protein names`),]
#MRP_Pstip[22,] #not a mitoribo subunit
#MRP_Pstip <- MRP_Pstip[1:21,]
#MRP_Pstip <- merge(MRP_Pstip[,c(1,3)], Pstip[,c(1,3:14)], by.x="Entry", by.y="Accession", all.x=T)
# only 21 MRP subunits found, remove this and pull from Scer orthologs instead
#remove(MRP_Pstip)

#pull in ortholog list:
orthologs_anyCopy <- read_delim("../databases/anyCopyOrthologs_allSpecies.txt", 
                                "\t", escape_double = FALSE, trim_ws = TRUE)


MRP_Pstip <- as.data.frame(MRP_Scer[,1])
colnames(MRP_Pstip) <- "Saccharomyces_cerevisiaeS288c"

MRP_Pstip$Orthogroup <- "a"
MRP_Pstip$Scheffersomyces_stipitis_CBS6054 <- "a"

for(i in 1:nrow(MRP_Pstip)){
  gene_name_to_match <- paste(MRP_Pstip[i,1], "\\b", sep="")
  MRP_Pstip[i,c("Orthogroup", "Scheffersomyces_stipitis_CBS6054")] <- 
    as.data.frame(orthologs_anyCopy[grep(gene_name_to_match, orthologs_anyCopy$Saccharomyces_cerevisiaeS288c),
                                    c("Orthogroup", "Scheffersomyces_stipitis_CBS6054")])
  
}

MRP_Pstip <- as.data.frame(unique(MRP_Pstip$Scheffersomyces_stipitis_CBS6054))
MRP_Pstip[,1] <- as.character(MRP_Pstip[,1])
MRP_Pstip[c(2,20,54),1] #clean these up by hand
MRP_Pstip[2,1] <- "A3LXK5"
MRP_Pstip[20,1] <- "A3LQ56"
MRP_Pstip[54,1] <- "A3GGQ1"
MRP_Pstip[73,1] <- "A3GGQ2"

colnames(MRP_Pstip) <- "Entry"
MRP_Pstip <- merge(MRP_Pstip, Pstip, by.x="Entry", by.y="Accession", all.x=T)
#now can plot Fig 4B

# pull out ATP synthase and COX proteins
# for Scer, Spombe, Kmarx this is from Uniprot
# for Pstip this is done based on orthologs with Scer

ATP_synthase_Scer <- uniprot_Scer[grep("ATP synthase subunit",uniprot_Scer$`Protein names`),]
ATP_synthase_Spombe <- uniprot_Spombe[grep("ATP synthase subunit",uniprot_Spombe$`Protein names`),]
ATP_synthase_Kmarx <- uniprot_Kmarx[grep("ATP synthase subunit",uniprot_Kmarx$`Protein names`),]

#Pstip: pull out from Sce orthologs
ATP_synthase_Pstip <- as.data.frame(ATP_synthase_Scer[,1])
colnames(ATP_synthase_Pstip) <- "Saccharomyces_cerevisiaeS288c"
ATP_synthase_Pstip$Orthogroup <- "a"
ATP_synthase_Pstip$Scheffersomyces_stipitis_CBS6054 <- "a"
for(i in 1:nrow(ATP_synthase_Pstip)){
  gene_name_to_match <- paste(ATP_synthase_Pstip[i,1], "\\b", sep="")
  ATP_synthase_Pstip[i,c("Orthogroup", "Scheffersomyces_stipitis_CBS6054")] <- 
    as.data.frame(orthologs_anyCopy[grep(gene_name_to_match, orthologs_anyCopy$Saccharomyces_cerevisiaeS288c),
                                    c("Orthogroup", "Scheffersomyces_stipitis_CBS6054")])
  
}
colnames(ATP_synthase_Pstip)[3] <- "Entry"
ATP_synthase_Pstip <- merge(uniprot_Pstip, ATP_synthase_Pstip)
ATP_synthase_Pstip <- ATP_synthase_Pstip[,1:3]

ATP_synthase_Scer <- merge(ATP_synthase_Scer, Scer, by.x="Entry", by.y="Accession", all.x=T)
ATP_synthase_Spombe <- merge(ATP_synthase_Spombe, Spombe, by.x="Entry", by.y="Accession", all.x=T)
ATP_synthase_Kmarx <- merge(ATP_synthase_Kmarx, Kmarx, by.x="Entry", by.y="Accession", all.x=T)
ATP_synthase_Pstip <- merge(ATP_synthase_Pstip, Pstip, by.x="Entry", by.y="Accession", all.x=T)

cytochrome_c_Scer <- uniprot_Scer[grep("Cytochrome c oxidase subunit",uniprot_Scer$`Protein names`),]
cytochrome_c_Spombe <- uniprot_Spombe[grep("Cytochrome c oxidase subunit",uniprot_Spombe$`Protein names`),]
cytochrome_c_Kmarx <- uniprot_Kmarx[grep("Cytochrome c oxidase subunit",uniprot_Kmarx$`Protein names`),]

#Pstip: pull out from Sce orthologs
cytochrome_c_Pstip <- as.data.frame(cytochrome_c_Scer[,1])
colnames(cytochrome_c_Pstip) <- "Saccharomyces_cerevisiaeS288c"
cytochrome_c_Pstip$Orthogroup <- "a"
cytochrome_c_Pstip$Scheffersomyces_stipitis_CBS6054 <- "a"
for(i in 1:nrow(cytochrome_c_Pstip)){
  gene_name_to_match <- paste(cytochrome_c_Pstip[i,1], "\\b", sep="")
  cytochrome_c_Pstip[i,c("Orthogroup", "Scheffersomyces_stipitis_CBS6054")] <- 
    as.data.frame(orthologs_anyCopy[grep(gene_name_to_match, orthologs_anyCopy$Saccharomyces_cerevisiaeS288c),
                                    c("Orthogroup", "Scheffersomyces_stipitis_CBS6054")])
  
}
colnames(cytochrome_c_Pstip)[3] <- "Entry"
cytochrome_c_Pstip <- merge(uniprot_Pstip, cytochrome_c_Pstip)
cytochrome_c_Pstip <- cytochrome_c_Pstip[,1:3]

cytochrome_c_Scer <- merge(cytochrome_c_Scer, Scer, by.x="Entry", by.y="Accession", all.x=T)
cytochrome_c_Spombe <- merge(cytochrome_c_Spombe, Spombe, by.x="Entry", by.y="Accession", all.x=T)
cytochrome_c_Kmarx <- merge(cytochrome_c_Kmarx, Kmarx, by.x="Entry", by.y="Accession", all.x=T)
cytochrome_c_Pstip <- merge(cytochrome_c_Pstip, Pstip, by.x="Entry", by.y="Accession", all.x=T)



# take mean (molar) of ATP synthase, cytochome c oxidase, and mitoribo subunits, then multiply by MW
# unit is ng/ugProt
ATPsynthase.ng.ugProt <- c(colMeans(ATP_synthase_Scer[,8:10], na.rm=T), colMeans(ATP_synthase_Kmarx[,8:10], na.rm=T),
                           colMeans(ATP_synthase_Spombe[,8:10], na.rm=T), colMeans(ATP_synthase_Pstip[,8:10], na.rm=T) ) * (29.099+5.822+7.759) /1000
cytochrome.ng.ugProt <- c(colMeans(cytochrome_c_Scer[,8:10], na.rm=T), colMeans(cytochrome_c_Kmarx[,8:10], na.rm=T),
                          colMeans(cytochrome_c_Spombe[,8:10], na.rm=T), colMeans(cytochrome_c_Pstip[,8:10], na.rm=T) ) * (58.798+28.567+30.360+34.055) /1000
mitoribosome.ng.ugProt <- c(colMeans(MRP_Scer[,6:8], na.rm=T), colMeans(MRP_Kmarx[,6:8], na.rm=T),
                            colMeans(MRP_Spombe[,6:8], na.rm=T), colMeans(MRP_Pstip[,6:8], na.rm=T) ) * 47.123 / 100


# calculate mitoribosome efficiency

MRP_translated_ng_ugProt <- ATPsynthase.ng.ugProt + cytochrome.ng.ugProt + mitoribosome.ng.ugProt
names(MRP_translated_ng_ugProt) <- rep(c("ng/ugProt_1", "ng/ugProt_2", "ng/ugProt_3"),4)

nmolMitoribo.mgProt <- c(colMeans(MRP_Scer[,6:8], na.rm=T), colMeans(MRP_Kmarx[,6:8], na.rm=T),
                         colMeans(MRP_Spombe[,6:8], na.rm=T), colMeans(MRP_Pstip[,6:8], na.rm=T) )
mitoribo.efficiency <- MRP_translated_ng_ugProt * growth / nmolMitoribo.mgProt
# now can plot Fig 4C

#take sum (g/gProt) of ATP synthase and cytochrome c oxidase & turn into percentage (multiply by 100)
ATPsynth.cytoch.percent <- (c(colSums(ATP_synthase_Scer[,11:13], na.rm=T), 
                              colSums(ATP_synthase_Kmarx[,11:13], na.rm=T),
                              colSums(ATP_synthase_Spombe[,11:13], na.rm=T), 
                              colSums(ATP_synthase_Pstip[,11:13], na.rm=T) ) + c(colSums(cytochrome_c_Scer[,11:13], na.rm=T), 
                                                                                 colSums(cytochrome_c_Kmarx[,11:13], na.rm=T), 
                                                                                 colSums(cytochrome_c_Spombe[,11:13], na.rm=T), 
                                                                                 colSums(cytochrome_c_Pstip[,11:13], na.rm=T) ) )*100
# now can plot Fig 4D

# Create figures

# fig 4A
ribo_eff_Sce <- t(ribo.efficiency[1:3])
ribo_eff_Kma <- t(ribo.efficiency[4:6])
ribo_eff_Spo <- t(ribo.efficiency[7:9])
ribo_eff_Sstip <- t(ribo.efficiency[10:12])
dataToPlot <- rbind(ribo_eff_Sce,ribo_eff_Kma,ribo_eff_Spo,ribo_eff_Sstip)
dataToPlot$organism <- c('Sce','Kma','Spo','Sstip')
colnames(dataToPlot)[1:3] <- c('repl1','repl2','repl3')
# Calculate mean and sd
dataToPlot$mean <- rowMeans(dataToPlot[,1:3])
dataToPlot$sd <- apply(dataToPlot[,1:3], 1, sd)
dataToPlot$organism <- factor(dataToPlot$organism,levels=c('Sstip','Kma','Spo','Sce'))
df1 <- melt(dataToPlot[,1:4],id.vars = 'organism',measure.vars = c('repl1','repl2','repl3'))

# Plot figure
fig4A <- ggplot(dataToPlot,aes(x=organism,y=mean,fill=organism)) +
  geom_bar(position=position_dodge2(), stat="identity",color='black',alpha=0.85) +
  geom_errorbar(aes(ymin=(mean-sd),ymax=(mean+sd)),width=0.4,
                position=position_dodge(0.9)) +
  geom_jitter(data=df1,aes(x=organism,y=value,color=organism),
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),size=1) +
  scale_y_continuous(breaks = seq(0,4,1),limits = c(0,4),expand = c(0,0)) +
  scale_fill_manual(values=c('#878787','#1D71B8','#CBBBA0','#000000')) +
  scale_color_manual(values=c('#878787','#1D71B8','#CBBBA0','#000000')) +
  theme_classic() +
  theme(text = element_text(size = 8, color = 'black',face = 'bold'),
        legend.position = 'top',
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10,face = 'bold',color = 'black'),
        axis.text.x = element_text(size = 7,face = 'bold',color = 'black'),
        axis.text.y = element_text(size = 8,face = 'bold',color = 'black'),
        axis.ticks = element_line(color = 'black')) +
  ylab('mg protein / nmol ribosomes / h')
fig4A
# Save figure
#ggsave('../../results/fig4A.pdf',
 #      fig4A,width=6,height = 6,units = 'cm',device = 'pdf',dpi=300)

# fig 4B
# collect data, converting g/g protein to %
mitoRibo_Scer <- data.frame(t(colSums(MRP_Scer[,9:11],na.rm = TRUE)*100))
mitoRibo_Spombe <- data.frame(t(colSums(MRP_Spombe[,9:11],na.rm = TRUE)*100))
mitoRibo_Kmarx <- data.frame(t(colSums(MRP_Kmarx[,9:11],na.rm = TRUE)*100))
mitoRibo_Sstip <- data.frame(t(colSums(MRP_Pstip[,9:11],na.rm = TRUE)*100))
mitoRibo_combined <- rbind(mitoRibo_Scer,mitoRibo_Spombe,mitoRibo_Kmarx,mitoRibo_Sstip)
colnames(mitoRibo_combined) <- c('repl1','repl2','repl3')
mitoRibo_combined$organism <- c('Sce','Spo','Kma','Sstip')
mitoRibo_combined$organism <- factor(mitoRibo_combined$organism,levels=c('Sstip','Kma','Spo','Sce'))
# Calculate mean and sd
mitoRibo_combined$mean <- rowMeans(mitoRibo_combined[,1:3])
mitoRibo_combined$sd <- apply(mitoRibo_combined[,1:3],1,sd)
df1 <- melt(mitoRibo_combined[,1:4],id.vars = 'organism',measure.vars = c('repl1','repl2','repl3'))

# plot data
fig4B <- ggplot(mitoRibo_combined,aes(x=organism,y=mean,fill=organism)) +
  geom_bar(position=position_dodge2(), stat="identity",color='black',alpha=0.85) +
  geom_errorbar(aes(ymin=(mean-sd),ymax=(mean+sd)),width=0.4,
                position=position_dodge(0.9)) +
  geom_jitter(data=df1,aes(x=organism,y=value,color=organism),
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),size=1) +
  scale_y_continuous(breaks = seq(0,0.65,0.2),limits = c(0,0.65),expand = c(0,0)) +
  scale_fill_manual(values=c('#878787','#1D71B8','#CBBBA0','#000000')) +
  scale_color_manual(values=c('#878787','#1D71B8','#CBBBA0','#000000')) +
  theme_classic() +
  theme(text = element_text(size = 8, color = 'black',face = 'bold'),
        legend.position = 'top',
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10,face = 'bold',color = 'black'),
        axis.text.x = element_text(size = 7,face = 'bold',color = 'black'),
        axis.text.y = element_text(size = 8,face = 'bold',color = 'black'),
        axis.ticks = element_line(color = 'black')) +
  ylab('Total % (g/g) of detected mitoribosomal proteins (MRPs)')
fig4B
# Save figure
#ggsave('../../results/fig4B.pdf',
 #      fig4B,width=6,height = 6,units = 'cm',device = 'pdf',dpi=300)

# fig 4C
dataToPlot <- data.frame(rbind(t(mitoribo.efficiency[1:3]),
                               t(mitoribo.efficiency[4:6]),
                               t(mitoribo.efficiency[7:9]),
                               t(mitoribo.efficiency[10:12])))
colnames(dataToPlot) <- c('repl1','repl2','repl3')
dataToPlot$organism <- c('Sce','Kma','Spo','Sstip')
dataToPlot$organism <- factor(dataToPlot$organism,levels=c('Sstip','Kma','Spo','Sce'))
# calculate mean and sd
dataToPlot$mean <- rowMeans(dataToPlot[,1:3])
dataToPlot$sd <- apply(dataToPlot[,1:3],1,sd)
df1 <- melt(dataToPlot[,1:4],id.vars = 'organism',measure.vars = c('repl1','repl2','repl3'))

# plot data
fig4C <- ggplot(dataToPlot,aes(x=organism,y=mean,fill=organism)) +
  geom_bar(position=position_dodge2(), stat="identity",color='black',alpha=0.85) +
  geom_errorbar(aes(ymin=(mean-sd),ymax=(mean+sd)),width=0.4,
                position=position_dodge(0.9)) +
  geom_jitter(data=df1,aes(x=organism,y=value,color=organism),
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),size=1) +
  scale_y_continuous(breaks = seq(0,1.5,0.5),limits = c(0,2),expand = c(0,0)) +
  scale_fill_manual(values=c('#878787','#1D71B8','#CBBBA0','#000000')) +
  scale_color_manual(values=c('#878787','#1D71B8','#CBBBA0','#000000')) +
  theme_classic() +
  theme(text = element_text(size = 8, color = 'black',face = 'bold'),
        legend.position = 'top',
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10,face = 'bold',color = 'black'),
        axis.text.x = element_text(size = 7,face = 'bold',color = 'black'),
        axis.text.y = element_text(size = 8,face = 'bold',color = 'black'),
        axis.ticks = element_line(color = 'black')) +
  ylab('mg protein (mitoribosome-translated) / nmol mitoribosomes / h')
fig4C
# Save figure
#ggsave('../../results/fig4C.pdf',
 #      fig4C,width=6,height = 6,units = 'cm',device = 'pdf',dpi=300)

# fig 4D
dataToPlot <- data.frame(rbind(t(ATPsynth.cytoch.percent[1:3]),
                               t(ATPsynth.cytoch.percent[4:6]),
                               t(ATPsynth.cytoch.percent[7:9]),
                               t(ATPsynth.cytoch.percent[10:12])))
colnames(dataToPlot) <- c('repl1','repl2','repl3')
dataToPlot$organism <- c('Sce','Kma','Spo','Sstip')
dataToPlot$organism <- factor(dataToPlot$organism,levels=c('Sstip','Kma','Spo','Sce'))
# calculate mean and sd
dataToPlot$mean <- rowMeans(dataToPlot[,1:3])
dataToPlot$sd <- apply(dataToPlot[,1:3],1,sd)
df1 <- melt(dataToPlot[,1:4],id.vars = 'organism',measure.vars = c('repl1','repl2','repl3'))

# plot data
fig4D <- ggplot(dataToPlot,aes(x=organism,y=mean,fill=organism)) +
  geom_bar(position=position_dodge2(), stat="identity",color='black',alpha=0.85) +
  geom_errorbar(aes(ymin=(mean-sd),ymax=(mean+sd)),width=0.4,
                position=position_dodge(0.9)) +
  geom_jitter(data=df1,aes(x=organism,y=value,color=organism),
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),size=1) +
  scale_y_continuous(breaks = seq(0,4,1),limits = c(0,3.5),expand = c(0,0)) +
  scale_fill_manual(values=c('#878787','#1D71B8','#CBBBA0','#000000')) +
  scale_color_manual(values=c('#878787','#1D71B8','#CBBBA0','#000000')) +
  theme_classic() +
  theme(text = element_text(size = 8, color = 'black',face = 'bold'),
        legend.position = 'top',
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10,face = 'bold',color = 'black'),
        axis.text.x = element_text(size = 7,face = 'bold',color = 'black'),
        axis.text.y = element_text(size = 8,face = 'bold',color = 'black'),
        axis.ticks = element_line(color = 'black')) +
  ylab('Total % (g/g) of ATP synthase and cytochrome c oxidase')
fig4D
# Save figure
#ggsave('../../results/fig4D.pdf',
 #      fig4D,width=6,height = 6,units = 'cm',device = 'pdf',dpi=300)



