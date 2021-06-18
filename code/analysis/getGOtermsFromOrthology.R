# Get GO term annotation for all species based on GO annotation from S. cerevisiae and gene orthology

#if (!requireNamespace("BiocManager", quietly = TRUE))
 # install.packages("BiocManager")

#BiocManager::install("biomaRt")

# load package biomaRt
library(biomaRt)

# read uniprot information for all species
setwd('../../data/databases/')
species <- c('Scer','Kmarx','Pstip','Spombe')
uniprot_Scer <- read.delim('uniprot_Scer.txt')
uniprot_Kmarx <- read.delim('uniprot_Kmarx.txt')
uniprot_Pstip <- read.delim('uniprot_Pstip.txt')
uniprot_Spombe <- read.delim('uniprot_Spombe.txt')

# load data for S. cerevisiae
bm <- useMart("ensembl")
bm <- useDataset("scerevisiae_gene_ensembl", mart=bm)

# Get ensembl gene ids and GO term information
EG2GO <- getBM(mart=bm, attributes=c('ensembl_gene_id','uniprotswissprot','go_id','name_1006','namespace_1003'))

# Remove entries with no GO term found, or information missing
EG2GO <- EG2GO[EG2GO$go_id != '' & EG2GO$name_1006 != '' & EG2GO$namespace_1003 != 'go',]
rownames(EG2GO) <- NULL

# map genes to uniprot IDs to retrieve missing IDs
for(i in 1:nrow(EG2GO)){
  id <- EG2GO$uniprotswissprot[i]
  if(id == ''){
    geneName <- EG2GO$ensembl_gene_id[i]
    uniprotID <- uniprot_Scer$Entry[grep(geneName,uniprot_Scer$Gene.names...ordered.locus..)]
    EG2GO$uniprotswissprot[i] <- uniprotID
  }
}

# load files with gene orthology (try mapping with single copy orthologs at first)
singleCopyOG <- read.delim('singleCopyOrthologs_allSpecies.txt')
# rearrange columns
singleCopyOG <- singleCopyOG[,c(1,3,2,4,5)]


# Retrieve GO term annotation through orthology  
mapFromSingleCopyOGtoGO <- function(uniprotData,singleCopyOG,organism,GOtermRef){
  GOtermAnnotation <- data.frame(matrix(ncol = 3, nrow = 0))
  for(id in uniprotData$Entry){
    colIndex <- grep(organism,colnames(singleCopyOG))
    rowIndex <- grep(id,singleCopyOG[,colIndex])
    if(length(rowIndex)>0){
      cerevisiaeID <- singleCopyOG[rowIndex,grep('cerevisiae',colnames(singleCopyOG))]
      GOterm_ind <- grep(cerevisiaeID,GOtermRef$uniprotswissprot)
      GOtermInfo <- GOtermRef[GOterm_ind,]
      #GOtermInfo <- GOtermInfo[which(GOtermInfo$namespace_1003 == 'biological_process'),]
      GOtermInfo <- data.frame(rep(id,nrow(GOtermInfo)),GOtermInfo[,3:4])
      colnames(GOtermInfo) <- c('UniprotID','GOterm','GOname')
      rownames(GOtermInfo) <- NULL
      GOtermAnnotation <- rbind(GOtermAnnotation,GOtermInfo)
    }
  }
  return(GOtermAnnotation)
}

# Retrieve GO term annotation through single copy orthologs
GOtermData_Kmarx <- mapFromSingleCopyOGtoGO(uniprot_Kmarx,singleCopyOG,'marxianus',EG2GO)
GOtermData_Pstip <- mapFromSingleCopyOGtoGO(uniprot_Pstip,singleCopyOG,'stipitis',EG2GO)
GOtermData_Spombe <- mapFromSingleCopyOGtoGO(uniprot_Spombe,singleCopyOG,'pombe',EG2GO)

# Retrieve GO term annotation through any copy orthologs

# Read file with any copy ortholog data
anyCopyOG <- read.delim('anyCopyOrthologs_allSpecies.txt')

mapFromAnyCopyOGtoGO <- function(uniprotData,anyCopyOG,organism,GOtermRef,IDsToExclude=''){
  GOtermData <- data.frame(matrix(nrow = 0,ncol = 3))
  for(id in uniprotData$Entry){
    if(length(grep(id,IDsToExclude)) == 0){
      colIndex <- grep(organism,colnames(anyCopyOG))
      rowIndex <- grep(id,anyCopyOG[,colIndex])
      if(length(rowIndex)>0){
        cerevisiaeIDs <- anyCopyOG[rowIndex,grep('cerevisiae',colnames(anyCopyOG))]
        cerevisiaeIDs <- strsplit(cerevisiaeIDs,' ')
        cerevisiaeIDs <- cerevisiaeIDs[[1]]
        collectedGOterms <- data.frame(matrix(nrow = 0,ncol = 2))
        for(entry in cerevisiaeIDs){
          GOind <- grep(entry,GOtermRef$uniprotswissprot)
          GOterms <- GOtermRef[GOind,]
          #GOterms <- GOterms[which(GOterms$namespace_1003 == 'biological_process'),]
          GOterms <- GOterms[,3:4]
          colnames(GOterms) <- c('GOterm','GOname')
          rownames(GOterms) <- NULL
          collectedGOterms <- rbind(collectedGOterms,GOterms)
        }
        collectedGOterms <- collectedGOterms[!duplicated(collectedGOterms$GOname),]
        collectedGOterms <- cbind(rep(id,nrow(collectedGOterms)),collectedGOterms)
        colnames(collectedGOterms)[1] <- 'UniprotID'
        GOtermData <- rbind(GOtermData,collectedGOterms)
      }
    }
  }
  return(GOtermData)
}

GOtermData_Kmarx_anyCopyOG <- mapFromAnyCopyOGtoGO(uniprot_Kmarx,anyCopyOG,'marxianus',EG2GO,singleCopyOG$Kluyveromyces_marxianus_DMKU3_1042)
GOtermData_Pstip_anyCopyOG <- mapFromAnyCopyOGtoGO(uniprot_Pstip,anyCopyOG,'stipitis',EG2GO,singleCopyOG$Scheffersomyces_stipitis_CBS6054)
GOtermData_Spombe_anyCopyOG <- mapFromAnyCopyOGtoGO(uniprot_Spombe,anyCopyOG,'pombe',EG2GO,singleCopyOG$Schizosaccharomyces_pombe_972)
# Combine GOannotation from singel and any copy orthologs dataset
GOtermAnnotation_Kmarx <- rbind(GOtermData_Kmarx,GOtermData_Kmarx_anyCopyOG)
GOtermAnnotation_Pstip <- rbind(GOtermData_Pstip,GOtermData_Pstip_anyCopyOG)
GOtermAnnotation_Spombe <- rbind(GOtermData_Spombe,GOtermData_Spombe_anyCopyOG)
GOtermAnnotation_Scer <- EG2GO[,c(2,3,4,1)]
colnames(GOtermAnnotation_Scer) <- c('UniprotID','GOterm','GOname','GeneName')

# Map IDs to gene names
mapIDtoGeneName <- function(dataToMap,uniprotData,organism){
  geneShortList <- c()
  geneNameList <- c()
  for(id in dataToMap$UniprotID){
    if(organism == 'cerevisiae'){
      geneShort <- uniprotData$Gene.names...primary..[grep(id,uniprotData$Entry)]
      geneShortList <- c(geneShortList,geneShort)
      }else{
        geneName <- uniprotData$Gene.names...ORF..[grep(id,uniprotData$Entry)]
        geneShort <- uniprotData$Gene.names...primary..[grep(id,uniprotData$Entry)]
        geneNameList <- c(geneNameList,geneName)
        geneShortList <- c(geneShortList,geneShort)
        }
  }
  if(!is.null(geneNameList)){
    dataToMap$GeneName <- geneNameList
  }
  dataToMap$GeneShort <- geneShortList
  return(dataToMap)
}

GOtermAnnotation_Scer <- mapIDtoGeneName(GOtermAnnotation_Scer,uniprot_Scer,'cerevisiae')
GOtermAnnotation_Kmarx <- mapIDtoGeneName(GOtermAnnotation_Kmarx,uniprot_Kmarx,'marxianus')
GOtermAnnotation_Pstip <- mapIDtoGeneName(GOtermAnnotation_Pstip,uniprot_Pstip,'stipitis')
GOtermAnnotation_Spombe <- mapIDtoGeneName(GOtermAnnotation_Spombe,uniprot_Spombe,'pombe')

# write GO term annotation to files
write.table(GOtermAnnotation_Scer,'../GOterms/Sce_GOterms.txt',sep = '\t',row.names = FALSE)
write.table(GOtermAnnotation_Kmarx,'../GOterms/Kma_GOterms.txt',sep = '\t',row.names = FALSE)
write.table(GOtermAnnotation_Pstip,'../GOterms/Pst_GOterms.txt',sep = '\t',row.names = FALSE)
write.table(GOtermAnnotation_Spombe,'../GOterms/Spo_GOterms.txt',sep = '\t',row.names = FALSE)

