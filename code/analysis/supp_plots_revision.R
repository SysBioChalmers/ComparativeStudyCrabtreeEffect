## Load required libraries
library(ggplot2)
library(RColorBrewer)

setwd('../../data/GOterms/')
filename <- 'HXT_and_glycolytic_proteins.csv'
data <- read.csv(filename,sep = ';',stringsAsFactors = F)
data$Organism <- data$Species
data$Species <- factor(data$Species, levels = data$Species)

# Create bar plot for hexose transporters
Hxts <- ggplot(data,aes(x=Species,y=HXT)) +
  geom_bar(stat="identity",color="black",position=position_dodge()) +
  theme_classic() +
  theme(plot.title = element_text(size = 10,face = 'bold',color = 'black',hjust = 0.5),
        axis.ticks = element_line(size = 1.1,color = 'black'),
        axis.text.x = element_text(size = 8,face = 'bold',color = 'black',angle = 90),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 8,face = 'bold',color = 'black'),
        axis.title.y = element_text(size = 8,face = 'bold',color = 'black'),
        panel.border = element_rect(size = 1.1, color = 'black',fill = NA)) +
  scale_y_continuous(expand=c(0,0),limits = c(0,23)) +
  ylab("Number of hexose transporters (HXTs)")

Hxts
ggsave('Hxts.pdf',device = 'pdf',dpi = 300,width = 15,height = 15,units = 'cm')

# Prepare dataset for glycolytic proteins
Species <- c()
values <- c()
nSpecies <- length(data$Species)
for(i in 1:nSpecies){
  Species <- cbind(Species,t(c(rep(data$Organism[i],ncol(data)-4))))
  values <- rbind(values,t(data[i,3:12]))
}
rxn <- c(rep(colnames(data)[3:12],nSpecies))
data_stacked_barplot <- data.frame(t(Species),values,rxn)
row.names(data_stacked_barplot) <- NULL
colnames(data_stacked_barplot) <- c('Species','value','rxn')
data_stacked_barplot$Species <- factor(data_stacked_barplot$Species,levels = data$Organism)
data_stacked_barplot$rxn <- factor(data_stacked_barplot$rxn,levels = rev(colnames(data)[3:12]))

# Create a stacked bar plot for glycolytic proteins
glycolytic_rxns <- ggplot(data_stacked_barplot, aes(fill=rxn, y=value, x=Species)) + 
  geom_bar(stat="identity") +
  theme_classic() +
  theme(plot.title = element_text(size = 10,face = 'bold',color = 'black',hjust = 0.5),
        axis.ticks = element_line(size = 1.1,color = 'black'),
        axis.text.x = element_text(size = 8,face = 'bold',color = 'black',angle = 90),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 8,face = 'bold',color = 'black'),
        axis.title.y = element_text(size = 8,face = 'bold',color = 'black'),
        panel.border = element_rect(size = 1.1, color = 'black',fill = NA),
        legend.text = element_text(size = 8,face = 'bold',color = 'black')) +
  scale_y_continuous(expand=c(0,0),limits = c(0,22),breaks = seq(0,22,1)) +
  scale_fill_brewer(palette = 'PRGn') +
  ylab("Number of proteins")
glycolytic_rxns
ggsave('Glycolytic_proteins.pdf',device = 'pdf',dpi = 300,width = 15,height = 15,
       units = 'cm')

