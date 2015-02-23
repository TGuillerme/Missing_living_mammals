#Loading the functions
setwd("~/PhD/Projects/Missing_living_mammals/Analysis")
source("functions.R")
load.functions(test=FALSE) #Set test=FALSE to speed up the loading
library(picante)
library(xtable)

#Loading the tree
full_trees<-read.nexus("../Data/Trees/FritzTree.rs200k.100trees.tre")
#Selecting one random tree
set.seed(1) ; one_tree<-full_trees[[sample(1:100, 1)]]

#Loading the taxonomic reference list
#Loading the taxonomic reference list
WR_list<-read.csv("../Data/Taxon_References/WilsonReederMSW.csv", header=T, stringsAsFactors=F)
#Changing ARTIODACTYLA and CETACEA to CETARTIODACTYLA
WR_list$Order[which(WR_list$Order == "ARTIODACTYLA")]<-"CETARTIODACTYLA"
WR_list$Order[which(WR_list$Order == "CETACEA")]<-"CETARTIODACTYLA"

#Loading the list of present taxa
#source("Extracting_living_taxa.R")
load("../Data/List_of_matching_taxa.Rda")

#Fixing the number of characters for truncated matrices
extraction_table[which(extraction_table$Matrix == "GS2007-R.nex"), 3]<-88
extraction_table[which(extraction_table$Matrix == "GS2010-VJD.nex"), 3]<-157
extraction_table[which(extraction_table$Matrix == "GS2012-GM.nex"), 3]<-45
extraction_table[which(extraction_table$Matrix == "GS2012-SM.nex"), 3]<-157

#Creating sub tables with a minimum number of characters
#1 Character
extraction_table1 <- extraction_table[which(as.numeric(extraction_table$Characters) >= 1),]
#75 Characters
extraction_table75<- extraction_table[which(as.numeric(extraction_table$Characters) >= 75),]
#150 Characters
extraction_table150<-extraction_table[which(as.numeric(extraction_table$Characters) >= 150),]
#300 Characters
extraction_table300<-extraction_table[which(as.numeric(extraction_table$Characters) >= 300),]
#list
extraction_list=list(extraction_table1,extraction_table75,extraction_table150,extraction_table300)
results_names=c("results_1.Rda","results_75.Rda","results_150.Rda","results_300.Rda")

#Select the character threshold
extraction_table<-extraction_list[[1]]

#Selecting the living taxa (for a minimal number of characters)
living_taxa<-extraction_table[which(extraction_table$Living == TRUE),]
living_taxa_list<-unique(living_taxa$Taxa)

######################
#Loading the results
######################

#Loading the results
load("results_150.Rda")

############################################
#Table with the number of taxa containing morphological data
############################################

#Reformat table (digit + only four first columns)
table_to_print<-results[,c(1:4)]
#Rounding to two digits only
table_to_print[,4]<-round(as.numeric(table_to_print[,4]), digit=2)
#Select the rows with less than 25% data = >75% missing data
low_data<-which(table_to_print[,4]<25)
#make the whole table as.character
for (column in 1:ncol(table_to_print)) {
    table_to_print[,column]<-as.character(table_to_print[,column])
}
#Lower case
table_to_print[,1]<-capwords(table_to_print[,1], strict=TRUE)
#Mark the low_data rows to be highlighted in the LaTeX table (BOLD)
for (column in 1:ncol(table_to_print)) {
    table_to_print[low_data,column]<-paste('BOLD',table_to_print[low_data,column], sep="")
}
#Fixing the column names
colnames(table_to_print)<-c("Order", "Taxonomic level", "Fraction of OTUs", "Percentage of OTUs")
#Saving the table
print(xtable(table_to_print), file="../Manuscript/morpho_taxa_proportion.tex", include.rownames=FALSE, tabular.environment="longtable", floating=FALSE, sanitize.text.function=bold.cells)

############################################
#Table with the data structure for each taxa
############################################

#Selecting the right columns
metric="PD"
metric_col<-grep(metric, names(results))
table_to_print<-results[,c(1:4,metric_col)]
#Removing the NAs
table_to_print<-table_to_print[-which(is.na(table_to_print[,5])),]
#Selecting the rows with significant difference from null
signif_data<-which(table_to_print[,6]<0.05)
#make the whole table as.character
for (column in 1:ncol(table_to_print)) {
    table_to_print[,column]<-as.character(table_to_print[,column])
}
#Lower case
table_to_print[,1]<-capwords(table_to_print[,1], strict=TRUE)
#Mark the low_data rows to be highlighted in the LaTeX table (BOLD)
for (column in 1:ncol(table_to_print)) {
    table_to_print[signif_data, column]<-paste('BOLD',table_to_print[signif_data, column], sep="")
}
#Fixing the column names
colnames(table_to_print)<-c("Order", "Taxonomic level", "Fraction of OTUs", "Percentage of OTUs", "PD", "p-value")
#Saving the table
print(xtable(table_to_print), file="../Manuscript/data_structure.tex", include.rownames=FALSE, tabular.environment="longtable", floating=FALSE, sanitize.text.function=bold.cells)

############################################
#Figure (phylogeny example)
############################################

#Setting the colour scheme (data, no data)
pdf("../Manuscript/example_coverage.pdf", with=16.6, height=8)
#quartz(width = 16.6, height = 8) #A4 landscape
#quartz(width = 8.3, height = 5.8) #A5 landscape
op<-par(mfrow=c(1,2), oma=c(0,1,0,1))
plotA<-plot.results(order="CETARTIODACTYLA", taxa=living_taxa_list, col_branch=c("red", "grey"), reference=WR_list, verbose=TRUE)
text(x=(plotA$x.lim[1]-plotA$x.lim[1]*0.05),y=(plotA$y.lim[2]-plotA$y.lim[2]*0.05),"A",cex=3)
plotB<-plot.results(order="CARNIVORA", taxa=living_taxa_list, col_branch=c("red", "grey"), reference=WR_list, verbose=TRUE)
text(x=(plotB$x.lim[1]-plotB$x.lim[1]*0.05),y=(plotB$y.lim[2]-plotB$y.lim[2]*0.05),"B",cex=3)
par(op)
dev.off()
