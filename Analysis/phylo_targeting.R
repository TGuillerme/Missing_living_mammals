#Function for preparing the data for the phylo targeting analysis on http://phylotargeting.fas.harvard.edu/Run_step1.php

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
extraction_table1<-extraction_table[which(as.numeric(extraction_table$Characters) >= 1),]
#75 Characters
extraction_table75<-extraction_table[which(as.numeric(extraction_table$Characters) >= 75),]
#94 Characters
extraction_table100<-extraction_table[which(as.numeric(extraction_table$Characters) >= 100),]
#138 Characters
extraction_table125<-extraction_table[which(as.numeric(extraction_table$Characters) >= 125),]
#150 Characters
extraction_table150<-extraction_table[which(as.numeric(extraction_table$Characters) >= 150),]
#300 Characters
extraction_table300<-extraction_table[which(as.numeric(extraction_table$Characters) >= 300),]
#list
extraction_list=list(extraction_table1,extraction_table75, extraction_table100, extraction_table125, extraction_table150,extraction_table300)
results_names=c("results_1.Rda","results_75.Rda","results_100.Rda","results_125.Rda","results_150.Rda","results_300.Rda")

#Select the character threshold
extraction_table<-extraction_list[[3]]

#Selecting the living taxa (for a minimal number of characters)
living_taxa<-extraction_table[which(extraction_table$Living == TRUE),]
living_taxa_list<-unique(living_taxa$Taxa)

######################
#Loading the results
######################

#Loading the results
load("results_100.Rda")

######################
#Building the nexus file
######################

#Improving data coverage for the Carnivora
order="CARNIVORA"
#Extract the tree and the taxa
tree<-extract.order(order, one_tree, WR_list, verbose=TRUE)
taxa<-extract.order(order, living_taxa_list, WR_list, verbose=TRUE)
#Make the taxa binomial
Sub_reference<-subset(WR_list[which(WR_list == order),])
taxa_binomial_non<-taxa.binomial(taxa, Sub_reference)
taxa_binomial<-taxa_binomial_non[[1]]
taxa_binomial<-taxa_binomial$Taxa
taxa_nonbinom<-taxa_binomial_non[[2]]

#Species level
structure<-data.structure(taxa_binomial, tree)
structure<-t(structure) #transpose
structure[,1]<-paste(as.numeric(structure[,1]), as.numeric(structure[,2]), sep="")
structure<-structure[,-2]
#Create the nexus file (tree only)
write.nexus(tree, file="../Data/PhyloTargeting-CARNIVORA-sp.nex")

#Create the dummy csv data file
write.table(structure, file="../Data/PhyloTargeting-CARNIVORA-sp.tmp", sep="\t", quote=FALSE, col.names=FALSE)

#Add the content of the tmp file to the nexus one with the proper header and bracketing (SHELL)
system('echo "" >> ../Data/PhyloTargeting-CARNIVORA-sp.nex')
system('echo "BEGIN CHARACTERS;" >> ../Data/PhyloTargeting-CARNIVORA-sp.nex')
system('echo "    TITLE data_structure;" >> ../Data/PhyloTargeting-CARNIVORA-sp.nex')
system('echo "    DIMENSIONS NCHAR=2;" >> ../Data/PhyloTargeting-CARNIVORA-sp.nex')
system('echo "    FORMAT DATATYPE = STANDARD SYMBOLS = \\"0 1\\";" >> ../Data/PhyloTargeting-CARNIVORA-sp.nex')
system('echo "    CHARSTATELABELS" >> ../Data/PhyloTargeting-CARNIVORA-sp.nex')
system('echo "         1 sampled, 2 complete;" >> ../Data/PhyloTargeting-CARNIVORA-sp.nex')
system('echo "    MATRIX" >> ../Data/PhyloTargeting-CARNIVORA-sp.nex')
system('cat ../Data/PhyloTargeting-CARNIVORA-sp.tmp >> ../Data/PhyloTargeting-CARNIVORA-sp.nex')
system('echo "    ;" >> ../Data/PhyloTargeting-CARNIVORA-sp.nex')
system('echo "END;" >> ../Data/PhyloTargeting-CARNIVORA-sp.nex')
system('rm ../Data/PhyloTargeting-CARNIVORA-sp.tmp')

