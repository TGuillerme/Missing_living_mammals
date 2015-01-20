#Draft analyses table

#Ones the living taxa are extracted from all the matrices, count the proportion of available taxa per order (species levels / genus level / family level) + Faith's Phylogenetic Diversity
#

library(picante)

#Test (example) on primates
setwd("~/PhD/Projects/Missing_living_mammals/Analysis")
source("functions.R")
load.functions(test=TRUE) #Set test=FALSE to speed up the loading

#Read Fritz tree
full_trees<-read.nexus("../Data/Trees/FritzTree.rs200k.100trees.tre")
set.seed(1) ; one_tree<-full_trees[[sample(1:100, 1)]]
#Taxonomic reference
WR_list<-read.csv("../Data/Taxon_References/WilsonReederMSW.csv", header=T, stringsAsFactors=F)

#Isolate only the primates
primatesMRCA<-getMRCA(one_tree, c("Homo_sapiens", "Loris_tardigradus"))
primates_tree<-extract.clade(one_tree, node=primatesMRCA)
#list of all taxa
primates_list<-primates_tree$tip.label

#Sub lists of "present taxa"
#Random
rand_primates<-sample(primates_list, 80)

#Only lorises and platirhines (clustered)
LorisMRCA<-getMRCA(primates_tree, c("Galago_zanzibaricus", "Loris_tardigradus"))
loris_tree<-extract.clade(primates_tree, LorisMRCA)
PlathMRCA<-getMRCA(primates_tree, c("Aotus_azarae", "Ateles_belzebuth"))
plath_tree<-extract.clade(primates_tree, PlathMRCA)
clustered_primates<-c(sample(loris_tree$tip.label, 20), sample(plath_tree$tip.label, 60))

#Test the different metrics (random vs. clustered)
#Creating the data set with presence / absence
#data_structure<-matrix(ncol=length(primates_list), nrow=2, data=0)
#colnames(data_structure)<-primates_list
#rownames(data_structure)<-c("random", "clustered")
#data_structure[1,match(rand_primates, primates_list)]<-1
#data_structure[2,match(clustered_primates, primates_list)]<-1

#Visualisation
#par(mfrow = c(1, 2))
#for (i in row.names(data_structure)) {
#    plot(primates_tree, show.tip.label = FALSE, main = i)
#    tiplabels(tip = which(primates_tree$tip.label %in% names(which(data_structure[i, ] > 0))), pch = 19, cex = 1)
#}



#Empty results table
result_table<-data.frame("Order"=rep("Primates",3), "Taxonomic level"=c("Family", "Genus", "Species"), "Number of OTUs"=rep(NA,3), "Percentage of OTUs"=rep(NA,3), "MPD(NRI)"=rep(NA,3), "MPD(NullDist)"=rep(NA,3), "p_value"=rep(NA,3), "MTD(NRI)"=rep(NA,3), "MTD(NullDist)"=rep(NA,3), "p_value"=rep(NA,3))

#SPECIES LEVEL ANALYSIS
data_structure<-data.structure(clustered_primates, primates_tree) #??? CHECK IF THE CALCULATIONS ARE CORRECT
result_table[3,-c(1,2)]<-ses.fun(data_structure, primates_tree)

#HIGHER CLADE ANALYSIS
#Genus
genus_level_data<-higher.clade(clustered_primates, primates_tree, taxonomic.level="Genus", reference=WR_list)[[1]]
genus_level_tree<-higher.clade(primates_tree$tip.label, primates_tree, taxonomic.level="Genus", reference=WR_list)[[2]]
genus_data_structure<-data.structure(genus_level_data, genus_level_tree)
result_table[2,-c(1,2)]<-ses.fun(genus_data_structure, genus_level_tree)

#Family
family_level_data<-higher.clade(clustered_primates, primates_tree, taxonomic.level="Family", reference=WR_list)[[1]]
family_level_tree<-higher.clade(primates_tree$tip.label, primates_tree, taxonomic.level="Family", reference=WR_list)[[2]]
family_data_structure<-data.structure(family_level_data, family_level_tree)
result_table[1,-c(1,2)]<-ses.fun(family_data_structure, family_level_tree)

