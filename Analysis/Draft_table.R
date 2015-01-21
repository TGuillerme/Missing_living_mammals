#Draft analyses table

#Ones the living taxa are extracted from all the matrices, count the proportion of available taxa per order (species levels / genus level / family level) + Faith's Phylogenetic Diversity
#

library(picante)

#Test (example) on primates
setwd("~/PhD/Projects/Missing_living_mammals/Analysis")
source("functions.R")
load.functions(test=FALSE) #Set test=FALSE to speed up the loading

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

#EXAMPLE OF ANALYSIS FOR CLUSTERED PRIMATES

#Only lorises and platirhines (clustered)
LorisMRCA<-getMRCA(primates_tree, c("Galago_zanzibaricus", "Loris_tardigradus"))
loris_tree<-extract.clade(primates_tree, LorisMRCA)
PlathMRCA<-getMRCA(primates_tree, c("Aotus_azarae", "Ateles_belzebuth"))
plath_tree<-extract.clade(primates_tree, PlathMRCA)
clustered_primates<-c(sample(loris_tree$tip.label, 20), sample(plath_tree$tip.label, 60))


clustered_results<-order.structure("Primates (clustered)", clustered_primates, primates_tree, WR_list)

#Empty results table
result_table_clustered<-data.frame("Order"=rep("Primates (cluster)",3), "Taxonomic level"=c("Family", "Genus", "Species"), "Number of OTUs"=rep(NA,3), "Percentage of OTUs"=rep(NA,3), "Relative PD"=rep(NA,3), "NRI"=rep(NA,3), "MPD p_value"=rep(NA,3), "NTI"=rep(NA,3), "MNTD p_value"=rep(NA,3))

#SPECIES LEVEL ANALYSIS
data_structure<-data.structure(clustered_primates, primates_tree)
result_table_clustered[3,-c(1,2)]<-community.structure(data_structure, primates_tree)

#HIGHER CLADE ANALYSIS
#Genus
genus_level_data<-higher.clade(clustered_primates, primates_tree, taxonomic.level="Genus", reference=WR_list)[[1]]
genus_level_tree<-higher.clade(primates_tree$tip.label, primates_tree, taxonomic.level="Genus", reference=WR_list)[[2]]
genus_data_structure<-data.structure(genus_level_data, genus_level_tree)
result_table_clustered[2,-c(1,2)]<-community.structure(genus_data_structure, genus_level_tree)

#Family
family_level_data<-higher.clade(clustered_primates, primates_tree, taxonomic.level="Family", reference=WR_list)[[1]]
family_level_tree<-higher.clade(primates_tree$tip.label, primates_tree, taxonomic.level="Family", reference=WR_list)[[2]]
family_data_structure<-data.structure(family_level_data, family_level_tree)
result_table_clustered[1,-c(1,2)]<-community.structure(family_data_structure, family_level_tree)

#EXAMPLE OF ANALYSIS FOR RANDOM PRIMATES

#Sub lists of "present taxa"
#Random
rand_primates<-sample(primates_list, 80)

#Empty results table
result_table_random<-data.frame("Order"=rep("Primates (random)",3), "Taxonomic level"=c("Family", "Genus", "Species"), "Number of OTUs"=rep(NA,3), "Percentage of OTUs"=rep(NA,3), "Relative PD"=rep(NA,3), "NRI"=rep(NA,3), "MPD p_value"=rep(NA,3), "NTI"=rep(NA,3), "MNTD p_value"=rep(NA,3))

#SPECIES LEVEL ANALYSIS
data_structure<-data.structure(rand_primates, primates_tree)
result_table_random[3,-c(1,2)]<-community.structure(data_structure, primates_tree)

#HIGHER CLADE ANALYSIS
#Genus
genus_level_data<-higher.clade(rand_primates, primates_tree, taxonomic.level="Genus", reference=WR_list)[[1]]
genus_level_tree<-higher.clade(primates_tree$tip.label, primates_tree, taxonomic.level="Genus", reference=WR_list)[[2]]
genus_data_structure<-data.structure(genus_level_data, genus_level_tree)
result_table_random[2,-c(1,2)]<-community.structure(genus_data_structure, genus_level_tree)

#Family
family_level_data<-higher.clade(rand_primates, primates_tree, taxonomic.level="Family", reference=WR_list)[[1]]
family_level_tree<-higher.clade(primates_tree$tip.label, primates_tree, taxonomic.level="Family", reference=WR_list)[[2]]
family_data_structure<-data.structure(family_level_data, family_level_tree)
result_table_random[1,-c(1,2)]<-community.structure(family_data_structure, family_level_tree)

#COMBINING THE RESULTS
results_table<-rbind(result_table_random, result_table_clustered)