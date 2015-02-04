#Loading the functions
setwd("~/PhD/Projects/Missing_living_mammals/Analysis")
source("functions.R")
load.functions(test=FALSE) #Set test=FALSE to speed up the loading
library(picante)

#Loading the tree
full_trees<-read.nexus("../Data/Trees/FritzTree.rs200k.100trees.tre")
set.seed(1) ; one_tree<-full_trees[[sample(1:100, 1)]]

#Loading the taxonomic reference list
WR_list<-read.csv("../Data/Taxon_References/WilsonReederMSW.csv", header=T, stringsAsFactors=F)

#Loading the list of present taxa
#source("Extracting_living_taxa.R")
load("../Data/")

#Sort by number of characters?

#Extract orders
tree<-extract.order("PRIMATES", one_tree, WR_list)
taxa<-extract.order("PRIMATES", living_taxa, WR_list)
#


#Calculate the data structure in the order
clustered_results<-order.structure("Primates (clustered)", clustered_primates, primates_tree, WR_list, verbose=TRUE)
