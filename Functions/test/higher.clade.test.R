#Test (example) on primates
setwd("~/PhD/Projects/Missing_living_mammals/Analysis")

#Read Fritz tree
full_trees<-read.nexus("../Data/Trees/FritzTree.rs200k.100trees.tre")
set.seed(1) ; one_tree<-full_trees[[sample(1:100, 1)]]

primatesMRCA<-getMRCA(one_tree, c("Homo_sapiens", "Loris_tardigradus"))
primates_tree<-extract.clade(one_tree, node=primatesMRCA)
#list of all taxa
primates_list<-primates_tree$tip.label

#Fix data loading
setwd("~/PhD/Projects/Missing_living_mammals/Functions")

#Loading the functions
sourceDir <- function(path, trace = TRUE, ...) {
    for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
    if(trace) cat(nm,":")
        source(file.path(path, nm), ...)
        if(trace) cat("\n")
    }
}
sourceDir(".")

#Testing
library(testthat)
setwd("test/")

WR_list<-read.csv("../../Data/Taxon_References/WilsonReederMSW.csv", header=T, stringsAsFactors=F)


#Data
reference<-WR_list
tree<-primates_tree
species<-primates_list
taxonomic.level="Genus"