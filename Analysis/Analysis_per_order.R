#Loading the functions
setwd("~/PhD/Projects/Missing_living_mammals/Analysis")
source("functions.R")
load.functions(test=FALSE) #Set test=FALSE to speed up the loading
library(picante)

#Loading the tree
full_trees<-read.nexus("../Data/Trees/FritzTree.rs200k.100trees.tre")
#Selecting one random tree
set.seed(1) ; one_tree<-full_trees[[sample(1:100, 1)]]

#Loading the taxonomic reference list
WR_list<-read.csv("../Data/Taxon_References/WilsonReederMSW.csv", header=T, stringsAsFactors=F)

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
#50 Characters
extraction_table50<-extraction_table[which(as.numeric(extraction_table$Characters) >= 50),]
#100 Characters
extraction_table100<-extraction_table[which(as.numeric(extraction_table$Characters) >= 100),]
#150 Characters
extraction_table150<-extraction_table[which(as.numeric(extraction_table$Characters) >= 150),]
#200 Characters
extraction_table200<-extraction_table[which(as.numeric(extraction_table$Characters) >= 200),]

#Select the number of characters
extraction_table<-extraction_table200

#Selecting the living taxa (for a minimal number of characters)
living_taxa<-extraction_table[which(extraction_table$Living == TRUE),]
living_taxa_list<-unique(living_taxa$Taxa)

#################################
#
#   TO DO: Sort by number of characters?
#
#################################



#List of orders
orders<-unique(WR_list$Order)


#Calculate the structure of the data for the first order to initiate the loop

#First order
#Extract orders
tree<-extract.order(orders[1], one_tree, WR_list, verbose=TRUE)
taxa<-extract.order(orders[1], living_taxa_list, WR_list, verbose=TRUE)
#Replace taxa by a vector if there was any taxonomic correction
if(class(taxa) != "character") {
    taxa<-taxa$Taxa
}

#Calculate the data structure in the order
results<-order.structure(orders[1], taxa, tree, WR_list, verbose=TRUE)



#Calculate the structure of the data for the other orders through a loop


#Loop through the rest
try(
for (order in 2:length(orders)) {
    #verbose
    message(paste("\n",orders[order], " analysis.\n", sep=""), appendLF=FALSE)

    #Extract orders
    tree<-extract.order(orders[order], one_tree, WR_list, verbose=TRUE)
    #If tree is null (monospecific order!) then abort the rest
    if(is.null(tree)) {
        results_tmp<-data.frame("Order"=rep(orders[order],3), "Taxonomic level"=c("Species", "Genus", "Family"), "Number of OTUs"=rep("1/1"), "Percentage of OTUs"=rep(1), "Relative PD"=rep(NA,3), "NRI"=rep(NA,3), "MPD p_value"=rep(NA,3), "NTI"=rep(NA,3), "MNTD p_value"=rep(NA,3))
    } else {
        taxa<-extract.order(orders[order], living_taxa_list, WR_list, verbose=TRUE)
        #Replace taxa by a vector if there was any taxonomic correction
        if(class(taxa) != "character") {
            taxa<-taxa$Taxa
        }

        #Calculate the data structure in the order
        try(
        results_tmp<-order.structure(orders[order], taxa, tree, WR_list, verbose=TRUE, runs=1000)
        ,silent=TRUE)
    }
    #Combine the temporary results
    results<-rbind(results, results_tmp)
}
,silent=TRUE)

write(results, file="results_200.Rda")