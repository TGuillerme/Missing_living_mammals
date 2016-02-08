#Loading the functions

source("functions.R")
load.functions(test=FALSE) #Set test=FALSE to speed up the loading
library(picante)

#Loading the tree
mam_tree<-read.nexus("../Data/Trees/mammalST_MSW05_best_chrono.tre")

#Loading the taxonomic reference list
WR_list<-read.csv("../Data/Taxon_References/WilsonReederMSW.csv", header=T, stringsAsFactors=F)
#Changing ARTIODACTYLA and CETACEA to CETARTIODACTYLA
WR_list$Order[which(WR_list$Order == "ARTIODACTYLA")]<-"CETARTIODACTYLA"
WR_list$Order[which(WR_list$Order == "CETACEA")]<-"CETARTIODACTYLA"

#Loading the list of present taxa
#source("Extracting_living_taxa.R")
load("../Data/List_of_matching_taxa/List_of_matching_taxa.Rda")

#Fixing the number of characters for truncated matrices
extraction_table[which(extraction_table$Matrix == "GS2007-R.nex"), 3]<-88
extraction_table[which(extraction_table$Matrix == "GS2010-VJD.nex"), 3]<-157
extraction_table[which(extraction_table$Matrix == "GS2012-GM.nex"), 3]<-45
extraction_table[which(extraction_table$Matrix == "GS2012-SM.nex"), 3]<-157
extraction_table[which(extraction_table$Matrix == "GS1997-V.nex"), 3]<-40
extraction_table[which(extraction_table$Matrix == "GS2013-MS.nex"), 3]<-48
extraction_table[which(extraction_table$Matrix == "GS2014-MG.nex"), 3]<-40
extraction_table[which(extraction_table$Matrix == "GS2014-MaL.nex"), 3]<-98
extraction_table[which(extraction_table$Matrix == "GS2015-RR.nex"), 3]<-149
extraction_table[which(extraction_table$Matrix == "GS2014-B.nex"), 3]<-27
extraction_table[which(extraction_table$Matrix == "GS2014-Y.nex"), 3]<-42

#Creating sub tables with a minimum number of characters
#1 Character
extraction_table1<-extraction_table[which(as.numeric(extraction_table$Characters) >= 1),]

#100 Characters
extraction_table100<-extraction_table[which(as.numeric(extraction_table$Characters) >= 100),]

#Combine the sub tables
extraction_list=list(extraction_table1)#, extraction_table100)
results_names=c("results_1.Rda")#,"results_100.Rda")

#Print the analysis by series of characters length

#for (character_threshold in 1:length(extraction_list)) {
for (character_threshold in 1:1) {
    #Select the number of characters
    extraction_table<-extraction_list[[character_threshold]]

    #Selecting the living taxa (for a minimal number of characters)
    living_taxa<-extraction_table[which(extraction_table$Living == TRUE),]
    living_taxa_list<-unique(living_taxa$Taxa)

    #List of orders
    orders<-unique(WR_list$Order)


    #Calculate the structure of the data for the first order to initiate the loop

    #Analysis for all mammals
    message(paste("\n","Mammalia (Class)", " analysis.\n", sep=""), appendLF=FALSE)
    tree<-extract.order("Mammalia (Class)", mam_tree, WR_list, verbose=TRUE)
    taxa<-extract.order("Mammalia (Class)", living_taxa_list, WR_list, verbose=TRUE)

    if(class(taxa) != "character") {
        taxa<-taxa$Taxa
    }

    #Calculate the data structure in the order
    results<-order.structure("Mammalia (Class)", taxa, tree, WR_list, metric=c("NTI", "NRI"), verbose=TRUE)

    #Calculate the structure of the data for the other orders through a loop

    #Loop through the rest
    for (order in 1:length(orders)) {
        #verbose
        message(paste("\n",orders[order], " analysis.\n", sep=""), appendLF=FALSE)

        #Extract orders
        tree<-extract.order(orders[order], mam_tree, WR_list, verbose=TRUE)
        #If tree is null (monospecific order!) then abort the rest
        if(is.null(tree)) {
            results_tmp<-matrix(nrow=3, ncol=ncol(results), data=NA)
            results_tmp[,1]<-rep(orders[order],3)
            results_tmp[,2]<-c("family", "genus", "species")
            results_tmp[,3]<-rep("1/1",3)
            results_tmp[,4]<-rep("100",3)
            results_tmp<-as.data.frame(results_tmp)
            names(results_tmp)<-names(results)
        } else {
            taxa<-extract.order(orders[order], living_taxa_list, WR_list, verbose=TRUE)
            #Replace taxa by a vector if there was any taxonomic correction
            if(class(taxa) != "character") {
                taxa<-taxa$Taxa
            }

            #Calculate the data structure in the order
            results_tmp<-order.structure(orders[order], taxa, tree, WR_list, verbose=TRUE, runs=1000)
        }
        #Combine the temporary results
        results<-rbind(results, results_tmp)
    }

    save(results, file=results_names[character_threshold])
}