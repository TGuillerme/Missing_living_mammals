#Rarefaction analysis

#Loading the functions
setwd("~/PhD/Projects/Missing_living_mammals/Analysis")
source("functions.R")
load.functions(test=FALSE) #Set test=FALSE to speed up the loading
library(picante)

#Loading the list of present taxa
#source("Extracting_living_taxa.R")
load("../Data/List_of_matching_taxa.Rda")

#Fixing the number of characters for truncated matrices
extraction_table[which(extraction_table$Matrix == "GS2007-R.nex"), 3]<-88
extraction_table[which(extraction_table$Matrix == "GS2010-VJD.nex"), 3]<-157
extraction_table[which(extraction_table$Matrix == "GS2012-GM.nex"), 3]<-45
extraction_table[which(extraction_table$Matrix == "GS2012-SM.nex"), 3]<-

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

#Selecting the living taxa per data source
#Morphobank
morphobank1<-living_taxa[grep("MB", living_taxa$Matrix),]
morphobank2<-living_taxa[grep("MP", living_taxa$Matrix),]
morphobank<-rbind(morphobank1, morphobank2)

#Ross Mounce
mounce<-living_taxa[grep("MM", living_taxa$Matrix),]

#Graeme Lloyd
lloyd<-living_taxa[grep("GL", living_taxa$Matrix),]

#Google Scholar
google<-living_taxa[grep("GS", living_taxa$Matrix),]

#Combine the repositories data
repositories<-rbind(morphobank, mounce, lloyd)


#Extracting the number of living taxa from the repositories
rep<-length(unique(repositories$Taxa))
#Cumulative taxa addition per google search containing additional living taxa
#Empty vector containing the number of taxa added
taxa_added<-vector()
for (mat in 1:length(levels(as.factor(google$Matrix)))) {
    #Extracting a matrix
    taxa_list<-google$Taxa[which(google$Matrix == levels(as.factor(google$Matrix))[mat])]
    #Counting taxa that are not present in the repositories list already
    taxa_unique<-length(which(is.na(match(taxa_list, repositories$Taxa))))
    #Adding the number of taxa added to the taxa existing in the repositories
    taxa_added[mat]<-taxa_unique
}
taxa_added<-sort(taxa_added, decreasing=TRUE)
taxa_add<-cumsum(taxa_added)

#Select the repo sp and then add the google sources one by one.

x<-c(10,11,12,13,14)
barplot(taxa_add, xlab="Google scholar searchs", ylab="Number of unique living OTUs", names.arg=c(seq(1:length(taxa_add))))
