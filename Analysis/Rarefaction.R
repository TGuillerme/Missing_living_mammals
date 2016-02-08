#Rarefaction analysis

#Loading the functions
source("functions.R")
load.functions(test=FALSE) #Set test=FALSE to speed up the loading
library(picante)

#Loading the list of present taxa
#source("Extracting_living_taxa.R")
load("../Data/List_of_matching_taxa/List_of_matching_taxa.Rda")

#Fixing the number of characters for truncated matrices
extraction_table[which(extraction_table$Matrix == "GS2007-R.nex"), 3]<-88
extraction_table[which(extraction_table$Matrix == "GS2010-VJD.nex"), 3]<-157
extraction_table[which(extraction_table$Matrix == "GS2012-GM.nex"), 3]<-45
extraction_table[which(extraction_table$Matrix == "GS2012-SM.nex"), 3]<-157

#Select the living taxa within the repositories sources (morphobank, Graeme Lloyd and Ross Mounce)
#Living taxa
living_taxa<-extraction_table[which(extraction_table$Living == TRUE),]

#Morphobank
morphobank1<-living_taxa[grep("MB", living_taxa$Matrix),]
morphobank2<-living_taxa[grep("MP", living_taxa$Matrix),]
morphobank<-rbind(morphobank1, morphobank2)

#Ross Mounce
mounce<-living_taxa[grep("MM", living_taxa$Matrix),]

#Graeme Lloyd
lloyd<-living_taxa[grep("GL", living_taxa$Matrix),]

#Combine the repositories data
repositories<-rbind(morphobank, mounce, lloyd)
#List of living taxa from the repositories
repo_taxa<-unique(repositories$Taxa)


#Loading the google search results
Google_search<-read.csv("../Data/Search/GSresults.csv")
#Setting the matrix column as character
Google_search$Counts<-as.character(Google_search$Matrix)

#Replacing the non google results
Google_search[grep("GS", Google_search$Counts, invert=TRUE),4]<-0
#Select the GS matrices
GS_matrices<-which(Google_search$Counts != 0)
#Within each matrix, select the number of extra taxa

#Selecting the living taxa in each Google scholar matrices
taxa_present<-list()
for (mat in 1:length(GS_matrices)) {
    #Selecting the living taxa from the matrix
    taxa_present[[mat]]<-living_taxa$Taxa[grep(Google_search$Counts[GS_matrices[mat]], living_taxa$Matrix)]
    names(taxa_present)[mat]<-Google_search$Counts[GS_matrices[mat]]
}

#Counting the number of living taxa per matrix
number_of_living<-unlist(lapply(taxa_present, length))

#Making a sub_table with the results
sub_results<-data.frame("Matrix"=Google_search$Matrix[GS_matrices], "Google_search"=GS_matrices, "N_Living"=number_of_living, "Extra_taxa"=rep(0, length(GS_matrices)))
sub_results<-sub_results[order(sub_results$N_Living, decreasing=TRUE),]

#Per matrix, count the number of extra taxa that are not present in the repository AND in the former matrices (additive count)
for (mat in 1:nrow(sub_results)) {
    #Select non matching taxa
    new_taxa<-taxa_present[[grep(sub_results$Matrix[mat], names(taxa_present))[1]]][which(is.na(match(taxa_present[[grep(sub_results$Matrix[mat], names(taxa_present))[1]]], repo_taxa)))]

    #WARNING: this selects anything regardless typos or taxonomic level. Include a check?

    #Counting the new taxa
    sub_results[mat,4]<-length(new_taxa)

    #Adding these new taxa to the repositories list
    if (length(new_taxa) > 0) {
        repo_taxa<-c(repo_taxa, new_taxa)
    }
}

#Replacing the Extra taxa in the Google_search results
for (mat in 1:nrow(sub_results)) {
    Google_search$Counts[sub_results$Google_search[mat]]<-sub_results$Extra_taxa[mat]
}

Additional_data<-sort(as.numeric(Google_search$Counts), decreasing=TRUE)
pdf("../Manuscript/Supplementary/Supp_figure_google_searches.pdf")
op<-par(bty="l")
plot(cumsum(Additional_data), type="l", ylab="Number of additional OTUs", xlab="Google Scholar papers")
par(op)
dev.off()