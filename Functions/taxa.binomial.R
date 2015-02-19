#Transforms a list of taxa into a fully binomial list of taxa
##########################
#SYNTAX:
#<taxa> a vector of taxa with morphological data
#<reference> the list of taxonomical references
##########################
#guillert(at)tcd.ie - 19/02/2015
##########################

taxa.binomial<-function(taxa,Sub_reference) {
    #Binomial taxa
    taxa_binomial<-taxa[grep("_", taxa)]
    #Non Binomial taxa
    taxa_nonbinom<-taxa[grep("_", taxa, invert=TRUE)]

    #If any non_binomial taxa
    if(length(taxa_nonbinom) >= 1) {
        #Checking through the nonbinomial taxa for monospecific taxa
        taxa_monospecific<-rep(NA, length(taxa_nonbinom))
        #Looping through the list to select the monospecific taxa (if any)
        for (taxon in 1:length(taxa_nonbinom)) {    
            #Selecting the species names
            monospecific_test<-Sub_reference$Species[which(Sub_reference == taxa_nonbinom[taxon], arr.ind=TRUE)[,1]]
            #Removing the blank names
            monospecific_test<-monospecific_test[-which(monospecific_test=="")]
            #Is there only one unique species name?
            if(length(unique(monospecific_test))==1) {
                #Selecting the genus name
                gen_name<-unique(Sub_reference$Genus[which(Sub_reference == taxa_nonbinom[taxon], arr.ind=TRUE)[,1]])
                #Removing the blank names if more than one name
                if(length(gen_name) > 1) {
                    gen_name<-gen_name[-which(gen_name=="")]
                }
                taxa_monospecific[taxon]<-paste(gen_name, unique(monospecific_test), sep="_")
            } else {
                taxa_monospecific[taxon]<-NA
            }
        }
        
        #Cleaning list if exists and adding it to the existing list
        if(!is.null(taxa_monospecific)) {
            #Removing the monospecific taxa from the non-binomial names list
            taxa_nonbinom<-taxa_nonbinom[-which(!is.na(taxa_monospecific))]
            #Removing NAs
            taxa_monospecific<-taxa_monospecific[-which(is.na(taxa_monospecific))]
            #Removing multiples
            taxa_monospecific<-unique(taxa_monospecific)
            #Adding to the species list
            taxa_binomial<-unique(c(taxa_binomial, taxa_monospecific))
        }

    }
    return(taxa_binomial)
}