#Extract all the taxa from a given order for a vector or a tree
##########################
#SYNTAX:
#<order> the name of the order
#<taxa> a vector of species with morphological data or a tree
#<reference> the list of taxonomical references
##########################
#guillert(at)tcd.ie - 21/01/2015
##########################

extract.order<-function(order, taxa, reference) {
    
    #SANTIZING
    #order
    check.class(order, "character", " must be a the name of the order to extract.")

    #reference
    check.class(reference, "data.frame", " must be a taxonomic reference list.")
    if(length(which(reference == order)) == 0) {
        stop("Order not found in provided reference list.")
    }

    #taxa
    if(class(taxa) == "character") {
        taxa_is_tree <- FALSE
    } else {
        if(class(taxa) == "phylo") {
            taxa_is_tree <- TRUE
        } else {
            stop("taxa must be either a vector of species names or a phylogenetic tree.")
        }
    }

    #SELECT ONLY THE TAXA FROM THE GIVEN ORDER
    Sub_reference<-reference[which(reference == order),]

    #If taxa is a vector
    if(taxa_is_tree == FALSE) {
        #Keep only the elements in taxa that do match with the Sub_reference
        #Empty vector with the taxa to keep
        taxa_to_keep<-rep(NA, length(taxa))
        #Empty data frame for eventual taxonomic changes
        taxa_change<-data.frame("Original"=NA, "Match"=NA)

        for(taxon in 1:length(taxa)) {
            #Is the taxon in the Sub_reference?
            test<-which(Sub_reference == taxa[taxon])

            if(length(test)!=0) {
                #The taxon is in the Sub_refence
                match_taxon<-TRUE

            } else {
                #Is it a binomial name?
                if(grep("_", taxa[taxon])==1) {
                    #Split the binomial name
                    binomial<-check.split(taxa[taxon], Sub_reference)

                    #No match
                    if(all(is.na(binomial))) {
                        match_taxon<-FALSE
                    }

                    #Full match
                    if(all(!is.na(binomial))) {
                        match_taxon<-TRUE
                    }

                    #Partial match
                    if(!is.na(binomial[1]) & is.na(binomial[2])) {
                        #Only genera name matches
                        taxon_change<-c(taxa[taxon], binomial[1])
                        message(paste(taxa[taxon], " did not match but ", binomial[1], " did.", sep=""))
                        taxa_change<-rbind(taxa_change, taxon_change)
                        }

                    }

                }

                #Replace the logical vector for which taxa to keep
                taxa_to_keep[taxon]<-match_taxon
            }
    }

    if(taxa_is_tree == TRUE) {
        #Select the list of taxa present in the tree
        full_taxa<-taxa$tip.label

        #Select the taxa present in the tree present in the order only
        taxa_to_keep<-rep(NA, length(full_taxa))

        for(taxon in 1:length(full_taxa)) {
            #Is the taxon in the Sub_reference?
            test<-which(Sub_reference == full_taxa[taxon])
            
            if(length(test)!=0) {
                #The taxon is in the Sub_refence
                match_taxon<-TRUE

            } else {
                #Is it a binomial name?
                if(grep("_", full_taxa[taxon])==1) {
                    #Split the binomial name
                    binomial<-check.split(full_taxa[taxon], Sub_reference)

                    #No match
                    if(all(is.na(binomial))) {
                        match_taxon<-FALSE
                    }

                    #Full match
                    if(all(!is.na(binomial))) {
                        match_taxon<-TRUE
                    }

                }
            }

            #Replace the logical vector for which taxa to keep
            taxa_to_keep[taxon]<-match_taxon
        }

        #Select the order MRCA containing all these taxa
        order_MRCA<-getMRCA(taxa, full_taxa[taxa_to_keep])

        #Prune the tree
        order_tree<-extract.clade(taxa, node=order_MRCA)
    }

    #Preparing the output
    if(taxa_is_tree == FALSE) {
        #If taxa_change is more than one row (first=NA)
        if(nrow(taxa_change) > 1) {
            taxa_change<-taxa_change[-1,]
            output<-list("Order"=taxa[taxa_to_keep],"Taxa changes"=taxa_change)
        } else {
            output<-taxa[taxa_to_keep]
        }
    } else {
        output<-order_tree
    }

    return(output)

#END
}