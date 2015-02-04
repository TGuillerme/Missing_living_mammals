#Extract all the taxa from a given order for a vector or a tree
##########################
#SYNTAX:
#<order> the name of the order
#<taxa> a vector of species with morphological data or a tree
#<reference> the list of taxonomical references
#<verbose> whether to be verbose or not
##########################
#guillert(at)tcd.ie - 04/02/2015
##########################

extract.order<-function(order, taxa, reference, verbose=FALSE) {
    
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

    #verbose
    check.class(verbose, "logical", " must be logical.")


    #SELECT ONLY THE TAXA FROM THE GIVEN ORDER
    Sub_reference<-reference[which(reference == order),]


    #If the order is monospecific, abort
    mono_test<-unique(Sub_reference$Species) ; mono_test<-mono_test[-which(mono_test=="")]
    if(length(mono_test) < 2) {
        message("Order is monospecific!")
        output<-NULL
    } else {

        #If taxa is a vector
        if(taxa_is_tree == FALSE) {
            #Keep only the elements in taxa that do match with the Sub_reference
            #Empty vector with the taxa to keep
            taxa_to_keep<-rep(NA, length(taxa))
            #Empty data frame for eventual taxonomic changes
            taxa_change<-data.frame("Original"=NA, "Match"=NA)

            #Verbose
            if(verbose == TRUE) {
                message("Scanning the reference list:", appendLF=FALSE)
            }        


            for(taxon in 1:length(taxa)) {
                #Is the taxon in the Sub_reference?
                test<-which(Sub_reference == taxa[taxon])

                if(length(test)!=0) {
                    #The taxon is in the Sub_refence
                    match_taxon<-TRUE

                } else {
                    #Is it a binomial name?
                    if(length(grep("_", taxa[taxon]))==1) {
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
                            warning(paste("\n", taxa[taxon], " did not match but ", binomial[1], " did.", sep=""))
                            taxa_change<-rbind(taxa_change, taxon_change)
                            }

                        } else {
                            #Name is not matching nor binomial
                            match_taxon<-FALSE
                        }

                    }

                    #Replace the logical vector for which taxa to keep
                    taxa_to_keep[taxon]<-match_taxon

                    #Verbose
                    if(verbose == TRUE) {
                        message(".", appendLF=FALSE)
                    }
                }
                #Verbose
                if(verbose == TRUE) {
                    message("Done.\n", appendLF=FALSE)
                }
        }

        if(taxa_is_tree == TRUE) {
            #Select the list of taxa present in the tree
            full_taxa<-taxa$tip.label

            #Select the taxa present in the tree present in the order only
            taxa_to_keep<-rep(NA, length(full_taxa))

            #Verbose
            if(verbose == TRUE) {
                message("Scanning the reference list:", appendLF=FALSE)
            }   

            for(taxon in 1:length(full_taxa)) {
                #Is the taxon in the Sub_reference?
                test<-which(Sub_reference == full_taxa[taxon])
                
                if(length(test)!=0) {
                    #The taxon is in the Sub_refence
                    match_taxon<-TRUE

                } else {
                    #Is it a binomial name?
                    if(length(grep("_", full_taxa[taxon]))==1) {
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

                #Verbose
                if(verbose == TRUE) {
                    message(".", appendLF=FALSE)
                }

            }
        
            #Verbose
            if(verbose == TRUE) {
                message("Done.\n", appendLF=FALSE)
            }

            #Select the order MRCA containing all these taxa
            order_MRCA<-getMRCA(taxa, full_taxa[taxa_to_keep])

            #Prune the tree
            order_tree<-extract.clade(taxa, node=order_MRCA)

        }

        #Preparing the output
        if(taxa_is_tree == FALSE) {
            #Kept taxa
            taxa_kept<-taxa[taxa_to_keep]

            #If taxa_change is more than one row (first=NA)
            if(nrow(taxa_change) > 1) {
                taxa_change<-taxa_change[-1,]

                #Replacing the original names in the table by the changes
                for (taxon in 1:nrow(taxa_change)) {
                    taxa_kept[match(taxa_change$Original[taxon], taxa_kept)]<-taxa_change$Match[taxon]
                }

                output<-list("Taxa"=taxa_kept,"Changes"=taxa_change)
            } else {
                output<-taxa_kept
            }
        } else {
            output<-order_tree
        }
    }
    return(output)

#END
}