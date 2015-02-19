#Analysis of the distribution of taxa per order with morphological data
##########################
#SYNTAX:
#<order> the name of the order
#<taxa> a vector of taxa with morphological data
#<tree> the full tree of the order
#<reference> the list of taxonomical references
#<metric> a list of metrics
#<runs> the number of randomization runs
#<null> the null model. See ?ses.pd for more details (default = "taxa.labels")
#<verbose> whether to be verbose or not (default=FALSE)
#<...> any optional arguments to be passed to ses.mpd and ses.mntd
##########################
#guillert(at)tcd.ie - 04/02/2015
##########################

order.structure<-function(order, taxa, tree, reference, metric=c("PD", "NRI", "NTI"), runs=1000, null="taxa.labels", verbose=FALSE, ...) {
    #LIBRARIES
    require(picante)

    #SANITIZING
    #order
    check.class(order, "character", " must be a single string of characters.")
    check.length(order, 1, " must be a single string of characters.")

    #taxa
    check.class(taxa, "character", " must be a vector of characters.")
    check.length(taxa, 1, " must be a vector of characters.", errorif=TRUE)

    #tree
    check.class(tree, "phylo", " must be a phylogenetic tree.")

    #reference
    check.class(reference, "data.frame", " must be a 'data.frame' object.")

    #runs
    check.class(runs, "numeric", " must be a single numerical value.")
    check.length(runs, 1, " must be a single numerical value.")

    #null
    check.class(null, "character", " must be a null model from the ses.mpd and ses.mntd functions.")
    check.length(null, 1, " must be a null model from the ses.mpd and ses.mntd functions.")

    #ANALYSING THE DISTRIBUTION OF SPECIES PER ORDER WITH MORPHOLOGICAL DATA
    #Empty table
    results<-data.frame("Order"=rep(order,3), "Taxonomic level"=c("Family", "Genus", "Species"), "Number of OTUs"=rep(NA,3), "Percentage of OTUs"=rep(NA,3))
    #Adding two columns per metric
    if(any(metric == "PD")) {
        results$PD<-rep(NA,3)
        results$PD_p<-rep(NA,3)
    }

    if(any(metric == "NRI")) {
        results$NRI<-rep(NA,3)
        results$NRI_p<-rep(NA,3)
    }

   if(any(metric == "NTI")) {
        results$NTI<-rep(NA,3)
        results$NTI_p<-rep(NA,3)
    }

    #Sort the list of taxa per species/genera/family
    #Reference list for the order
    Sub_reference<-subset(reference[which(reference == order),])

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

    #Test if the order is monogeneric or monofamilial
    #Monogeneric?
    mono_test<-unique(Sub_reference$Genus) ; mono_test<-mono_test[-which(mono_test=="")]
    if(length(mono_test) < 2) {
        monogeneric<-TRUE
        if(verbose == TRUE) {
            message("Order is monogeneric.\n", appendLF=FALSE)
        }
    } else {
        monogeneric<-FALSE
    }

    #Monofamilial?
    mono_test<-unique(Sub_reference$Family) ; mono_test<-mono_test[-which(mono_test=="")]
    if(length(mono_test) < 2) {
        monofamilial<-TRUE
        if(verbose == TRUE) {
            message("Order is monofamilial.\n", appendLF=FALSE)
        }
    } else {
        monofamilial<-FALSE
    }

    #If there is only one sampled taxa, throw NAs
    if(length(taxa_binomial) == 1) {
        if(verbose == TRUE) {
            message("Only one taxa has been sampled.\n", appendLF=FALSE)
        }
        
        #Calculate the number of OTUS at the species level
        OTUs<-paste(1, Ntip(tree), sep="/")
        #OTU presence / absence ratio at the species level
        OTUs_ratio<-1/Ntip(tree)*100
        results[3,c(3,4)]<-c(OTUs, round(OTUs_ratio, digit=2))

        if(monogeneric == FALSE) { 
            #Calculate the number of OTUS at the genus level
            OTUs<-paste(1, Ntip(higher.clade(tree$tip.label, tree, taxonomic.level="Genus", reference=Sub_reference)[[2]]), sep="/")
            #OTU presence / absence ratio at the genus level
            OTUs_ratio<-1/Ntip(higher.clade(tree$tip.label, tree, taxonomic.level="Genus", reference=Sub_reference)[[2]])*100
            results[2,c(3,4)]<-c(OTUs, round(OTUs_ratio, digit=2))
        } else {
            results[2,c(3,4)]<-c("1/1", 100)
        }

        if(monofamilial == FALSE) {
            #Calculate the number of OTUS at the family level
            OTUs<-paste(1, Ntip(higher.clade(tree$tip.label, tree, taxonomic.level="Family", reference=Sub_reference)[[2]]), sep="/")
            #OTU presence / absence ratio at the family level
            OTUs_ratio<-1/Ntip(higher.clade(tree$tip.label, tree, taxonomic.level="Family", reference=Sub_reference)[[2]])*100
            results[1,c(3,4)]<-c(OTUs, round(OTUs_ratio, digit=2))
        } else {
            results[1,c(3,4)]<-c("1/1", 100)
        }

    } else {

        #CALCULATE THE DATA STRUCTURE

        #Species level

        #data structure
        data_structure<-data.structure(taxa_binomial, tree)
        if(verbose == TRUE) {
            message("Calculating the community structure at species level: ...", appendLF=FALSE)
        }
        #Don't calculate the data structure if all taxa are sampled
        if(all(data_structure[1,] == 1)) {
            results[3, c(3,4)]<-community.structure(data_structure, tree, metric="PD", runs=2, null)[1:2]
        } else {
            #Calculate the community structure for the different metrics
            results[3,-c(1,2)]<-community.structure(data_structure, tree, metric, runs, null, ...)
        }
        if(verbose == TRUE) {
            message("Done\n", appendLF=FALSE)
        } 

        #Genus level

        #Check if the Order is not monogeneric
        if(monogeneric == TRUE) {
            results[2,c(3,4)]<-c("1/1", 100)
        } else {
            #Tree
            genus_level_tree<-higher.clade(tree$tip.label, tree, taxonomic.level="Genus", reference=Sub_reference)[[2]]
            #Data for binomial names
            genus_level_data<-higher.clade(taxa_binomial, tree, taxonomic.level="Genus", reference=Sub_reference)[[1]]
            #Check through the data if their are taxa that are present in the nonbionmial data set
            add_nonbinom<-taxa_nonbinom[which(is.na(match(taxa_nonbinom, genus_level_data)))]
            #Adding the new genera
            genus_level_data<-c(genus_level_data, add_nonbinom[-which(is.na(match(add_nonbinom,reference$Genus)))])

            #Data structure
            genus_data_structure<-data.structure(genus_level_data, genus_level_tree)
            if(verbose == TRUE) {
                message("Calculating the community structure at genus level: ...", appendLF=FALSE)
            }
            #Don't calculate the data structure if all taxa are sampled
            if(all(genus_data_structure[1,] == 1)) {
                results[2, c(3,4)]<-community.structure(genus_data_structure, genus_level_tree, metric="PD", runs=2, null)[1:2]
            } else {
                #Calculate the community structure for the different metrics
                results[2,-c(1,2)]<-community.structure(genus_data_structure, genus_level_tree, metric, runs, null, ...)
            }
            if(verbose == TRUE) {
                message("Done\n", appendLF=FALSE)
            }
        }

        #Family level

        #Check if the Order is not monofamilial
        if(monofamilial == TRUE) {
            results[1,c(3,4)]<-c("1/1", 100)
        } else {    
            #Tree
            family_level_tree<-higher.clade(tree$tip.label, tree, taxonomic.level="Family", reference=reference)[[2]]
            #Data
            family_level_data<-higher.clade(taxa_binomial, tree, taxonomic.level="Family", reference=reference)[[1]]
            family_data_structure<-data.structure(family_level_data, family_level_tree)
            if(verbose == TRUE) {
                message("Calculating the community structure at family level: ...", appendLF=FALSE)
            }    
            #Don't calculate the data structure if all taxa are sampled
            if(all(family_data_structure[1,] == 1)) {
                results[1, c(3,4)]<-community.structure(family_data_structure, family_level_tree, metric="PD", runs=2, null)[1:2]
            } else {
                #Calculate the community structure for the different metrics
                results[1,-c(1,2)]<-community.structure(family_data_structure, family_level_tree, metric, runs, null, ...)
            }
            if(verbose == TRUE) {
                message("Done.\n", appendLF=FALSE)
            }
        }
    
    }
    #Output
    return(results)
}