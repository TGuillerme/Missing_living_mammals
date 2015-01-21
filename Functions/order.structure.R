#Analysis of the distribution of species per order with morphological data
##########################
#SYNTAX:
#<order> the name of the order
#<species> a vector of species with morphological data
#<tree> the full tree of the order
#<reference> the list of taxonomical references
#<runs> the number of randomization runs
#<null> the null model. See ?ses.pd for more details (default = "taxa.labels")
#<verbose> whether to be verbose or not (default=FALSE)
#<...> any optional arguments to be passed to ses.mpd and ses.mntd
##########################
#guillert(at)tcd.ie - 21/01/2015
##########################

order.structure<-function(order, species, tree, reference, runs=1000, null="taxa.labels", verbose=FALSE, ...) {
    #LIBRARIES
    require(picante)

    #SANITIZING
    #order
    check.class(order, "character", " must be a single string of characters.")
    check.length(order, 1, " must be a single string of characters.")

    #species
    check.class(species, "character", " must be a vector of characters.")
    check.length(species, 1, " must be a vector of characters.", errorif=TRUE)

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
    results<-data.frame("Order"=rep(order,3), "Taxonomic level"=c("Family", "Genus", "Species"), "Number of OTUs"=rep(NA,3), "Percentage of OTUs"=rep(NA,3), "Relative PD"=rep(NA,3), "NRI"=rep(NA,3), "MPD p_value"=rep(NA,3), "NTI"=rep(NA,3), "MNTD p_value"=rep(NA,3))

    #Species level
    data_structure<-data.structure(species, tree)
    if(verbose == TRUE) {
        message("Calculating the community structure at species level: ...", appendLF=FALSE)
    }
    results[3,-c(1,2)]<-community.structure(data_structure, tree, runs, null, ...)
    if(verbose == TRUE) {
        message("Done\n", appendLF=FALSE)
    }

    #Genus level
    genus_level_data<-higher.clade(species, tree, taxonomic.level="Genus", reference=reference)[[1]]
    genus_level_tree<-higher.clade(tree$tip.label, tree, taxonomic.level="Genus", reference=reference)[[2]]
    genus_data_structure<-data.structure(genus_level_data, genus_level_tree)
    if(verbose == TRUE) {
        message("Calculating the community structure at genus level: ...", appendLF=FALSE)
    }    
    results[2,-c(1,2)]<-community.structure(genus_data_structure, genus_level_tree, runs, null, ...)
    if(verbose == TRUE) {
        message("Done\n", appendLF=FALSE)
    }

    #Family level
    family_level_data<-higher.clade(species, tree, taxonomic.level="Family", reference=reference)[[1]]
    family_level_tree<-higher.clade(tree$tip.label, tree, taxonomic.level="Family", reference=reference)[[2]]
    family_data_structure<-data.structure(family_level_data, family_level_tree)
    if(verbose == TRUE) {
        message("Calculating the community structure at genus level: ...", appendLF=FALSE)
    }    
    results[1,-c(1,2)]<-community.structure(family_data_structure, family_level_tree, runs, null, ...)
    if(verbose == TRUE) {
        message("Done\n", appendLF=FALSE)
    }

    #Output
    return(results)

}