#Matching list of taxa and tree for higher clade based on Wilson Reeder's MSW.
##########################
#SYNTAX:
#<species> being a list of binomial names
#<tree> being a species trees
#<taxonomic.level> which taxonomic level, should match with the reference file (case sensitive)
#<reference> taxonomic reference list
##########################
#guillert(at)tcd.ie - 16/01/2015
##########################


#FUNCTIONS
#functions for matching one name of the matrix
#Extracting the first match

higher.clade<-function(species, tree, taxonomic.level, reference) {
    #SANITIZING
    #species
    check.class(species, "character", " must be a vector of binomial names.")
    check.length(species, 1, " must be a vector of binomial names.", errorif=TRUE)

    #tree
    check.class(tree, "phylo", " must be a 'phylo' object.")

    #reference
    check.class(reference, "data.frame", " must be a 'data.frame' object.")

    #taxonomic level (must be a column name in reference)
    check.class(taxonomic.level, "character", " must be matching with one of the column names of the given reference.")
    if(length(grep(taxonomic.level, colnames(reference)))==0){
        stop("taxonomic.level must be matching with one of the column names of the given reference.")
    } else {
        taxonomic_level<-grep(taxonomic.level, colnames(reference))
    }

    #Get the sublist of species containing only the genus names
    binom.split<-function(taxon) {
        return(strsplit(taxon, split="_")[[1]][1])
    }

    unique_genus<-unique(unlist(lapply(as.list(species), binom.split)))

    #Subset of reference containing only the genus names
    sub_ref<-reference[match(unique_genus, reference$Genus),]

    #Extract the taxonomic level elements
    pruned_species<-sub_ref[match(unique(sub_ref[,taxonomic_level]), sub_ref[,taxonomic_level]),8]

    #Remove the extra species from the tree
    binom.make<-function(taxon, tree) {
        return(grep(taxon, tree$tip.label, value=TRUE)[1])
    }
    species_to_keep<-unlist(lapply(as.list(pruned_species), binom.make, tree=tree))
    species_to_remove<-tree$tip.label[-c(match(species_to_keep, tree$tip.label))]
    pruned_tree<-drop.tip(tree, species_to_remove)
    return(list(pruned_species, pruned_tree))
    #End
}