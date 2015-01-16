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



    return(list(pruned_species, pruned_tree))
    #End
}