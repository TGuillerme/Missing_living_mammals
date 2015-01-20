#Making a data table with the structure of presences/absences of species for further community structure analysis using picante
##########################
#SYNTAX:
#<species> being a list of species of interest
#<tree> being a species trees
#<plot> logical, whether to plot the data structure or not
#<...> any optional arguments to be passed to plot
##########################
#guillert(at)tcd.ie - 20/01/2015
##########################

data.structure<-function(species, tree, plot=FALSE, ...){
    #SANTIZING
    #species
    check.class(species, "character", " must be a vector of a subset of taxa names present in the tree.")

    #tree
    check.class(tree, "phylo", " must be a phylogenetic tree.")

    #plot
    check.class(plot, "logical", " must be logical.")

    #MAKING THE DATA STRUCTURE TABLE
    #list of all taxa
    species_list<-tree$tip.label

    #Creating the data structure tabke
    data_structure<-matrix(ncol=length(species_list), nrow=1, data=0)

    #column names
    colnames(data_structure)<-species_list

    #filling the presence absences
    data_structure[1, match(species, species_list)]<-1

    #Plotting (optional)
    if(plot == TRUE) {
        plot(tree, show.tip.label = FALSE, ...)
        tiplabels(tip = which(primates_tree$tip.label %in% names(which(data_structure[1, ] == 1))), pch = 19, cex = 1)
    }

    #Output
    return(data_structure)

}
