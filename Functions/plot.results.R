plot.results<-function(order, col_branch, reference, verbose=FALSE) {
    #Extracting the tree
    order_tree<-extract.order(order, one_tree, reference, verbose)
    #Extracting the taxa with morphological data
    order_taxa<-extract.order(order, living_taxa_list, reference, verbose)
    #Saving only the taxa (discard the change log)
    if(class(order_taxa) != "character") {
        order_taxa<-order_taxa$Taxa
    }
    #Substracting the reference list for the order only
    Sub_reference<-subset(reference[which(reference == order),])
    #Making the list of taxa binomial
    order_taxa<-taxa.binomial(order_taxa, Sub_reference)[[1]]
    #Selecting the right edges
    taxa_edges<-which.edge(order_tree, order_taxa)
    edge_colors<-rep(col_branch[2], nrow(order_tree$edge))
    edge_colors[taxa_edges]<-col_branch[1]

    #Making dummy tip states
    order_tree$tip.state<-data.frame(row.names=order_tree$tip.label, A=rep(0,Ntip(order_tree)))

    #Getting the list of family names for each taxa in the tree
    family_class<-NULL
    for (taxon in 1:Ntip(order_tree)) {
        #Selecting the genera name
        genus<-strsplit(order_tree$tip.label[taxon],"_")[[1]][1]
        #Selecting the corresponding family name in the reference list
        family_class[[taxon]]<-Sub_reference$Family[grep(genus, Sub_reference$Genus)[1]]
    }
    #Plot the phylogeny (radial)
    trait.plot(order_tree, order_tree$tip.state, cols=list(A=("white")), class=family_class, font=1, cex.lab=1, edge.color=edge_colors)
}