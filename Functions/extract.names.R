#Matching the names of a nexus matrix with a TaxonReference list and saving them in a table
##########################
#SYNTAX:
#<matrix> being a nexus format matrix
#<references> being a list of references to be checked chronologically
##########################
#guillert(at)tcd.ie - 13/01/2015
##########################

extract.names<-function(matrix, references) {

    #DEBUG
    matrix<-read.nexus.file("Matrices/Gl1988-GH.nex")
    tree_list<-read.table("Taxon_References/FritzTree_taxa.txt", header=F, stringsAsFactors=F)
    WR_list<-read.csv("Taxon_References/WilsonReederMSW.csv", header=T, stringsAsFactors=F)
    Paleo_list<-read.csv("Taxon_References/pdb_taxa.csv", header=T, stringsAsFactors=F)
    Ref_list<-list(tree_list, WR_list, Paleo_list)
    names(Ref_list)<-c("tree", "WR", "Paleo")

    #requirements
    require(ape)

    #Sanitizing
    #Sanitizing functions
    source("sanitizing.R")
    
    #matrix
    check.class(matrix, "list", " must be a nexus matrix.\n Use read.nexus.data().")
    #tip labels
    check.class(class(names(matrix)), "character", " must be a nexus matrix.\n Use read.nexus.data().")
    #number of characters
    n_characters<-length(matrix[[1]])
    #list of taxa
    taxa_list<-names(matrix)

    #references
    check.class(Ref_list, "list", " must be a list of data.frames.")
    #checking each element
    for(n in 1:length(Ref_list)) {
        check.class(Ref_list[[n]], "data.frame", " must be a list of data.frames.")
    }
    #the list must have names
    if(is.null(names(Ref_list))) {
        stop("The Reference list must have names provided.")
    }

    #FUNCTION
    #functions for matching one name of the matrix

    #Extracting the first match
    find.match<-function(matrix_taxon, Ref_list_element) {
        options(warn=-1)
        return(grep(matrix_taxon, as.character(unlist(Ref_list_element)), ignore.case=TRUE, value=TRUE))
    }

    #Extracting the number of the column containing is.match
    level.column<-function(is_match, Ref_list_element) {
        level_column<-NA
        for(n in 1:ncol(Ref_list_element)) {
        level_column[n]<-match(is_match, as.character(Ref_list_element[,n]))[1]
        }
        return(which(!is.na(level_column)))
    }

    #Try matching by separating genus and species
    check.split<-function(matrix_taxon, Ref_list_element) {
        if(grep("_", matrix_taxon)==1) {
            #Genus names
            genus<-strsplit(matrix_taxon, split="_")[[1]][1]
            #Species name
            species<-strsplit(matrix_taxon, split="_")[[1]][2]

            #Finding the genus name
            genus_match<-find.match(genus, Ref_list_element)

            #Finding the species name
            species_match<-grep(species, Ref_list_element[genus_match,])

        } else {
            #impossible to split
            genus_match<-NA
            species_match<-NA
        }

        return(c(genus_match, species_match))
    }

    match.taxon<-function(matrix_taxon, Ref_list) {
        #Checking through the first element of the list (Fritz tree)
        is_match<-find.match(matrix_taxon, Ref_list[[1]])[1]

        if(!is.na(is_match)) {
            #Saving the information
            is_living<-TRUE
            is_source<-names(Ref_list)[1]
            is_level<-"Species"

        } else {
            #Checking through the second element of the list (Wilson Reeder list)
            is_match<-find.match(matrix_taxon, Ref_list[[2]])[1]
            
            if(!is.na(is_match)) {
                #Saving the information
                is_living<-TRUE
                is_source<-names(Ref_list)[2]
                is_level<-names(Ref_list[[2]])[(level.column(is_match, Ref_list[[2]]))]
            
            } else {

                #Try to split the genera/species
                genus_species<-check.split(matrix_taxon, Ref_list[[2]])

                #If genus or genus_species is present, save the information
                if(!is.na(genus_species[1])) {
                    #Saving the information
                    is_living<-TRUE
                    is_source<-names(Ref_list)[2]

                    #Checking the level
                    if(!is.na(genus_species[2])) {
                        is_level<-"Species"
                    } else {
                        is_level<-"Genus"
                    }

                } else {
                    #Checking through the third element of the list (Paleo DB)
                    is_match<-find.match(matrix_taxon, Ref_list[[3]])[1]

                    if(!is.na(is_match)) {
                        #Saving the information
                        is_living<-FALSE
                        is_source<-names(Ref_list)[3]

                        #Checking the level
                        if(level.column(is_match, Ref_list[[3]])==1) {
                            #Present in the first column
                            is_level<-"Genus"
                        } else {
                            #Present in the second column
                            is_level<-"Species"
                        }

                    } else {
                        #Try to split the genera/species
                        genus_species<-check.split(matrix_taxon, Ref_list[[3]])[1]

                        #If genus or genus_species is present, save the information
                        if(!is.na(genus_species[1])) {
                            #Saving the information
                            is_living<-FALSE
                            is_source<-names(Ref_list)[3]

                            #Checking the level
                            if(!is.na(genus_species[2])) {
                                is_level<-"Species"
                            } else {
                                is_level<-"Genus"
                            }

                        } else {
                            #No match at all!
                            is_living<-NA
                            is_source<-NA
                            is_level<-NA
                        } 
                    }
                }
            }
        }

        #Output
        return(c(is_living, is_source, is_level))
    }

    #DEBUG
    matrix_taxon<-"Bob_the_builder" ; match.taxon(matrix_taxon, Ref_list) #Should be NAs !!
    matrix_taxon<-"Tachyglossus_aculeatus" ; match.taxon(matrix_taxon, Ref_list)
    matrix_taxon<-"T._aculeatus" ; match.taxon(matrix_taxon, Ref_list) #Should be NAs !!
    matrix_taxon<-"Suidae" ; match.taxon(matrix_taxon, Ref_list)
    matrix_taxon<-"Amphiperatherium_maximum" ; match.taxon(matrix_taxon, Ref_list)
    matrix_taxon<-"Amphiperatherium" ; match.taxon(matrix_taxon, Ref_list)



#End
}

#Output table

#matrix, taxon, characters, living, source, level
#GS254, Homo_sapiens, 160, TRUE, tree, species