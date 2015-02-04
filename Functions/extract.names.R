#Matching the names of a nexus matrix with a TaxonReference list and saving them in a table
##########################
#SYNTAX:
#<matrix> being a nexus format matrix
#<references> being a list of references to be checked chronologically
##########################
#guillert(at)tcd.ie - 13/01/2015
##########################


#FUNCTIONS
#functions for matching one name of the matrix
#Extracting the first match
find.match<-function(matrix_taxon, Ref_list_element) {
    options(warn=-1)
    to_match<-grep(matrix_taxon, as.character(unlist(Ref_list_element)), ignore.case=FALSE, value=TRUE)
    if(length(to_match) == 0) {
        match_found<-FALSE
    } else {
        #Check if the number of matching characters in to_match is the same length as the ones in matrix_taxon (i.e. full match)
        match_length<-nchar(matrix_taxon)

        #If any grep match has the same number of characters
        if(any(nchar(to_match) == match_length)) {
            match_found<-TRUE
        } else {
            match_found<-FALSE
        }
    }
    return(match_found)
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
    if(length(grep("_", matrix_taxon))==1) {
        #Genus names
        genus<-strsplit(matrix_taxon, split="_")[[1]][1]
        #Species name
        species<-strsplit(matrix_taxon, split="_")[[1]][2]

        #Is the genus name matching?
        genus_match<-find.match(genus, Ref_list_element)

        if(genus_match == TRUE) {
            #Find the genus name
            genus_match<-grep(genus, as.character(unlist(Ref_list_element)), ignore.case=TRUE, value=TRUE)
            #Selecting the first one that has the same number of characters
            genus_match<-genus_match[which(nchar(genus_match) == nchar(genus))[1]]

            #Is the species name matching?
            species_match<-find.match(genus, Ref_list_element)
            if(species_match == TRUE) {
                #Find the species name
                species_match<-grep(species, as.character(unlist(Ref_list_element)), ignore.case=TRUE, value=TRUE)
                #Selecting the first one that has the same number of characters
                species_match<-species_match[which(nchar(species_match) == nchar(species))[1]]
                
            } else {
                species_match<-NA
            }

        } else {
            #No matching
            genus_match<-NA
            species_match<-NA
        }

    } else {
        #impossible to split
        genus_match<-NA
        species_match<-NA
    }

    return(c(genus_match, species_match))
}


match.taxon<-function(matrix_taxon, Ref_list) {
    #Checking through the first element of the list (Fritz tree)
    match_found<-find.match(matrix_taxon, Ref_list[[1]])

    if(match_found==TRUE) {
        #Saving the information
        is_living<-TRUE
        is_source<-names(Ref_list)[1]
        is_level<-"Species"

    } else {
        #Checking through the second element of the list (Wilson Reeder list)
        match_found<-find.match(matrix_taxon, Ref_list[[2]])
            
        if(match_found==TRUE) {
            #Saving the information
            is_living<-TRUE
            is_source<-names(Ref_list)[2]
            is_match<-grep(matrix_taxon, as.character(unlist(Ref_list[[2]])), ignore.case=TRUE, value=TRUE)[1]

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
                match_found<-find.match(matrix_taxon, Ref_list[[3]])

                if(match_found==TRUE) {
                    #Saving the information
                    is_living<-FALSE
                    is_source<-names(Ref_list)[3]
                    is_match<-grep(matrix_taxon, as.character(unlist(Ref_list[[3]])), ignore.case=TRUE, value=TRUE)[1]

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
                    genus_species<-check.split(matrix_taxon, Ref_list[[3]])

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

extract.names<-function(matrix, references, verbose=FALSE) {
    #requirements
    require(ape)
    
    #matrix
    check.class(matrix, "list", " must be a nexus matrix.\n Use read.nexus.data().")
    #tip labels
    check.class(class(names(matrix)), "character", " must be a nexus matrix.\n Use read.nexus.data().")
    #number of characters
    n_characters<-length(matrix[[1]])
    #list of taxa
    taxa_list<-names(matrix)

    #references
    check.class(references, "list", " must be a list of data.frames.")
    #checking each element
    for(n in 1:length(references)) {
        check.class(references[[n]], "data.frame", " must be a list of data.frames.")
    }
    #the list must have names
    if(is.null(names(references))) {
        stop("The Reference list must have names provided.")
    }

    #verbose
    check.class(verbose, "logical", " must be logical.")

    #EXTRACTING NAMES

    #Building the empty data frame.
    output<-data.frame("Taxa"=taxa_list, "Characters"=rep(n_characters, length(taxa_list)), "Living"=rep(NA, length(taxa_list)), "Source"=rep(NA, length(taxa_list)), "Tax.level"=rep(NA, length(taxa_list)))

    for (taxon in 1: length(taxa_list)) {
        match_results<-match.taxon(taxa_list[taxon], references)
        output$Living[taxon]<-match_results[1]
        output$Source[taxon]<-match_results[2]
        output$Tax.level[taxon]<-match_results[3]
        if(verbose==TRUE) {
            message('.', appendLF=FALSE)
        }
    }

    return(output)
}

#Extra checking function to match the names measured as NAs with Wilson Reeder's list 
check.NA<-function(NA_taxa, Reference) {
    #Check if the name contains a "_" or not
    if(length(grep("_", NA_taxa))!=0) {
        NA_taxa<-strsplit(NA_taxa, split="_")[[1]][1]
    }

    #Check if the name contains a "." or not
    if(length(grep(".", NA_taxa))!=0) {
        NA_taxa<-strsplit(NA_taxa, split="_")[[1]][1]
    }

    #Check if it correspond to any column name in WR
    column_names<-names(Reference)[match(NA_taxa, Reference)]

    #Replacing column_names by NA if empty or by the first column if not NA
    if(length(column_names) == 0) {
        column_names<-NA
    } else {
        if(length(column_names) > 1) {
            column_names<-column_names[1]
        }
    }
    return(column_names)
}
