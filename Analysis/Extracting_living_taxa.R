#Extracting the living taxa from all the matrices
setwd("~/PhD/Projects/Missing_living_mammals/Analysis")

#Load functions
message("Loading/testing functions:", appendLF=FALSE)
source("functions.R")
load.functions(test=FALSE) #Set test=FALSE to speed up the loading
message("Done\n", appendLF=FALSE)

#List of matrices
matrices_list<-list.files("../Data/Matrices/")

#transform verbose warnings in read.nexus.data in messages
#Renaming the function (backup)
read.nex.dat<-read.nexus.data
#Substitute element 19 (stop message if wrong number of taxa)
body(read.nex.dat)[[19]]<-substitute(
    if (tot.ntax != 0) {
        message(paste("ntax:", ntax, "differ from actual number of taxa in file?", sep=" "), appendLF=FALSE)
        stop("nexus parser did not read names correctly (tot.ntax!=0)")
    }
)
#Substitute element 20 (stop message if wrong number of characters)
body(read.nex.dat)[[20]]<-substitute(
    for (i in 1:length(Obj)) {
        if (length(Obj[[i]]) != nchar) {
            message(paste(names(Obj[i]), "has", length(Obj[[i]]), "characters", sep=" "), appendLF=FALSE)
            stop("nchar differ from sequence length (length(Obj[[i]])!=nchar)")
        }
    }
)

#List of references
#Fritz tree
tree_list<-read.table("../Data/Taxon_References/FritzTree_taxa.txt", header=F, stringsAsFactors=F)
#Wilson Reeder's Mammals Species of the World
WR_list<-read.csv("../Data/Taxon_References/WilsonReederMSW.csv", header=T, stringsAsFactors=F)
#Palaeobiology database
Paleo_list<-read.csv("../Data/Taxon_References/pdb_taxa.csv", header=T, stringsAsFactors=F)
#Combining the references
Reference_list<-list(tree_list, WR_list, Paleo_list)
names(Reference_list)<-c("Fritz Tree", "Wilson & Reeder", "Paleobio database")

#Building the empty table
extraction_table<-data.frame("Matrix"=NA,"Taxa"=NA, "Characters"=NA, "Living"=NA,"Source"=NA,"Tax.level"=NA)

#Build the table for each matrix
for(matrix in 1:length(matrices_list)){
    #trial mode (continue if one matrix is bugged)
    #be verbose!
    message(paste(matrices_list[matrix], ":",sep=""), appendLF=FALSE)
    #Try to load the matrix
    current_matrix<-NULL
    try(current_matrix<-read.nex.dat(paste("../Data/Matrices/", matrices_list[matrix], sep="")), silent=TRUE)

    if(!is.null(current_matrix)) {
        #If successfully, proceed:
        #Extracting the taxa
        extraction_table_tmp<-extract.names(current_matrix, Reference_list, verbose=TRUE)
        #Saving the matrix name
        extraction_table_tmp$Matrix<-rep(matrices_list[matrix], length(current_matrix))
        #Binding with the previous table
        extraction_table<-rbind(extraction_table, extraction_table_tmp)
        message("Done\n", appendLF=FALSE)
    } else {
        #Else, save the matrix name and continue
        extraction_table_tmp<-data.frame("Matrix"=matrices_list[matrix],"Taxa"="FAILURE", "Characters"="FAILURE", "Living"="FAILURE","Source"="FAILURE","Tax.level"="FAILURE")
        extraction_table<-rbind(extraction_table, extraction_table_tmp)
        message("... FAILED\n", appendLF=FALSE)
    }
}

#Removing the first (dummy) line of the extraction_table
extraction_table<-extraction_table[-1,]

#Save the results as .rda
save(extraction_table, file="../Data/List_of_matching_taxa.Rda")

#Save the result as csv
write.csv(extraction_table, file="../Data/List_of_matrching_taxa.csv")
