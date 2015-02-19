#Loading the functions
setwd("~/PhD/Projects/Missing_living_mammals/Analysis")
source("functions.R")
load.functions(test=FALSE) #Set test=FALSE to speed up the loading
library(picante)

#Loading the tree
full_trees<-read.nexus("../Data/Trees/FritzTree.rs200k.100trees.tre")
#Selecting one random tree
set.seed(1) ; one_tree<-full_trees[[sample(1:100, 1)]]

#Loading the taxonomic reference list
WR_list<-read.csv("../Data/Taxon_References/WilsonReederMSW.csv", header=T, stringsAsFactors=F)

#Loading the list of present taxa
#source("Extracting_living_taxa.R")
load("../Data/List_of_matching_taxa.Rda")

#Loading the results
load("results_1.Rda")

#xtable
library(xtable)
#Reformat table (digit + only four first columns)
table_to_print<-results[,c(1:4)]
#Rounding to two digits only
table_to_print[,4]<-round(as.numeric(table_to_print[,4]), digit=2)
#Select the rows with less than 25% data
low_data<-which(table_to_print[,4]<25)
#make the whole table as.character
for (column in 1:ncol(table_to_print)) {
    table_to_print[,column]<-as.character(table_to_print[,column])
}
#Mark the low_data rows to be highlighted in the LaTeX table (BOLD)
for (column in 1:ncol(table_to_print)) {
    table_to_print[low_data,column]<-paste('BOLD',table_to_print[low_data,column])
}

#Bold rows for xtable
bold.cells<-function(x) gsub('BOLD(.*)',paste('\\\\textbf{\\1','}'),x)

#Saving the table
print(xtable(table_to_print), file="../Manuscript/Tables/results_1-table.tex", include.rownames=FALSE, tabular.environment="longtable", floating=FALSE, sanitize.text.function=bold.cells)
