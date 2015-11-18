#Analyses of the distribution of living taxa with morphological data

#Function input
warning("Set the directory to '/Missing_living_mammals/Analysis'")
#setwd("~/PhD/Projects/Missing_living_mammals/Analysis")
source("functions.R")
load.functions(test=FALSE) #Set test=FALSE to speed up the loading

#Example data input
#Read Fritz tree
one_tree<-read.nexus("../Data/Trees/mammalST_MSW05_best_chrono.tre")
#Taxonomic reference
WR_list<-read.csv("../Data/Taxon_References/WilsonReederMSW.csv", header=T, stringsAsFactors=F)

#Example on primates
#Isolate the clade
primatesMRCA<-getMRCA(one_tree, c("Homo_sapiens", "Loris_tardigradus"))
primates_tree<-extract.clade(one_tree, node=primatesMRCA)
#list of all taxa
primates_list<-primates_tree$tip.label

#Building three sampling data sets
#Clustered one (biased)
#Only lorises and platirhines (clustered)
LorisMRCA<-getMRCA(primates_tree, c("Galago_zanzibaricus", "Loris_tardigradus")) ; loris_tree<-extract.clade(primates_tree, LorisMRCA)
PlathMRCA<-getMRCA(primates_tree, c("Aotus_azarae", "Ateles_belzebuth")) ; plath_tree<-extract.clade(primates_tree, PlathMRCA)
clustered_primates<-c(sample(loris_tree$tip.label, 20), sample(plath_tree$tip.label, 60))

#Random primates
random_primates<-sample(primates_list, 80)

#Even primates
even_primates<-primates_list[seq(from=1, to=351, by=5)]

#Visualizing the different data distributions
par(mfrow=c(1,3))
silent<-data.structure(clustered_primates, primates_tree, plot=TRUE, main="Clustered")
silent<-data.structure(random_primates, primates_tree, plot=TRUE, main="Random")
silent<-data.structure(even_primates, primates_tree, plot=TRUE, main="Even")

#Calculating the distribution of present taxa
structure_clustered<-order.structure("Primates (clustered)", clustered_primates, primates_tree, reference=WR_list, runs=100, null="taxa.labels", verbose=TRUE)
structure_random<-order.structure("Primates (random)", random_primates, primates_tree, reference=WR_list, runs=100, null="taxa.labels", verbose=TRUE)
structure_even<-order.structure("Primates (even)", even_primates, primates_tree, reference=WR_list, runs=100, null="taxa.labels", verbose=TRUE)

#Printing the results
rbind(structure_clustered, structure_random, structure_even)