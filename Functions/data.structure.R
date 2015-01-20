#Making a data table with the structure of presences/absences of species for further community structure analysis using picante
##########################
#SYNTAX:
#<species> being a list of species of interest
#<tree> being a species trees
#<taxonomic.level> which taxonomic level, should match with the reference file (case sensitive - default="Species")
#<reference> taxonomic reference list
##########################
#guillert(at)tcd.ie - 20/01/2015
##########################

#Test
#Tree
#Read Fritz tree
full_trees<-read.nexus("../../Data/Trees/FritzTree.rs200k.100trees.tre")
set.seed(1) ; one_tree<-full_trees[[sample(1:100, 1)]]
#Isolate primates
primatesMRCA<-getMRCA(one_tree, c("Homo_sapiens", "Loris_tardigradus"))

tree<-extract.clade(one_tree, node=primatesMRCA)

#Species list
LorisMRCA<-getMRCA(primates_tree, c("Galago_zanzibaricus", "Loris_tardigradus"))
loris_tree<-extract.clade(primates_tree, LorisMRCA)
PlathMRCA<-getMRCA(primates_tree, c("Aotus_azarae", "Ateles_belzebuth"))
plath_tree<-extract.clade(primates_tree, PlathMRCA)

species<-c(sample(loris_tree$tip.label, 20), sample(plath_tree$tip.label, 60))

#taxonomic.level
taxonomic.level="Species"

#Taxonomic reference
WR_list<-read.csv("../../Data/Taxon_References/WilsonReederMSW.csv", header=T, stringsAsFactors=F)

reference<-WR_list


#FUNCTIONS


#MAKING THE DATA STRUCTURE TABLE
data.structure<-function(species, tree, taxonomic.level, reference){}

#list of all taxa
primates_list<-primates_tree$tip.label

#Creating the data set with presence / absence
data_structure<-matrix(ncol=length(primates_list), nrow=2, data=0)
colnames(data_structure)<-primates_list
rownames(data_structure)<-c("random", "clustered")
data_structure[1,match(rand_primates, primates_list)]<-1
data_structure[2,match(clustered_primates, primates_list)]<-1
