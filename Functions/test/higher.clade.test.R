#Test (example) on primates
#setwd("~/PhD/Projects/Missing_living_mammals/Functions")

#Read Fritz tree
full_trees<-read.nexus("../../Data/Trees/FritzTree.rs200k.100trees.tre")
set.seed(1) ; one_tree<-full_trees[[sample(1:100, 1)]]

primatesMRCA<-getMRCA(one_tree, c("Homo_sapiens", "Loris_tardigradus"))
primates_tree<-extract.clade(one_tree, node=primatesMRCA)
#list of all taxa
primates_list<-primates_tree$tip.label
#Reference list
WR_list<-read.csv("../../Data/Taxon_References/WilsonReederMSW.csv", header=T, stringsAsFactors=F)

#Testing

#Data
reference<-WR_list
tree<-primates_tree
species<-primates_list
taxonomic.level="Genus"

#Test
expect_error(higher.clade(species, tree, taxonomic.level="Bob_the_builder", reference)) ; message('.', appendLF=FALSE)
#FAMILY LEVEL
res_fam<-higher.clade(species, tree, taxonomic.level="Family", reference)
#Output class
expect_is(res_fam, "list") ; message('.', appendLF=FALSE)
expect_equal(length(res_fam), 3) ; message('.', appendLF=FALSE)
#OTUs list
expect_equal(names(res_fam[1]), "OTUs") ; message('.', appendLF=FALSE)
expect_equal(length(res_fam[[1]]), 15) ; message('.', appendLF=FALSE)
expect_is(res_fam[[1]], "character") ; message('.', appendLF=FALSE)
#Tree
expect_equal(names(res_fam[2]), "tree") ; message('.', appendLF=FALSE)
expect_equal(length(res_fam[[2]]), 4) ; message('.', appendLF=FALSE)
expect_is(res_fam[[2]], "phylo") ; message('.', appendLF=FALSE)
#Bionomial list
expect_equal(names(res_fam[3]), "Binomial") ; message('.', appendLF=FALSE)
expect_equal(length(res_fam[[3]]), 15) ; message('.', appendLF=FALSE)
expect_is(res_fam[[3]], "character") ; message('.', appendLF=FALSE)

#GENUS LEVEL
res_gen<-higher.clade(species, tree, taxonomic.level="Genus", reference)
#Output class
expect_is(res_gen, "list") ; message('.', appendLF=FALSE)
expect_equal(length(res_gen), 3) ; message('.', appendLF=FALSE)
#OTUs list
expect_equal(names(res_gen[1]), "OTUs") ; message('.', appendLF=FALSE)
expect_equal(length(res_gen[[1]]), 68) ; message('.', appendLF=FALSE)
expect_is(res_gen[[1]], "character") ; message('.', appendLF=FALSE)
#Tree
expect_equal(names(res_gen[2]), "tree") ; message('.', appendLF=FALSE)
expect_equal(length(res_gen[[2]]), 4) ; message('.', appendLF=FALSE)
expect_is(res_gen[[2]], "phylo") ; message('.', appendLF=FALSE)
#Bionomial list
expect_equal(names(res_gen[3]), "Binomial") ; message('.', appendLF=FALSE)
expect_equal(length(res_gen[[3]]), 68) ; message('.', appendLF=FALSE)
expect_is(res_gen[[3]], "character") ; message('.', appendLF=FALSE)
