#Test (example) on primates
#setwd("~/PhD/Projects/Missing_living_mammals/Functions")


#Test
#Tree
#Read Fritz tree
one_tree<-read.nexus("../../Data/Trees/FritzTree.rs200k.1tree.tre")
#Isolate primates
primatesMRCA<-getMRCA(one_tree, c("Homo_sapiens", "Loris_tardigradus"))

primates_tree<-extract.clade(one_tree, node=primatesMRCA)

#Species list
LorisMRCA<-getMRCA(primates_tree, c("Galago_zanzibaricus", "Loris_tardigradus"))
loris_tree<-extract.clade(primates_tree, LorisMRCA)
PlathMRCA<-getMRCA(primates_tree, c("Aotus_azarae", "Ateles_belzebuth"))
plath_tree<-extract.clade(primates_tree, PlathMRCA)

species<-c(sample(loris_tree$tip.label, 20), sample(plath_tree$tip.label, 60))

#testing
expect_error(data.structure("Bob", "The_builder")) ; message('.', appendLF=FALSE)
data_structure<-data.structure(species, tree)
expect_is(data_structure, "matrix") ; message('.', appendLF=FALSE)
expect_equal(nrow(data_structure), 1) ; message('.', appendLF=FALSE)
expect_equal(ncol(data_structure), Ntip(tree)) ; message('.', appendLF=FALSE)

