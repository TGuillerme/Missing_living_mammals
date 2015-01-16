#Draft analyses table

#Ones the living taxa are extracted from all the matrices, count the proportion of available taxa per order (species levels / genus level / family level) + Faith's Phylogenetic Diversity
#

#Test (example) on primates
setwd("~/PhD/Projects/Missing_living_mammals/Analysis")

#Read Fritz tree
full_trees<-read.nexus("../Data/Trees/FritzTree.rs200k.100trees.tre")
set.seed(1) ; one_tree<-full_trees[[sample(1:100, 1)]]


#library(picante)
#psv
#pd
#ses.pd

data(phylocom)
names(phylocom)
phy <- phylocom$phylo
comm <- phylocom$sample
traits <- phylocom$traits

#Remove extra taxa
prunedphy <- prune.sample(comm, phy)

#Communities
par(mfrow = c(2, 3))
for (i in row.names(comm)) {
    plot(prunedphy, show.tip.label = FALSE, main = i)
    tiplabels(tip = which(prunedphy$tip.label %in% names(which(comm[i, ] > 0))), pch = 19, cex = 2)
}

#Traits
par(mfrow = c(2, 2))
for (i in names(traits)) {
    plot(phy, show.tip.label = FALSE, main = i)
    tiplabels(pch = 22, col = traits[, i] + 1, bg = traits[, i] + 1, cex = 1.5)
}

#Faith's PD
pd.result <- pd(comm, phy, include.root = TRUE)
pd.result #low when clustered, high when even or random

#Mean pairwise distance or mean nearest taxon distance
phydist <- cophenetic(phy)
#Mean pairwise distance
ses.mpd.result <- ses.mpd(comm, phydist, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 1000)
ses.mpd.result
ses.mpd.result$mpd.obs.z*-1
#Mean nearest taxon distance
ses.mntd.result <- ses.mntd(comm, phydist, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 1000)
ses.mntd.result
#ntaxa Number of taxa in community
#mpd.obs Observed mpd in community
#mpd.rand.mean Mean mpd in null communities
#mpd.rand.sd Standard deviation of mpd in null communities
#mpd.obs.rank Rank of observed mpd vs. null communities
#mpd.obs.z Standardized effect size of mpd vs. null communities (equivalent to -NRI)
#mpd.obs.p P-value (quantile) of observed mpd vs. null communities (= mpd.obs.rank/ runs + 1)
#runs Number of randomizations

#Isolate only the primates
primatesMRCA<-getMRCA(one_tree, c("Homo_sapiens", "Loris_tardigradus"))
primates_tree<-extract.clade(one_tree, node=primatesMRCA)
#list of all taxa
primates_list<-primates_tree$tip.label

#Sub lists of "present taxa"
#Random
rand_primates<-sample(primates_list, 80)

#Only lorises and platirhines (clustered)
LorisMRCA<-getMRCA(primates_tree, c("Galago_zanzibaricus", "Loris_tardigradus"))
loris_tree<-extract.clade(primates_tree, LorisMRCA)
PlathMRCA<-getMRCA(primates_tree, c("Aotus_azarae", "Ateles_belzebuth"))
plath_tree<-extract.clade(primates_tree, PlathMRCA)
clustered_primates<-c(sample(loris_tree$tip.label, 20), sample(plath_tree$tip.label, 60))

#Test the different metrics (random vs. clustered)
#Creating the data set with presence / absence
data_structure<-matrix(ncol=length(primates_list), nrow=2, data=0)
colnames(data_structure)<-primates_list
rownames(data_structure)<-c("random", "clustered")
data_structure[1,match(rand_primates, primates_list)]<-1
data_structure[2,match(clustered_primates, primates_list)]<-1

#Visualisation
par(mfrow = c(1, 2))
for (i in row.names(data_structure)) {
    plot(primates_tree, show.tip.label = FALSE, main = i)
    tiplabels(tip = which(primates_tree$tip.label %in% names(which(data_structure[i, ] > 0))), pch = 19, cex = 1)
}

#Faith's pd
pd_results<-pd(data_structure, primates_tree)
#Corrected Faith's pd
corpd_resu<-pd_results[,1]/pd_results[,2]

#Mean pairwise distance
distance_matrix<-cophenetic(primates_tree)
mpd_result<-ses.mpd(data_structure, distance_matrix, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 1000)
mtd_result<-ses.mntd(data_structure, distance_matrix, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 1000)
MPD_NRI<-mpd_result[,6]*-1 #NRI
MPD_Null_distance<-abs(mpd_result[,3]-mpd_result[,2]) #Absolute distance from expected mpd (null)
MPD_p_value<-mpd_result[,7] #Significance of the absolute distance

MTD_NRI<-mtd_result[,6]*-1 #NRI
MTD_Null_distance<-abs(mtd_result[,3]-mtd_result[,2]) #Absolute distance from expected mtd (null)
MTD_p_value<-mtd_result[,7] #Significance of the absolute distance


#Empty results table
result_table<-data.frame("Order"=rep("Primates",3), "Taxonomic level"=c("Family", "Genus", "Species"), "Number of OTUs"=rep(NA,3), "Percentage of OTUs"=rep(NA,3), "MPD(NRI)"=rep(NA,3), "MPD(NullDist)"=rep(NA,3), "p_value"=rep(NA,3), "MTD(NRI)"=rep(NA,3), "MTD(NullDist)"=rep(NA,3), "p_value"=rep(NA,3))

#Fill in the table
OTUs<-paste(length(clustered_primates), Ntip(primates_tree), sep="/")
OTUs_ratio<-length(clustered_primates)/Ntip(primates_tree)*100
MPD_NRI<-mpd_result[,6]*-1
MPD_Null_distance<-abs(mpd_result[,3]-mpd_result[,2])
MPD_p_value<-mpd_result[,7]
MTD_NRI<-mtd_result[,6]*-1
MTD_Null_distance<-abs(mtd_result[,3]-mtd_result[,2])
MTD_p_value<-mtd_result[,7] 
primates_species<-c(OTUs_ratio, MPD_NRI[2], MPD_Null_distance[2], MPD_p_value[2], MTD_NRI[2], MTD_Null_distance[2], MTD_p_value[2])

result_table[3,3]<-OTUs ; result_table[3,4:10]<-round(primates_species, digit=3)

#Add function for scoring per genus and per family