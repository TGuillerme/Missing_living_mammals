#Loading the functions
#Be sure to be in the right repository
source("functions.R")
load.functions(test=FALSE) #Set test=FALSE to speed up the loading
library(picante)
library(xtable)

#Loading the tree
mam_tree<-read.nexus("../Data/Trees/mammalST_MSW05_best_chrono.tre")
#Resolve polytomies
mam_tree<-multi2di(mam_tree)
#Selecting one random tree
#set.seed(1) ; mam_tree<-mam_tree[[sample(1:100, 1)]]

#Loading the taxonomic reference list
#Loading the taxonomic reference list
WR_list<-read.csv("../Data/Taxon_References/WilsonReederMSW.csv", header=T, stringsAsFactors=F)
#Changing ARTIODACTYLA and CETACEA to CETARTIODACTYLA
WR_list$Order[which(WR_list$Order == "ARTIODACTYLA")]<-"CETARTIODACTYLA"
WR_list$Order[which(WR_list$Order == "CETACEA")]<-"CETARTIODACTYLA"

#Loading the list of present taxa
#source("Extracting_living_taxa.R")
load("../Data/List_of_matching_taxa.Rda")

#Fixing the number of characters for truncated matrices
extraction_table[which(extraction_table$Matrix == "GS2007-R.nex"), 3]<-88
extraction_table[which(extraction_table$Matrix == "GS2010-VJD.nex"), 3]<-157
extraction_table[which(extraction_table$Matrix == "GS2012-GM.nex"), 3]<-45
extraction_table[which(extraction_table$Matrix == "GS2012-SM.nex"), 3]<-157
extraction_table[which(extraction_table$Matrix == "GS1997-V.nex"), 3]<-40
extraction_table[which(extraction_table$Matrix == "GS2013-MS.nex"), 3]<-48
extraction_table[which(extraction_table$Matrix == "GS2014-MG.nex"), 3]<-40
extraction_table[which(extraction_table$Matrix == "GS2014-MaL.nex"), 3]<-98
extraction_table[which(extraction_table$Matrix == "GS2015-RR.nex"), 3]<-149
extraction_table[which(extraction_table$Matrix == "GS2014-B.nex"), 3]<-27
extraction_table[which(extraction_table$Matrix == "GS2014-Y.nex"), 3]<-42

#Creating sub tables with a minimum number of characters
#1 Character
extraction_table1<-extraction_table[which(as.numeric(extraction_table$Characters) >= 1),]

#100 Characters
extraction_table100<-extraction_table[which(as.numeric(extraction_table$Characters) >= 100),]

#Combine the sub tables
extraction_list=list(extraction_table1, extraction_table100)
results_names=c("results_1.Rda","results_100.Rda")

#Select the character threshold
extraction_table<-extraction_list[[2]]

#Selecting the living taxa (for a minimal number of characters)
living_taxa<-extraction_table[which(extraction_table$Living == TRUE),]
living_taxa_list<-unique(living_taxa$Taxa)

######################
#Loading the results
######################

#Loading the results
load("results_100.Rda")


############################################
#Fable with number of taxa with morpho data+data structure
############################################
caption<-"Number of taxa with available cladistic data for mammalian orders at three taxonomic levels. The coverage represents the proportion of taxa with available morphological data. The left vertical bar represents 25\\% of available data (``low'' coverage if \\textless 25\\%); The right vertical bar represents 75\\% of available data (``high'' coverage if \\textgreater 75\\%). When the Net Relatedness Index (NRI) and the Nearest Taxon Index (NTI) are negative, taxa are more phylogenetically dispersed than expected by chance; when NRI or NTI are positive, taxa are more phylogenetically clustered by expected by chance. Significant NRI or NTI are highlighted in bold. One star (*) represents a p-value between 0.05 and 0.005; two starts between 0.005 and 0.0005 and three stars a p-value less than 0.0005."
table.result(results, metric=c("NRI","NTI"), threshold=c(25,75), save.path="../Manuscript/", file.save="Table_results", caption=caption, environement="longtable")

#Bugged!
#caption<-"Number of taxa with available cladistic data for mammalian orders at three taxonomic levels. The coverage represents the proportion of taxa with available morphological data. The left vertical bar represents 25\\% of available data (``low'' coverage if \\textless 25\\%); The right vertical bar represents 75\\% of available data (``high'' coverage if \\textgreater 75\\%). When the Net Relatedness Index (NRI) and the Nearest Taxon Index (NTI) are negative, taxa are more phylogenetically dispersed than expected by chance; when NRI or NTI are positive, taxa are more phylogenetically clustered by expected by chance. Significant NRI or NTI are highlighted in bold. One star (*) represents a p-value between 0.05 and 0.005; two starts between 0.005 and 0.0005 and three stars a p-value less than 0.0005."
#table.result(results, metric=c("PD","NRI","NTI"), threshold=c(25,75), save.path="../Manuscript/Supplementary/results_PD_100c/", file.save="Table_results_PD", caption=caption, environement="longtable")


############################################
#Figure (phylogeny example)
############################################

#Setting the colour scheme (data, no data)
pdf("../Manuscript/example_coverage.pdf", width=12.45, height=6)
#quartz(width = 16.6, height = 8) #A4 landscape
#quartz(width = 8.3, height = 5.8) #A5 landscape
op<-par(mfrow=c(1,2), oma=c(0,1,0,1))
plotA<-plot.results(order="PRIMATES", tree=mam_tree, taxa=living_taxa_list, col_branch=c("blue", "grey"), reference=WR_list, verbose=TRUE)
text(x=(plotA$x.lim[1]-plotA$x.lim[1]*0.05),y=(plotA$y.lim[2]-plotA$y.lim[2]*0.05),"A",cex=3)
plotB<-plot.results(order="CARNIVORA", tree=mam_tree, taxa=living_taxa_list, col_branch=c("blue", "grey"), reference=WR_list, verbose=TRUE)
text(x=(plotB$x.lim[1]-plotB$x.lim[1]*0.05),y=(plotB$y.lim[2]-plotB$y.lim[2]*0.05),"B",cex=3)
par(op)
dev.off()

#Isolating the order names
orders<-as.character(unique(results[,1]))
#Isolating the number of OTUs per order
OTUs_number<-results[which(as.character(results[,2])=="species"),3]
OTUs_number<-strsplit(OTUs_number, split="/")
OTUS_number_val<-as.numeric(unlist(OTUs_number))
OTUS_total<-OTUS_number_val[-seq(from=1, to=length(OTUS_number_val), by=2)]
#Selecting the orders with more than 20 OTUs
supp_order<-orders[which(as.numeric(OTUS_total) >= 20)]

#Supplementary figures (orders with more than 20 taxa)
for (order in 11:length(supp_order)) {
    pdf(paste("../Manuscript/Supplementary/Supp_figure_", supp_order[[order]], ".pdf", sep=""))
    op<-par(mfrow=c(1,1), oma=c(0,1,0,1))
    plot.results(order=supp_order[[order]], tree=mam_tree, taxa=living_taxa_list, col_branch=c("blue", "grey"), reference=WR_list, verbose=TRUE)
    par(op)
    dev.off()
}


#Supplementary (number of characters per matrix)
#Threshold=1 character
char<-sort(unique(as.numeric(extraction_table1$Characters)))
pdf("../Manuscript/Supplementary/Supp_figure_1MorphoCharperMatrices.pdf")
hist(char, breaks=100, xlab="Number of morphological characters", ylab="Matrices", main="")
dev.off()



############################################
#Table with the data structure for each taxa
############################################
caption<-"Data structure for the orders with OTUs without morphological data per taxonomic level. When the Nearest Taxon Index (NTI) is negative, the OTUs are more dispersed than expected by chance (random); when the NTI is positive, the OTUs are more clustered by expected by chance. The p-value indicates the significance in difference from the null model (random)."
summary.results(results, metric="NTI", save.path="../Manuscript/Supplementary/",file.save="Supp_data_structureNTI", caption=caption, environement="longtable")
caption<-"Data structure for the orders with OTUs without morphological data per taxonomic level. When the Faith's Phylogenetic Distance (PD). The p-value indicates the significance in difference from the null model (random)."
summary.results(results, metric="PD", save.path="../Manuscript/Supplementary/",file.save="Supp_data_structurePD", caption=caption, environement="longtable")

#Results for threshold 1
load("results_1.Rda")
#Supplementary tables (wrapper)
caption<-"Proportion of available OTUs with morphological data per order and per taxonomic level (Character threshold = 1). We highlighted in bold the orders that have more than 75\\% of missing data for each taxonomic level. Note that it is possible that more data is available at a higher taxonomic level (Genus $>$ Species) since if the species name for an OTU was not or miss specified, we still counted the OTU for higher taxonomic level analysis."
summary.results(results, metric="proportion", save.path="../Manuscript/Supplementary/",file.save="morpho_taxa_proportion_threshold1", caption=caption, environement="longtable")
caption<-"Data structure for the orders with OTUs without morphological data per taxonomic level (Character threshold = 1). When the Net Relatedness Index (NRI) is negative, the OTUs are more dispersed than expected by chance (random); when the NRI is positive, the OTUs are more clustered by expected by chance. The p-value indicates the significance in difference from the null model (random)."
summary.results(results, metric="NRI", save.path="../Manuscript/Supplementary/",file.save="Supp_data_structureNRI_threshold1", caption=caption, environement="longtable")
caption<-"Data structure for the orders with OTUs without morphological data per taxonomic level (Character threshold = 1). When the Nearest Taxon Index (NTI) is negative, the OTUs are more dispersed than expected by chance (random); when the NTI is positive, the OTUs are more clustered by expected by chance. The p-value indicates the significance in difference from the null model (random)."
summary.results(results, metric="NTI", save.path="../Manuscript/Supplementary/",file.save="Supp_data_structureNTI_threshold1", caption=caption, environement="longtable")
caption<-"Data structure for the orders with OTUs without morphological data per taxonomic level (Character threshold = 1). When the Faith's Phylogenetic Distance (PD). The p-value indicates the significance in difference from the null model (random)."
summary.results(results, metric="PD", save.path="../Manuscript/Supplementary/",file.save="Supp_data_structurePD_threshold1", caption=caption, environement="longtable")




############################################
#Fable with number of taxa with morpho data+data structure (no threshold)
############################################
#Loading the results
load("results_1.Rda")

caption<-"Number of taxa with available cladistic data for mammalian orders at three taxonomic levels. The coverage represents the proportion of taxa with available morphological data. The left vertical bar represents 25\\% of available data (``low'' coverage if \\textless 25\\%); The right vertical bar represents 75\\% of available data (``high'' coverage if \\textgreater 75\\%). When the Net Relatedness Index (NRI) and the Nearest Taxon Index (NTI) are negative, taxa are more phylogenetically dispersed than expected by chance; when NRI or NTI are positive, taxa are more phylogenetically clustered by expected by chance. Significant NRI or NTI are highlighted in bold. One star (*) represents a p-value between 0.05 and 0.005; two starts between 0.005 and 0.0005 and three stars a p-value less than 0.0005."
table.result(results, metric=c("NRI","NTI"), threshold=c(25,75), save.path="../Manuscript/Supplementary/Results_1c", file.save="Table_results", caption=caption, environement="longtable")
