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
load("../Data/List_of_matching_taxa/List_of_matching_taxa.Rda")

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
extraction_table<-extraction_list[[1]]

#Selecting the living taxa (for a minimal number of characters)
living_taxa<-extraction_table[which(extraction_table$Living == TRUE),]
living_taxa_list<-unique(living_taxa$Taxa)

#Loading the results
load("results_1.Rda")

#Remove the whole class data (first 3 rows)
results <- results[-c(1:3),]

######################
#Sampling effort
######################
# Looking at the correlation between the number of OTUs with data and the total number of OTUs at each taxonomic level

sampled_families <- log(as.numeric(sapply(strsplit(results[which(results[,2] == "family"),3], split="/"), "[", 1)))
total_families <- log(as.numeric(sapply(strsplit(results[which(results[,2] == "family"),3], split="/"), "[", 2)))
sampled_genera <- log(as.numeric(sapply(strsplit(results[which(results[,2] == "genus"),3], split="/"), "[", 1)))
total_genera <- log(as.numeric(sapply(strsplit(results[which(results[,2] == "genus"),3], split="/"), "[", 2)))
sampled_species <- log(as.numeric(sapply(strsplit(results[which(results[,2] == "species"),3], split="/"), "[", 1)))
if(any(is.infinite(sampled_species))) {sampled_species[which(is.infinite(sampled_species))] <- 0} #removing 0s for loging
total_species <- log(as.numeric(sapply(strsplit(results[which(results[,2] == "species"),3], split="/"), "[", 2)))



pdf("../Manuscript/Supplementary/Supp_figure_sampling_effort.pdf", width = 8, height = 3.6)
#quartz(width = 8, height = 3.6)
par(mfrow = (c(1,3)), bty = "n")

total <- total_families
sampled <- sampled_families
model <- summary(lm(sampled + 0 ~ total))
plot(total, sampled, pch=19, col="grey", xlim = c(0,ceiling(max(total))), ylim = c(0,ceiling(max(total))),
#plot(total, sampled, pch=19, col="grey", xlim = c(0,8), ylim = c(0,8),
    ylab = "sampled families (log)", xlab = "total families (log)", las = 1)
abline(model$coefficients[[1]], model$coefficients[[2]])
text(ceiling(max(total)) - 0.85*ceiling(max(total)), ceiling(max(total)) - 0.035*ceiling(max(total)), "A", cex=2)
text(ceiling(max(total)) - 0.85*ceiling(max(total)), ceiling(max(total)) - 0.1*ceiling(max(total)),
    bquote(italic(R)^2 == .(format(model$adj.r.squared, digits = 3))) )
text(ceiling(max(total)) - 0.85*ceiling(max(total)), ceiling(max(total)) - 0.15*ceiling(max(total)),
    bquote(italic(b) == .(format(model$coefficients[[2]], digits = 3))))

total <- total_genera
sampled <- sampled_genera
model <- summary(lm(sampled + 0 ~ total))
plot(total, sampled, pch=19, col="grey", xlim = c(0,ceiling(max(total))), ylim = c(0,ceiling(max(total))),
#plot(total, sampled, pch=19, col="grey", xlim = c(0,8), ylim = c(0,8),
    ylab = "sampled genera (log)", xlab = "total genera (log)", las = 1)
abline(model$coefficients[[1]], model$coefficients[[2]])
text(ceiling(max(total)) - 0.85*ceiling(max(total)), ceiling(max(total)) - 0.035*ceiling(max(total)), "B", cex=2)
text(ceiling(max(total)) - 0.85*ceiling(max(total)), ceiling(max(total)) - 0.1*ceiling(max(total)),
    bquote(italic(R)^2 == .(format(model$adj.r.squared, digits = 3))) )
text(ceiling(max(total)) - 0.85*ceiling(max(total)), ceiling(max(total)) - 0.15*ceiling(max(total)),
    bquote(italic(b) == .(format(model$coefficients[[2]], digits = 3))))


total <- total_species
sampled <- sampled_species
model <- summary(lm(sampled + 0 ~ total))
plot(total, sampled, pch=19, col="grey", xlim = c(0,ceiling(max(total))), ylim = c(0,ceiling(max(total))),
#plot(total, sampled, pch=19, col="grey", xlim = c(0,8), ylim = c(0,8),
    ylab = "sampled species (log)", xlab = "total species (log)", las = 1)
abline(model$coefficients[[1]], model$coefficients[[2]])
text(ceiling(max(total)) - 0.85*ceiling(max(total)), ceiling(max(total)) - 0.035*ceiling(max(total)), "C", cex=2)
text(ceiling(max(total)) - 0.85*ceiling(max(total)), ceiling(max(total)) - 0.1*ceiling(max(total)),
    bquote(italic(R)^2 == .(format(model$adj.r.squared, digits = 3))) )
text(ceiling(max(total)) - 0.85*ceiling(max(total)), ceiling(max(total)) - 0.15*ceiling(max(total)),
    bquote(italic(b) == .(format(model$coefficients[[2]], digits = 3))))
dev.off()