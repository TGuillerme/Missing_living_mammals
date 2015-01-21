#Community structure analyses
##########################
#SYNTAX:
#<data_structure> the presence/absence list of OTUs
#<tree> the full tree
#<runs> the number of randomization runs
#<null> the null model
##########################
#guillert(at)tcd.ie - 20/01/2015
##########################

community.structure<-function(data_structure, tree, runs=1000, null) {

    library(picante)

    #CALCULATE THE DIFFERENT METRICS
    #Faith's pd
    faith_pd<-pd(data_structure, tree)

    #Corrected Faith's pd
    cor_faith_pd<-faith_pd[,1]/faith_pd[,2]

    #Cophenetic Distances
    distance_matrix<-cophenetic(tree)

    #Mean pairwise distance
    mpd_result<-ses.mpd(data_structure, distance_matrix, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 1000)

    #Mean pairwise distance #NTI
    MPD_NRI<-mpd_result[,6]*-1

    #Absolute distance from expected mpd (null)
    MPD_Null_distance<-abs(mpd_result[,3]-mpd_result[,2])

    #Significance of the absolute distance
    MPD_p_value<-mpd_result[,7]

    #Mean nearest taxon distance
    mtd_result<-ses.mntd(data_structure, distance_matrix, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 1000)

    #Mean nearest taxon distance #NRI
    MTD_NRI<-mtd_result[,6]*-1

    #Absolute distance from expected mtd (null)
    MTD_Null_distance<-abs(mtd_result[,3]-mtd_result[,2])

    #Significance of the absolute distance
    MTD_p_value<-mtd_result[,7]

    #OTU presence / absence
    OTUs<-paste(length(which(data_structure[1,]==1)), Ntip(tree), sep="/")

    #OTU presence / absence ratio
    OTUs_ratio<-length(which(data_structure[1,]==1))/Ntip(tree)*100

    #MPD Null distance
    MPD_Null_distance<-abs(mpd_result[,3]-mpd_result[,2])

    #MPD null distance p-value
    MPD_p_value<-mpd_result[,7]

    #MTD Null distance
    MTD_Null_distance<-abs(mtd_result[,3]-mtd_result[,2])

    #MTD null distance p-value
    MTD_p_value<-mtd_result[,7]

    #output
    results<-c(OTUs_ratio, MPD_NRI[1], MPD_Null_distance[2], MPD_p_value[2], MTD_NRI[2], MTD_Null_distance[2], MTD_p_value[2])
    results<-c(OTUs, round(results, digit=3))
    return(results)

}
