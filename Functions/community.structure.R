#Community structure analyses
##########################
#SYNTAX:
#<data_structure> the presence/absence list of OTUs
#<tree> the full tree
#<runs> the number of randomization runs
#<null> the null model. See ?ses.pd for more details (default = "taxa.labels")
#<...> any optional arguments to be passed to ses.mpd and ses.mntd
##########################
#guillert(at)tcd.ie - 20/01/2015
##########################

community.structure<-function(data_structure, tree, runs=1000, null="taxa.labels", ...) {

    library(picante)

    #CALCULATE THE DIFFERENT METRICS
    #Phylogenetic species variability. psv=1 when species are unrelated (star phylogeny) and psv=0 when species are more related. 
    #PSV<-psv(data_structure, tree)[1,1]

    #Phylogenetic species richness psr=n when sp are unrelated and psr=0 when are more related
    #PSR<-psr(data_structure, tree)[1,1]

    #Phylogenetic species evenness. pse=1 when species are evenely spread (good), and 0 when clumped
    #PSE<-pse(data_structure, tree)[1,1]

    #Phylogenetic species clustering. psc=1 when species are less clustered and 0 when clumped
    #PSC<-psc(data_structure, tree)[1,1]

    #Phylogenetic diversity (sum of the total phylogenetic branch length)
    PD<-pd(data_structure, tree)
    #Relative PD = sum of branch length in sample / sum of branch length in the tree (when rel_PD=1, all the species are represented in the tree, when rel_PD=0, no species are present)
    rel_PD<-PD[1,1]/PD[2,1]

    #ses_pd<-ses.pd(data_structure, primates_tree, runs=1000, null.model=null)
    #ntaxa   Number of taxa in community
    #pd.obs  Observed PD in community
    #pd.rand.mean    Mean PD in null communities
    #pd.rand.sd  Standard deviation of PD in null communities
    #pd.obs.rank Rank of observed PD vs. null communities
    #pd.obs.z    Standardized effect size of PD vs. null communities (= (pd.obs - pd.rand.mean) / pd.rand.sd)
    #pd.obs.p    P-value (quantile) of observed PD vs. null communities (= mpd.obs.rank / runs + 1)
    #runs    Number of randomizations
    #PD<-ses_pd$pd.obs.z[1]#

    #Mean pairwise distance
    mpd_result<-ses.mpd(data_structure, cophenetic(tree), null.model=null, runs=runs, ...)

    #Nearest relatedness index; NRI=0 -> random ; NRI<0 -> over-dispersion (good) ; NRI>0 -> clustering (bad)
    NRI<-mpd_result$mpd.obs.z[1]*-1

    #Significance of difference from NULL
    MPD_p_value<-mpd_result$mpd.obs.p[1]

    #Mean nearest taxon distance
    mtnd_result<-ses.mntd(data_structure, cophenetic(tree), null.model=null, runs=runs, ...)

    #Nearest taxon index; NTI=0 -> random ; NTI<0 -> over-dispersion (good) ; NTI>0 -> clustering (bad)
    NTI<-mtnd_result$mntd.obs.z[1]*-1

    #Significance of difference from NULL
    MNTD_p_value<-mtnd_result$mntd.obs.p[1]

    #OTU presence / absence
    OTUs<-paste(length(which(data_structure[1,]==1)), Ntip(tree), sep="/")

    #OTU presence / absence ratio
    OTUs_ratio<-length(which(data_structure[1,]==1))/Ntip(tree)*100

    #output
    results<-c(OTUs_ratio, rel_PD, NRI, MPD_p_value, NTI, MNTD_p_value)
    results<-c(OTUs, round(results, digit=3))
    return(results)

}
