#Making a hybrid table/figure with the taxonomic coverage result
##########################
#SYNTAX:
#<results> a result table from previous functions
#<order> which order to plot
#<metric> which metric to print
#<thresholds> threshold values (default=c(25,75))
##########################
#guillert(at)tcd.ie - 22/05/2015
##########################

table.result<-function(results, order, metric, thresholds=c(25,75)) {
    #SANITIZING

    #results
    check.class(results, 'data.frame')
    #names
    if(length( match(names(results), c("Order", "Taxonomic.level", "Number.of.OTUs", "Percentage.of.OTUs"))[-which(is.na(match(names(results), c("Order", "Taxonomic.level", "Number.of.OTUs", "Percentage.of.OTUs"))))] ) != 4) {
        stop("Wrong headers names.")
    }

    #order
    check.class(order, 'character')

    #metric
    check.class(metric, 'character')
    if(length(grep(metric, names(results))) == 0) {
        stop(paste(metric, "not found in the data.frame."))
    }
    if(length(grep(paste(metric, "p", sep="_"), names(results))) == 0) {
        stop(paste(metric, "p-value not found in the data.frame."))
    }

    #thresholds
    check.class(thresholds, 'numeric')

    #CREATING THE TABLE/FIGURE
    #Isolating the orders
    if(length(order != 1)) {
        results_tmp<-results[which(results$Order == order[1]),]
        for(ord in 2:length(order)) {
            results_tmp<-rbind(results_tmp, results[which(results$Order == order[ord]),])
        }
    } else {
        results_tmp<-results[which(results$Order == order),]
    }

    #barplot
    #percentages<-table(as.numeric(results_tmp$Percentage.of.OTUs))
    #names(percentages) <- rep(c("Family", "Genera", "Species"), length(order))
    #percentages[]<-as.numeric(results_tmp$Percentage.of.OTUs)
    

    barplot(rev(as.numeric(results_tmp$Percentage.of.OTUs)), horiz=TRUE, xaxt="n", yaxt="n", main="percentage of OTUs")
    for(thresh in 1:length(thresholds)) {
        abline(v=thresholds[thresh], lty=2)
    }
    axis(side=3)

#    axis(side=2, at=1:(length(order)*3), labels=rev(rep(c("Family", "Genera", "Species"), 3)), las=2)
    axis(side=2)
    axis(side=2, at=1:(length(order)*3), labels=rev(rep(c("family", "genus", "species"), 3)), las=2)




            plot(1,1 , xlab="", ylab="", ylim = c(1 - xspc, ncol(mulTree.mcmc) + xspc), xlim = coeff.lim, type = "n", xaxt = "n", yaxt="n", bty = "n", ...)
            #coeff.estimates (is x)
            axis(side = 3)
            #terms (is y)
            axis(side = 2, at = 1:ncol(mulTree.mcmc), labels = rev(terms), las=2) #reverse the terms to go from top to bottom

}