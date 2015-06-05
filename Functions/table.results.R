#Making a hybrid table/figure with the taxonomic coverage result
##########################
#SYNTAX:
#<results> a result table from previous functions
#<metric> which metric to print
#<thresholds> threshold values (default=c(25,75))
#<save.path> where to save the table
#<file.save> name (or chain name) for saving the table
#<caption> table caption (if split option is called, the caption will only appear for the first table)
#<environement> table environement
##########################
#environement(at)tcd.ie - 05/06/2015
##########################

table.result<-function(results, metric, threshold=c(25,75), save.path, file.save, caption, environement) {
    
    #DEBUG
    #warning("Debug mode is ON (table.result)")
    #results -> must be a table
    #metric=c("NRI", "NTI") #-> can be a list (or 0) but must be present in the results header
    #threshold=c(25,75) #-> can be any values (or NULL) between 0 and 100
    #split=c(24,30,30) #-> can be a list of values (sub-tables) or NULL
    #save.path="../Manuscript/" #-> must be text
    #file.save="Tab_test" #-> must be text
    #caption="Here's my caption" #-> must be text
    #environement="longtable" #if split != null set to table. Else must be text

    #Selecting the metric columns
    metric_col<-NULL
    #Select the columns of metrics
    for(met in 1:length(metric)) {
        metric_col[met]<-grep(metric[met], names(results))
    }
    #Select the columns of p_values
    for(met in 1:length(metric)) {
        metric_col[length(metric)+met]<-grep(metric[met], names(results))+1
    }

    #Creating the sub table
    table_to_print<-results[,c(1:4,metric_col)]
    #make the text part of the table as.character
    table_to_print[,c(1,2)]<-apply(table_to_print[,c(1,2)], 2, as.character)
    #make the numeric part as numeric
    table_to_print[,c(4:(4+2*length(metric)))]<-apply(table_to_print[,c(4:(4+2*length(metric)))], 2, as.numeric)
    #Change cases
    table_to_print[,1]<-capwords(table_to_print[,1], strict=TRUE)
    table_to_print[,2]<-tolower(table_to_print[,2])
    #Order alphabetically
    table_to_print<-table_to_print[order(table_to_print[,1]),]

    #Marking the significant differences with stars
    significant_storage<-NULL
    for(met in 1:length(metric)) {
        #Isolating the metrics (and rounding)
        metrics<-round(table_to_print[,4+met], digit=2)
        #Isolating the pvalues
        p_values<-table_to_print[,4+met+2]

        #replacing the metrics with stars if pvalue > 0.05
        #One star
        one_star<-which(p_values <= 0.05)
        significant_storage<-c(significant_storage, one_star)
        if(length(one_star) != 0) {
            metrics[one_star]<-paste(metrics[one_star], "*", sep="")
        }
        #Two stars
        two_star<-which(p_values <= 0.005)
        significant_storage<-c(significant_storage, two_star)
        if(length(one_star) != 0) {
            metrics[two_star]<-paste(metrics[two_star], "*", sep="")
        }
        #Three stars
        three_star<-which(p_values <= 0.0005)
        significant_storage<-c(significant_storage, three_star)
        if(length(three_star) != 0) {
            metrics[three_star]<-paste(metrics[three_star], "*", sep="")
        }

        #Replace NA by blank (space)
        metrics[which(is.na(metrics))]<-" "

        #Replacing the results in the final tables
        table_to_print[,4+met]<-metrics
    }

    #Removing the p_value columns
    table_to_print<-table_to_print[,1:(length(metric)+4)]

    #Mark the significant rows to be highlighted in the LaTeX table (BOLD)
    #Saving the percentage column for later
    percentages<-table_to_print[,4]
    #Selecting the significant data (rows)
    for (column in 1:ncol(table_to_print)) {
        table_to_print[unique(significant_storage), column]<-paste('BOLD',table_to_print[unique(significant_storage), column], sep="")
    }
    bold.cells<-function(x) gsub('BOLD(.*)',paste('\\\\textbf{\\1','}',sep=""),x)
    #refilling the percentage columne
    table_to_print[,4]<-percentages


    #Creating the folder for storing the graphical elements
    new_folder<-paste(save.path, "Table_figures/", sep="")
    #If folder does not exist or is empty, create it
    if(length(list.files(new_folder)) == 0) {
        system(paste("mkdir", new_folder))
    } else {
    #Else remove it and create a new one
        system(paste("rm -R", new_folder))
        system(paste("mkdir", new_folder))
    }

    #Make sure the percentage column is numeric

    #Creating all the barplots
    for(row in 1:nrow(table_to_print)) {
        #pdf settings
        pdf(paste(new_folder, "bar", row , ".pdf", sep=""), width=4, height=1)
        #margin settings
        par(mar=c(0,1,0,1))
        #par plot
        barplot(table_to_print[ row ,4], horiz=TRUE, xlim=c(1,100), xaxt="n")
        #threshold lines
        for(line in 1:length(threshold)) {
            abline(v=threshold[line], lty=3,lwd=1)
        }
        dev.off()
    }

    #Replacing percentage column by the graphs links
    bar_template<-"\\includegraphics[width=0.20\\linewidth, height=0.05\\linewidth]{Table_figures/"
    for(row in 1:nrow(table_to_print)) {
        table_to_print[ row ,4]<-paste(bar_template, "bar", row ,".pdf}", sep="")
    }


    #Fixing the column names
    colnames(table_to_print)<-c("Order", "Taxonomic level", "Proportion of taxa", "Coverage", c(metric))
    #Saving the table
    table<-xtable(table_to_print)
    caption(table)<-caption
    label(table)<-file.save
    align(table)<-c("l", "l", "L{1.8cm}", "C{2cm}", "l", "c", "c")
    print(table, file=paste(save.path,file.save, ".tex", sep=""), include.rownames=FALSE, tabular.environment=environement, floating=FALSE, sanitize.text.function=bold.cells, caption=caption, caption.placement="top")

    cat("Add the following to your document header:\n
\\newcolumntype{L}[1]{>{\\raggedright\\let\\newline\\\\\\arraybackslash\\hspace{0pt}}m{#1}}
\\newcolumntype{C}[1]{>{\\centering\\let\\newline\\\\\\arraybackslash\\hspace{0pt}}m{#1}}
\\newcolumntype{R}[1]{>{\\raggedleft\\let\\newline\\\\\\arraybackslash\\hspace{0pt}}m{#1}}")

}