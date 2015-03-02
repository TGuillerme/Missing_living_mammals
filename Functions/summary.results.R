summary.results<-function(results, metric, save.path, file.save, caption, environement) {

    if(metric == "proportion") {
        ############################################
        #Table with the number of taxa containing morphological data
        ############################################

        #Reformat table (digit + only four first columns)
        table_to_print<-results[,c(1:4)]
        #Rounding to two digits only
        table_to_print[,4]<-round(as.numeric(table_to_print[,4]), digit=2)
        #Select the rows with less than 25% data = >75% missing data
        low_data<-which(table_to_print[,4]<25)
        #make the whole table as.character
        for (column in 1:ncol(table_to_print)) {
            table_to_print[,column]<-as.character(table_to_print[,column])
        }
        #Lower case
        table_to_print[,1]<-capwords(table_to_print[,1], strict=TRUE)
        #Mark the low_data rows to be highlighted in the LaTeX table (BOLD)
        for (column in 1:ncol(table_to_print)) {
            table_to_print[low_data,column]<-paste('BOLD',table_to_print[low_data,column], sep="")
        }
        #Fixing the column names
        colnames(table_to_print)<-c("Order", "Taxonomic level", "Fraction of OTUs", "Percentage of OTUs")
        #Saving the table
        table<-xtable(table_to_print)
        caption(table)<-caption
        label(table)<-file.save
        print(table, file=paste(save.path,file.save, ".tex",sep=""), include.rownames=FALSE, tabular.environment=environement, floating=FALSE, sanitize.text.function=bold.cells, caption.placement="top")

    } else {

        ############################################
        #Table with the data structure for each taxa
        ############################################

        #Selecting the right columns
        metric_col<-grep(metric, names(results))
        table_to_print<-results[,c(1:4,metric_col)]
        #Removing the NAs
        table_to_print<-table_to_print[-which(is.na(table_to_print[,5])),]
        #Selecting the rows with significant difference from null
        signif_data<-which(table_to_print[,6]<0.05)
        #make the whole table as.character
        for (column in 1:ncol(table_to_print)) {
            table_to_print[,column]<-as.character(table_to_print[,column])
        }
        #Lower case
        table_to_print[,1]<-capwords(table_to_print[,1], strict=TRUE)
        #Mark the low_data rows to be highlighted in the LaTeX table (BOLD)
        for (column in 1:ncol(table_to_print)) {
            table_to_print[signif_data, column]<-paste('BOLD',table_to_print[signif_data, column], sep="")
        }
        #Removing column Fraction
        table_to_print<-table_to_print[,-3]
        #Fixing the column names
        colnames(table_to_print)<-c("Order", "Taxonomic level", "Percentage of OTUs", metric, "p-value")
        #Saving the table
        table<-xtable(table_to_print)
        caption(table)<-caption
        label(table)<-file.save
        print(table, file=paste(save.path,file.save, ".tex", sep=""), include.rownames=FALSE, tabular.environment=environement, floating=FALSE, sanitize.text.function=bold.cells, caption=caption, caption.placement="top")
    }
}
