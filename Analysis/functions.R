#sourceDir function (from man source)

load.functions<-function(test=FALSE) {
    sourceDir <- function(path, trace = TRUE, ...) {
        for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
        if(trace) cat(nm,":")
            source(file.path(path, nm), ...)
            if(trace) cat("\n")
        }
    }

    setwd("../Functions")
    sourceDir(path=".")
    setwd("../")
    if(test==TRUE) {
        library(testthat)
        setwd("Functions/test")
        sourceDir(path=".")
        setwd("../../")
    }
    setwd("Analysis/")
}