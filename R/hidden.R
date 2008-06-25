
.onAttach <- function(...) {
     cat("##\n## Dirichlet Process Package (DPpackage)\n")
     cat("## Copyright (C) 2006-2007-2008, Alejandro Jara\n")
     cat("## Biostatistical Centre \n")
     cat("## Catholic University of Leuven \n")
     cat("##\n## Support provided by the KUL-PUC bilateral\n")
     cat("## (Belgium-Chile) grant BIL05/03 and the NIH/NCI\n")
     cat("## R01CA75981 grant \n##\n")
}




.onUnload <- function(libpath) {
    library.dynam.unload("DPpackage", libpath)
}


"print.anovaPsCP"<-
function (x, digits = max(getOption("digits") - 2, 3), signif.stars = getOption("show.signif.stars"), 
    ...) 
{
    if (!is.null(heading <- attr(x, "heading"))) 
        cat(heading, sep = "\n")

    nc <- dim(x)[2]
    if (is.null(cn <- colnames(x))) 
        stop("'anova' object must have colnames")

    has.P <- substr(cn[nc], 1, 3) == "PsC"

    printCoefmat(x, digits = digits, signif.stars = signif.stars, 
        has.Pvalue = has.P, P.values = has.P,eps.Pvalue=0.01)
    invisible(x)
}

