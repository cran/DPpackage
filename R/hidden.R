
.onAttach <- function(...) {
     cat("##\n## Dirichlet Process Package (DPpackage)\n")
     cat("## Copyright (C) 2006, Alejandro Jara\n")
     cat("## Biostatistical Centre \n")
     cat("## Catholic University of Leuven \n")
     cat("##\n## Support provided by the KUL-PUC bilateral\n")
     cat("## (Belgium-Chile) grant BIL05/03\n##\n")
}

.onUnload <- function(libpath) {
    library.dynam.unload("DPpackage", libpath)
}

