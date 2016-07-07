# project count matrix factor
# sourcing all R files

project.sources = function(subpath="pkg") {

    ## library
    library(fields)

    ########################################################################
    ##### R
    ########################################################################


    ## path to R src
    src.path = paste0(WORKINGDIR, "/", subpath, "/R/")

    ## list of R files
    file.list = system(paste0("cd ", src.path, " && git ls-files | grep \"\\\\.R\""), intern=TRUE)

    ## sourcing all R files
    path.list = paste0(src.path, file.list)
    sapply(path.list, source)


}

project.sources()
