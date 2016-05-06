# project hierarchical clustering
# sourcing all files (R and Cpp)

project.sources = function(subpath="pkg") {

    ## library


    ########################################################################
    ##### Cpp
    ########################################################################

    library(Rcpp)

    ## path to Cpp src
    src.path = paste0(WORKINGDIR, "/", subpath, "/src/")

    ## list of R files
    file.list = system(paste0("cd ", src.path, " && git ls-files | grep \"\\\\.cpp\""), intern=TRUE)

    file.list = file.list[file.list %in% c("gamPoisFactor_wrapper.cpp")]

    print(file.list)

    ## sourcing all R files
    path.list = paste0(src.path, file.list)
    sapply(path.list, sourceCpp)



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
