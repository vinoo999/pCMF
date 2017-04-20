# project count matrix factor
# sourcing all files (R and Cpp)

project.sources = function(subpath="pkg") {

    ## library
    library(fields)

    ## remove *.o and *.so files ?
    #system(paste0("bash ", WORKINGDIR, "/admin/rm_library_files.sh"))

    ########################################################################
    ##### Cpp
    ########################################################################

    library(Rcpp)

    ## path to Cpp src
    src.path = paste0(WORKINGDIR, "/", subpath, "/src/")

    ## list of R files
    file.list = system(paste0("cd ", src.path, " && git ls-files | grep \"\\\\.cpp\""), intern=TRUE)

    fileToCompile = c("wrapper_varEM_GaP.cpp",
                      "wrapper_varEM_sparse_GaP.cpp",
                      "wrapper_varEM_sparse_ZI_GaP.cpp",
                      "wrapper_varEM_ZI_GaP.cpp",
                      "wrapper_variational_GaP.cpp",
                      "wrapper_variational_sparse_GaP.cpp",
                      "wrapper_variational_ZI_GaP.cpp"
                      )

    # fileToCompile = c("wrapper_varEM_sparse_ZI_GaP.cpp")

    file.list = file.list[file.list %in% fileToCompile]

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
