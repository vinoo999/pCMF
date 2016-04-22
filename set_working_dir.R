#### set working dir
PROJECT="countMatrixFactor"
HOME=Sys.getenv("HOME")
PROJECTDIR=paste0("source_code", "/", PROJECT)
WORKINGDIR=paste0(HOME, "/", PROJECTDIR)
setwd(WORKINGDIR)
LOCLIB=NULL
myLib="installDir"
# FIGUREDIR="figures"
# if(!dir.exists(FIGUREDIR)) {
#     dir.create(FIGUREDIR, recursive=TRUE, mode="0755")
# }
