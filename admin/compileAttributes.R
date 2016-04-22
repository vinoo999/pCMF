#### compile attributes

library(Rcpp)

## PATH
args = commandArgs(trailingOnly=TRUE);
path = args

## dir
pkgDir = paste0(path, "/pkg")

## compile
compileAttributes(pkgdir=pkgDir)