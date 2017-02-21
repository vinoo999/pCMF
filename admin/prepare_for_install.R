#!/usr/bin/Rscript

R_LIBS <- Sys.getenv("R_LIBS")


install.packages(c("Rcpp",
                   "RcppEigen",
                   "roxygen2",
                   "BH",
                   "fields"),
                 lib=R_LIBS)
