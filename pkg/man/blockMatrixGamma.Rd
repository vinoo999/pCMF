% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/blockMatrixGamma.R
\name{blockMatrixGamma}
\alias{blockMatrixGamma}
\title{blockMatrix}
\usage{
blockMatrixGamma(nrow, ncol, nRowBlock, nColBlock, signalBlock = NULL)
}
\arguments{
\item{nrow}{number of rows}

\item{ncol}{number of columns}

\item{nRowBlock}{number of blocks in rows and cols}

\item{nColBlock}{number of blocks in rows and cols}

\item{signalBlock}{matrix nrow x ncol of parameters used to generate
the matrix entries within each block}
}
\value{
list containing the following
\item{mat}{block matrix}
\item{signalBlock}{see input parameters}
\item{idRows}{vector indicating the division of rows into blocks}
\item{idRows}{vector indicating the division of columns into blocks}
}
\description{
Generation by block with gamma distribution,
each block has specific parameter values
}
\details{
The matrix contains \code{nblock} blocks in rows and columns,
the entries of each block are generated from a Gamma distribution
with specific parameters for each block. In fact, it uses
distribution of the form Gamma(1,alpha), i.e. Exponential(alpha).
The parameter alpha is set for each block.
}
\author{
Ghislain Durif, \email{gd.dev@libertymail.net}
}
\keyword{internal}
