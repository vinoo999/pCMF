### Copyright 2016-05 Ghislain DURIF
###
### This file is part of the `countMatrixFactor' library for R and related languages.
### It is made available under the terms of the GNU General Public
### License, version 2, or at your option, any later version,
### incorporated herein by reference.
###
### This program is distributed in the hope that it will be
### useful, but WITHOUT ANY WARRANTY; without even the implied
### warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
### PURPOSE.  See the GNU General Public License for more
### details.
###
### You should have received a copy of the GNU General Public
### License along with this program; if not, write to the Free
### Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
### MA 02111-1307, USA

#' @title blockMatrix
#'
#' @description
#' Generation by block with gamma distribution,
#' each block has specific parameter values
#'
#' @details
#' The matrix contains \code{nblock} blocks in rows and columns,
#' the entries of each block are generated from a Gamma distribution
#' with specific parameters for each block. In fact, it uses
#' distribution of the form Gamma(1,alpha), i.e. Exponential(alpha).
#' The parameter alpha is set for each block.
#'
#' @author
#' Ghislain Durif, \email{ghislain.durif@univ-lyon1.fr}
#'
#' @param nrow number of rows
#' @param ncol number of columns
#' @param nRowBlock number of blocks in rows and cols
#' @param nColBlock number of blocks in rows and cols
#' @param signalBlock matrix nrow x ncol of parameters used to generate
#' the matrix entries within each block
#'
#' @return list containing the following
#' \item{mat}{block matrix}
#' \item{signalBlock}{see input parameters}
#' \item{idRows}{vector indicating the division of rows into blocks}
#' \item{idRows}{vector indicating the division of columns into blocks}
#'
#' @export
blockMatrixGamma = function(nrow, ncol, nRowBlock, nColBlock, signalBlock=NULL) {

    ###### verification on input paramters
    if((ncol < nColBlock) || (nrow < nRowBlock)) {
        stop("message from blockMatrix: more blocks than columns or rows")
    }

    if(is.null(signalBlock)) {
        if(nRowBlock != nColBlock) {
            stop("message from blockMatrix: consider supplying signalBlock in input")
        } else {
            signalBlock=diag(2.9, nRowBlock)+0.1
        }
    } else {
        if((nrow(signalBlock) != nRowBlock) || (ncol(signalBlock) != nColBlock)) {
            stop("message from blockMatrix: matrix signalBlock of wrong dimensions")
        }
    }

    mat = matrix(NA, nrow=nrow, ncol=ncol)

    ###### significant blocks

    # identification of the blocks in rows and columns
    id.rowblock = 1:nRowBlock
    id.colblock = 1:nColBlock

    # assignation to a blocks for each row and columns
    idRows = sort(rep(id.rowblock, length=nrow))
    idCols = sort(rep(id.colblock, length=ncol))

    ###### construction of the matrix (by blocks)
    for(rowBlock in id.rowblock) {
        for(colBlock in id.colblock) {
            rowsBlock = (1:nrow)[idRows ==  rowBlock]
            nrowsBlock = length(rowsBlock)
            colsBlock = (1:ncol)[idCols ==  colBlock]
            ncolsBlock = length(colsBlock)
            blockSize = nrowsBlock * ncolsBlock

            mat[rowsBlock, colsBlock] = matrix(rexp(blockSize, rate=1/signalBlock[rowBlock, colBlock]), nrow=nrowsBlock, ncol=ncolsBlock) #matrix(signalBlock[rowBlock, colBlock], nrow=nrowsBlock, ncol=ncolsBlock)
        }
    }

    ###### return
    return(list(mat=mat, signalBlock=signalBlock, idRows=idRows, idCols=idCols))

}