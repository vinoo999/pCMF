### Copyright 2016-06 Ghislain DURIF
###
### This file is part of the `pCMF' library for R and related languages.
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

#' @title matrixHeatmap
#'
#' @description
#' plot the entries of a matrix as a heatmap
#'
#' @details
#' hello world
#'
#' @author
#' Ghislain Durif, \email{gd.dev@libertymail.net}
#'
#' @importFrom fields image.plot
#'
#' @param mat the matrix to plot
#'
#' @export
matrixHeatmap <- function(mat, ...) {
    image.plot(t(apply(mat, 2, rev)), xaxt="n", yaxt="n", ...)
}
