### Copyright 2017-07 Ghislain DURIF
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


#' @title graphU
#'
#' @description
#' Graphical representation of the individuals coordinates (matrix U)
#' in the Gamma-Poisson Factor model
#'
#' @details
#' Coordinates of individuals in the Gamma-Poisson factor model in a 2-D
#' space, along two of the axes in the dimensional sub-space
#' in the Euclidean geometry, corresponding to the
#' matrix U, or in the geometry related to the Gamma distribution in the
#' exponential family (log), corresponding to matrix logU
#'
#' Graphical representation based on ggplot2
#'
#' see pCMF function output
#'
#' @author
#' Ghislain Durif, \email{gd.dev@libertymail.net}
#'
#' @seealso \code{\link{pCMF}}
#'
#' @useDynLib pCMF
#'
#' @param model a Gamma-Poisson factor model output by the
#' function \code{\link{pCMF}}
#' @param axes a vector of 2 indexes corresponding to the 2 axes to represent
#' @param labels a vector of indidividual labels
#' @param log_representation boolean, indicating if the representation is in
#' the natural geometry associated to the Gamma distribution (log) or in the
#' Euclidean space, default is TRUE
#' @param edit_theme boolean, indicating if the ggplot2 standard should be
#' edited or not
#' @param graph boolean, indicating if the graph should be drawn or not
#'
#' @return the matrix U of individuals coordinates in the lower
#' dimensional sub-space
#'
#' @export
graphU <- function(model, axes, labels=NULL, log_representation=TRUE,
                   edit_theme=TRUE, graph=TRUE) {

    return(matrixPlot(mat=getU(model, log_representation=log_representation),
                      axes=axes, labels=labels,
                      log_representation=log_representation,
                      edit_theme=edit_theme, graph=graph))
}


#' @title graphV
#'
#' @description
#' Graphical representation of the variables contributions (matrix V)
#' in the Gamma-Poisson Factor model
#'
#' @details
#' Coordinates of the variables in the Gamma-Poisson factor model in a 2-D
#' space, along two of the axes in the dimensional sub-space
#' in the Euclidean geometry, corresponding to the
#' matrix V, or in the geometry related to the Gamma distribution in the
#' exponential family (log), corresponding to matrix logV
#'
#' Graphical representation based on ggplot2
#'
#' see pCMF function output
#'
#' @author
#' Ghislain Durif, \email{gd.dev@libertymail.net}
#'
#' @seealso \code{\link{pCMF}}
#'
#' @useDynLib pCMF
#'
#' @param model a Gamma-Poisson factor model output by the
#' function \code{\link{pCMF}}
#' @param axes a vector of 2 indexes corresponding to the 2 axes to represent
#' @param labels a vector of indidividual labels
#' @param log_representation boolean, indicating if the representation is in
#' the natural geometry associated to the Gamma distribution (log) or in the
#' Euclidean space, default is TRUE
#' @param edit_theme boolean, indicating if the ggplot2 standard should be
#' edited or not
#' @param graph boolean, indicating if the graph should be drawn or not
#'
#' @return the matrix U of individuals coordinates in the lower
#' dimensional sub-space
#'
#' @export
graphU <- function(model, axes, labels=NULL, log_representation=TRUE,
                   edit_theme=TRUE, graph=TRUE) {

    return(matrixPlot(mat=getV(model, log_representation=log_representation),
                      axes=axes, labels=labels,
                      log_representation=log_representation,
                      edit_theme=edit_theme, graph=graph))
}


#' @title matrixPlot
#' @keywords internal
#'
#' @description
#' Graphical representation of the variables contributions (matrix V)
#' in the Gamma-Poisson Factor model
#'
#' @details
#' Coordinates of the variables in the Gamma-Poisson factor model in a 2-D
#' space, along two of the axes in the dimensional sub-space
#' in the Euclidean geometry, corresponding to the
#' matrix V, or in the geometry related to the Gamma distribution in the
#' exponential family (log), corresponding to matrix logV
#'
#' Graphical representation based on ggplot2
#'
#' see pCMF function output
#'
#' @author
#' Ghislain Durif, \email{gd.dev@libertymail.net}
#'
#' @seealso \code{\link{pCMF}}
#'
#' @useDynLib pCMF
#'
#' @importFrom ggplot2 ggplot geom_point theme element_line element_blank
#'
#' @param model a Gamma-Poisson factor model output by the
#' function \code{\link{pCMF}}
#' @param axes a vector of 2 indexes corresponding to the 2 axes to represent
#' @param labels a vector of indidividual labels
#' @param log_representation boolean, indicating if the representation is in
#' the natural geometry associated to the Gamma distribution (log) or in the
#' Euclidean space, default is TRUE
#' @param edit_theme boolean, indicating if the ggplot2 standard should be
#' edited or not
#' @param graph boolean, indicating if the graph should be drawn or not
#'
#' @return the matrix U of individuals coordinates in the lower
#' dimensional sub-space
#'
#' @export
matrixPlot <- function(mat, axes=c(1:2), labels=NULL,
                       log_representation=TRUE,
                       edit_theme=TRUE, graph=TRUE) {

    ## check input
    Kmax <- max(axes)
    if(Kmax <= ncol(mat)) {
        stop("'axes' argument is not compatible with 'mat' dimension")
    }
    if(!is.null(labels)) {
        if(length(labels) != nrow(mat)) {
            stop("'labels' argument length is not compatible with 'mat' dimension")
        }
    }

    ## format the data
    dataToPlot <- data.frame(comp1=mat[,axes[1]], comp2=mat[,axes[1]])
    if(!is.null(labels)) {
        dataToPlot$labels <- labels
    }

    ## graph representation
    g <- ggplot(dataToPlot, aes(x=comp1, y=comp2, color=labels))
    g <- g + geom_point()
    if(edit_theme) {
        g <- g + theme(legend.text=element_text(size=14),
                        legend.title=element_text(size=14),
                        axis.text.x=element_text(size=10),
                        axis.title.x=element_text(size=14),
                        axis.title.y=element_text(size=14),
                        axis.text.y=element_text(size=10),
                        strip.text.x=element_text(size=14),
                        strip.text.y=element_text(size=14),
                        plot.title=element_text(size=14))
        g <- g + theme(panel.background=element_rect(fill="white", colour="black"),
                        panel.grid.major=element_line(color="grey90"),
                        panel.grid.minor=element_line(color="grey90"),
                        strip.background=element_blank())
    }
    ## plot graph ?
    if(graph) {
        g
    }
    ## output
    return(g)
}

