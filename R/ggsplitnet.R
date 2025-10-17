#' @method fortify networx
#' @export
fortify.networx <- function(model, data, layout = "unrooted", ladderize = FALSE,
                            right = FALSE, mrsd = NULL, as.Date = FALSE, ...) {
    nTips <- length(model$tip.label)
    label <-  rep(NA_character_, nrow(model$edge))
    isTip <- logical(nrow(model$edge))  ## edge leading to tip
    if (!is.null(model$translate)) {
        ind <- match(model$translate$node, model$edge[, 2])
        label[ind] <- model$translate$label
    } else {
        ind <- match(seq_len(nTips), model$edge[, 2])
        label[ind] <- model$tip.label
    }
    isTip[ind] <- TRUE
    if(is.null(model$edge.length)) model$edge.length <- rep(1, nrow(model$edge))
    df <- data.frame(node = model$edge[, 2], parent = model$edge[, 1],
                    branch.length = model$edge.length,
                    split = model$splitIndex, label = label, isTip = isTip)
    if (!is.null(model$.plot)) coord <- model$.plot$vertices
    else coord <- phangorn::coords(model, dim = "equal_angle")
    df <- cbind(df, x = coord[df$node, 1], y = coord[df$node, 2],
                xend = coord[df$parent, 1], yend = coord[df$parent, 2])
    angle <- atan2(df$y - df$yend, df$x - df$xend) * 360/(2 * pi)
    angle[angle < 0] <- angle[angle < 0] + 360
    df <- cbind(df, angle = angle)
    df
}


#' drawing phylogenetic tree from phylo object
#'
#'
#' @title ggsplitnet
#' @param tr a networx object
#' @param mapping aes mapping
#' @param layout so far only 'slanted' is supported.
#' @param mrsd most recent sampling date
#' @param as.Date logical whether using Date class in time tree
#' @param yscale y scale
#' @param yscale_mapping yscale mapping for category variable
#' @param ladderize logical
#' @param right logical
#' @param branch.length variable for scaling branch, if 'none' draw cladogram
#' @param ndigits number of digits to round numerical annotation variable
#' @param ... additional parameter
#' @return tree
#' @seealso \code{\link[ggtree]{ggtree}}, \code{\link[phangorn]{networx}},
#' \code{\link[phangorn]{consensusNet}}, \code{\link[phangorn]{neighborNet}}
#' @references Schliep, K., Potts, A. J., Morrison, D. A. and Grimm, G. W.
#' (2017), Intertwining phylogenetic trees and networks.
#' \emph{Methods Ecol Evol}. \bold{8}, 1212--1220. doi:10.1111/2041-210X.12760
#'
#' Dress, A.W.M. and Huson, D.H. (2004) Constructing Splits Graphs
#' \emph{IEEE/ACM Transactions on Computational Biology and Bioinformatics
#' (TCBB)}, \bold{1(3)}, 109--115
#'
#' Bagci, C., Bryant, D., Cetinkaya, B. and Huson, D.H. (2021), Microbial
#' Phylogenetic Context Using Phylogenetic Outlines. \emph{Genome Biology and
#' Evolution}. Volume 13. Issue 9. evab213
#'
#' @importFrom utils modifyList
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 annotate
#' @importFrom ggplot2 scale_x_reverse
#' @importFrom ggplot2 ylim
#' @importFrom ggplot2 coord_flip
#' @importFrom ggplot2 coord_polar
#' @importFrom rlang sym
#' @importFrom ggtree theme_tree
#' @importFrom methods is
#' @author Klaus Schliep
#' @examples
#' data(yeast, package='phangorn')
#' dm <- phangorn::dist.ml(yeast)
#' nnet <- phangorn::neighborNet(dm)
#' ggsplitnet(nnet) + geom_tiplab2()
#' @export
ggsplitnet <- function(tr, mapping = NULL, layout = "slanted",
        mrsd = NULL, as.Date = FALSE, yscale = "none", yscale_mapping = NULL,
        ladderize = FALSE, right = FALSE, branch.length = "branch.length",
        ndigits = NULL, ...) {
    if(!is(tr, "networx")) stop("tr must be of class 'networx'")
    layout <- match.arg(layout, c("slanted"))
    if (is.null(mapping)) {
      mapping <- aes(!!sym("x"), !!sym("y"))
    } else {
        mapping <- modifyList(aes(!!sym("x"), !!sym("y")), mapping)
    }
    p <- suppressWarnings(ggplot(tr, mapping = mapping, layout = layout,
          mrsd = mrsd, as.Date = as.Date, yscale = yscale,
          yscale_mapping = yscale_mapping, ladderize = ladderize,
          right = right, branch.length = branch.length, ndigits = ndigits, ...))
    p <- p + geom_splitnet(layout = layout, ...)
    p <- p + theme_tree()
    class(p) <- c("ggtree", class(p))
    return(p)
}





##' add splitnet layer
##'
##'
##' @title geom_splitnet
##' @param layout one of 'rectangular', 'slanted', 'circular', 'radial' or
##' 'unrooted'
##' @param ... additional parameter
##' @return splitnet layer
##' @examples
##' data(yeast, package='phangorn')
##' dm <- phangorn::dist.ml(yeast)
##' nnet <- phangorn::neighborNet(dm)
##' ggplot(nnet, aes(x, y))  + geom_splitnet() + theme_tree()
##' @importFrom ggplot2 geom_segment
##' @importFrom ggplot2 aes
##' @export
##' @author Klaus Schliep
geom_splitnet <- function(layout = "slanted", ...) {
    x <- y <- xend <- yend <- parent <- NULL
    lineend <- "round"
    if (layout == "rectangular" || layout == "fan" || layout == "circular") {
        list(geom_segment(aes(x = x, xend = xend, y = y, yend = y),
            lineend = lineend, ...),
            geom_segment(aes(x = xend, xend = xend, y = y, yend = yend),
            lineend = lineend, ...))
    } else if (layout == "slanted" || layout == "radial" ||
            layout == "unrooted") {
        geom_segment(aes(x = x, xend = xend, y = y, yend = yend),
                    lineend = lineend, ...)
    }
}
