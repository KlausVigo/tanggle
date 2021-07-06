#' @importFrom magrittr %>%
#' @export
magrittr::"%>%"

#' @method fortify evonet
#' @importFrom ggplot2 fortify
#' @export
fortify.evonet <- function(model, data, layout = "rectangular",
        ladderize = FALSE, right = FALSE, mrsd = NULL, as.Date = FALSE, ...) {
    class(model) <- "phylo"
    # ggtree:::fortify.phylo
    df <- fortify(model, ladderize = ladderize)

    hybridEdge <- logical(nrow(df))
    hybridEdge[grep("#", df$label)] <- TRUE
    df <- cbind(df, hybridEdge = hybridEdge)

    reticulation <- model$reticulation
    df.ret <- df[reticulation[, 1], , drop = FALSE]
    # df.ret <- df[reticulation[,2], , drop=FALSE]
    df.ret[, c("node", "parent")] <- reticulation
    df.ret[, "hybridEdge"] <- TRUE
    df <- rbind(df, df.ret)
    df
}


#' drawing phylogenetic tree from phylo object
#'
#'
#' @title ggevonet
#' @param tr phylo object
#' @param mapping aes mapping
#' @param layout one of 'rectangular', 'slanted'
#' @param open.angle open angle, only for 'fan' layout
#' @param mrsd most recent sampling date
#' @param as.Date logical whether using Date class in time tree
#' @param yscale y scale
#' @param yscale_mapping yscale mapping for category variable
#' @param ladderize logical
#' @param right logical
#' @param branch.length variable for scaling branch, if 'none' draw cladogram
#' @param ndigits number of digits to round numerical annotation variable
#' @param min_crossing logical, rotate clades to minimize crossings
#' @param ... additional parameter
#' @return tree
#' @seealso \code{\link[ape]{evonet}}, \code{\link[ggtree]{ggtree}}
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 annotate
#' @importFrom ggplot2 scale_x_reverse
#' @importFrom ggplot2 ylim
#' @importFrom ggplot2 coord_flip
#' @importFrom ggplot2 coord_polar
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 aes_
#' @importFrom ggplot2 aes_string
#' @importFrom ggtree geom_tree2
#' @importFrom ggtree theme_tree
#' @importFrom phangorn coords
#' @importFrom ape read.evonet
#' @author Klaus Schliep
#' @examples
#' (enet <- ape::read.evonet(text='((a:2,(b:1)#H1:1):1,(#H1,c:1):2);'))
#' ggevonet(enet) + geom_tiplab()
#' @export
ggevonet <- function(tr, mapping = NULL, layout = "slanted", open.angle = 0,
        mrsd = NULL, as.Date = FALSE, yscale = "none", yscale_mapping = NULL,
        ladderize = FALSE, right = FALSE, branch.length = "branch.length",
        ndigits = NULL, min_crossing = TRUE, ...) {
    layout <- match.arg(layout, c("rectangular", "slanted"))
    # , 'fan', 'circular', 'radial', 'unrooted', 'equal_angle', 'daylight'

    if (is.null(tr$edge.length)) {
        nh <- node.depth.evonet(tr)
        tr$edge.length <- nh[tr$edge[, 1]] - nh[tr$edge[, 2]]
    }
    if (min_crossing) {
        tr <- minimize_overlap(tr)
    }
    if (yscale != "none") {
        layout <- "slanted"
    }
    if (is.null(mapping)) {
        mapping <- aes_(~x, ~y)
    } else {
        mapping <- modifyList(aes_(~x, ~y), mapping)
    }
    mapping <- modifyList(aes_string(linetype = "hybridEdge"), mapping)
    p <- ggplot(tr, mapping = mapping, layout = layout, mrsd = mrsd,
                as.Date = as.Date, yscale = yscale,
                yscale_mapping = yscale_mapping, ladderize = ladderize,
                right = right, branch.length = branch.length,
                ndigits = ndigits, ...)
    p <- p + geom_tree2(layout = layout, ...)
    p <- p + theme_tree(legend.position = "none")
    class(p) <- c("ggtree", class(p))
    return(p)
}


#' @title minimize_overlap
#' reduces reticulation lines crossing over in plots
#' @param x Tree of class 'evonet'
#' @return  A Tree with rotated nodes of class 'evonet'
#' @importFrom ape node.height rotate Ntip
#' @importFrom phangorn Ancestors
#' @author L. Francisco Henao Diaz
#' @examples
#' fishnet <- ape::read.evonet(text='(Xalvarezi,Xmayae,((Xsignum,((Xmonticolus,
#' (Xclemenciae_F2,#H25)),(((((((((Xgordoni,Xmeyeri),Xcouchianus),Xvariatus),
#' Xevelynae),(Xxiphidium,#H24)),Xmilleri),Xandersi),Xmaculatus),(((Xmontezumae,
#' (Xcortezi,(Xbirchmanni_GARC,Xmalinche_CHIC2))),((Xnigrensis,Xmultilineatus),
#' (Xpygmaeus,Xcontinens))))#H24))),(Xhellerii)#H25));')
#' fishnet$edge.length <- NULL
#' new_tre <- minimize_overlap(fishnet)
#'
#' par(mfrow=c(1,2))
#' ggevonet(fishnet, min_crossing = FALSE)
#' ggevonet(new_tre)
#'
#' net2 <- ape::read.evonet(text='(15,(1,((14,(#H1,(((12,13),(11,#H3)),(7,
#'     ((10)#H3,(8,9)))))),((((2,3))#H2,(6,(5,(#H2,4)))))#H1)));')
#' # Cui et al. 2013 Evol.
#' new_net2 <- minimize_overlap(net2)
#' ggevonet(net2, min_crossing = FALSE)
#' ggevonet(new_net2)
#'
#' @export
minimize_overlap <- function(x) {
    if (!inherits(x, "evonet"))
        stop("x should be an 'evonet' class")
    n_iter <- round(x$Nnode * 3/4)
    # r_hist <- c()
    for (j in seq_len(n_iter)) {
        h <- node.height(x)
        best_r <- sum(abs(h[x$reticulation[, 1]] - h[x$reticulation[, 2]]))
        best_c <- -1
        # r_hist[j] <- best_r
        nodes2rot <- intersect(sort(unique(unlist(Ancestors(x,
                c(x$reticulation))))), which(tabulate(x$edge[, 1]) > 1))
        for (i in seq_along(nodes2rot)) {
            nh <- node.height(ape::rotate(x, nodes2rot[i]))
            best_nr <- sum(abs(nh[x$reticulation[, 1]] -
                                nh[x$reticulation[, 2]]))
            if (best_nr < best_r) {
                best_c <- nodes2rot[i]
                best_r <- best_nr
            }
        }
        if (best_c > 0)
            x <- ape::rotate(x, best_c) else (break)()
    }
    x
}

#' These functions return the depths or heights of nodes and tips.
#'
#' @title Depth of Nodes
#' @param x an object of class 'evonet'
#' @param \dots Further arguments passed to or from other methods.
#' @return a vector with the depth of the nodes
#' @seealso \code{\link[ape]{node.depth}}
#' @examples
#' z <- ape::read.evonet(text = '((1,((2,(3,(4)Y#H1)g)e,
#' (((Y#H1, 5)h,6)f)X#H2)c)a,((X#H2,7)d,8)b)r;')
#' nd <- node.depth.evonet(z)
#' z$edge.length <- nd[z$edge[,1]] - nd[z$edge[,2]]
#' ggevonet(z)
#'
#' @importFrom phangorn getRoot Ancestors Descendants
#' @export
node.depth.evonet <- function(x, ...) {
    x <- ape::reorder.phylo(x)
    root <- getRoot(x)
    max_nodes <- max(x$edge)
    nTip <- Ntip(x)
    desc <- Descendants(x, seq_len(max_nodes), "children")
    anc <- Ancestors(x)
    pa <- vector("list", max_nodes)
    ind <- which(x$edge[, 2] > Ntip(x))
    pa[x$edge[ind, 2]] <- x$edge[ind, 1]
    for (i in seq_len(nrow(x$reticulation))) {
        pa[[x$reticulation[i, 2]]] <- sort(c(pa[[x$reticulation[i, 2]]],
                                            x$reticulation[i, 1]))
    }
    ind <- which(lengths(pa) > 0)
    depth <- numeric(max_nodes)
    depth[root] <- 1
    done <- logical(max_nodes)
    done[root] <- TRUE
    candidates <- desc[[root]]
    candidates <- candidates[candidates > nTip]
    d <- 1
    while (length(candidates) > 0) {
        active <- vapply(candidates, function(x) all(done[pa[[x]]]), FALSE)
        tmp <- which(active)[1]  #sample(active,1)
        candidates <- c(candidates, desc[[candidates[tmp]]])
        candidates <- candidates[candidates > nTip]
        d <- d + 1
        done[candidates[tmp]] <- TRUE
        depth[candidates[tmp]] <- d
        candidates <- candidates[-tmp]
    }
    depth <- d + 2 - depth
    depth[seq_len(nTip)] <- 1
    depth
}

