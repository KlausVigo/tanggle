
#' Swapping the minor edges of an evonet object
#'
#'
#' @title swap_hybrid_minor
#' @param x evonet object
#' @param hybrid_nodes a vector of hybrid nodes to have their minor edges swapped
#' @param node_times an optional argument with node times
#' @return network
#' @examples
#' (enet <- ape::read.evonet(text='((a:2,(b:1)#H1:1):1,(#H1,c:1):2);'))
#' ggevonet(enet) + geom_tiplab()
#' swapped_enet<-swap_hybrid_minor(enet,6)
#' ggevonet(swapped_enet) + geom_tiplab()
#' @importFrom ape node.depth.edgelength
#' @export


swap_hybrid_minor <- function(x, hybrid_nodes, node_times=NULL){
  if(is.null(node_times)){
    node_times <- node.depth.edgelength(x)
  }

  if(any(!(hybrid_nodes %in% x$reticulation[,2]))){
    stop('You goofed. hybrid_nodes only takes hybrid nodes but was given a tree node. Look at x$reticulation[,2] for hybrid nodes')
  }

  for(hyb_nd in hybrid_nodes){
    ##get the edges that point to the hybrid node of interest
    edge_row <- which(x$edge[,2]==hyb_nd) ##get the edge in the 'edge' matrix
    ret_row <- which(x$reticulation[,2]==hyb_nd) ##get the edge in the 'reticulation' matrix

    ##Get the edges
    edge_edge <- x$edge[edge_row,]
    ret_edge <- x$reticulation[ret_row,]

    ##swap edges
    x$edge[edge_row,] <- ret_edge
    x$reticulation[ret_row,] <- edge_edge

    ##Update edge lengths
    x$edge.length[edge_row]<-node_times[ret_edge[2]]-node_times[ret_edge[1]]

    ##update inheritance probabilities if present
    if(!is.null(x$inheritance)){
      x$inheritance[ret_row]<-1-x$inheritance[ret_row]
    }
  }
  return(x)
}
