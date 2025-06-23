
#' Swapping the minor edges of an evonet object
#'
#'
#' @title swap_hybrid_minor
#' @param  a evonet object
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


swap_hybrid_minor <-function(net,hybrid_nodes,node_times=NULL){
  if(is.null(node_times)){
    node_times<-node.depth.edgelength(net)
  }
  
  if(any(!(hybrid_nodes %in% net$reticulation[,2]))){
    error('You goofed. hybrid_nodes only takes hybrid nodes but was given a tree node. Look at net$reticulation[,2] for hybrid nodes')
  }
  
  for(hyb_nd in hybrid_nodes){
    ##get the edges that point to the hybrid node of interest
    edge_row<-which(net$edge[,2]==hyb_nd) ##get the edge in the 'edge' matrix
    ret_row <-which(net$reticulation[,2]==hyb_nd) ##get the edge in the 'reticulation' matrix
    
    ##Get the edges
    edge_edge<-net$edge[edge_row,]
    ret_edge<-net$reticulation[ret_row,]
    
    ##swap edges
    net$edge[edge_row,] <- ret_edge
    net$reticulation[ret_row,]<-edge_edge
    
    ##Update edge lengths
    net$edge.length[edge_row]<-node_times[ret_edge[2]]-node_times[ret_edge[1]]
    
    ##update inheritance probabilities if present
    if(!is.null(net$inheritance)){
      net$inheritance[ret_row]<-1-net$inheritance[ret_row]
    }
  }
  return(net)
}
