# swap_hybrid_minor

Swapping the minor edges of an evonet object

## Usage

``` r
swap_hybrid_minor(x, hybrid_nodes, node_times = NULL)
```

## Arguments

- x:

  evonet object

- hybrid_nodes:

  a vector of hybrid nodes to have their minor edges swapped

- node_times:

  an optional argument with node times

## Value

network

## Examples

``` r
(enet <- ape::read.evonet(text='((a:2,(b:1)#H1:1):1,(#H1,c:1):2);'))
#> 
#>     Evolutionary network with 1 reticulation
#> 
#>                --- Base tree ---
#> Phylogenetic tree with 3 tips and 4 internal nodes.
#> 
#> Tip labels:
#>   a, b, c
#> Node labels:
#>   , , #H1, 
#> 
#> Rooted; includes branch length(s).
ggevonet(enet) + geom_tiplab()

swapped_enet<-swap_hybrid_minor(enet,6)
ggevonet(swapped_enet) + geom_tiplab()
```
