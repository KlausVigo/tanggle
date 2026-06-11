# Depth of Nodes

These functions return the depths or heights of nodes and tips.

## Usage

``` r
node_depth_evonet(x, ...)
```

## Arguments

- x:

  an object of class 'evonet'

- ...:

  Further arguments passed to or from other methods.

## Value

a vector with the depth of the nodes

## See also

[`node.depth`](https://rdrr.io/pkg/ape/man/node.depth.html)

## Examples

``` r
z <- ape::read.evonet(text = '((1,((2,(3,(4)Y#H1)g)e,
(((Y#H1, 5)h,6)f)X#H2)c)a,((X#H2,7)d,8)b)r;')
nd <- node_depth_evonet(z)
z$edge.length <- nd[z$edge[,1]] - nd[z$edge[,2]]
ggevonet(z)

```
