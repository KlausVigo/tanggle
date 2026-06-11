# minimize_overlap reduces reticulation lines crossing over in plots

minimize_overlap reduces reticulation lines crossing over in plots

## Usage

``` r
minimize_overlap(x)
```

## Arguments

- x:

  Tree of class 'evonet'

## Value

A Tree with rotated nodes of class 'evonet'

## Author

L. Francisco Henao Diaz

## Examples

``` r
fishnet <- ape::read.evonet(text='(Xalvarezi,Xmayae,((Xsignum,((Xmonticolus,
(Xclemenciae_F2,#H25)),(((((((((Xgordoni,Xmeyeri),Xcouchianus),Xvariatus),
Xevelynae),(Xxiphidium,#H24)),Xmilleri),Xandersi),Xmaculatus),(((Xmontezumae,
(Xcortezi,(Xbirchmanni_GARC,Xmalinche_CHIC2))),((Xnigrensis,Xmultilineatus),
(Xpygmaeus,Xcontinens))))#H24))),(Xhellerii)#H25));')
fishnet$edge.length <- NULL
new_tre <- minimize_overlap(fishnet)

par(mfrow=c(1,2))
ggevonet(fishnet, min_crossing = FALSE)

ggevonet(new_tre)


net2 <- ape::read.evonet(text='(15,(1,((14,(#H1,(((12,13),(11,#H3)),(7,
    ((10)#H3,(8,9)))))),((((2,3))#H2,(6,(5,(#H2,4)))))#H1)));')
# Cui et al. 2013 Evol.
new_net2 <- minimize_overlap(net2)
ggevonet(net2, min_crossing = FALSE)

ggevonet(new_net2)

```
