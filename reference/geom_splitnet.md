# geom_splitnet

add splitnet layer

## Usage

``` r
geom_splitnet(layout = "slanted", ...)
```

## Arguments

- layout:

  one of 'rectangular', 'slanted', 'circular', 'radial' or 'unrooted'

- ...:

  additional parameter

## Value

splitnet layer

## Author

Klaus Schliep

## Examples

``` r
data(yeast, package='phangorn')
dm <- phangorn::dist.ml(yeast)
nnet <- phangorn::neighborNet(dm)
ggplot(nnet, aes(x, y))  + geom_splitnet() + theme_tree()
```
