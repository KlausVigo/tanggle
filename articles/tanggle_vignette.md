# \*\*\*tanggle\*\*\*: Visualization of phylogenetic networks in a \*ggplot2\* framework

## Introduction

Here we present a vignette for the R package ***tanggle***, and provide
an overview of its functions and their usage. ***Tanggle*** extends the
*ggtree* R package (Yu et al. 2017) to allow for the visualization of
several types of phylogenetic networks using the *ggplot2* (Wickham
2016) syntax. More specifically, *tanggle* contains functions to allow
the user to effectively plot: (1) split (i.e. implicit) networks
(unrooted, undirected) and (2) explicit networks (rooted, directed) with
reticulations. It offers an alternative to the plot functions already
available in *ape* (Paradis and Schliep 2018) and *phangorn* (Schliep
2011).

## List of functions

| Function name | Brief description |
|:---|:---|
| `geom_splitnet` | Adds a *splitnet* layer to a ggplot, to combine visualising data and the network |
| `ggevonet` | Plots an explicit network from a *phylo* object |
| `ggsplitnet` | Plots an implicit network from a *phylo* object |
| `minimize_overlap` | Reduces the number of reticulation lines crossing over in the plot |
| `node_depth_evonet` | Returns the depths or heights of nodes and tips in the phylogenetic network |

## Getting started

Install the package from Bioconductor directly:

``` r

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("tanggle")
```

Or install the development version of the package from
[Github](https://github.com/KlausVigo/tanggle).

``` r

if (!requireNamespace("remotes", quietly=TRUE))
    install.packages("remotes")
remotes::install_github("KlausVigo/tanggle")
```

If you need to install *ggtree* from github:

``` r

remotes::install_github("YuLab-SMU/ggtree")
```

And load all the libraries:

``` r

library(tanggle)
library(phangorn)
library(ggtree)
```

------------------------------------------------------------------------

## Split Networks

Split networks are data-display objects which allow for the definition
of 2 (or more) options for non-compatible splits. Split networks are
most often used to visualize consensus networks (Holland et al. 2004) or
neighbor-nets (Bryant and Moulton 2004). This can be done either by
using the `consensusNet` or `neighbor-net` functions in *phangorn*
(Schliep 2011) or by importing nexus files from SplitsTree (Huson and
Bryant 2006).

### Data Types

*tanggle* accepts three forms of input data for split networks. The
following input options all generate a *networx* object for plotting.

- Nexus file created with SplitsTree (Huson and Bryant 2006) and read
  with the `read.nexus.network` function in *phangorn* (Schliep 2011).

- Read in a split network in nexus format:

``` r

fdir <- system.file("extdata/trees", package = "phangorn")
Nnet <- phangorn::read.nexus.networx(file.path(fdir,"woodmouse.nxs"))
```

2.  A collection of gene trees (e.g.~from RAxML (Stamatakis 2014)) in
    one of the following formats:
    - Nexus file read with the function `read.nexus`
    - Text file in Newick format (one gene tree per line) read with the
      function `read.tree` A consensus split network is then computed
      using the function `consensusNet` in *phangorn* (Schliep 2011).

- Sequences in nexus, fasta or phylip format, read with the function
  `read.phyDat` in *phangorn* (Schliep 2011) or the function `read.dna`
  in *ape* (Paradis and Schliep 2018). Distances matrices are then
  computed for specific models of evolution using the function `dist.ml`
  in *phangorn* (Schliep 2011) or `dist.dna` in *ape* (Paradis and
  Schliep 2018). From the distance matrix, a split network is
  reconstructed using the function `neighborNet` in *phangorn* (Schliep
  2011). ***Optional***: branch lengths may be estimated using the
  function `splitsNetworks` in *phangorn* (Schliep 2011).

### Plotting a Split Network:

We can plot the network with the default options:

``` r

p <- ggsplitnet(Nnet) + geom_tiplab2()
p
```

![](tanggle_vignette_files/figure-html/plot_splitnet-1.png)

When we can set the limits for the x and y axis so that the labels are
readable.

``` r

p <- p + xlim(-0.019, .003) + ylim(-.01,.012) 
p
```

![](tanggle_vignette_files/figure-html/nicer_plot-1.png)

You can rename tip labels. Here we changed the names to species from 1
to 15:

``` r

Nnet$translate$label <- seq_along(Nnet$tip.label)
```

We can include the tip labels with `geom_tiplab2`, and customize some of
the options. For example, here the tip labels are in blue and both in
bold and italics, and we show the internal nodes in green:

``` r

ggsplitnet(Nnet) + geom_tiplab2(col = "blue", font = 4, hjust = -0.15) + 
    geom_nodepoint(col = "green", size = 0.25)
```

![](tanggle_vignette_files/figure-html/color_plot-1.png)

Nodes can also be annotated with `geom_point`.

``` r

ggsplitnet(Nnet) + geom_point(aes(shape=isTip, color=isTip), size=2)
```

![](tanggle_vignette_files/figure-html/geom_point-1.png)

### Plotting Explicit Networks

The function `ggevonet` plots explicit networks (phylogenetic trees with
reticulations). A recent addition to *ape* (Paradis and Schliep 2018)
made it possible to read in trees in extended newick format (Cardona et
al. 2008).

Read in an explicit network (example from Fig. 2 in Cardona et
al. 2008):

``` r

z <- read.evonet(text = "((1,((2,(3,(4)Y#H1)g)e,(((Y#H1,5)h,6)f)X#H2)c)a,
                            ((X#H2,7)d,8)b)r;")
```

Plot an explicit network:

``` r

ggevonet(z, layout = "rectangular") + geom_tiplab() + geom_nodelab()
```

![](tanggle_vignette_files/figure-html/plot_evonet-1.png)

``` r

p <- ggevonet(z, layout = "slanted") + geom_tiplab() + geom_nodelab()
p + geom_tiplab(size=3, color="purple")
```

![](tanggle_vignette_files/figure-html/plot_evonet-2.png)

``` r

p + geom_nodepoint(color="#b5e521", alpha=1/4, size=10)
```

![](tanggle_vignette_files/figure-html/plot_evonet-3.png)

## Summary

This vignette illustrates all the functions in the R package
***tanggle***, and provides some examples on how to plot both explicit
and implicit networks. The split network plots should take most of the
functions compatible with unrooted trees in ggtree. The layout options
for explicit network plots are *rectangular* or *slanted*.

## Session info

``` r

sessionInfo()
#> R version 4.6.0 (2026-04-24)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.4 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] phangorn_2.12.1  ape_5.8-1        tanggle_1.19.2   ggtree_4.2.0    
#> [5] ggplot2_4.0.3    BiocStyle_2.40.0
#> 
#> loaded via a namespace (and not attached):
#>  [1] fastmatch_1.1-8         gtable_0.3.6            xfun_0.58              
#>  [4] bslib_0.11.0            htmlwidgets_1.6.4       lattice_0.22-9         
#>  [7] quadprog_1.5-8          vctrs_0.7.3             tools_4.6.0            
#> [10] generics_0.1.4          yulab.utils_0.2.4       parallel_4.6.0         
#> [13] tibble_3.3.1            pkgconfig_2.0.3         Matrix_1.7-5           
#> [16] ggplotify_0.1.3         RColorBrewer_1.1-3      S7_0.2.2               
#> [19] desc_1.4.3              lifecycle_1.0.5         compiler_4.6.0         
#> [22] farver_2.1.2            treeio_1.36.1           textshaping_1.0.5      
#> [25] codetools_0.2-20        ggfun_0.2.0             fontquiver_0.2.1       
#> [28] fontLiberation_0.1.0    htmltools_0.5.9         sass_0.4.10            
#> [31] yaml_2.3.12             lazyeval_0.2.3          pillar_1.11.1          
#> [34] pkgdown_2.2.0           jquerylib_0.1.4         tidyr_1.3.2            
#> [37] MASS_7.3-65             cachem_1.1.0            nlme_3.1-169           
#> [40] fontBitstreamVera_0.1.1 tidyselect_1.2.1        aplot_0.2.9            
#> [43] digest_0.6.39           dplyr_1.2.1             purrr_1.2.2            
#> [46] bookdown_0.46           labeling_0.4.3          fastmap_1.2.0          
#> [49] grid_4.6.0              cli_3.6.6               magrittr_2.0.5         
#> [52] patchwork_1.3.2         withr_3.0.2             gdtools_0.5.1          
#> [55] scales_1.4.0            rappdirs_0.3.4          rmarkdown_2.31         
#> [58] igraph_2.3.2            otel_0.2.0              ragg_1.5.2             
#> [61] evaluate_1.0.5          knitr_1.51              gridGraphics_0.5-1     
#> [64] rlang_1.2.0             ggiraph_0.9.6           Rcpp_1.1.1-1.1         
#> [67] glue_1.8.1              tidytree_0.4.7          BiocManager_1.30.27    
#> [70] jsonlite_2.0.0          R6_2.6.1                systemfonts_1.3.2      
#> [73] fs_2.1.0
```

## References

Bryant, David, and Vincent Moulton. 2004. “Neighbor-Net: An
Agglomerative Method for the Construction of Phylogenetic Networks.”
*Molecular Biology and Evolution* 21 (2): 255–65.
<https://doi.org/10.1093/molbev/msh018>.

Cardona, Gabriel, Francesc Rosselló, and Gabriel Valiente. 2008.
“Extended Newick: It Is Time for a Standard Representation of
Phylogenetic Networks.” *BMC Bioinformatics* 9 (1): 532.
<https://doi.org/10.1186/1471-2105-9-532>.

Holland, Barbara R., Katharina T. Huber, Vincent Moulton, and Peter J.
Lockhart. 2004. “Using Consensus Networks to Visualize Contradictory
Evidence for Species Phylogeny.” *Molecular Biology and Evolution* 21
(7): 1459–61. <https://doi.org/10.1093/molbev/msh145>.

Huson, D. H., and D. Bryant. 2006. “Application of Phylogenetic Networks
in Evolutionary Studies.” *Molecular Biology and Evolution* 23 (2):
254–67.

Paradis, Emmanuel, and Klaus Schliep. 2018. “Ape 5.0: An Environment for
Modern Phylogenetics and Evolutionary Analyses in r.” *Bioinformatics*
35 (3): 526–28.

Schliep, Klaus Peter. 2011. “Phangorn: Phylogenetic Analysis in R.”
*Bioinformatics* 27 (4): 592–93.
<https://doi.org/10.1093/bioinformatics/btq706>.

Stamatakis, A. 2014. “RAxML Version 8: A Tool for Phylogenetic Analysis
and Post-Analysis of Large Phylogenies.” *Bioinformatics* 30 (9):
1312–13.

Wickham, Hadley. 2016. *Ggplot2: Elegant Graphics for Data Analysis*.
Springer-Verlag New York. <https://ggplot2.tidyverse.org>.

Yu, Guangchuang, David Smith, Huachen Zhu, Yi Guan, and Tommy Tsan-Yuk
Lam. 2017. “Ggtree: An r Package for Visualization and Annotation of
Phylogenetic Trees with Their Covariates and Other Associated Data.”
*Methods in Ecology and Evolution* 8: 28–36.
<https://doi.org/10.1111/2041-210X.12628>.
