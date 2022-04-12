---
title: 'tanggle: an R package for the visualization of phylogenetic networks'
tags:
  - R
  - phylogenetics
  - networks
  - ggtree
  - ggplot2
authors:
  - name: Klaus Schliep^[Corresponding author]
    affiliation: 1
  - name: Marta Vidal-Garcia
    affiliation: 2
  - name: Leann Biancani
    affiliation: 3
  - name: L. Francisco Henao-Diaz
    affiliation: 4
  - name: Eren Ada
    affiliation: 5
  - name: Claudia Solis-Lemus^[Corresponding author]
    affiliation: 6
affiliations:
 - name: Institute of Environmental Biotechnology, Graz University of Technology, Viena
   index: 1
 - name: Department of Cell Biology and Anatomy, University of Calgary, Canada
   index: 2
 - name: University of Maryland College Park, Smithsonian National Museum of Natural History, USA
   index: 3
 - name: Department of Zoology and Biodiversity Research Center, University of British Columbia, Canada
   index: 4
 - name: Department of Biological Sciences, University of Rhode Island, USA
   index: 5
 - name: Wisconsin Institute for Discovery and Department of Plant Pathology, University of Wisconsin-Madison, USA
   index: 6
date: 12 April 2022
bibliography: tanggle.bib

---

# Summary

Phylogenetic trees depict evolutionary relationships among taxa. However, they are strictly bifurcating structures that do not take into account several evolutionary events such as horizontal gene transfer, hybridization, or introgression. Recently, there has been an increased production of phylogenetic network methods. Yet, there is limited visualization software to plot the phylogenetic networks. Here, we present the R package `tanggle`, a visualization package for phylogenetic networks. Our package extends the widely used
visualization package `ggtree`.
Our R package allows a variety of input data from DNA sequences to extended Newick format and it builds on the flexibility of `ggplot2` to manipulate colors and other plot characteristics. In
addition, our package allows for the inclusion of images and mapped morphological and geographical characteristics on the network. 
With the increased (and timely) pressure for reproducible and open-source research,
`tanggle` works towards script-based publication figures, as
opposed to figures created with design software such as Adobe Illustrator. 

# Statement of need

Recent years have seen an explosion of phylogenetic network inference methods
[@Than2008; @Huson2010; @Yu2012; @Gruenewald2013; @Yang2014; @Yu2014; @Solis-Lemus2016; @Wen2016; @2017netsymposium-blischak; @2017netsymposium-degnan; @Wen2018; @Zhang2018] which can be divided into two main classes: methods to
reconstruct implicit (or split) networks, and methods to reconstruct
explicit networks. Implicit networks are data-displayed objects that
include two (or more) options for every uncertain split and these options can be resolved by collapsing parallel edges. For example, in
Figure \autoref{fig:nettypesR}, there are two possible ways to resolve
the uncertain split: 1) by considering the longer edges of the rectangle
to be the correct split and collapsing the shorter edges (making the two spiders on the bottom sister
species) or 2) by considering the shorter edges of the rectangle to be
the correct split and collapsing the longer edges (making the two spiders on the right sister species).

Explicit networks, on the contrary, assign a biological mechanism to
every internal node, and thus, they are easier to interpret. For
example, in Figure \autoref{fig:nettypesL}, the hybrid node (where the
green arrow is pointing) corresponds to some reticulation event, and the
tree nodes correspond to speciation events. 
Hybrid edges are parametrized by an inheritance
probability ($\gamma$) that can provide information on the proportion
of genes that were transferred through the hybrid edge (green arrow).
For example, an estimated $\hat{\gamma} \approx 0.5$ could indicate a
hybridization event, whereas smaller estimated
$\hat{\gamma} \approx 0.01$ could refer to horizontal gene transfer.


![Explicit network. Internal nodes represent biological processes (e.g. tree nodes correspond to speciation events and hybrid nodes correspond to reticulation events). In addition, hybrid edges (green arrow) have a numerical parameter (inheritance probability) that represents the proportion of genes transferred through the arrow (here 17\%).\label{fig:nettypesL}](figures/explicit.pdf)

![Split network. Split configurations are represented as parallel lines which represent discordance in the input data (sequences or gene trees). It is not possible to know the specific source of the discordance (e.g. estimation error, incomplete lineage sorting (ILS) or gene flow).\label{fig:nettypesR}](figures/implicit.pdf)


Alongside the development of methods to reconstruct phylogenetic networks, multiple plotting options have appeared (see table below). Here, we introduce `tanggle`, a novel R package that extends `ggtree` [@ggtree] to plot split and explicit networks.
The package `ggtree` itself builds on the R package `ggplot2` [@wickham2011ggplot2] which implemented ideas from the influential book the "grammar of graphics" [@wilkinson1999] and as such, it provides a whole philosophy for composing sophisticated graphics. 


Our package `tanggle` builds on other existing popular R packages such as `ape` [@Paradis2019] and `phangorn` [@Schliep2010] using the functions and the data structures defined therein. The data structure for explicit networks (class `evonet`) is defined in the package `ape` along with functions to import and export these data as extended Newick format (see Appendix and Applications to real data). The class `evonet` extends the class `phylo`, which is the main data structure to store phylogenetic trees in R.  
Similarly, the data structure for split networks (`networx`) also extending the class `phylo` is defined in `phangorn`. `phangorn` provides functions to import or export networks in nexus format and different algorithms to infer split networks (e.g. consensus networks [@Holland2004], NeighborNet [@Bryant2004]).
Methods to plot implicit and explicit networks have been already available using base R graphics in the `ape` and `phangorn` packages. Base R graphics are often a good choice to create a plot, but it might not be as easy to extend the graphics as with `ggplot2` without writing a new function. Furthermore `phangorn` allows to plot interactive 3D split networks making use of the package `rgl` [@rgl].



| R class (S3) | function name | network |  plotting system | package |
| :---: | :---: | :---: | :---: | :---: |
| phylo | plot.phylo, plotTree | phylogenetic tree | base R | ape, phytools [@phytools] |
| phylo | ggtree | phylogenetic tree | ggplot2 | ggtree |
| evonet | plot.evonet  | explicit | base R | ape |
| evonet | ggevonet     | explicit | ggplot2 | tanggle |
| networx | ggsplitnet    | implicit | ggplot2 | tanggle |
| networx | plot.networx | implicit | base R, rgl | phangorn | 
| haploNet | plot.haploNet | haplo type | base R  | pegas [@Paradis2010] | 

# The `tanggle` package

The `tanggle` package extends the already widely used
visualization package `ggtree} [@ggtree] to implicit and
explicit phylogenetic networks. 

## Functions in `tanggle`

Currently we provide the following five functions in the `tanggle` package. Our package has two main functions: 
`ggsplitnet` that plots split networks, and `ggevonet` that plots explicit networks.
There are three further functions, which are mainly called internally by these two functions:  

1. `geom_splitnet` provides a geom object to represents a data representation for plotting split networks. It is called by `ggsplitnet(net)`, which is a shortcut for `ggplot(net) + geom_splitnet`.
2. `node_depth_evonet` can be used to assign edge length to an explicit networks. The function is called by `ggevonet` if no edge lengths are provided. 
3. `minimize_overlap` reduces reticulation lines crossing over in explicit networks. This function is by default called in the function `ggevonet` by the argument `min_crossing`.

Klaus 2do: “More focus in the paper on visualization strategies, and more examples exploring the nuts and bolts of tanggle (and for extension, of ggtree and ggplot2) to straightforwardly improve and automating the production of phylogenetic network figures in R”


## Installation steps

Claudia 2do: add the installation steps per reviews

## Basic examples

### Reconstructing and plotting a split network


There are three types of input data that our package takes to reconstruct
a split network:

1. Nexus file created with SplitsTree [@Huson1998] which is read with the `read.nexus.networx` function in the `phangorn` package [@Schliep2010].
2. Collection of gene trees in a nexus file, or in a text file with one row per gene tree in Newick format. These trees are read with the functions `read.tree` or `read.nexus`. A consensus split network is then computed using the `consensusNet` function in the `phangorn` package which implements the algorithm in @Holland2004.
3. Sequences in nexus, fasta or phylip format read with the `read.phyDat` function in `phangorn` package or the `read.dna` function in `ape` package [@Paradis2019]. Distance matrices are computed for specific models of evolution using the function `dist.ml` from `phangorn` package or the function `dist.dna` from `ape` package. From the distance matrix, a split network is reconstructed using the `neighborNet` function in the `phangorn` package which implements the algorithm in `Bryant2004`. In this case, there is the option to estimate branch lengths as well (which is not usually the case for split networks) with the function `splitsNetworks` in the `phangorn` package.


After any of these three steps, we have a `networkx` object that
represents a split network by two matrices: 1) standard edge matrix with
two columns for the two nodes connecting every edge (row), and 2) split
vector with integer identifier for every edge. Each entry in the split
vector represents the split for that edge. In split networks, there are
more edges than splits as two (or more) edges could represent the same
split (see Figure \autoref{fig:nettypesR} where parallel edges represent the
same split).

The `ggsplitnet` function takes a `networkx` object as
input and produces a plot. In Figure \autoref{splitplot}, we show the
split network corresponding to the yeast dataset in the
`phangorn` package produces with the code below:

```{r,eval=FALSE}
library(tanggle)
data(yeast, package="phangorn")
dm <- phangorn::dist.ml(yeast)
nnet <- phangorn::neighborNet(dm)
p <- ggsplitnet(nnet) + geom_tiplab2(size=5)
```

![Split network for yeast data.\label{splitplot}](figures/yeast.pdf)

### Plotting explicit networks

Algorithms to reconstruct explicit networks are more computationally
intensive than those to reconstruct split networks. Thus, to the best of
our knowledge, there are no R functions that would estimate an explicit
network from sequences or from gene trees. To obtain an explicit
network, users need to use existing tools like PhyloNet
[@Yu2014; @Wen2016; @Wen2018], BEAST [@Zhang2018] or SNaQ
[@Solis-Lemus2016; @phylonetworks]. These methods take sequences or gene trees as
input and estimate an explicit network which is represented in extended
Newick format (description in the Appendix). The explicit network is then read with the `ape`
function `read.evonet` and an `evonet` object is created.

The `evonet` object has two matrices: 1) standard edge matrix
with two columns for the two nodes connecting every edge (row), and 2)
hybrid edges matrix with two columns for the nodes connected by the
hybrid edge (row). In Figure \autoref{enet}, we show an explicit
network for a toy example in extended Newick format produced with the code to the left of the figure. Hybrid edges are drawn as dotted lines by default to distinguish them from the tree
edges.

```{r, eval=FALSE}
library(ape)
net <- "((a:2,(b:1)#H1:1):1,(#H1,c:1):2);"
enet <- read.evonet(text=net)
p2 <- ggevonet(enet) + geom_tiplab(size=5)
```

![Explicit network from extended Newick format.\label{enet}](figures/ex-explicit.pdf)

## Advanced examples

Marta 2do: “More elaborated toy examples, with details about the image formats required/permitted, how they should be named and stored (directory, folders), which are the resolutions recommended, etc. An example where you include images of graphs (e.g., pie charts or any other type) generated within R, or examples where you combine R-generated graphs with external images (such as the combinations of frogs and pie-charts in Fig. 4) would be fine”. We can use the ggtree tutorial as guideline for some potential examples


# Applications to real data

## Split network for the _Neobatrachus_ frog genus

_Neobatrachus_ is a genus of endemic Australian burrowing frogs, widely distributed across mainland Australia. Even though they can be identified based on external morphology [@Mahony1986; @Roberts1991], all _Neobatrachus_ species display relatively similar body shape patterns, mostly due to selective pressures imposed by the arid and semi-arid environments they inhabit and their burrowing behaviour [@Vidal-Garcia2014; @Roberts2010]. 
This genus has a particularly interesting and complex evolutionary story, as it comprises six diploid (_N. albipes_, _N. fulvus_, _N. pelobatoides_, _N. pictus_, _N. sutor_, _N. wilsmorei_; $2n = 24$) and three tetraploid (_N. aquilonius_, _N. kunapalari_, _N. sudellae_; $4n = 48$) species [@Roberts2010; @Novikova2020]. The authors in @Novikova2020 reconstructed their evolutionary history by analysing nucleotide sequence data for 439 loci in 87 individuals covering all nine species, and observed non-bifurcating evolutionary relationships. They suggest that this complex evolutionary history could be due to rapid speciation events (particularly in diploid species) or shared genetic variation among species due to incomplete lineage sorting [@Novikova2020]. Here, we used the set of 349 gene trees generated in @Novikova2020 to estimate and plot the consensus split network for the _Neobatrachus_ species with `tanggle` in Figure \autoref{fig:frognet} (see the Appendix for the code to produce the plot).

![Split network for the _Neobatrachus_ genus frogs. Frog pictures obtained from @Novikova2020, with permission from the authors, and added to the plot using `tanggle`.\label{fig:frognet}](figures/frogs_amended.pdf)

## Explicit network for the _Mysticeti_ whales genus

Baleen whales (_Mysticeti_) represent the largest animals on Earth.
Their evolutionary history is complex given the lack of physical
barriers in the ocean, and the fact that they hybridize frequently.
@Arnason2018 reconstructed the network-like
evolutionary history from a variety of methods using 34,192 genome
fragments, each with 20 kbp long, and estimated gene trees with RAxML
[@Stamatakis2014] assuming a GTR nodel with gamma-distributed rate
variation with invariable sites. Here, we used SNaQ
[@Solis-Lemus2016] to reconstruct the whales network with two
hybridization events, one of which was already reported in the original
paper (from blue whale to gray whale). SNaQ produces the network in
extended Newick format, which we can read into R and then plot with
`tanggle` (Figure \autoref{fig:whalesnet}). Details on the reconstruction of the network as well as the plotting code can be found in the Appendix.

![Explicit network for whales in @Arnason2018. Hybrid edges are represented as dotted lines. Whale images and squared taxon names are automatically added using `tanggle` functionalities. Silhouettes for each whale species were obtained from PhyloPic (http://phylopic.org), and coloured blue for visualization purposes.\label{fig:whalesnet}](figures/whales.pdf)

# Discussion

Evidence of reticulate evolution across multiple domains in the Tree of
Life have caused great attention to the area of phylogenetic networks from new methods to reconstruct networks to network as well as extensive research on the evolution of traits
on phylogenetic networks via comparative methods
[@2017netsymposium-bastide; @2017netsymposium-jhwueng].

Our R package `tanggle` is not the first software for the
visualization of networks. Dendroscope [@Huson2007] and SplitsTree
[@Huson1998] are the pioneers of network visualization, with
multiple functionalities being added constantly. In R, the widely used
`ape` package [@Paradis2019] has extended the base `plot`
function to visualize explicit networks, but it does not have functionalities still to
plot implicit networks yet. In Julia, `PhyloPlots` is the
accompanying package of `PhyloNetworks` [@phylonetworks]
which is able to plot rooted explicit networks, but no implicit networks
or unrooted networks.

The `tanggle` package is complementing the `ape` plotting
function in two main ways: 1) by allowing to plot implicit networks, and 2)
by building on the flexibility of `ggtree` [@ggtree] (which in
turns builds on the flexibility of `ggplot2` [@ggplot2]) to
manipulate the structure, nodes and edges in the network. With the
increased (and timely) pressure for reproducible research, limiting the
use of design software like Adobe Illustrator to produce publication
figures is paramount. The R package `ggtree` and now
`tanggle` work towards script-based publication figures that can
incorporate colors, images and mapped morphological and geographical
characteristics. Furthermore, `ggtree` and thus `tanggle`
are based on underlying data frames to represent the phylogenies, which
allow a lot of flexibility on mapped traits, colors and linetypes.



## Author contributions
KS created the R package, programmed all the main functions and functionalities, and compiled the R package. MVG, LB, FHD, EA and CSL aided in further developing the R package by contributing code and documentation. MVG and CSL analyzed the frog data. MVG, LB and CSL wrote the vignette. CSL analyzed the whales data. CSL wrote the initial draft of the manuscript, and KS and MVG contributed to writing by editing and reviewing subsequent drafts of the manuscript. All authors approved the final manuscript draft for submission.


## Acknowledgements
The R package `tanggle` was partly developed during the Nantucket DevelopR workshop on the field station of UMass Boston. KS and the workshop were supported by NSF DBI-1759940 to Liam Revell and NSF DBI-2113424. This work was also supported by Department of Energy DE-SC0021016 (CSL) and the NSF DEB-2144367 (CSL). MVG was supported by an ACHRI Postdoctoral Fellowship and an Alberta Innovates Postdoctoral Fellowship in Health Innovation. We are grateful to I.G. Brennan for providing the \emph{Neobatrachus} gene trees, and for the original frog pictures from M. Mahony and S. Mahony. The authors declare no conflicts of interest.

## Data availability
Data was not generated as part of this research. The frogs and whales datasets are publicly available or available by request per the original publications. All code is publicly available on [Bioconductor](https://bioconductor.org/packages/tanggle) and the development version on [GitHub](https://github.com/KlausVigo/tanggle).

# Appendix

## Extended Newick format

Traditionally, phylogenetic trees are written in Newick format. To
represent explicit network objects, we use the extended Newick format
[@Morin2006; @Than2008; @Cardona2008b; @Cardona2008c]. This format uses
the concept of minor hybrid edges (edges with $\gamma < 0.5$) and
major hybrid edges ($\gamma > 0.5$). By default, we detach the minor
hybrid edge at each hybrid node to write the extended Newick description
of a network as we would for a tree with the repeated label of the
hybrid node (`\#H1' in Figure \autoref{fig:netNewick}). This description can
include edge information formatted as $:length:support:\gamma$.

![Extended Newick format for networks. The network is written as a tree with two nodes having the same label \#H1: (((A,(B)\#H1),(C,\#H1)),D);\label{fig:netNewick}](figures/extendedNewick.pdf)


For example, the parenthetical format of the network in Figure
\autoref{fig:netNewick} can include $\gamma$ values:
`(((A,(B)\#H1:::0.8),(C,\#H1:::0.2)),D);`,
which are written after three colons, because there is no information
about branch length nor support for the hybrid edges. Other internal
edges (tree edges) have information about branch lengths, which follow
just one colon. These tree edges do not have information about
$\gamma$, because all tree edges have $\gamma=1$.


## Plotting code for the split network for the _Neobatrachus_ frog genus

The collection of gene trees is stored in the `Neobatrachus.trees` file, and there is a `figures` folder containing the images of the frogs as PNG files with the same filenames as the taxon names in the network. For example, the frog image for `N_kunapalari` should be in the file `N_kunapalari.png` inside the `figures` folder. 

A `consensusNet` [@Holland2004] is a generalization of a (majority) consensus tree. Instead of only representing splits (taxon bipartitions) occurring in at least 50\% of the trees in a bootstrap or MCMC sample one can use a lower threshold (in the example below 15\%) and explore competing splits. 

```{r, eval=FALSE}
library(tanggle)
library(phangorn)
trees = read.tree("Neobatrachus.trees")
frognet = consensusNet(trees, 0.15)
p = ggsplitnet(frognet) +
geom_tiplab(geom=’label’, offset=1.5, hjust=1.2) +
geom_tiplab(aes(image=paste0(’figures/’,label, ’.png’)), geom="image",offset=0)
```

## Reconstruction and plotting of the explicit network for the _Mysticeti_

The 34,192 gene trees had already been estimated by @Arnason2018 and stored in the file `gene-trees.tre`. Information on how to obtain the gene trees can be found in the original publication. 
Here, We first ran ASTRAL (v5.6.2) to estimate a starting tree for SNaQ:
```
java -jar astral.5.6.2.jar -i gene-trees.tre -o whales.tre
```

Then, we read the starting tree `whales.tre` and the gene trees `gene-trees.tre` in Julia, and load the `PhyloNetworks` [@phylonetworks] package to estimate an explicit network using `SNaQ` [@Solis-Lemus2016]. More information about this inference process can be found in the `PhyloNetworks` [wiki page](\url{https://github.com/crsl4/PhyloNetworks.jl/wiki).

```{julia, eval=FALSE}
using PhyloNetworks
cf = readTrees2CF("gene-trees.tre")
net0 = readTopology("whales.tre")
net1 = snaq!(net0, cf, hmax=1, filename="snaq1", runs=10)
net2 = snaq!(net1, cf, hmax=2, filename="snaq2", runs=10)
```

The extended Newick format of the network object `net2` is stored in the `snaq2.out` file, which can later be read into R for plotting:

```{r, eval=FALSE}
library(tanggle)
library(ggimage)
library(viridis)
whalesnet = read.evonet("snaq2.out")
p = ggevonet(whalesnet,layout ="slanted") +
xlim(0,8) + ylim(0,15) +
geom_tiplab(geom=’label’, offset=0.5) +
geom_tiplab(aes(image=paste0(’figures/’,label, ’.png’)), geom="image", offset=0)
```

As for the frogs data, there is a `figures` folder containing the images of the frogs as PNG files with the same filenames as the taxon names in the network.

# References