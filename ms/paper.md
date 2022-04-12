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
`[@Than2008; @Huson2010; @Yu2012; @Gruenewald2013; @Yang2014; @Yu2014; @Solis-Lemus2016; @Wen2016; @2017netsymposium-blischak; @2017netsymposium-degnan; @Wen2018; @Zhang2018]` which can be divided into two main classes: methods to
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


Alongside the development of methods to reconstruct phylogenetic networks, multiple plotting options have appeared (see table below). Here, we introduce `tanggle`, a novel R package that extends `ggtree` `[@ggtree]` to plot split and explicit networks.
The package `ggtree` itself builds on the R package `ggplot2` `[@wickham2011ggplot2]` which implemented ideas from the influential book the "grammar of graphics" `[@wilkinson1999]` and as such, it provides a whole philosophy for composing sophisticated graphics. 


Our package `tanggle` builds on other existing popular R packages such as `ape` `[@Paradis2019]` and `phangorn` `[@Schliep2010]` using the functions and the data structures defined therein. The data structure for explicit networks (class `evonet`) is defined in the package `ape` along with functions to import and export these data as extended Newick format (see Appendix and Applications to real data). The class `evonet` extends the class `phylo`, which is the main data structure to store phylogenetic trees in R.  
Similarly, the data structure for split networks (`networx`) also extending the class `phylo` is defined in `phangorn`. `phangorn` provides functions to import or export networks in nexus format and different algorithms to infer split networks (e.g. consensus networks `[@Holland2004]`, NeighborNet `[@Bryant2004]`).
Methods to plot implicit and explicit networks have been already available using base R graphics in the `ape` and `phangorn` packages. Base R graphics are often a good choice to create a plot, but it might not be as easy to extend the graphics as with `ggplot2` without writing a new function. Furthermore `phangorn` allows to plot interactive 3D split networks making use of the package `rgl` `[@rgl]`.



| R class (S3) | function name | network |  plotting system | package |
| :---: | :---: | :---: | :---: | :---: |
| phylo | plot.phylo, plotTree | phylogenetic tree | base R | ape, phytools `[@phytools]`|
| phylo | ggtree | phylogenetic tree | ggplot2 | ggtree |
| evonet | plot.evonet  | explicit | base R | ape |
| evonet | ggevonet     | explicit | ggplot2 | tanggle |
| networx | ggsplitnet    | implicit | ggplot2 | tanggle |
| networx | plot.networx | implicit | base R, rgl | phangorn | 
| haploNet | plot.haploNet | haplo type | base R  | pegas `[@Paradis2010]` | 

