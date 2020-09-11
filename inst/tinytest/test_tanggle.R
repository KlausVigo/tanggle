library(tanggle)
library(ape)

enet <- ape::read.evonet(text="((a:2,(b:1)#H1:1):1,(#H1,c:1):2);")
# ggevonet(enet) + geom_tiplab()

expect_true(inherits(enet, "evonet"))
min_net <- minimize_overlap(enet)
expect_true(inherits(min_net, "evonet"))

nd <- node.depth.evonet(enet)

fe <- tanggle:::fortify.evonet(enet)
expect_equal(class(fe), "data.frame")


data(yeast, package="phangorn")
dm <- phangorn::dist.ml(yeast)
nnet <- phangorn::neighborNet(dm)
# ggsplitnet(nnet) + geom_tiplab2()

fspl <- tanggle:::fortify.networx(nnet)
expect_equal(class(fspl), "data.frame")

