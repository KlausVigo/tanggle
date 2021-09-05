library(tanggle)
library(ape)

enet <- ape::read.evonet(text="((a:2,(b:1)#H1:1):1,(#H1,c:1):2);")

expect_true(inherits(enet, "evonet"))
min_net <- minimize_overlap(enet)
expect_true(inherits(min_net, "evonet"))

nd <- node_depth_evonet(enet)

fe <- tanggle:::fortify.evonet(enet)
expect_true(is(fe, "data.frame"))


data(yeast, package="phangorn")
dm <- phangorn::dist.ml(yeast)
nnet <- phangorn::neighborNet(dm)

fspl <- tanggle:::fortify.networx(nnet)
expect_equal(class(fspl), "data.frame")

