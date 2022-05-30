## Installation
#install.packages("devtools")
#devtools::install_github("jjustison/SiPhyNetwork")

library(SiPhyNetwork)
set.seed(530246)
networks <- sim.bdh.taxa.ssa(n = 15, ## number of tips
                             numbsim = 100, ## number of networks to simulate
                             lambda = 0.9, ## speciation rate
                             mu = 0.0,## extinction rate
                             nu = 0.04,  ## the larger, the more complex the networks
                             hybprops = c(1, 1, 1), ## type of hybridization
                             hyb.inher.fxn = make.beta.draw(1, 1), 
                             frac = 1,
                             mrca = FALSE, 
                             complete = TRUE, 
                             stochsampling = FALSE, 
                             hyb.rate.fxn = NULL, 
                             trait.model = NULL)
plot(networks[[1]])
