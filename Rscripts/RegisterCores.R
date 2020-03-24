# register cores
library(foreach)
library(doMC)

.minCores <- 3
.perc <- 0.5

.cores <- detectCores()
.cores <- max(floor(.cores*.perc), .minCores)
.cores <- min(.cores, detectCores())
registerDoMC(cores = .cores)
