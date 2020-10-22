
# Wrapper function which source all others and the main computatinal complete.master.model.epidemic.metrics.combined.features.cov script 
# which is the master model

# Define directory

work.dir <- "/home/dniyukuri/lustre/benchmark_master_model" # on CHPC


# work.dir <- "/home/david/benchmark_master_model" # on laptop


setwd(paste0(work.dir))

# work.dir <- "/user/data/gent/vsc400/vsc40070/phylo/" # on cluster


pacman::p_load(snow, parallel, RSimpactCyan, RSimpactHelper, ape, Rsamtools)


wrapper.benchmark.master.model <- function(inputvector = inputvector){
  
  
  work.dir <- "/home/dniyukuri/lustre/benchmark_master_model" # on CHPC
  
  
  # work.dir <- "/home/david/benchmark_master_model" # on laptop
  
  
  setwd(paste0(work.dir))
  
  library(RSimpactCyan)
  library(RSimpactHelper)
  library(Rcpp)
  library(ape)
  library(expoTree)
  library(data.table)
  # library(readr)
  library(phangorn)
  library(lme4)
  library(nlme)
  library(dplyr)
  library(adephylo)
  library(treedater)
  library(geiger)
  library(picante)
  library(igraph)
  library(phyloTop)
  library(phytools)
  library(Rsamtools)
  library(robustbase)
  library(intergraph)
  library(lubridate)
  library(tidyr)
  
  
  # source("/home/dniyukuri/lustre/benchmark_master_model/advanced.transmission.network.builder.R")
  
  source("/home/dniyukuri/lustre/benchmark_master_model/needed.functions.RSimpactHelp.R")
  
  source("/home/dniyukuri/lustre/benchmark_master_model/complete.master.epidemic.metrics.R")
  
  source("/home/dniyukuri/lustre/benchmark_master_model/compute.summary.statistics.hiv.epi.behavior.R")
  
  source("/home/dniyukuri/lustre/benchmark_master_model/compute.summary.statistics.phylo.MCAR.R")
  
  # source("/home/dniyukuri/lustre/benchmark_master_model/compute.summary.statistics.phylo.MAR.R")
  
  source("/home/dniyukuri/lustre/benchmark_master_model/complete.master.model.epidemic.metrics.combined.features.cov.R")
  
  # 
  # 
  
  # source("/home/david/benchmark_master_model/advanced.transmission.network.builder.R") # ok
  
  # source("/home/david/benchmark_master_model/needed.functions.RSimpactHelp.R") # ok
  # 
  # source("/home/david/benchmark_master_model/complete.master.epidemic.metrics.R") # ok
  # 
  # source("/home/david/benchmark_master_model/compute.summary.statistics.hiv.epi.behavior.R") # ok
  # 
  # source("/home/david/benchmark_master_model/compute.summary.statistics.phylo.MCAR.R") # ok
  # 
  # # source("/home/david/benchmark_master_model/compute.summary.statistics.phylo.MAR.R") # ok
  # 
  # source("/home/david/benchmark_master_model/complete.master.model.epidemic.metrics.combined.features.cov.R") # ok
  
  
  # 
  
  results.f <- tryCatch(complete.master.model.epidemic.metrics.combined.features.cov(inputvector = inputvector),
                        error=function(e) return(rep(NA, 126))) # 69 MM, 39 Class, and 18 Phylo
  
                        
                        
  # results.f <- complete.master.model.epidemic.metrics.combined.features.cov(inputvector = inputvector)
  
  
  return(results.f)
  
  
}



reps <- 1120


inputvector <- c(-0.52, -0.05, 2, 10, 5, 0.25, -0.3, -0.1,
                 -1, -90, 0.5, 0.048, -0.14) 



inputmatrix <- matrix(rep(inputvector, reps), byrow = TRUE, nrow = reps)


epi.mm.stats <- simpact.parallel(model = wrapper.benchmark.master.model,
                                 actual.input.matrix = inputmatrix,
                                 seed_count = 1,
                                 n_cluster = 56)

write.csv(epi.mm.stats, file = "results.benchmark.epi.mm.stats_seed_1.csv")


# width = 16, height = 12, units = "cm"
# width = 25, height = 15, units = "cm"

