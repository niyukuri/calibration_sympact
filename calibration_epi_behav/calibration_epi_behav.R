
# Loading libraries

library(RSimpactCyan)
library(RSimpactHelper)
library(magrittr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(EasyABC)
library(robustbase)
library(snow)
library(parallel)


# Define directory

work.dir <-  "/home/dniyukuri/lustre/calibration_epi_behav"  # on CHPC

# work.dir <- "/home/david/benchmark_master_model/calibration_epi_behav" # on laptop


setwd(paste0(work.dir))


# df <- read.csv("/home/david/benchmark_master_model/calibration_epi_behav/results.benchmark.epi.mm.stats_seed_1.csv")


df <- read.csv("/home/dniyukuri/lustre/calibration_epi_behav/results.benchmark.epi.mm.stats_seed_1.csv")


benchmark_features_hiv_behavior_data <- dplyr::select(df, 
                                                      
                                                      Pop.growthrate,
                                                      
                                                      prev.m.15.24,                   
                                                      prev.m.25.29,  prev.m.30.34, prev.m.35.39,                   
                                                      prev.m.40.44,  prev.m.45.49, 
                                                      
                                                      prev.w.15.24,                   
                                                      prev.w.25.29,  prev.w.30.34, prev.w.35.39,                   
                                                      prev.w.40.44,  prev.w.45.49, 
                                                      
                                                      incid.m.15.24,                  
                                                      incid.m.25.29, incid.m.30.34, incid.m.35.39,                  
                                                      incid.m.40.44, incid.m.45.49, 
                                                      
                                                      incid.w.15.24,                  
                                                      incid.w.25.29, incid.w.30.34, incid.w.35.39,                  
                                                      incid.w.40.44, incid.w.45.49,
                                                      
                                                      pp.cp.6months.male.rels,
                                                      
                                                      relsperpersonperyear, agegap.mean, agegap.med, agegap.sd,
                                                      ART.33, ART.34, ART.35, ART.36, ART.37, ART.38, ART.39, ART.40, 
                                                      VL.suppr.)



benchmark_features_hiv_behavior_data <- na.omit(benchmark_features_hiv_behavior_data)

mean.classic.features <- sapply(benchmark_features_hiv_behavior_data, mean)


# I. one simpact model calibrated to classic summary statistics (full)
######################################################################


# Priors

inputvector <- c(-0.52, -0.05, 2, 10, 5, 0.25, -0.3, -0.1,
                 -1, -90, 0.5, 0.048, -0.14) 

inputvector_1_2 <- (inputvector)/2

in_upp <- inputvector + inputvector_1_2

in_low <- inputvector - inputvector_1_2



simpact_prior <-list( c("unif", in_low[1], in_upp[1]), # dissolution.alpha_0 = - 0.52
                      c("unif", in_low[2], in_upp[2]), # dissolution.alpha_4 = -0.05
                      
                      c("unif", in_low[3], in_upp[3]), # formation.hazard.agegapry.baseline = 2
                      c("unif", in_low[4], in_upp[4]), # person.agegap.man.dist.normal.mu and ~.woman.~ = 10
                      c("unif", in_low[5], in_upp[5]), # person.agegap.man.dist.normal.sigma and ~.woman.~ = 5
                      c("unif", in_low[6], in_upp[6]), # formation.hazard.agegapry.gap_agescale_man and ~_woman = 0.25
                      c("unif", in_low[7], in_upp[7]), # formation.hazard.agegapry.numrel_man and ~._woman = -0.3
                      c("unif", in_low[8], in_upp[8]), # formation.hazard.agegapry.numrel_diff = -0.1
                      
                      c("unif", in_low[9], in_upp[9]), # hivtransmission.param.a = -1
                      c("unif", in_low[10], in_upp[10]), # hivtransmission.param.b = -90
                      c("unif", in_low[11], in_upp[11]), # hivtransmission.param.c = 0.5
                      c("unif", in_low[12], in_upp[12]), # hivtransmission.param.f1 = 0.048
                      c("unif", in_low[13], in_upp[13]) # hivtransmission.param.f2 = -0.14
                      
                      
)



# Default model --------------------

wrapper.simpact.default <- function(inputvector = inputvector){
  
  
  
  work.dir <-  "/home/dniyukuri/lustre/calibration_epi_behav"  # on CHPC
  
  # work.dir <- "/home/david/benchmark_master_model" # on laptop
  
  source("/home/dniyukuri/lustre/calibration_epi_behav/needed.functions.RSimpactHelp.R")
  source("/home/dniyukuri/lustre/calibration_epi_behav/complete.master.epidemic.metrics.R")
  
  # source("/home/david/benchmark_master_model/calibration_epi_behav/needed.functions.RSimpactHelp.R")
  # source("/home/david/benchmark_master_model/calibration_epi_behav/complete.master.epidemic.metrics.R")
  
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
  
  
  # default mode itself
  
  
  simpact.default <- function(inputvector = inputvector){
    
    
    
    work.dir <-  "/home/dniyukuri/lustre/calibration_epi_behav"  # on CHPC
    
    # work.dir <- "/home/david/benchmark_master_model/" # on laptop
    
    source("/home/dniyukuri/lustre/calibration_epi_behav/needed.functions.RSimpactHelp.R")
    source("/home/dniyukuri/lustre/calibration_epi_behav/complete.master.epidemic.metrics.R")
    
    # source("/home/david/benchmark_master_model/needed.functions.RSimpactHelp.R")
    # source("/home/david/benchmark_master_model/complete.master.epidemic.metrics.R")
    
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
    
    
    
    ## Run Simpact for specific parameter combination
    
    age.distr <- agedistr.creator(shape = 5, scale = 65)
    #
    cfg.list <- input.params.creator(population.eyecap.fraction = 0.2,
                                     population.simtime = 40, 
                                     population.nummen = 5000, 
                                     population.numwomen = 5000,
                                     hivseed.time = 10, 
                                     hivseed.type = "amount",
                                     hivseed.amount = 10, 
                                     hivseed.age.min = 20,
                                     hivseed.age.max = 50,
                                     formation.hazard.agegapry.meanage = -0.025,
                                     debut.debutage = 15
    )
    
    # # Assumption of nature of sexual network
    # #########################################
    #
    cfg.list["population.msm"] = "no"
    
    
    # # Sexual behaviour
    # ###################
    #
    
    seedid <- inputvector[1]
    
    cfg.list["dissolution.alpha_0"] <- inputvector[2] # [1] # -0.52 c("unif", -1, 0)
    cfg.list["dissolution.alpha_4"] <- inputvector [3] # [2] # -0.05 c("unif", -0.5, 0)
    cfg.list["formation.hazard.agegapry.baseline"] <- inputvector[4] # [3] # 2 c("unif", 1, 3)
    cfg.list["person.agegap.man.dist.normal.mu"] <- inputvector[5] # [4] # 0 c("unif", -0.5, 0.5)
    cfg.list["person.agegap.woman.dist.normal.mu"] <- inputvector[5] # [4] # 0
    cfg.list["person.agegap.man.dist.normal.sigma"] <- inputvector[6] # [5] # 3 c("unif", 2, 4)
    cfg.list["person.agegap.woman.dist.normal.sigma"] <- inputvector[6] # [5] # 3 
    cfg.list["formation.hazard.agegapry.gap_agescale_man"] <- inputvector[7] # [6] # 0.25 c("unif", 0, 1)
    cfg.list["formation.hazard.agegapry.gap_agescale_woman"] <- inputvector[7] # [6] # 0.25
    cfg.list["formation.hazard.agegapry.numrel_man"] <- inputvector[8] # [7] # -0.3 c("unif", -1, 0)
    cfg.list["formation.hazard.agegapry.numrel_woman"] <- inputvector[8] # [7] # -0.3
    cfg.list["formation.hazard.agegapry.numrel_diff"] <- inputvector[9] # [8] # -0.1 c("unif", -0.9, 0)
    
    
    # # HIV transmission
    # ###################
    #
    
    cfg.list["hivtransmission.param.a"] <- inputvector[10] # [10] # -1 c("unif", -2, 0)
    cfg.list["hivtransmission.param.b"] <- inputvector[11] # [11] # -90 c("unif", -100, -80)
    cfg.list["hivtransmission.param.c"] <- inputvector[12] # [12] # 0.5 c("unif", 0, 1)
    cfg.list["hivtransmission.param.f1"] <- inputvector[13] # [13] # 0.04879016 c("unif", 0, 0.5)
    cfg.list["hivtransmission.param.f2"] <- inputvector[14] # [14] # -0.1386294 c("unif", -0.5, 0)
    
    # Disease progression > may be remove in parameter to estimates
    
    cfg.list["person.vsp.toacute.x"] <- 5 # inputvector[15] # [15] # 5 c("unif", 3, 7)
    cfg.list["person.vsp.toaids.x"] <- 7 # inputvector[16] # [16] # 7 c("unif", 5, 9)
    cfg.list["person.vsp.tofinalaids.x"] <- 12 # inputvector[17] # [17] # 12 c("unif", 10, 14)
    
    
    #
    # # Demographic
    # ##############
    #
    
    cfg.list["conception.alpha_base"] <- -2.7 # inputvector[18] # [18] # -2.7 c("unif", -3.5, -1.7)
    
    
    
    # # Assumptions to avoid negative branch lengths
    # ###############################################
    # # + sampling == start ART
    # # when someone start ART, he/she is sampled and becomes non-infectious
    
    cfg.list["monitoring.fraction.log_viralload"] <- 0 # very important for transmission tree and sequence simulation
    
    # Note: If treatment is started, the person’s set-point viral load value will be lowered according 
    # to the setting in monitoring.fraction.log_viralload, if we set it to 0, it means the person is no more infectious
    
    #
    # ## Add-ons
    #
    ### BEGIN Add-on
    cfg.list["mortality.aids.survtime.C"] <- 65
    cfg.list["mortality.aids.survtime.k"] <- -0.2
    cfg.list["dropout.interval.dist.type"] <- "uniform"
    cfg.list["dropout.interval.dist.uniform.min"] <- 1000
    cfg.list["dropout.interval.dist.uniform.max"] <- 2000
    
    cfg.list["person.survtime.logoffset.dist.type"] <- "normal"
    cfg.list["person.survtime.logoffset.dist.normal.mu"] <- 0
    cfg.list["person.survtime.logoffset.dist.normal.sigma"] <- 0.1
    
    cfg.list["person.agegap.man.dist.type"] <- "normal" #fixed
    #cfg.list["person.agegap.man.dist.fixed.value"] <- -6
    cfg.list["person.agegap.woman.dist.type"] <- "normal" #"fixed"
    #cfg.list["person.agegap.woman.dist.fixed.value"] <- -6
    
    
    cfg.list["person.eagerness.man.dist.gamma.a"] <- 0.23 # 0.23
    cfg.list["person.eagerness.woman.dist.gamma.a"] <- 0.23 # 0.23
    cfg.list["person.eagerness.man.dist.gamma.b"] <- 45 # 45
    cfg.list["person.eagerness.woman.dist.gamma.b"] <- 45 # 45
    
    
    
    cfg.list["monitoring.cd4.threshold"] <- 0 # 0 means nobody qualifies for ART
    cfg.list["person.art.accept.threshold.dist.fixed.value"] <- 0.75
    cfg.list["diagnosis.baseline"] <- -99999 # -2
    
    
    #### END Add-ons
    
    
    # # ART intervention
    # ###################
    
    
    # Let's introduce ART, and evaluate whether the HIV prevalence drops less  rapidly
    
    # Introducing ART
    art.intro <- list()
    art.intro["time"] <- 23    # ~2000
    art.intro["diagnosis.baseline"] <- -2
    art.intro["monitoring.cd4.threshold"] <- 100
    
    art.intro1 <- list()
    art.intro1["time"] <- 25     # ~2002
    art.intro1["diagnosis.baseline"] <- -1.8
    art.intro1["monitoring.cd4.threshold"] <- 150
    
    art.intro2 <- list()
    art.intro2["time"] <- 28     # ~2005
    art.intro2["diagnosis.baseline"] <- -1.5
    art.intro2["monitoring.cd4.threshold"] <- 200
    
    art.intro3 <- list()
    art.intro3["time"] <- 33     # ~2010
    art.intro3["diagnosis.baseline"] <- -1
    art.intro3["monitoring.cd4.threshold"] <- 350
    
    art.intro4 <- list()
    art.intro4["time"] <- 36     # ~2013
    art.intro4["monitoring.cd4.threshold"] <- 500
    
    art.intro5 <- list()
    art.intro5["time"] <- 39     # ~2016
    art.intro5["monitoring.cd4.threshold"] <- 700 # 
    
    
    interventionlist <- list(art.intro, art.intro1, art.intro2, art.intro3, art.intro4, art.intro5)
    
    
    intervention <- interventionlist
    
    # Events
    cfg.list["population.maxevents"] <- as.numeric(cfg.list["population.simtime"][1]) * as.numeric(cfg.list["population.nummen"][1]) * 3
    
    # Avoid overlaping in same directory
    
    #creating subfolder with unique name for each simulation
    generate.filename <- function(how.long){
      
      rn <- sample(1:100,1)
      t <- as.numeric(Sys.time())
      set.seed((t - floor(t)) * 1e8)
      chars <- c(letters, LETTERS)
      sub.dir.sim.id <-  paste0(sample(chars,how.long), collapse = "")
      
      noise.sample1 <- sample(8:15,1, replace = TRUE)
      sub.dir.sim.id.ext <- paste0(sample(chars,noise.sample1), collapse = "")
      noise.sample <- sample(1:1000,1)
      noise.sample2 <- sample(8:17,1, replace = TRUE)
      sub.dir.sim.id <- paste0(sub.dir.sim.id.ext,
                               paste0(sample(chars,noise.sample2), collapse = ""),noise.sample, rn)
      
      return(sub.dir.sim.id)
    }
    
    
    
    sub.dir.rename <- paste0(work.dir,"/temp/",generate.filename(10))
    
    
    
    
    # Running Simpact 
    
    
    results <- tryCatch(simpact.run(configParams = cfg.list,
                                    destDir = sub.dir.rename,
                                    agedist = age.distr,
                                    seed = seedid,
                                    intervention = intervention),
                        error = simpact.errFunction)
    
    
    
    
    if (length(results) == 0){
      
      epidemic.metrics <- rep(NA, 69)
    } else {
      if (as.numeric(results["eventsexecuted"]) >= (as.numeric(cfg.list["population.maxevents"]) - 1)) {
        
        epidemic.metrics <- rep(NA, 69)
        
      } else {
        
        datalist.agemix <- readthedata(results)
        
        epidemic.metrics <- complete.master.epidemic.metrics(datalist = datalist.agemix)
        
      }
    }
    
    
    return(epidemic.metrics)
    
    unlink(paste0(sub.dir.rename), recursive = TRUE)
    
    
  }
  
  
  
  results.epidemic.metrics <- tryCatch(simpact.default(inputvector = inputvector),
                                       error=function(e) return(rep(NA, 69)))
  
  return(results.epidemic.metrics)
  
  
}






# Calibration with classic features ---------------------------------------



# 
wrapper.simpact4ABC.epi.behav <- function(inputvector = inputvector){
  
  
  work.dir <-  "/home/dniyukuri/lustre/calibration_epi_behav"  # on CHPC
  
  # work.dir <- "/home/david/benchmark_master_model/calibration_epi_behav" # on laptop
  
  # source("/home/david/benchmark_master_model/calibration_epi_behav/needed.functions.RSimpactHelp.R")
  # source("/home/david/benchmark_master_model/calibration_epi_behav/compute.summary.statistics.hiv.epi.behavior.R")
  
  source("/home/dniyukuri/lustre/calibration_epi_behav/needed.functions.RSimpactHelp.R")
  source("/home/dniyukuri/lustre/calibration_epi_behav/compute.summary.statistics.hiv.epi.behavior.R")
  
  
  
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
  
  
  
  simpact4ABC.epi.behav <- function(inputvector = inputvector){
    

    work.dir <-  "/home/dniyukuri/lustre/calibration_epi_behav"  # on CHPC
    
    # work.dir <- "/home/david/benchmark_master_model/calibration_epi_behav" # on laptop
    
    # source("/home/david/benchmark_master_model/calibration_epi_behav/needed.functions.RSimpactHelp.R")
    # source("/home/david/benchmark_master_model/calibration_epi_behav/compute.summary.statistics.hiv.epi.behavior.R")
    
    source("/home/dniyukuri/lustre/calibration_epi_behav/needed.functions.RSimpactHelp.R")
    source("/home/dniyukuri/lustre/calibration_epi_behav/compute.summary.statistics.hiv.epi.behavior.R")
    
    
    
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
    
    
    
    ## Run Simpact for specific parameter combination
    ## Run Simpact for specific parameter combination
    
    age.distr <- agedistr.creator(shape = 5, scale = 65)
    #
    cfg.list <- input.params.creator(population.eyecap.fraction = 0.2,
                                     population.simtime = 40, 
                                     population.nummen = 5000, 
                                     population.numwomen = 5000,
                                     hivseed.time = 10, 
                                     hivseed.type = "amount",
                                     hivseed.amount = 10, 
                                     hivseed.age.min = 20,
                                     hivseed.age.max = 50,
                                     formation.hazard.agegapry.meanage = -0.025,
                                     debut.debutage = 15
    )
    
    # # Assumption of nature of sexual network
    # #########################################
    #
    cfg.list["population.msm"] = "no"
    
    # # Sexual behaviour
    # ###################
    #
    
    seedid <- inputvector[1]
    
    cfg.list["dissolution.alpha_0"] <- inputvector[2] # [1] # -0.52 c("unif", -1, 0)
    cfg.list["dissolution.alpha_4"] <- inputvector [3] # [2] # -0.05 c("unif", -0.5, 0)
    cfg.list["formation.hazard.agegapry.baseline"] <- inputvector[4] # [3] # 2 c("unif", 1, 3)
    cfg.list["person.agegap.man.dist.normal.mu"] <- inputvector[5] # [4] # 0 c("unif", -0.5, 0.5)
    cfg.list["person.agegap.woman.dist.normal.mu"] <- inputvector[5] # [4] # 0
    cfg.list["person.agegap.man.dist.normal.sigma"] <- inputvector[6] # [5] # 3 c("unif", 2, 4)
    cfg.list["person.agegap.woman.dist.normal.sigma"] <- inputvector[6] # [5] # 3 
    cfg.list["formation.hazard.agegapry.gap_agescale_man"] <- inputvector[7] # [6] # 0.25 c("unif", 0, 1)
    cfg.list["formation.hazard.agegapry.gap_agescale_woman"] <- inputvector[7] # [6] # 0.25
    cfg.list["formation.hazard.agegapry.numrel_man"] <- inputvector[8] # [7] # -0.3 c("unif", -1, 0)
    cfg.list["formation.hazard.agegapry.numrel_woman"] <- inputvector[8] # [7] # -0.3
    cfg.list["formation.hazard.agegapry.numrel_diff"] <- inputvector[9] # [8] # -0.1 c("unif", -0.9, 0)
    
    
    # # HIV transmission
    # ###################
    #
    
    cfg.list["hivtransmission.param.a"] <- inputvector[10] # [10] # -1 c("unif", -2, 0)
    cfg.list["hivtransmission.param.b"] <- inputvector[11] # [11] # -90 c("unif", -100, -80)
    cfg.list["hivtransmission.param.c"] <- inputvector[12] # [12] # 0.5 c("unif", 0, 1)
    cfg.list["hivtransmission.param.f1"] <- inputvector[13] # [13] # 0.04879016 c("unif", 0, 0.5)
    cfg.list["hivtransmission.param.f2"] <- inputvector[14] # [14] # -0.1386294 c("unif", -0.5, 0)
    
    # Disease progression > may be remove in parameter to estimates
    
    cfg.list["person.vsp.toacute.x"] <- 5 # inputvector[15] # [15] # 5 c("unif", 3, 7)
    cfg.list["person.vsp.toaids.x"] <- 7 # inputvector[16] # [16] # 7 c("unif", 5, 9)
    cfg.list["person.vsp.tofinalaids.x"] <- 12 # inputvector[17] # [17] # 12 c("unif", 10, 14)
    
    
    #
    # # Demographic
    # ##############
    #
    
    cfg.list["conception.alpha_base"] <- -2.7 # inputvector[18] # [18] # -2.7 c("unif", -3.5, -1.7)
    
    
    
    # # Assumptions to avoid negative branch lengths
    # ###############################################
    # # + sampling == start ART
    # # when someone start ART, he/she is sampled and becomes non-infectious
    
    cfg.list["monitoring.fraction.log_viralload"] <- 0 # very important for transmission tree and sequence simulation
    
    # Note: If treatment is started, the person’s set-point viral load value will be lowered according 
    # to the setting in monitoring.fraction.log_viralload, if we set it to 0, it means the person is no more infectious
    
    #
    # ## Add-ons
    #
    ### BEGIN Add-on
    cfg.list["mortality.aids.survtime.C"] <- 65
    cfg.list["mortality.aids.survtime.k"] <- -0.2
    cfg.list["dropout.interval.dist.type"] <- "uniform"
    cfg.list["dropout.interval.dist.uniform.min"] <- 1000
    cfg.list["dropout.interval.dist.uniform.max"] <- 2000
    
    cfg.list["person.survtime.logoffset.dist.type"] <- "normal"
    cfg.list["person.survtime.logoffset.dist.normal.mu"] <- 0
    cfg.list["person.survtime.logoffset.dist.normal.sigma"] <- 0.1
    
    cfg.list["person.agegap.man.dist.type"] <- "normal" #fixed
    #cfg.list["person.agegap.man.dist.fixed.value"] <- -6
    cfg.list["person.agegap.woman.dist.type"] <- "normal" #"fixed"
    #cfg.list["person.agegap.woman.dist.fixed.value"] <- -6
    
    
    cfg.list["person.eagerness.man.dist.gamma.a"] <- 0.23 # 0.23
    cfg.list["person.eagerness.woman.dist.gamma.a"] <- 0.23 # 0.23
    cfg.list["person.eagerness.man.dist.gamma.b"] <- 45 # 45
    cfg.list["person.eagerness.woman.dist.gamma.b"] <- 45 # 45
    
    
    
    cfg.list["monitoring.cd4.threshold"] <- 0 # 0 means nobody qualifies for ART
    cfg.list["person.art.accept.threshold.dist.fixed.value"] <- 0.75
    cfg.list["diagnosis.baseline"] <- -99999 # -2
    
    
    #### END Add-ons
    
    
    # # ART intervention
    # ###################
    
    
    # Let's introduce ART, and evaluate whether the HIV prevalence drops less  rapidly
    
    # Introducing ART
    art.intro <- list()
    art.intro["time"] <- 23    # ~2000
    art.intro["diagnosis.baseline"] <- -2
    art.intro["monitoring.cd4.threshold"] <- 100
    
    art.intro1 <- list()
    art.intro1["time"] <- 25     # ~2002
    art.intro1["diagnosis.baseline"] <- -1.8
    art.intro1["monitoring.cd4.threshold"] <- 150
    
    art.intro2 <- list()
    art.intro2["time"] <- 28     # ~2005
    art.intro2["diagnosis.baseline"] <- -1.5
    art.intro2["monitoring.cd4.threshold"] <- 200
    
    art.intro3 <- list()
    art.intro3["time"] <- 33     # ~2010
    art.intro3["diagnosis.baseline"] <- -1
    art.intro3["monitoring.cd4.threshold"] <- 350
    
    art.intro4 <- list()
    art.intro4["time"] <- 36     # ~2013
    art.intro4["monitoring.cd4.threshold"] <- 500
    
    art.intro5 <- list()
    art.intro5["time"] <- 39     # ~2016
    art.intro5["monitoring.cd4.threshold"] <- 700 # 
    
    
    interventionlist <- list(art.intro, art.intro1, art.intro2, art.intro3, art.intro4, art.intro5)
    
    intervention <- interventionlist
    
    # Events
    cfg.list["population.maxevents"] <- as.numeric(cfg.list["population.simtime"][1]) * as.numeric(cfg.list["population.nummen"][1]) * 3
    
    # Avoid overlaping in same directory
    
    #creating subfolder with unique name for each simulation
    generate.filename <- function(how.long){
      
      rn <- sample(1:100,1)
      t <- as.numeric(Sys.time())
      set.seed((t - floor(t)) * 1e8)
      chars <- c(letters, LETTERS)
      sub.dir.sim.id <-  paste0(sample(chars,how.long), collapse = "")
      
      noise.sample1 <- sample(8:15,1, replace = TRUE)
      sub.dir.sim.id.ext <- paste0(sample(chars,noise.sample1), collapse = "")
      noise.sample <- sample(1:1000,1)
      noise.sample2 <- sample(8:17,1, replace = TRUE)
      sub.dir.sim.id <- paste0(sub.dir.sim.id.ext,
                               paste0(sample(chars,noise.sample2), collapse = ""),noise.sample, rn)
      
      return(sub.dir.sim.id)
    }
    
    
    
    sub.dir.rename <- paste0(work.dir,"/temp/",generate.filename(10))
    
    
    
    
    # Running Simpact 
    
    
    results <- tryCatch(simpact.run(configParams = cfg.list,
                                    destDir = sub.dir.rename,
                                    agedist = age.distr,
                                    seed = seedid,
                                    intervention = intervention),
                        error = simpact.errFunction)
    
    
    
    
    if (length(results) == 0){
      
      summ.stats <- rep(NA, 39)
      
    } else {
      
      if (as.numeric(results["eventsexecuted"]) >= (as.numeric(cfg.list["population.maxevents"]) - 1)) {
        
        summ.stats <- rep(NA, 39)
        
      } else {
        
        datalist.agemix <- readthedata(results)
        
        summ.stats <- compute.summary.statistics.hiv.epi.behavior(datalist = datalist.agemix,
                                                                  timewindow = c(35, 40))
        
      }
    }
    
    
    return(summ.stats)
    
    
    unlink(paste0(sub.dir.rename), recursive = TRUE)
    
    
    
  }
  
  
  
  results.summ.stats <- tryCatch(simpact4ABC.epi.behav(inputvector = inputvector),
                                 error=function(e) return(rep(NA, 39)))
  
  return(results.summ.stats)
  
  
}




# Calibration -------------------


# source("/home/dniyukuri/lustre/calibration_epi_behav/calibration.ABC.R")

# source("/home/david/benchmark_master_model/calibration_epi_behav/calibration.ABC.R")



# being calibration -----------

model.sim <- wrapper.simpact4ABC.epi.behav
sum_stat_obs <- mean.classic.features
simpact_prior <- simpact_prior
design.points <- 3752
alpha <- 0.30
seed.val <- 777
n_cores <- 56
method <- "rejection"




# Calibration with ABC approach

# (i). Having plausible ranges for the parameters, sample the parameter spaces by latin hypercube several times
# (ii). Run the default model and compute the summary statistics
# (iii). Use different ABC-based methods to fit the model, this will give parameters estimates and associated summary statistics



simpact_prior <- simpact_prior

min.v <- vector()
max.v <- vector()

for( i in 1:length(simpact_prior)){
  
  min.i <- as.numeric(simpact_prior[[i]][[2]])
  max.i <- as.numeric(simpact_prior[[i]][[3]]) 
  
  if(min.i < max.i){
    
    min.v <- c(min.v, as.numeric(simpact_prior[[i]][[2]]))
    max.v <- c(max.v, as.numeric(simpact_prior[[i]][[3]])) 
    
  }else{
    
    min.v <- c(min.v, as.numeric(simpact_prior[[i]][[3]]))
    max.v <- c(max.v, as.numeric(simpact_prior[[i]][[2]]))  
    
  }
  
}


variables <- length(simpact_prior)

set.seed(seed.val)

rlhs <- lhs::randomLHS(design.points, variables)


lhs.df <- rlhs


for (j in 1:ncol(lhs.df)){
  
  min.var <- min.v[j]
  max.var <- max.v[j]
  
  for(k in 1:nrow(lhs.df)){
    
    lhs.df[k,j] <- qunif(lhs.df[k,j], min = min.var, max = max.var)
    
  }
  
}



par.sim <- lhs.df # results of (i): parameter matrix

# par.sim <- par.sim[, -c(ncol(lhs.df))]

# (ii)

stat.sim.cdf <- simpact.parallel(model = model.sim, # simpact4ABC.epi.behav,
                                 actual.input.matrix = par.sim,
                                 seed_count = seed.val,
                                 n_cluster = 56)

saveRDS(stat.sim.cdf, file = "raw_simulation_calibration_hiv_epi_behavior.RDS")



# -------- ABC calibration



stat.sim <- stat.sim.cdf[,1:length(sum_stat_obs)] # results of (ii): summary statistics matrix obtained from simulations done with parameter matrix



# (iii)

# Condition: nrow(par.sim) == nrow(stat.sim)


# (iv) removing NA

# NA.stat.sim <- na.omit(stat.sim)

index.na.fun <- function(stat.sim = stat.sim){
  index.vec <- vector()
  for(i in 1:nrow(stat.sim)){
    row.i <- stat.sim[i,] 
    if(NA%in%row.i){
      index.vec <- c(index.vec, i)
    }
  }
  return(index.vec)
}

na.index <- index.na.fun(stat.sim = stat.sim)

if(length(na.index)<1){
  
  na.stat.sim <- stat.sim
  na.par.sim <- par.sim
  
}else{
  
  na.stat.sim <- stat.sim[-c(na.index),]
  na.par.sim <- par.sim[-c(na.index), ]
  
}




abc_res <- abc::abc(target=sum_stat_obs, 
                    param=na.par.sim, 
                    sumstat=na.stat.sim, 
                    tol=alpha, 
                    method = paste0(method)) 


saveRDS(abc_res, file = "abc_res_calibration_hiv_epi_behavior.RDS")


# --------

# Run calibrated models after calibration ---------------------

cdf <- abc_res

cal.val <- cdf$unadj.values


stats.parameters <- function(input = input){
  
  input <- na.omit(input)
  
  colm.stats.cal.va <- NA
  
  for(i in 1:ncol(input)){
    
    input.i <- input[,i]
    
    quantirles.v <- quantile(input.i, probs = seq(0, 1, 0.25))
    
    quantirles.25.50.75 <- as.numeric(quantirles.v)[2:4]
    
    
    min.v <- min(input.i)
    max.v <- max(input.i)
    mean.v <- mean(input.i)
    median.v <- median(input.i)
    sd.v <- sd(input.i)
    
    param.stats <- c(quantirles.25.50.75, min.v, max.v, mean.v, median.v, sd.v)
    
    
    colm.stats.cal.va <- rbind(colm.stats.cal.va, param.stats)
    
  }
  
  
  colm.stats.cal.va <- colm.stats.cal.va[-1,]
  
  colnames(colm.stats.cal.va) <- c("Q1", "Q2", "Q3", "min", "max", "mean", "median", "sd")
  
  return(colm.stats.cal.va)
  
}


stats.params.table <- stats.parameters(input = cal.val)



saveRDS(stats.params.table, file = "params.table_raw_simulation_calibration_hiv_epi_behavior.RDS")





# 1st with all selected parameters


inputmatrix.params <- cal.val


epi.metrics_params <- simpact.parallel(model = wrapper.simpact.default,
                                       actual.input.matrix = inputmatrix.params,
                                       seed_count = 1,
                                       n_cluster = 56)


saveRDS(epi.metrics_params, file = "epi.metrics_after_calibration_run_params_hiv_epi_behavior.RDS")




# 2nd with median of selected params


median.select_parms <- colMedians(as.matrix(cal.val), na.rm = TRUE)

reps <- 1120

inputmatrix <- matrix(rep(median.select_parms, reps), byrow = TRUE, nrow = reps)

epi.metrics_median_params <- simpact.parallel(model = wrapper.simpact.default,
                                              actual.input.matrix = inputmatrix,
                                              seed_count = 1,
                                              n_cluster = 56)


saveRDS(epi.metrics_median_params, file = "epi.metrics_after_calibration_run_median_params_hiv_epi_behavior.RDS")



