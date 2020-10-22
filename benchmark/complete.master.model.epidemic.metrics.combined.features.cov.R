
# Master model which simulates the epidemic and compute benchmark values for 
# HIV transmission network determinants and summary statistics (targets during calibration) in different scenarios


# The output is a vector of values of

# 1. Transmission network determinants (characteristics)

# - age mixing statistics
# - mean, median, and standard deviation of onward transmissions
# - temporal trend of incidence


# 2. Epidmiological, demographical, sexual behavioural, and interventions realted summary statistics

# 3. Phylogenetic summary statistics in MCAR (35:95, by 5) scenarios, each scenario returns measurements which are described
# in compute.summary.statistics.phylo.MCAR scripts):

# Missing Completly at Random has 13 scenarios




complete.master.model.epidemic.metrics.combined.features.cov <- function(inputvector = inputvector){
  
  
  # 
  # source("/home/dniyukuri/lustre/benchmark_master_model/advanced.transmission.network.builder.R")
  
  
  source("/home/dniyukuri/lustre/benchmark_master_model/needed.functions.RSimpactHelp.R")
  
  source("/home/dniyukuri/lustre/benchmark_master_model/complete.master.epidemic.metrics.R")
  
  source("/home/dniyukuri/lustre/benchmark_master_model/compute.summary.statistics.hiv.epi.behavior.R")
  
  source("/home/dniyukuri/lustre/benchmark_master_model/compute.summary.statistics.phylo.MCAR.R")
  
  # source("/home/dniyukuri/lustre/benchmark_master_model/compute.summary.statistics.phylo.MAR.R")
  # 
  # 
  
  work.dir <- "/home/dniyukuri/lustre/benchmark_master_model" # on CHPC
  
  
  # work.dir <- "/home/david/benchmark_master_model" # on laptop
  
  
  # source("/home/david/benchmark_master_model/advanced.transmission.network.builder.R")
  # 
  # source("/home/david/benchmark_master_model/needed.functions.RSimpactHelp.R")
  # 
  # source("/home/david/benchmark_master_model/complete.master.epidemic.metrics.R")
  # 
  # source("/home/david/benchmark_master_model/compute.summary.statistics.hiv.epi.behavior.R")
  # 
  # source("/home/david/benchmark_master_model/compute.summary.statistics.phylo.MCAR.R")
  
  # source("/home/david/benchmark_master_model/compute.summary.statistics.phylo.MAR.R")
  
  
  # work.dir <-  "/home/david/benchmark_master_model"
  
  
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
  
  
  
  ###########################################
  # Step 1: Setup and running simpact      #
  ###########################################
  
  
  
  
  
  
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
  #################
  
  results <- tryCatch(simpact.run(configParams = cfg.list,
                                  destDir = sub.dir.rename,
                                  agedist = age.distr,
                                  seed = seedid,
                                  intervention = intervention),
                      error = simpact.errFunction)
  
  
  
  datalist.agemix <- readthedata(results)
  
  
  # datalist.agemix <- readRDS("datalist.agemix.RDS")
  
  
  ###########################################
  # Step 2: Construct transmission networks #
  ###########################################
  
  
  simpact.trans.net <- transmission.network.builder(datalist = datalist.agemix, endpoint = 40)
  
  # simpact.trans.net.projection <- transmission.network.builder(datalist = datalist.agemix, endpoint = 45)
  
  
  
  net.size.vector <- vector() # i_th seed in the list of seeds
  
  for(i in 1:length(simpact.trans.net)){
    
    tree.n <- simpact.trans.net[[i]] # transmission network for i^th seed
    
    net.size.vector <- c(net.size.vector, nrow(as.data.frame(tree.n)))
    
  }
  
  
  big.index <- which(net.size.vector>=50)
  
  
  
  
  
  ###############################
  # Step 3: Sequence simulation #
  ###############################
  
  
  dirseqgen <- work.dir
  
  seeds.num <- inputvector[1]
  
  # Sequence simulation is done for at least a transmission network with 6 individuals
  # This means that limitTransmEvents equal at least 7
  
  sequence.simulation.seqgen.par(dir.seq = work.dir,
                                 sub.dir.rename = sub.dir.rename,
                                 simpact.trans.net = simpact.trans.net, 
                                 seq.gen.tool = "seq-gen",
                                 seeds.num = seeds.num,
                                 endpoint = 40,
                                 limitTransmEvents = 7, # no less than 7
                                 hiv.seq.file = "hiv.seq.C.pol.j.fasta",
                                 clust = TRUE) # hiv.seq.file lodged in work.dir
  
  
  
  # Transform the sequence format to be handled by ClusterPicker
  sequ.dna <- read.dna(file = paste0(sub.dir.rename,"/C.Epidemic_seed.seq.bis.sim.nwk.fasta"), format = "interleaved")
  write.dna(sequ.dna, file = paste0(sub.dir.rename,"/C.Epidemic.fas") , format = "fasta")
  
  
  #####################################################
  ### I. Compute transmission network characteristics #: 69
  #####################################################
  
  # source("/home/niyukuri/phylosimpact_simulation_studies_2018/stress_testing/stress_testing_final/complete.master.epic.metrics.R")
  
  
  epidemic.metrics <- tryCatch(complete.master.epidemic.metrics(datalist = datalist.agemix),
                               error=function(e) return(rep(NA, 69)))
  
  
  ##################################
  ### II. Compute classic features #: 39
  ################################## ??? change arguments
  
  # source("/home/niyukuri/phylosimpact_simulation_studies_2018/stress_testing/stress_testing_final/compute.summary.statistics.hiv.epi.behavior.R")
  
  epi.behav.stats <- tryCatch(compute.summary.statistics.hiv.epi.behavior(datalist = datalist.agemix,
                                                                          timewindow = c(35, 40)),
                              error=function(e) return(rep(NA, 39)))
  
  
  ########################################
  ## III. Compute phylogenetic features  #: 18
  ########################################
  
  
  MCAR.cov.100 <- tryCatch(compute.summary.statistics.phylo.MCAR(simpact.trans.net = simpact.trans.net,
                                                                 datalist.agemix = datalist.agemix,
                                                                 work.dir = work.dir,
                                                                 sub.dir.rename = sub.dir.rename,
                                                                 dirfasttree = work.dir,
                                                                 limitTransmEvents = 7,
                                                                 seq.cov = 100,
                                                                 age.group.15.25 = c(15,25),
                                                                 age.group.25.40 = c(25,40),
                                                                 age.group.40.50 = c(40,50),
                                                                 endpoint = 40,
                                                                 timewindow = c(35,40),
                                                                 cut.off = 7),
                           error=function(e) return(rep(NA, 18)))
  
  
  
  # MCAR.All <- c(MCAR.cov.35, MCAR.cov.40, MCAR.cov.45, MCAR.cov.50, MCAR.cov.55, MCAR.cov.60, MCAR.cov.65, MCAR.cov.70, 
  #               MCAR.cov.75, MCAR.cov.80, MCAR.cov.85, MCAR.cov.90, MCAR.cov.95, MCAR.cov.100)
  
  names.columns <- names(MCAR.cov.100)
  
  
  results.mcar <- as.numeric(MCAR.cov.100)
  
  
  names(results.mcar) <- c(paste0("cov.MCAR.",100,".",paste0(names.columns)))
  
  
  # Values
  
  # outputvector.values <- c(epidemic.metrics, epi.behav.stats, 
  #                          results.mcar, results.mar.a, 
  #                          results.mar.b, results.mar.c)
  
  
  outputvector.values <- c(epidemic.metrics, epi.behav.stats, 
                           results.mcar)
  
  
  
  return(outputvector.values)
  
  unlink(paste0(sub.dir.rename), recursive = TRUE)
  
  
  
  
}


