#' Computing epidmiological, demographical, sexual behavioural, and interventions realted summary statistics
#'

#' @param datalist.agemix Data list of simpact output produced by \code{\link{readthedata()}}
#' @param timewindow Time window in which the experience is carried out

#' @export


# The outputs of the function are:

# (i) Demographic
# - Population growth

# (ii) Epidemic

# - Prevalence
# - Incidence

# (iii) Sexual behaviour

# - point prevalence of conurrent relationships
# - relationship per person per year
# - mean, median and standard deviation of age gap


# (iv) Interventions

# - ART coverage
# - Viral load suppression


# "Pop.growthrate", 
# 
# Prevalence
# 
# Incidence
# 
# "pp.cp.6months.male.rels",
# 
# "relsperpersonperyear", 
# "agegap.mean", "agegap.med", "agegap.sd", 
# 
# "ART.cov.vector", "VL.suppression.fraction"



compute.summary.statistics.hiv.epi.behavior <- function(datalist = datalist.agemix,
                                                        timewindow = c(30, 40)){
  
  source("/home/dniyukuri/lustre/calibration_combined/needed.functions.RSimpactHelp.R")
  
  # source("/home/david/benchmark_master_model/needed.functions.RSimpactHelp.R")
  
  
  datalist.agemix <- datalist
  
  
  ########################################
  # I. Behavioural and epidemic features #
  ########################################
  
  
  # 1.2. Features from sexual and transmission network
  
  # 
  # 1.2.1. Demographic feature:
  
  #   (i) Population growth rate (pop.growth.calculator function)
  
  growthrate <- pop.growth.calculator(datalist = datalist.agemix,
                                      timewindow = timewindow) # c(0, datalist.agemix$itable$population.simtime[1])
  
  
  # 1.2.2. Prevalence
  
  
  # Women
  
  hiv.prev.lt25.women <- prevalence.calculator(datalist = datalist.agemix,
                                               agegroup = c(15, 25),
                                               timepoint = 40) %>%
    dplyr::select(pointprevalence) %>%
    dplyr::slice(2) %>%
    as.numeric()
  
  hiv.prev.25.30.women <- prevalence.calculator(datalist = datalist.agemix,
                                                agegroup = c(25, 30),
                                                timepoint = 40) %>%
    dplyr::select(pointprevalence) %>%
    dplyr::slice(2) %>%
    as.numeric()
  
  hiv.prev.30.35.women <- prevalence.calculator(datalist = datalist.agemix,
                                                agegroup = c(30, 35),
                                                timepoint = 40) %>%
    dplyr::select(pointprevalence) %>%
    dplyr::slice(2) %>%
    as.numeric()
  
  hiv.prev.35.40.women <- prevalence.calculator(datalist = datalist.agemix,
                                                agegroup = c(35, 40),
                                                timepoint = 40) %>%
    dplyr::select(pointprevalence) %>%
    dplyr::slice(2) %>%
    as.numeric()
  
  hiv.prev.40.45.women <- prevalence.calculator(datalist = datalist.agemix,
                                                agegroup = c(40, 45),
                                                timepoint = 40) %>%
    dplyr::select(pointprevalence) %>%
    dplyr::slice(2) %>%
    as.numeric()
  
  hiv.prev.45.50.women <- prevalence.calculator(datalist = datalist.agemix,
                                                agegroup = c(45, 50),
                                                timepoint = 40) %>%
    dplyr::select(pointprevalence) %>%
    dplyr::slice(2) %>%
    as.numeric()
  
  
  
  # Men
  
  hiv.prev.lt25.men <- prevalence.calculator(datalist = datalist.agemix,
                                             agegroup = c(15, 25),
                                             timepoint = 40) %>%
    dplyr::select(pointprevalence) %>%
    dplyr::slice(1) %>%
    as.numeric()
  
  hiv.prev.25.30.men <- prevalence.calculator(datalist = datalist.agemix,
                                              agegroup = c(25, 30),
                                              timepoint = 40) %>%
    dplyr::select(pointprevalence) %>%
    dplyr::slice(1) %>%
    as.numeric()
  
  hiv.prev.30.35.men <- prevalence.calculator(datalist = datalist.agemix,
                                              agegroup = c(30, 35),
                                              timepoint = 40) %>%
    dplyr::select(pointprevalence) %>%
    dplyr::slice(1) %>%
    as.numeric()
  
  hiv.prev.35.40.men <- prevalence.calculator(datalist = datalist.agemix,
                                              agegroup = c(35, 40),
                                              timepoint = 40) %>%
    dplyr::select(pointprevalence) %>%
    dplyr::slice(1) %>%
    as.numeric()
  
  hiv.prev.40.45.men <- prevalence.calculator(datalist = datalist.agemix,
                                              agegroup = c(40, 45),
                                              timepoint = 40) %>%
    dplyr::select(pointprevalence) %>%
    dplyr::slice(1) %>%
    as.numeric()
  
  hiv.prev.45.50.men <- prevalence.calculator(datalist = datalist.agemix,
                                              agegroup = c(45, 50),
                                              timepoint = 40) %>%
    dplyr::select(pointprevalence) %>%
    dplyr::slice(1) %>%
    as.numeric()
  
  
  
  
  
  # agegroup <- c("15-24", "25-29", "30-34", "35-39", "40-44", "45-49")
  
  men.prev <- c(hiv.prev.lt25.men, hiv.prev.25.30.men, hiv.prev.30.35.men, hiv.prev.35.40.men, hiv.prev.40.45.men, hiv.prev.45.50.men)
  names(men.prev) <- c("prev.m.15.24", "prev.m.25.29", "prev.m.30.34", "prev.m.35.39", "prev.m.40.44", "prev.m.45.49")
  women.prev <- c(hiv.prev.lt25.women, hiv.prev.25.30.women, hiv.prev.30.35.women, hiv.prev.35.40.women, hiv.prev.40.45.women, hiv.prev.45.50.women)
  names(women.prev) <- c("prev.w.15.24", "prev.w.25.29", "prev.w.30.34", "prev.w.35.39", "prev.w.40.44", "prev.w.45.49")        
  
  
  
  
  
  # 1.2.3. Incidence
  
  
  # Women
  
  hiv.incid.lt25.women <- incidence.calculator(datalist = datalist.agemix,
                                               agegroup = c(15, 25),
                                               timewindow = c(39, 40),
                                               only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(2) %>%
    as.numeric()
  
  hiv.incid.25.30.women <- incidence.calculator(datalist = datalist.agemix,
                                                agegroup = c(25, 30),
                                                timewindow = c(39, 40),
                                                only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(2) %>%
    as.numeric()
  
  hiv.incid.30.35.women <- incidence.calculator(datalist = datalist.agemix,
                                                agegroup = c(30, 35),
                                                timewindow = c(39, 40),
                                                only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(2) %>%
    as.numeric()
  
  hiv.incid.35.40.women <- incidence.calculator(datalist = datalist.agemix,
                                                agegroup = c(35, 40),
                                                timewindow = c(39, 40),
                                                only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(2) %>%
    as.numeric()
  
  hiv.incid.40.45.women <- incidence.calculator(datalist = datalist.agemix,
                                                agegroup = c(40, 45),
                                                timewindow = c(39, 40),
                                                only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(2) %>%
    as.numeric()
  
  hiv.incid.45.50.women <- incidence.calculator(datalist = datalist.agemix,
                                                agegroup = c(45, 50),
                                                timewindow = c(39, 40),
                                                only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(2) %>%
    as.numeric()
  
  
  
  # Men
  
  hiv.incid.lt25.men <- incidence.calculator(datalist = datalist.agemix,
                                             agegroup = c(15, 25),
                                             timewindow = c(39, 40),
                                             only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(1) %>%
    as.numeric()
  
  hiv.incid.25.30.men <- incidence.calculator(datalist = datalist.agemix,
                                              agegroup = c(25, 30),
                                              timewindow = c(39, 40),
                                              only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(1) %>%
    as.numeric()
  
  hiv.incid.30.35.men <- incidence.calculator(datalist = datalist.agemix,
                                              agegroup = c(30, 35),
                                              timewindow = c(39, 40),
                                              only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(1) %>%
    as.numeric()
  
  hiv.incid.35.40.men <- incidence.calculator(datalist = datalist.agemix,
                                              agegroup = c(35, 40),
                                              timewindow = c(39, 40),
                                              only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(1) %>%
    as.numeric()
  
  hiv.incid.40.45.men <- incidence.calculator(datalist = datalist.agemix,
                                              agegroup = c(40, 45),
                                              timewindow = c(39, 40),
                                              only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(1) %>%
    as.numeric()
  
  hiv.incid.45.50.men <- incidence.calculator(datalist = datalist.agemix,
                                              agegroup = c(45, 50),
                                              timewindow = c(39, 40),
                                              only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(1) %>%
    as.numeric()
  
  
  
  
  
  # agegroup <- c("15-24", "25-29", "30-34", "35-39", "40-44", "45-49")
  
  men.incid <- c(hiv.incid.lt25.men, hiv.incid.25.30.men, hiv.incid.30.35.men, hiv.incid.35.40.men, hiv.incid.40.45.men, hiv.incid.45.50.men)
  names(men.incid) <- c("incid.m.15.24", "incid.m.25.29", "incid.m.30.34", "incid.m.35.39", "incid.m.40.44", "incid.m.45.49")
  women.incid <- c(hiv.incid.lt25.women, hiv.incid.25.30.women, hiv.incid.30.35.women, hiv.incid.35.40.women, hiv.incid.40.45.women, hiv.incid.45.50.women)
  names(women.incid) <- c("incid.w.15.24", "incid.w.25.29", "incid.w.30.34", "incid.w.35.39", "incid.w.40.44", "incid.w.45.49")        
  
  
  
  # Sexual behaviour
  
  
  #  Point 	prevalence of concurrency in the adult population
  
  # Concurrency point prevalence 6 months before a survey, among men
  
  pp.cp.6months.male.rels <- concurr.pointprev.calculator(datalist = datalist.agemix,
                                                          timepoint = timewindow[2] - 0.5) # @CHPC
  
  
  
  # pp.cp.6months.male.rels <- concurr.pointprev.calculator(datalist = datalist.agemix,
  #                                                         timepoint = timewindow[2] - 0.5) %>% @Local
  #   dplyr::select(concurr.pointprev) %>%
  #   dplyr::slice(1) %>%
  #   as.numeric()
  # 
  # 
  # pp.cp.6months.female.rels <- concurr.pointprev.calculator(datalist = datalist.agemix,
  #                                                           timepoint = timewindow[2] - 0.5) %>%
  #   dplyr::select(concurr.pointprev) %>%
  #   dplyr::slice(2) %>%
  #   as.numeric()
  # 
  
  
  # (ii) Relationship per person per year ??
  
  # relationships in timewindow
  
  dat.rels.df <- datalist.agemix$rtable
  
  dat.rels.df$FormTime
  
  window.dat.rels.df <- dplyr::filter(dat.rels.df, dat.rels.df$FormTime >= timewindow[1] & dat.rels.df$FormTime <= timewindow[2])
  
  # relsperpersonperyear <- nrow(datalist.agemix$rtable) / (nrow(datalist.agemix$ptable)/2) / (timewindow[2] - timewindow[1])
  
  relsperpersonperyear <- nrow(window.dat.rels.df) / (nrow(window.dat.rels.df)/2) / (timewindow[2] - timewindow[1])
  
  
  # (iv) SD age gap between couples
  
  # agegap.mean <- mean(datalist.agemix$rtable$AgeGap)
  # 
  # agegap.med <- median(datalist.agemix$rtable$AgeGap)
  # 
  # agegap.sd <- sd(datalist.agemix$rtable$AgeGap)
  
  
  agegap.mean <- mean(window.dat.rels.df$AgeGap)
  
  agegap.med <- median(window.dat.rels.df$AgeGap)
  
  agegap.sd <- sd(window.dat.rels.df$AgeGap)
  
  
  
  ####
  # ART coverage among adults 15+ years old from UNAIDS (2010 - 2017 estimates)
  ####
  ART.cov.eval.timepoints <- seq(from = 33,
                                 to = 40)
  
  ART.cov.vector <- rep(NA, length(ART.cov.eval.timepoints))
  
  for (art.cov.index in 1:length(ART.cov.vector)){
    ART.cov.vector[art.cov.index] <- sum(ART.coverage.calculator(datalist = datalist.agemix,
                                                                 agegroup = c(15, 50),
                                                                 timepoint = ART.cov.eval.timepoints[art.cov.index])$sum.onART) /
      sum(ART.coverage.calculator(datalist = datalist.agemix,
                                  agegroup = c(15, 50),
                                  timepoint = ART.cov.eval.timepoints[art.cov.index])$sum.cases)
  }
  
  names(ART.cov.vector) <- paste0("ART.", ART.cov.eval.timepoints)
  
  ####
  # VL suppression fraction (all ages in 2017 ~ >= 15 yo) 0.74
  ####
  VL.suppression.fraction <- VL.suppression.calculator(datalist = datalist.agemix,
                                                       agegroup = c(15, 50),
                                                       timepoint = timewindow[2],
                                                       vl.cutoff = 1000,
                                                       site="All") %>%
    dplyr::select(vl.suppr.frac) %>%
    dplyr::slice(3) %>%
    as.numeric()
  
  names(VL.suppression.fraction) <- "VL.suppr." 
  
  
  classic.features <-   c(exp(growthrate), 
                          
                          men.prev,
                          women.prev,
                          
                          men.incid, # exp(as.numeric(men.incid)), 
                          women.incid, # exp(as.numeric(women.incid)), 
                          
                          pp.cp.6months.male.rels, # pp.cp.6months.female.rels,
                          
                          relsperpersonperyear, 
                          agegap.mean, agegap.med, agegap.sd,
                          
                          ART.cov.vector, VL.suppression.fraction)
  
  
  
  
  classic.features.names <- c("Pop.growthrate", 
                              
                              names(men.prev),
                              names(women.prev),
                              names(men.incid),
                              names(women.incid),
                              
                              "pp.cp.6months.male.rels",
                              
                              "relsperpersonperyear", 
                              "agegap.mean", "agegap.med", "agegap.sd", 
                              
                              paste0(names(ART.cov.vector)), paste0(names(VL.suppression.fraction)))
  
  names(classic.features) <- classic.features.names 
  
  # classic.features.num <- as.numeric(classic.features)
  
  return(classic.features)
  
}




