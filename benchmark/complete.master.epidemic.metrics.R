#' Computing transmission network characteristics: age mixing patterns, time trend incidence, and 
#' onward transmission distribution

#' @param datalist.agemix Data list of simpact output produced by \code{\link{readthedata()}}

#' @export


# The outputs of the function are:

# - temporal trend of incidence
# - age mixing statistics
# - mean, median, and standard deviation of onward transmissions



complete.master.epidemic.metrics <- function(datalist = datalist.agemix){
  
  
  source("/home/dniyukuri/lustre/benchmark_master_model/needed.functions.RSimpactHelp.R")
  
  # source("/home/david/benchmark_master_model/needed.functions.RSimpactHelp.R")
  
  
  datalist.agemix <- datalist
  
  simpact.trans.net <- transmission.network.builder(datalist = datalist.agemix, endpoint = 40)
  
  
  ##################################################
  ### Compute transmission network characteristics #
  ##################################################
  
  
  
  # 1. Incidence trend #
  ######################
  
  # 35 - 36
  
  
  # Women
  
  hiv.incid.lt25.women.35.36 <- incidence.calculator(datalist = datalist.agemix,
                                                     agegroup = c(15, 25),
                                                     timewindow = c(35, 36),
                                                     only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(2) %>%
    as.numeric()
  
  hiv.incid.25.30.women.35.36 <- incidence.calculator(datalist = datalist.agemix,
                                                      agegroup = c(25, 30),
                                                      timewindow = c(35, 36),
                                                      only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(2) %>%
    as.numeric()
  
  hiv.incid.30.35.women.35.36 <- incidence.calculator(datalist = datalist.agemix,
                                                      agegroup = c(30, 35),
                                                      timewindow = c(35, 36),
                                                      only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(2) %>%
    as.numeric()
  
  hiv.incid.35.40.women.35.36 <- incidence.calculator(datalist = datalist.agemix,
                                                      agegroup = c(35, 40),
                                                      timewindow = c(35, 36),
                                                      only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(2) %>%
    as.numeric()
  
  hiv.incid.40.45.women.35.36 <- incidence.calculator(datalist = datalist.agemix,
                                                      agegroup = c(40, 45),
                                                      timewindow = c(35, 36),
                                                      only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(2) %>%
    as.numeric()
  
  hiv.incid.45.50.women.35.36 <- incidence.calculator(datalist = datalist.agemix,
                                                      agegroup = c(45, 50),
                                                      timewindow = c(35, 36),
                                                      only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(2) %>%
    as.numeric()
  
  
  
  # Men
  
  hiv.incid.lt25.men.35.36 <- incidence.calculator(datalist = datalist.agemix,
                                                   agegroup = c(15, 25),
                                                   timewindow = c(35, 36),
                                                   only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(1) %>%
    as.numeric()
  
  hiv.incid.25.30.men.35.36 <- incidence.calculator(datalist = datalist.agemix,
                                                    agegroup = c(25, 30),
                                                    timewindow = c(35, 36),
                                                    only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(1) %>%
    as.numeric()
  
  hiv.incid.30.35.men.35.36 <- incidence.calculator(datalist = datalist.agemix,
                                                    agegroup = c(30, 35),
                                                    timewindow = c(35, 36),
                                                    only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(1) %>%
    as.numeric()
  
  hiv.incid.35.40.men.35.36 <- incidence.calculator(datalist = datalist.agemix,
                                                    agegroup = c(35, 40),
                                                    timewindow = c(35, 36),
                                                    only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(1) %>%
    as.numeric()
  
  hiv.incid.40.45.men.35.36 <- incidence.calculator(datalist = datalist.agemix,
                                                    agegroup = c(40, 45),
                                                    timewindow = c(35, 36),
                                                    only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(1) %>%
    as.numeric()
  
  hiv.incid.45.50.men.35.36 <- incidence.calculator(datalist = datalist.agemix,
                                                    agegroup = c(45, 50),
                                                    timewindow = c(35, 36),
                                                    only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(1) %>%
    as.numeric()
  
  
  
  
  
  # agegroup <- c("15-24", "25-29", "30-34", "35-39", "40-44", "45-49")
  
  men.incid.35.36 <- c(hiv.incid.lt25.men.35.36, hiv.incid.25.30.men.35.36, hiv.incid.30.35.men.35.36, hiv.incid.35.40.men.35.36, hiv.incid.40.45.men.35.36, hiv.incid.45.50.men.35.36)
  names(men.incid.35.36) <- c("incid.35.36.m.15.24", "incid.35.36.m.25.29", "incid.35.36.m.30.34", "incid.35.36.m.35.39", "incid.35.36.m.40.44", "incid.35.36.m.45.49")
  women.incid.35.36 <- c(hiv.incid.lt25.women.35.36, hiv.incid.25.30.women.35.36, hiv.incid.30.35.women.35.36, hiv.incid.35.40.women.35.36, hiv.incid.40.45.women.35.36, hiv.incid.45.50.women.35.36)
  names(women.incid.35.36) <- c("incid.35.36.w.15.24", "incid.35.36.w.25.29", "incid.35.36.w.30.34", "incid.35.36.w.35.39", "incid.35.36.w.40.44", "incid.35.36.w.45.49")        
  
  
  # 36 - 37
  
  # Women
  
  hiv.incid.lt25.women.36.37 <- incidence.calculator(datalist = datalist.agemix,
                                                     agegroup = c(15, 25),
                                                     timewindow = c(36, 37),
                                                     only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(2) %>%
    as.numeric()
  
  hiv.incid.25.30.women.36.37 <- incidence.calculator(datalist = datalist.agemix,
                                                      agegroup = c(25, 30),
                                                      timewindow = c(36, 37),
                                                      only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(2) %>%
    as.numeric()
  
  hiv.incid.30.35.women.36.37 <- incidence.calculator(datalist = datalist.agemix,
                                                      agegroup = c(30, 35),
                                                      timewindow = c(36, 37),
                                                      only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(2) %>%
    as.numeric()
  
  hiv.incid.35.40.women.36.37 <- incidence.calculator(datalist = datalist.agemix,
                                                      agegroup = c(35, 40),
                                                      timewindow = c(36, 37),
                                                      only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(2) %>%
    as.numeric()
  
  hiv.incid.40.45.women.36.37 <- incidence.calculator(datalist = datalist.agemix,
                                                      agegroup = c(40, 45),
                                                      timewindow = c(36, 37),
                                                      only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(2) %>%
    as.numeric()
  
  hiv.incid.45.50.women.36.37 <- incidence.calculator(datalist = datalist.agemix,
                                                      agegroup = c(45, 50),
                                                      timewindow = c(36, 37),
                                                      only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(2) %>%
    as.numeric()
  
  
  
  # Men
  
  hiv.incid.lt25.men.36.37 <- incidence.calculator(datalist = datalist.agemix,
                                                   agegroup = c(15, 25),
                                                   timewindow = c(36, 37),
                                                   only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(1) %>%
    as.numeric()
  
  hiv.incid.25.30.men.36.37 <- incidence.calculator(datalist = datalist.agemix,
                                                    agegroup = c(25, 30),
                                                    timewindow = c(36, 37),
                                                    only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(1) %>%
    as.numeric()
  
  hiv.incid.30.35.men.36.37 <- incidence.calculator(datalist = datalist.agemix,
                                                    agegroup = c(30, 35),
                                                    timewindow = c(36, 37),
                                                    only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(1) %>%
    as.numeric()
  
  hiv.incid.35.40.men.36.37 <- incidence.calculator(datalist = datalist.agemix,
                                                    agegroup = c(35, 40),
                                                    timewindow = c(36, 37),
                                                    only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(1) %>%
    as.numeric()
  
  hiv.incid.40.45.men.36.37 <- incidence.calculator(datalist = datalist.agemix,
                                                    agegroup = c(40, 45),
                                                    timewindow = c(36, 37),
                                                    only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(1) %>%
    as.numeric()
  
  hiv.incid.45.50.men.36.37 <- incidence.calculator(datalist = datalist.agemix,
                                                    agegroup = c(45, 50),
                                                    timewindow = c(36, 37),
                                                    only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(1) %>%
    as.numeric()
  
  
  
  
  
  # agegroup <- c("15-24", "25-29", "30-34", "35-39", "40-44", "45-49")
  
  men.incid.36.37 <- c(hiv.incid.lt25.men.36.37, hiv.incid.25.30.men.36.37, hiv.incid.30.35.men.36.37, hiv.incid.35.40.men.36.37, hiv.incid.40.45.men.36.37, hiv.incid.45.50.men.36.37)
  names(men.incid.36.37) <- c("incid.36.37.m.15.24", "incid.36.37.m.25.29", "incid.36.37.m.30.34", "incid.36.37.m.35.39", "incid.36.37.m.40.44", "incid.36.37.m.45.49")
  women.incid.36.37 <- c(hiv.incid.lt25.women.36.37, hiv.incid.25.30.women.36.37, hiv.incid.30.35.women.36.37, hiv.incid.35.40.women.36.37, hiv.incid.40.45.women.36.37, hiv.incid.45.50.women.36.37)
  names(women.incid.36.37) <- c("incid.36.37.w.15.24", "incid.36.37.w.25.29", "incid.36.37.w.30.34", "incid.36.37.w.35.39", "incid.36.37.w.40.44", "incid.36.37.w.45.49")        
  
  
  # 37 - 38
  
  # Women
  
  hiv.incid.lt25.women.37.38 <- incidence.calculator(datalist = datalist.agemix,
                                                     agegroup = c(15, 25),
                                                     timewindow = c(37, 38),
                                                     only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(2) %>%
    as.numeric()
  
  hiv.incid.25.30.women.37.38 <- incidence.calculator(datalist = datalist.agemix,
                                                      agegroup = c(25, 30),
                                                      timewindow = c(37, 38),
                                                      only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(2) %>%
    as.numeric()
  
  hiv.incid.30.35.women.37.38 <- incidence.calculator(datalist = datalist.agemix,
                                                      agegroup = c(30, 35),
                                                      timewindow = c(37, 38),
                                                      only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(2) %>%
    as.numeric()
  
  hiv.incid.35.40.women.37.38 <- incidence.calculator(datalist = datalist.agemix,
                                                      agegroup = c(35, 40),
                                                      timewindow = c(37, 38),
                                                      only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(2) %>%
    as.numeric()
  
  hiv.incid.40.45.women.37.38 <- incidence.calculator(datalist = datalist.agemix,
                                                      agegroup = c(40, 45),
                                                      timewindow = c(37, 38),
                                                      only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(2) %>%
    as.numeric()
  
  hiv.incid.45.50.women.37.38 <- incidence.calculator(datalist = datalist.agemix,
                                                      agegroup = c(45, 50),
                                                      timewindow = c(37, 38),
                                                      only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(2) %>%
    as.numeric()
  
  
  
  # Men
  
  hiv.incid.lt25.men.37.38 <- incidence.calculator(datalist = datalist.agemix,
                                                   agegroup = c(15, 25),
                                                   timewindow = c(37, 38),
                                                   only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(1) %>%
    as.numeric()
  
  hiv.incid.25.30.men.37.38 <- incidence.calculator(datalist = datalist.agemix,
                                                    agegroup = c(25, 30),
                                                    timewindow = c(37, 38),
                                                    only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(1) %>%
    as.numeric()
  
  hiv.incid.30.35.men.37.38 <- incidence.calculator(datalist = datalist.agemix,
                                                    agegroup = c(30, 35),
                                                    timewindow = c(37, 38),
                                                    only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(1) %>%
    as.numeric()
  
  hiv.incid.35.40.men.37.38 <- incidence.calculator(datalist = datalist.agemix,
                                                    agegroup = c(35, 40),
                                                    timewindow = c(37, 38),
                                                    only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(1) %>%
    as.numeric()
  
  hiv.incid.40.45.men.37.38 <- incidence.calculator(datalist = datalist.agemix,
                                                    agegroup = c(40, 45),
                                                    timewindow = c(37, 38),
                                                    only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(1) %>%
    as.numeric()
  
  hiv.incid.45.50.men.37.38 <- incidence.calculator(datalist = datalist.agemix,
                                                    agegroup = c(45, 50),
                                                    timewindow = c(37, 38),
                                                    only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(1) %>%
    as.numeric()
  
  
  
  
  
  # agegroup <- c("15-24", "25-29", "30-34", "35-39", "40-44", "45-49")
  
  men.incid.37.38 <- c(hiv.incid.lt25.men.37.38, hiv.incid.25.30.men.37.38, hiv.incid.30.35.men.37.38, hiv.incid.35.40.men.37.38, hiv.incid.40.45.men.37.38, hiv.incid.45.50.men.37.38)
  names(men.incid.37.38) <- c("incid.37.38.m.15.24", "incid.37.38.m.25.29", "incid.37.38.m.30.34", "incid.37.38.m.35.39", "incid.37.38.m.40.44", "incid.37.38.m.45.49")
  women.incid.37.38 <- c(hiv.incid.lt25.women.37.38, hiv.incid.25.30.women.37.38, hiv.incid.30.35.women.37.38, hiv.incid.35.40.women.37.38, hiv.incid.40.45.women.37.38, hiv.incid.45.50.women.37.38)
  names(women.incid.37.38) <- c("incid.37.38.w.15.24", "incid.37.38.w.25.29", "incid.37.38.w.30.34", "incid.37.38.w.35.39", "incid.37.38.w.40.44", "incid.37.38.w.45.49")        
  
  
  
  # 38 - 39
  
  # Women
  
  hiv.incid.lt25.women.38.39 <- incidence.calculator(datalist = datalist.agemix,
                                                     agegroup = c(15, 25),
                                                     timewindow = c(38, 39),
                                                     only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(2) %>%
    as.numeric()
  
  hiv.incid.25.30.women.38.39 <- incidence.calculator(datalist = datalist.agemix,
                                                      agegroup = c(25, 30),
                                                      timewindow = c(38, 39),
                                                      only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(2) %>%
    as.numeric()
  
  hiv.incid.30.35.women.38.39 <- incidence.calculator(datalist = datalist.agemix,
                                                      agegroup = c(30, 35),
                                                      timewindow = c(38, 39),
                                                      only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(2) %>%
    as.numeric()
  
  hiv.incid.35.40.women.38.39 <- incidence.calculator(datalist = datalist.agemix,
                                                      agegroup = c(35, 40),
                                                      timewindow = c(38, 39),
                                                      only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(2) %>%
    as.numeric()
  
  hiv.incid.40.45.women.38.39 <- incidence.calculator(datalist = datalist.agemix,
                                                      agegroup = c(40, 45),
                                                      timewindow = c(38, 39),
                                                      only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(2) %>%
    as.numeric()
  
  hiv.incid.45.50.women.38.39 <- incidence.calculator(datalist = datalist.agemix,
                                                      agegroup = c(45, 50),
                                                      timewindow = c(38, 39),
                                                      only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(2) %>%
    as.numeric()
  
  
  
  # Men
  
  hiv.incid.lt25.men.38.39 <- incidence.calculator(datalist = datalist.agemix,
                                                   agegroup = c(15, 25),
                                                   timewindow = c(38, 39),
                                                   only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(1) %>%
    as.numeric()
  
  hiv.incid.25.30.men.38.39 <- incidence.calculator(datalist = datalist.agemix,
                                                    agegroup = c(25, 30),
                                                    timewindow = c(38, 39),
                                                    only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(1) %>%
    as.numeric()
  
  hiv.incid.30.35.men.38.39 <- incidence.calculator(datalist = datalist.agemix,
                                                    agegroup = c(30, 35),
                                                    timewindow = c(38, 39),
                                                    only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(1) %>%
    as.numeric()
  
  hiv.incid.35.40.men.38.39 <- incidence.calculator(datalist = datalist.agemix,
                                                    agegroup = c(35, 40),
                                                    timewindow = c(38, 39),
                                                    only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(1) %>%
    as.numeric()
  
  hiv.incid.40.45.men.38.39 <- incidence.calculator(datalist = datalist.agemix,
                                                    agegroup = c(40, 45),
                                                    timewindow = c(38, 39),
                                                    only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(1) %>%
    as.numeric()
  
  hiv.incid.45.50.men.38.39 <- incidence.calculator(datalist = datalist.agemix,
                                                    agegroup = c(45, 50),
                                                    timewindow = c(38, 39),
                                                    only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(1) %>%
    as.numeric()
  
  
  
  
  
  # agegroup <- c("15-24", "25-29", "30-34", "35-39", "40-44", "45-49")
  
  men.incid.38.39 <- c(hiv.incid.lt25.men.38.39, hiv.incid.25.30.men.38.39, hiv.incid.30.35.men.38.39, hiv.incid.35.40.men.38.39, hiv.incid.40.45.men.38.39, hiv.incid.45.50.men.38.39)
  names(men.incid.38.39) <- c("incid.38.39.m.15.24", "incid.38.39.m.25.29", "incid.38.39.m.30.34", "incid.38.39.m.35.39", "incid.38.39.m.40.44", "incid.38.39.m.45.49")
  women.incid.38.39 <- c(hiv.incid.lt25.women.38.39, hiv.incid.25.30.women.38.39, hiv.incid.30.35.women.38.39, hiv.incid.35.40.women.38.39, hiv.incid.40.45.women.38.39, hiv.incid.45.50.women.38.39)
  names(women.incid.38.39) <- c("incid.38.39.w.15.24", "incid.38.39.w.25.29", "incid.38.39.w.30.34", "incid.38.39.w.35.39", "incid.38.39.w.40.44", "incid.38.39.w.45.49")        
  
  
  
  # 39 - 40
  
  # Women
  
  hiv.incid.lt25.women.39.40 <- incidence.calculator(datalist = datalist.agemix,
                                                     agegroup = c(15, 25),
                                                     timewindow = c(39, 40),
                                                     only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(2) %>%
    as.numeric()
  
  hiv.incid.25.30.women.39.40 <- incidence.calculator(datalist = datalist.agemix,
                                                      agegroup = c(25, 30),
                                                      timewindow = c(39, 40),
                                                      only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(2) %>%
    as.numeric()
  
  hiv.incid.30.35.women.39.40 <- incidence.calculator(datalist = datalist.agemix,
                                                      agegroup = c(30, 35),
                                                      timewindow = c(39, 40),
                                                      only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(2) %>%
    as.numeric()
  
  hiv.incid.35.40.women.39.40 <- incidence.calculator(datalist = datalist.agemix,
                                                      agegroup = c(35, 40),
                                                      timewindow = c(39, 40),
                                                      only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(2) %>%
    as.numeric()
  
  hiv.incid.40.45.women.39.40 <- incidence.calculator(datalist = datalist.agemix,
                                                      agegroup = c(40, 45),
                                                      timewindow = c(39, 40),
                                                      only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(2) %>%
    as.numeric()
  
  hiv.incid.45.50.women.39.40 <- incidence.calculator(datalist = datalist.agemix,
                                                      agegroup = c(45, 50),
                                                      timewindow = c(39, 40),
                                                      only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(2) %>%
    as.numeric()
  
  
  
  # Men
  
  hiv.incid.lt25.men.39.40 <- incidence.calculator(datalist = datalist.agemix,
                                                   agegroup = c(15, 25),
                                                   timewindow = c(39, 40),
                                                   only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(1) %>%
    as.numeric()
  
  hiv.incid.25.30.men.39.40 <- incidence.calculator(datalist = datalist.agemix,
                                                    agegroup = c(25, 30),
                                                    timewindow = c(39, 40),
                                                    only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(1) %>%
    as.numeric()
  
  hiv.incid.30.35.men.39.40 <- incidence.calculator(datalist = datalist.agemix,
                                                    agegroup = c(30, 35),
                                                    timewindow = c(39, 40),
                                                    only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(1) %>%
    as.numeric()
  
  hiv.incid.35.40.men.39.40 <- incidence.calculator(datalist = datalist.agemix,
                                                    agegroup = c(35, 40),
                                                    timewindow = c(39, 40),
                                                    only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(1) %>%
    as.numeric()
  
  hiv.incid.40.45.men.39.40 <- incidence.calculator(datalist = datalist.agemix,
                                                    agegroup = c(40, 45),
                                                    timewindow = c(39, 40),
                                                    only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(1) %>%
    as.numeric()
  
  hiv.incid.45.50.men.39.40 <- incidence.calculator(datalist = datalist.agemix,
                                                    agegroup = c(45, 50),
                                                    timewindow = c(39, 40),
                                                    only.active = "No") %>%
    dplyr::select(incidence) %>%
    dplyr::slice(1) %>%
    as.numeric()
  
  
  
  
  
  # agegroup <- c("15-24", "25-29", "30-34", "35-39", "40-44", "45-49")
  
  men.incid.39.40 <- c(hiv.incid.lt25.men.39.40, hiv.incid.25.30.men.39.40, hiv.incid.30.35.men.39.40, hiv.incid.35.40.men.39.40, hiv.incid.40.45.men.39.40, hiv.incid.45.50.men.39.40)
  names(men.incid.39.40) <- c("incid.39.40.m.15.24", "incid.39.40.m.25.29", "incid.39.40.m.30.34", "incid.39.40.m.35.39", "incid.39.40.m.40.44", "incid.39.40.m.45.49")
  women.incid.39.40 <- c(hiv.incid.lt25.women.39.40, hiv.incid.25.30.women.39.40, hiv.incid.30.35.women.39.40, hiv.incid.35.40.women.39.40, hiv.incid.40.45.women.39.40, hiv.incid.45.50.women.39.40)
  names(women.incid.39.40) <- c("incid.39.40.w.15.24", "incid.39.40.w.25.29", "incid.39.40.w.30.34", "incid.39.40.w.35.39", "incid.39.40.w.40.44", "incid.39.40.w.45.49")        
  
  
  METRICS.incidence.men <- c(men.incid.35.36, men.incid.36.37, men.incid.37.38, men.incid.38.39, men.incid.39.40)
  
  METRICS.incidence.women <- c(women.incid.35.36, women.incid.36.37, women.incid.37.38, women.incid.38.39, women.incid.39.40)  
  
  # 2. Age mixing #
  ################
  
  # (i) Age mixing in relationships
  
  # 
  
  agemix.rels.df <- agemix.df.maker(datalist.agemix)
  
  # 
  agemix.model <- pattern.modeller(dataframe = agemix.rels.df,
                                   agegroup = c(15, 50),
                                   timepoint = 40, # datalist.agemix$itable$population.simtime[1],
                                   timewindow = 5)#1)#3)
  # 
  # # men.lme <- tryCatch(agemixing.lme.fitter(data = dplyr::filter(agemix.model[[1]], Gender =="male")),
  # #                     error = agemixing.lme.errFunction) # Returns an empty list if the lme model can't be fitted
  #
  # men.lmer <- ampmodel(data = dplyr::filter(agemix.model[[1]], Gender =="male"))
  
  data = dplyr::filter(agemix.model[[1]], Gender =="male")
  
  if( nrow(data) > length(unique(data$ID)) & length(unique(data$ID)) > 1 ){
    
    men.lmer <- lmer(pagerelform ~ agerelform0 + (1 | ID),
                     data = dplyr::filter(agemix.model[[1]], Gender =="male"),
                     REML = TRUE,
                     control=lmerControl(check.nobs.vs.nlev = "ignore",
                                         check.nobs.vs.rankZ = "ignore",
                                         check.nobs.vs.nRE="ignore"))
    
    bignumber <- NA # let's try if NA works (instead of 9999 for example)
    AAD.male <- ifelse(length(men.lmer) > 0, mean(dplyr::filter(agemix.model[[1]], Gender =="male")$AgeGap), bignumber)
    SDAD.male <- ifelse(length(men.lmer) > 0, sd(dplyr::filter(agemix.model[[1]], Gender =="male")$AgeGap), bignumber)
    #powerm <- ifelse(length(men.lme) > 0, as.numeric(attributes(men.lme$apVar)$Pars["varStruct.power"]), bignumber)
    slope.male <- ifelse(length(men.lmer) > 0, summary(men.lmer)$coefficients[2, 1], bignumber) #summary(men.lmer)$tTable[2, 1], bignumber)
    WSD.male <- ifelse(length(men.lmer) > 0, summary(men.lmer)$sigma, bignumber) #WVAD.base <- ifelse(length(men.lme) > 0, men.lme$sigma^2, bignumber)
    
    BSD.male <- ifelse(length(men.lmer) > 0, bvar(men.lmer), bignumber) # Bad name for the function because it actually extracts between subject standard deviation # BVAD <- ifelse(length(men.lmer) > 0, getVarCov(men.lme)[1,1], bignumber)
    
    intercept.male <- ifelse(length(men.lmer) > 0, summary(men.lmer)$coefficients[1,1] - 15, bignumber)
    
    # c(AAD.male, SDAD.male, slope.male, WSD.male, BSD.male, intercept.male)
    
    ## AAD: average age difference across all relationship
    ## VAD: variance of these age differences
    ## SDAD: standard deviation of age differences
    ## BSD: between-subject standard deviation of age differences
    
    mix.rels.dat <- c(AAD.male, SDAD.male, slope.male, WSD.male, BSD.male, intercept.male)
    
    names(mix.rels.dat) <- c("AAD.male", "SDAD.male", "slope.male", "WSD.male", "BSD.male", "intercept.male")
    
  }else{
    
    mix.rels.dat <- rep(NA, 6)
    
    names(mix.rels.dat) <- c("AAD.male", "SDAD.male", "slope.male", "WSD.male", "BSD.male", "intercept.male")
    
  }
  
  
  
  METRICS.LMEM.rels.age.mix <-  mix.rels.dat 
  
  
  
  
  # 3. Onward transmissions # starting ART interventions time
  ###########################
  
  transm.count <- onwardtransmissions.dat(datalist = datalist.agemix, 
                                          trans.network = simpact.trans.net,
                                          time.window=c(35, 40))
  
  
  
  
  transm.average <- mean(transm.count)
  
  transm.median <- median(transm.count) # add
  
  transm.sd <- sd(transm.count) # add
  
  METRICS.onwardtransmissions <- c(transm.average, transm.median, transm.sd)
  
  
  
  epi.Metrics <- c(as.numeric(METRICS.incidence.men),
                   as.numeric(METRICS.incidence.women),
                   
                   as.numeric(METRICS.LMEM.rels.age.mix),
                   
                   METRICS.onwardtransmissions)
  
  
  
  metric.names <- c(paste0("metr.",c("incid.35.36.m.15.24", "incid.35.36.m.25.29", "incid.35.36.m.30.34", "incid.35.36.m.35.39",
                                     "incid.35.36.m.40.44", "incid.35.36.m.45.49", "incid.36.37.m.15.24", "incid.36.37.m.25.29",
                                     "incid.36.37.m.30.34", "incid.36.37.m.35.39", "incid.36.37.m.40.44", "incid.36.37.m.45.49",
                                     "incid.37.38.m.15.24", "incid.37.38.m.25.29", "incid.37.38.m.30.34", "incid.37.38.m.35.39",
                                     "incid.37.38.m.40.44", "incid.37.38.m.45.49", "incid.38.39.m.15.24", "incid.38.39.m.25.29",
                                     "incid.38.39.m.30.34", "incid.38.39.m.35.39", "incid.38.39.m.40.44", "incid.38.39.m.45.49",
                                     "incid.39.40.m.15.24", "incid.39.40.m.25.29", "incid.39.40.m.30.34", "incid.39.40.m.35.39",
                                     "incid.39.40.m.40.44", "incid.39.40.m.45.49",
                                     
                                     "incid.35.36.w.15.24", "incid.35.36.w.25.29", "incid.35.36.w.30.34", "incid.35.36.w.35.39",
                                     "incid.35.36.w.40.44", "incid.35.36.w.45.49", "incid.36.37.w.15.24", "incid.36.37.w.25.29",
                                     "incid.36.37.w.30.34", "incid.36.37.w.35.39", "incid.36.37.w.40.44", "incid.36.37.w.45.49",
                                     "incid.37.38.w.15.24", "incid.37.38.w.25.29", "incid.37.38.w.30.34", "incid.37.38.w.35.39",
                                     "incid.37.38.w.40.44", "incid.37.38.w.45.49", "incid.38.39.w.15.24", "incid.38.39.w.25.29",
                                     "incid.38.39.w.30.34", "incid.38.39.w.35.39", "incid.38.39.w.40.44", "incid.38.39.w.45.49",
                                     "incid.39.40.w.15.24", "incid.39.40.w.25.29", "incid.39.40.w.30.34", "incid.39.40.w.35.39",
                                     "incid.39.40.w.40.44", "incid.39.40.w.45.49")),
                    
                    paste0("metr.",c("AAD.male", "SDAD.male", "slope.male", "WSD.male", "BSD.male", "intercept.male")), 
                    
                    "metr.transm.av",
                    "metr.transm.med", 
                    "metr.transm.sd") # c(transm.average, transm.median, transm.sd)
  
  
  names(epi.Metrics) <- metric.names
  
  return(epi.Metrics)
  
  
}

