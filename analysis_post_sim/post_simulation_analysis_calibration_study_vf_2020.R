

# Post simulation analysis for

# Inferring transmission network characteristics usingagent-based models calibrated 
# to epidemiological,behavioural, demographic and phylogeneticsummary statistics

set.seed(777)

# Packages ------------

library(tidyverse)
library(kableExtra)
library(robustbase)
library(reshape2)

# Functions --------------

quant.mean <- function(input){
  
  input <- na.omit(input)
  
  quantirles.v <- quantile(input, probs = seq(0, 1, 0.25))
  
  mean.i <- mean(input)
  
  quantirles.25.50.75 <- as.numeric(quantirles.v)[2:4]
  
  quantirles.25.mean.75 <- c(quantirles.25.50.75[1], mean.i, quantirles.25.50.75[3])
  
  # return(quantirles.25.50.75)
  
  return(quantirles.25.mean.75)
  
}


source("/home/david/benchmark_master_model/analysis.complete.master.epidemic.metrics.R")




# Read benchmark data  ------------------


benckmark.df <- read.csv("/home/david/benchmark_master_model/results.benchmark.epi.mm.stats_seed_1.csv")


epi.metric.df <- benckmark.df %>%
  select(contains("metr.")) 

epi.metric.df <- na.omit(epi.metric.df)

# The data is normal if the p-value is above 0.05.
shapiro.test(epi.metric.df$metr.transm.av)

name <- vector()
p <- vector()
for(i in 1:ncol(epi.metric.df)){
  
  name.i <- names(epi.metric.df)[i]
  p.i <- shapiro.test(epi.metric.df[,i])[[2]]
  name <- c(name, name.i)
  p <- c(p, p.i)
}

p_df <- data.frame(name, p)


# Prevalence at 40 ---------------
# extracted for epi-behav SS

hiv.prev <- benckmark.df %>% 
  select(starts_with("prev.")) 

prev.m.15.24 <- quant.mean(hiv.prev$prev.m.15.24)
prev.m.25.29 <- quant.mean(hiv.prev$prev.m.25.29)
prev.m.30.34 <- quant.mean(hiv.prev$prev.m.30.34)
prev.m.35.39 <- quant.mean(hiv.prev$prev.m.35.39)
prev.m.40.44 <- quant.mean(hiv.prev$prev.m.40.44)
prev.m.45.49 <- quant.mean(hiv.prev$prev.m.45.49)

prev.w.15.24 <- quant.mean(hiv.prev$prev.w.15.24)
prev.w.25.29 <- quant.mean(hiv.prev$prev.w.25.29)
prev.w.30.34 <- quant.mean(hiv.prev$prev.w.30.34)
prev.w.35.39 <- quant.mean(hiv.prev$prev.w.35.39)
prev.w.40.44 <- quant.mean(hiv.prev$prev.w.40.44)
prev.w.45.49 <- quant.mean(hiv.prev$prev.w.45.49)

val.prev.F <- c(prev.m.15.24[2], prev.m.25.29[2], prev.m.30.34[2], prev.m.35.39[2], prev.m.40.44[2], prev.m.45.49[2], 
                prev.w.15.24[2], prev.w.25.29[2], prev.w.30.34[2], prev.w.35.39[2], prev.w.40.44[2], prev.w.45.49[2])

val.prev.L <- c(prev.m.15.24[1], prev.m.25.29[1], prev.m.30.34[1], prev.m.35.39[1], prev.m.40.44[1], prev.m.45.49[1], 
                prev.w.15.24[1], prev.w.25.29[1], prev.w.30.34[1], prev.w.35.39[1], prev.w.40.44[1], prev.w.45.49[1])

val.prev.U <- c(prev.m.15.24[3], prev.m.25.29[3], prev.m.30.34[3], prev.m.35.39[3], prev.m.40.44[3], prev.m.45.49[3], 
                prev.w.15.24[3], prev.w.25.29[3], prev.w.30.34[3], prev.w.35.39[3], prev.w.40.44[3], prev.w.45.49[3])


par.gender <- c(rep("men", 6), rep("women", 6))
agegroup <- rep(c("15-24", "25-29", "30-34", "35-39", "40-44", "45-49"), 2)

val.prev <- data.frame(agegroup, val.prev.L, val.prev.F, val.prev.U, par.gender)

val.prev.tab <- val.prev

names(val.prev.tab) <- c("age_group", "lower.Q1", "median", "upper.Q3", "gender")

val.prev.tab %>% 
  kable() %>% 
  kable_styling("striped") 


write.csv(val.prev.tab, file = "/home/david/benchmark_master_model/results/000_master.model.table.prevalence.csv")


names(val.prev) <- c("age_group", "lower.Q1", "mean", "upper.Q3", "gender")


plot.prev.men.women <- ggplot(val.prev, aes(x=age_group, y=mean, colour=gender, group = gender)) + 
  # geom_errorbar(aes(ymin=lower.Q1, ymax=upper.Q3), width=.1) +
  geom_line(size=1) +
  geom_point(size=2) + 
  xlab("age group") + ylab("HIV prevalence")+
  theme(legend.position="bottom")


# print(plot.prev.men.women)

ggsave(filename = "000_master.model.prevalence.at.40.pdf",
       plot = plot.prev.men.women,
       path = "/home/david/benchmark_master_model/results",
       width = 16, height = 12, units = "cm")



# Incidence at 39:40 ------------------
# extracted for epi-behav SS

hiv.inc <- benckmark.df %>% 
  select(starts_with("incid.")) 


incid.m.15.24 <- quant.mean(hiv.inc$incid.m.15.24)
incid.m.25.29 <- quant.mean(hiv.inc$incid.m.25.29)
incid.m.30.34 <- quant.mean(hiv.inc$incid.m.30.34)
incid.m.35.39 <- quant.mean(hiv.inc$incid.m.35.39)
incid.m.40.44 <- quant.mean(hiv.inc$incid.m.40.44)
incid.m.45.49 <- quant.mean(hiv.inc$incid.m.45.49)

incid.w.15.24 <- quant.mean(hiv.inc$incid.w.15.24)
incid.w.25.29 <- quant.mean(hiv.inc$incid.w.25.29)
incid.w.30.34 <- quant.mean(hiv.inc$incid.w.30.34)
incid.w.35.39 <- quant.mean(hiv.inc$incid.w.35.39)
incid.w.40.44 <- quant.mean(hiv.inc$incid.w.40.44)
incid.w.45.49 <- quant.mean(hiv.inc$incid.w.45.49)


val.incid.F <- c(incid.m.15.24[2], incid.m.25.29[2], incid.m.30.34[2], incid.m.35.39[2], incid.m.40.44[2], incid.m.45.49[2], 
                 incid.w.15.24[2], incid.w.25.29[2], incid.w.30.34[2], incid.w.35.39[2], incid.w.40.44[2], incid.w.45.49[2])

val.incid.L <- c(incid.m.15.24[1], incid.m.25.29[1], incid.m.30.34[1], incid.m.35.39[1], incid.m.40.44[1], incid.m.45.49[1], 
                 incid.w.15.24[1], incid.w.25.29[1], incid.w.30.34[1], incid.w.35.39[1], incid.w.40.44[1], incid.w.45.49[1])

val.incid.U <- c(incid.m.15.24[3], incid.m.25.29[3], incid.m.30.34[3], incid.m.35.39[3], incid.m.40.44[3], incid.m.45.49[3], 
                 incid.w.15.24[3], incid.w.25.29[3], incid.w.30.34[3], incid.w.35.39[3], incid.w.40.44[3], incid.w.45.49[3])


par.gender <- c(rep("men", 6), rep("women", 6))
agegroup <- rep(c("15-24", "25-29", "30-34", "35-39", "40-44", "45-49"), 2)

val.incid <- data.frame(agegroup, val.incid.L, val.incid.F, val.incid.U, par.gender)

val.incid.tab <- val.incid

names(val.incid.tab) <- c("age_group", "lower.Q1", "mean", "upper.Q3", "gender")

val.incid.tab %>% 
  kable() %>% 
  kable_styling("striped") 


write.csv(val.incid.tab, file = "/home/david/benchmark_master_model/results/000_master.model.table..incidence.csv")


names(val.incid) <- c("age_group", "lower.Q1", "mean", "upper.Q3", "gender")

plot.incid.men.women <- ggplot(val.incid, aes(x=age_group, y=mean, colour=gender, group = gender)) + 
  # geom_errorbar(aes(ymin=lower.Q1, ymax=upper.Q3), width=.1) +
  geom_line(size=1) +
  geom_point(size=2) + 
  xlab("age group") + ylab("HIV incidence")+
  theme(legend.position="bottom")


# print(plot.incid.men.women)

ggsave(filename = "000_master.model.incidence.at.40.pdf",
       plot = plot.incid.men.women,
       path = "/home/david/benchmark_master_model/results",
       width = 16, height = 12, units = "cm")




# Master model predictions outputs -------------- 

# (i) Incidence trend
# (ii) Age mixing in partnerships 
# (iii) Onwards transmissions 



analysis.complete.master.epidemic.metrics(keyword = "000_master.model.",  # ------------------------------------
                                          epi.metric.df = epi.metric.df)





# HIV epidemiology and sexual behaviour summary statistics  ---------------


benchmark_features_hiv_behavior_data <- dplyr::select(benckmark.df, 
                                                      
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

mean.epi.behav.features <- sapply(benchmark_features_hiv_behavior_data, mean)

names_ss_epi_behav <- names(benchmark_features_hiv_behavior_data)



# Phylogenetic summary statistics  ---------------


benchmark_features_phylo_data <- dplyr::select(benckmark.df, 
                                               cov.MCAR.100.mean_nodeHeights,   
                                               cov.MCAR.100.median_nodeHeights, 
                                               cov.MCAR.100.sd_nodeHeights,     
                                               cov.MCAR.100.colless ,          
                                               cov.MCAR.100.sackin,             
                                               cov.MCAR.100.min_BL,             
                                               cov.MCAR.100.1Q_BL,              
                                               cov.MCAR.100.median_BL,         
                                               cov.MCAR.100.mean_BL,            
                                               cov.MCAR.100.3Q_BL,              
                                               cov.MCAR.100.max_BL,             
                                               cov.MCAR.100.mean.tipsDepths,   
                                               cov.MCAR.100.median.tipsDepths,  
                                               cov.MCAR.100.sd.tipsDepths,      
                                               cov.MCAR.100.mean.nodesDepths,   
                                               cov.MCAR.100.median.nodesDepths,
                                               cov.MCAR.100.sd.nodesDepths,     
                                               cov.MCAR.100.maxHeight)

names(benchmark_features_phylo_data) <- c("mean_nodeHeights", "median_nodeHeights", "sd_nodeHeights",
                                          "colless", "sackin",
                                          "min_BL", "1Q_BL", "median_BL", "mean_BL", "3Q_BL", "max_BL",
                                          "mean.tipsDepths", "median.tipsDepths", "sd.tipsDepths",
                                          "mean.nodesDepths", "median.nodesDepths", "sd.nodesDepths",
                                          "maxHeight")

benchmark_features_phylo_data <- na.omit(benchmark_features_phylo_data)

mean.phylo_MCAR_100.features <- sapply(benchmark_features_phylo_data, mean)

names_ss_phylo <- names(mean.phylo_MCAR_100.features)



# Input parameter vector ---------------



# dissolution.alpha_0 = - 0.52
# dissolution.alpha_4 = -0.05
# formation.hazard.agegapry.baseline = 2
# person.agegap.man.dist.normal.mu  = 10 # and ~.woman.~
# person.agegap.man.dist.normal.sigma = 5 # and ~.woman.~
# formation.hazard.agegapry.gap_agescale_man = 0.25 # and ~.woman.~
# formation.hazard.agegapry.numrel_man = -0.3 # and ~.woman.~
# formation.hazard.agegapry.numrel_diff = -0.1
# hivtransmission.param.a = -1
# hivtransmission.param.b = -90
# hivtransmission.param.c = 0.5
# hivtransmission.param.f1 = 0.048
# hivtransmission.param.f2 = -0.14


params.names <- c("dissolution.alpha_0", "dissolution.alpha_4", "formation.hazard.agegapry.baseline",
                  "person.agegap.man.dist.normal.mu", "formation.hazard.agegapry.gap_agescale_man",
                  "person.agegap.man.dist.normal.sigma", "formation.hazard.agegapry.numrel_man",
                  "formation.hazard.agegapry.numrel_diff",
                  
                  "hivtransmission.param.a", "hivtransmission.param.b",
                  "hivtransmission.param.c", "hivtransmission.param.f1",
                  "hivtransmission.param.f2")

inputvector <- c(-0.52, -0.05, 2, 10, 5, 0.25, -0.3, -0.1,
                 -1, -90, 0.5, 0.048, -0.14) 


names(inputvector) <- params.names





# Predictions of age-mixing, onward transmission, and temporal incidence --------------



# Calibration - 1 - hiv epi and sexual behaviour - outputs and post-calibration results -----------



# After generating parameters with Latin Hypercube

# Outputs after running simulations and compute summary statistics

raw_simulation_calibration_hiv_epi_behavior <- readRDS("/home/david/benchmark_master_model/calibration_epi_behav/raw_simulation_calibration_hiv_epi_behavior.RDS")



# Calibration outputs: given the summary statistics generated by the parameters 

# which ones are best fit

abc_res_calibration_hiv_epi_behavior <- readRDS("/home/david/benchmark_master_model/calibration_epi_behav/abc_res_calibration_hiv_epi_behavior.RDS")




# Summarise the selected parameters

params.table_raw_simulation_calibration_hiv_epi_behavior <- readRDS("/home/david/benchmark_master_model/calibration_epi_behav/params.table_raw_simulation_calibration_hiv_epi_behavior.RDS")


rownames(params.table_raw_simulation_calibration_hiv_epi_behavior) <- params.names

params.table_raw_simulation_calibration_hiv_epi_behavior <- round(params.table_raw_simulation_calibration_hiv_epi_behavior, digits = 5)

params.table_raw_simulation_calibration_hiv_epi_behavior %>% 
  kable() %>% 
  kable_styling("striped") 


# Outputs after running the default model with new parameters (after calibration) ----------

epi.metrics_after_calibration_run_params_hiv_epi_behavior <- readRDS("/home/david/benchmark_master_model/calibration_epi_behav/epi.metrics_after_calibration_run_params.RDS")

x <- as.data.frame(epi.metrics_after_calibration_run_params_hiv_epi_behavior)

y <- x[,-(dim(x)[2])]

epi.metrics_after_calibration_run_params_hiv_epi_behavior <- y


epi.metrics_after_calibration_run_params_hiv_epi_behavior <- na.omit(epi.metrics_after_calibration_run_params_hiv_epi_behavior)



# Default model after calibration - hiv epi and sexual behavior

# (i) Incidence trend
# (ii) Age mixing in partnerships 
# (iii) Onwards transmissions 



analysis.complete.master.epidemic.metrics(keyword = "111_classic.calibration.parms.",  # ------------------
                                          epi.metric.df = epi.metrics_after_calibration_run_params_hiv_epi_behavior)




# Outputs after running the default model with median values of new parameters (after calibration) ----------


# Outputs after running the default model with median value of new parameters

epi.metrics_after_calibration_run_median_params_epi_behavior <- readRDS("/home/david/benchmark_master_model/calibration_epi_behav/epi.metrics_after_calibration_run_median_params.RDS")

x <- as.data.frame(epi.metrics_after_calibration_run_median_params_epi_behavior)

y <- x[,-(dim(x)[2])]

epi.metrics_after_calibration_run_median_params_epi_behavior <- y


epi.metrics_after_calibration_run_median_params_epi_behavior <- na.omit(epi.metrics_after_calibration_run_median_params_epi_behavior)


analysis.complete.master.epidemic.metrics(keyword = "111_2_classic.calibration.median.parms.", # --------------------------
                                          epi.metric.df = epi.metrics_after_calibration_run_median_params_epi_behavior)








# Calibration - 2 - phylo 100 - outputs and post-calibration results -----------


# After generating parameters with Latin Hypercube


# Outputs after running simulaitons and compute summary statistics 


raw_simulation_calibration_phylo_100 <- readRDS("/home/david/benchmark_master_model/calibration_phylo/raw_simulation_calibration_phylo.RDS")



# Calibration outputs: given the summary statistics generated by the parameters 

# which ones are best fit

abc_res_calibration_phylo_100 <- readRDS("/home/david/benchmark_master_model/calibration_phylo/abc_res_calibration_phylo.RDS")


# Summarise the selected parameters

params.table_raw_simulation_calibration_phylo_100 <- readRDS("/home/david/benchmark_master_model/calibration_phylo/params.table_raw_simulation_calibration_phylo.RDS")


rownames(params.table_raw_simulation_calibration_phylo_100) <- params.names


params.table_raw_simulation_calibration_phylo_100 <- round(params.table_raw_simulation_calibration_phylo_100, digits = 5)


params.table_raw_simulation_calibration_phylo_100 %>% 
  kable() %>% 
  kable_styling("striped") 



# Outputs after running the default model with new parameters (after calibration) ----------


epi.metrics_after_calibration_run_params_phylo_100 <- readRDS("/home/david/benchmark_master_model/calibration_phylo/epi.metrics_after_calibration_run_params_phylo.RDS")

x <- as.data.frame(epi.metrics_after_calibration_run_params_phylo_100)

y <- x[,-(dim(x)[2])]

epi.metrics_after_calibration_run_params_phylo_100 <- y

epi.metrics_after_calibration_run_params_phylo_100 <- na.omit(epi.metrics_after_calibration_run_params_phylo_100)



# Default model after calibration  - phylo 70

# (i) Incidence trend
# (ii) Age mixing in partnerships 
# (iii) Onwards transmissions 



analysis.complete.master.epidemic.metrics(keyword = "222_phylo.100.calibration.params.",  # -----------------------------
                                          epi.metric.df = epi.metrics_after_calibration_run_params_phylo_100)



# Outputs after running the default model with median values of new parameters (after calibration) ----------


epi.metrics_after_calibration_run_median_params_phylo <- readRDS("/home/david/benchmark_master_model/calibration_phylo/epi.metrics_after_calibration_run_median_params_phylo.RDS")

x <- as.data.frame(epi.metrics_after_calibration_run_median_params_phylo)

y <- x[,-(dim(x)[2])]


epi.metrics_after_calibration_run_median_params_phylo <- y

epi.metrics_after_calibration_run_median_params_phylo <- na.omit(epi.metrics_after_calibration_run_median_params_phylo)



analysis.complete.master.epidemic.metrics(keyword = "222_2_phylo.100.calibration.median.params.",  # -----------------------------
                                          epi.metric.df = epi.metrics_after_calibration_run_median_params_phylo)







# Calibration - 3 - combined phylo  70 and hiv epi and sexual behaviour - outputs and post-calibration results -----------




# After generating parameters with Latin Hypercube

# Outputs after running simulations and compute summary statistics

raw_simulation_calibration_combined_phylo_100 <- readRDS("/home/david/benchmark_master_model/calibration_combined/raw_simulation_calibration_combined_phylo.RDS")



# Calibration outputs: given the summary statistics generated by the parameters 

# which ones are best fit

abc_res_calibration_combined_phylo_100 <- readRDS("/home/david/benchmark_master_model/calibration_combined/abc_res_calibration_combined_phylo.RDS")



# Summarise the selected parameters

params.table_raw_simulation_calibration_combined_phylo_100 <- readRDS("/home/david/benchmark_master_model/calibration_combined/params.table_raw_simulation_calibration_combined_phylo.RDS")


rownames(params.table_raw_simulation_calibration_combined_phylo_100) <- params.names


params.table_raw_simulation_calibration_combined_phylo_100 <- round(params.table_raw_simulation_calibration_combined_phylo_100, digits = 5)

params.table_raw_simulation_calibration_combined_phylo_100 %>% 
  kable() %>% 
  kable_styling("striped") 


# Outputs after running the default model with new parameters (after calibration) ----------


epi.metrics_after_calibration_run_params_combined <- readRDS("/home/david/benchmark_master_model/calibration_combined/epi.metrics_after_calibration_run_params_combined_phylo.RDS")


x <- as.data.frame(epi.metrics_after_calibration_run_params_combined)

y <- x[,-(dim(x)[2])]

epi.metrics_after_calibration_run_params_combined <- y


epi.metrics_after_calibration_run_params_combined <- na.omit(epi.metrics_after_calibration_run_params_combined)


# Default model after calibration  - combined with phylo 70 

# (i) Incidence trend
# (ii) Age mixing in partnerships 
# (iii) Onwards transmissions 



analysis.complete.master.epidemic.metrics(keyword = "333_combined.phylo.100.calibration.parms.", # --------------------
                                          epi.metric.df = epi.metrics_after_calibration_run_params_combined)



# Outputs after running the default model with median values of new parameters (after calibration) ----------


epi.metrics_after_calibration_run_median_params_combined <- readRDS("/home/david/benchmark_master_model/calibration_combined/epi.metrics_after_calibration_run_median_params_combined_phylo.RDS")

x <- as.data.frame(epi.metrics_after_calibration_run_median_params_combined)

y <- x[,-(dim(x)[2])]

epi.metrics_after_calibration_run_median_params_combined <- y

epi.metrics_after_calibration_run_median_params_combined <- na.omit(epi.metrics_after_calibration_run_median_params_combined)




analysis.complete.master.epidemic.metrics(keyword = "333_2_combined.phylo.100.calibration.median.parms.", # --------------------
                                          epi.metric.df = epi.metrics_after_calibration_run_median_params_combined)



# Produce tables -------------


# Table 1: Selected parameters -------------


# 
# params.table_raw_simulation_calibration_hiv_epi_behavior
# 
# params.table_raw_simulation_calibration_phylo_100
# 
# params.table_raw_simulation_calibration_combined_phylo_100



params.table_hiv_epi_behavior <- paste(round(params.table_raw_simulation_calibration_hiv_epi_behavior[,7], digits = 3), "(", # column 7 is median values
                                       round(params.table_raw_simulation_calibration_hiv_epi_behavior[,1], digits = 3), ",",
                                       round(params.table_raw_simulation_calibration_hiv_epi_behavior[,3], digits = 3), ")")


params.table_phylo_100 <- paste(round(params.table_raw_simulation_calibration_phylo_100[,7], digits = 3), "(", 
                                round(params.table_raw_simulation_calibration_phylo_100[,1], digits = 3), ",",
                                round(params.table_raw_simulation_calibration_phylo_100[,3], digits = 3), ")")


params.table_combined_phylo_100 <- paste(round(params.table_raw_simulation_calibration_combined_phylo_100[,7], digits = 3), "(", 
                                         round(params.table_raw_simulation_calibration_combined_phylo_100[,1], digits = 3), ",",
                                         round(params.table_raw_simulation_calibration_combined_phylo_100[,3], digits = 3), ")")




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



params.table_master <- inputvector


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


params_espaces <- paste("(", min.v, ",",  max.v,")")


params_all <- data.frame(params.table_hiv_epi_behavior, params.table_phylo_100, 
                         params.table_combined_phylo_100, params.table_master, 
                         params_espaces)


names(params_all) <- c("params_calibration_HIV_Epi_Behavior",
                       "params_calibration_Phylo",
                       "params_calibration_combined",
                       "params_master_model", "parameter space")

# names(params.table_master) <- rownames(params.table_raw_simulation_calibration_hiv_epi_behavior)

rownames(params_all) <- params.names


write.csv(params_all, file = "/home/david/benchmark_master_model/results/CALIBRATION_SELECTED_PARAMS_ALL.csv")





# Table 2: Age mixing patterns in sexual partnership -------------


age_mix_master <- read.csv("/home/david/benchmark_master_model/results/000_master.model.age.mixing.pop.csv")


age_mix_calib_HIV_Epi_Behavior <- read.csv("/home/david/benchmark_master_model/results/111_classic.calibration.parms.age.mixing.pop.csv")

age_mix_calib_HIV_Epi_Behavior_med <- read.csv("/home/david/benchmark_master_model/results/111_2_classic.calibration.median.parms.age.mixing.pop.csv")


age_mix_calib_phylo <- read.csv("/home/david/benchmark_master_model/results/222_phylo.100.calibration.params.age.mixing.pop.csv")

age_mix_calib_phylo_med <- read.csv("/home/david/benchmark_master_model/results/222_2_phylo.100.calibration.median.params.age.mixing.pop.csv")


age_mix_calib_comb <- read.csv("/home/david/benchmark_master_model/results/333_combined.phylo.100.calibration.parms.age.mixing.pop.csv")

age_mix_calib_comb_med <- read.csv("/home/david/benchmark_master_model/results/333_2_combined.phylo.100.calibration.median.parms.age.mixing.pop.csv")



v.age_mix_master <- age_mix_master$mean
v.age_mix_calib_HIV_Epi_Behavior <- age_mix_calib_HIV_Epi_Behavior$mean
v.age_mix_calib_HIV_Epi_Behavior_med <- age_mix_calib_HIV_Epi_Behavior_med$mean
v.age_mix_calib_phylo <- age_mix_calib_phylo$mean
v.age_mix_calib_phylo_med <- age_mix_calib_phylo_med$mean
v.age_mix_calib_comb <- age_mix_calib_comb$mean
v.age_mix_calib_comb_med <- age_mix_calib_comb_med$mean


age_mixing_df <- data.frame(v.age_mix_master, 
                            v.age_mix_calib_HIV_Epi_Behavior,
                            v.age_mix_calib_HIV_Epi_Behavior_med,
                            v.age_mix_calib_phylo,
                            v.age_mix_calib_phylo_med,
                            v.age_mix_calib_comb,
                            v.age_mix_calib_comb_med)

names(age_mixing_df) <- c("master model", "calibration 1", "calibration 1*", 
                          "calibration 2", "calibration 2*",
                          "calibration 3", "calibration 3*")

rownames(age_mixing_df) <- paste(age_mix_calib_comb_med$param)


write.csv(age_mixing_df, file = "/home/david/benchmark_master_model/results/PREDIC_AGE_MIXING_ALL.csv")



# Table 3 : Onward transmissions -------------


onwardtransmi_master <- read.csv("/home/david/benchmark_master_model/results/000_master.model.onwards.transmission.csv")



onwardtransmi_calib_HIV_Epi_Behavior <- read.csv("/home/david/benchmark_master_model/results/111_classic.calibration.parms.onwards.transmission.csv")

onwardtransmi_calib_HIV_Epi_Behavior_med <- read.csv("/home/david/benchmark_master_model/results/111_2_classic.calibration.median.parms.onwards.transmission.csv")


onwardtransmi_calib_phylo <- read.csv("/home/david/benchmark_master_model/results/222_phylo.100.calibration.params.onwards.transmission.csv")

onwardtransmi_calib_phylo_med <- read.csv("/home/david/benchmark_master_model/results/222_2_phylo.100.calibration.median.params.onwards.transmission.csv")


onwardtransmi_calib_comb <- read.csv("/home/david/benchmark_master_model/results/333_combined.phylo.100.calibration.parms.onwards.transmission.csv")

onwardtransmi_calib_comb_med <- read.csv("/home/david/benchmark_master_model/results/333_2_combined.phylo.100.calibration.median.parms.onwards.transmission.csv")


v.onwardtransmi_master <- onwardtransmi_master$mean
v.onwardtransmi_calib_HIV_Epi_Behavior <- onwardtransmi_calib_HIV_Epi_Behavior$mean
v.onwardtransmi_calib_HIV_Epi_Behavior_med <- onwardtransmi_calib_HIV_Epi_Behavior_med$mean
v.onwardtransmi_calib_phylo <- onwardtransmi_calib_phylo$mean
v.onwardtransmi_calib_phylo_med <- onwardtransmi_calib_phylo_med$mean
v.onwardtransmi_calib_comb <- onwardtransmi_calib_comb$mean
v.onwardtransmi_calib_comb_med <- onwardtransmi_calib_comb_med$mean


onwardtransmi_df <- data.frame(v.onwardtransmi_master, 
                               v.onwardtransmi_calib_HIV_Epi_Behavior,
                               v.onwardtransmi_calib_HIV_Epi_Behavior_med,
                               v.onwardtransmi_calib_phylo,
                               v.onwardtransmi_calib_phylo_med,
                               v.onwardtransmi_calib_comb,
                               v.onwardtransmi_calib_comb_med)

names(onwardtransmi_df) <- c("master model", "calibration 1", "caibration 1*", 
                             "calibration 2", "calibration 2*",
                             "calibration 3", "calibration 3*")

rownames(onwardtransmi_df) <- paste(onwardtransmi_calib_comb_med$param)


write.csv(onwardtransmi_df, file = "/home/david/benchmark_master_model/results/PREDIC_ONWARD_TRANSMISSION_ALL.csv")




# Table 4: Temporal trend of incidence ---------------



# Master model

incid.master.35.36 <- read.csv("/home/david/benchmark_master_model/results/000_master.model.table.val.inc.35.36.csv")
incid.master.36.37 <- read.csv("/home/david/benchmark_master_model/results/000_master.model.table.val.inc.36.37.csv")
incid.master.37.38 <- read.csv("/home/david/benchmark_master_model/results/000_master.model.table.val.inc.37.38.csv")
incid.master.38.39 <- read.csv("/home/david/benchmark_master_model/results/000_master.model.table.val.inc.38.39.csv")
incid.master.39.40 <- read.csv("/home/david/benchmark_master_model/results/000_master.model.table.val.inc.39.40.csv")

incid.master <- rbind(incid.master.35.36, incid.master.36.37, 
                      incid.master.37.38, incid.master.38.39, 
                      incid.master.39.40)

incid.master$incid_year <- c(rep("35-36", nrow(incid.master.35.36)),
                             rep("36-37", nrow(incid.master.36.37)),
                             rep("37-38", nrow(incid.master.37.38)),
                             rep("38-39", nrow(incid.master.38.39)),
                             rep("39-40", nrow(incid.master.39.40)))

incid.master$Calibration <- rep("Benchmark", length(incid.master$incid_year))


# ( 1) With all selected parameters values

# Epi-bio-behavioural calibration

incid.classic.35.36 <- read.csv("/home/david/benchmark_master_model/results/111_classic.calibration.parms.table.val.inc.35.36.csv")
incid.classic.36.37 <- read.csv("/home/david/benchmark_master_model/results/111_classic.calibration.parms.table.val.inc.36.37.csv")
incid.classic.37.38 <- read.csv("/home/david/benchmark_master_model/results/111_classic.calibration.parms.table.val.inc.37.38.csv")
incid.classic.38.39 <- read.csv("/home/david/benchmark_master_model/results/111_classic.calibration.parms.table.val.inc.38.39.csv")
incid.classic.39.40 <- read.csv("/home/david/benchmark_master_model/results/111_classic.calibration.parms.table.val.inc.39.40.csv")

incid.classics <- rbind(incid.classic.35.36, incid.classic.36.37, 
                        incid.classic.37.38, incid.classic.38.39, 
                        incid.classic.39.40)

incid.classics$incid_year <- c(rep("35-36", nrow(incid.classic.35.36)),
                               rep("36-37", nrow(incid.classic.36.37)),
                               rep("37-38", nrow(incid.classic.37.38)),
                               rep("38-39", nrow(incid.classic.38.39)),
                               rep("39-40", nrow(incid.classic.39.40)))

incid.classics$Calibration <- rep("HIV_Epi_Behavior", length(incid.classics$incid_year))



# Phylo calibration

incid.phylo.35.36 <- read.csv("/home/david/benchmark_master_model/results/222_phylo.100.calibration.params.table.val.inc.35.36.csv")
incid.phylo.36.37 <- read.csv("/home/david/benchmark_master_model/results/222_phylo.100.calibration.params.table.val.inc.36.37.csv")
incid.phylo.37.38 <- read.csv("/home/david/benchmark_master_model/results/222_phylo.100.calibration.params.table.val.inc.37.38.csv")
incid.phylo.38.39 <- read.csv("/home/david/benchmark_master_model/results/222_phylo.100.calibration.params.table.val.inc.38.39.csv")
incid.phylo.39.40 <- read.csv("/home/david/benchmark_master_model/results/222_phylo.100.calibration.params.table.val.inc.39.40.csv")

incid.phylos <- rbind(incid.phylo.35.36, incid.phylo.36.37, 
                      incid.phylo.37.38, incid.phylo.38.39, 
                      incid.phylo.39.40)

incid.phylos$incid_year <- c(rep("35-36", nrow(incid.phylo.35.36)),
                             rep("36-37", nrow(incid.phylo.36.37)),
                             rep("37-38", nrow(incid.phylo.37.38)),
                             rep("38-39", nrow(incid.phylo.38.39)),
                             rep("39-40", nrow(incid.phylo.39.40)))

incid.phylos$Calibration <- rep("PhyloTree", length(incid.phylos$incid_year))



# Combined calibration

incid.combined.35.36 <- read.csv("/home/david/benchmark_master_model/results/333_combined.phylo.100.calibration.parms.table.val.inc.35.36.csv")
incid.combined.36.37 <- read.csv("/home/david/benchmark_master_model/results/333_combined.phylo.100.calibration.parms.table.val.inc.36.37.csv")
incid.combined.37.38 <- read.csv("/home/david/benchmark_master_model/results/333_combined.phylo.100.calibration.parms.table.val.inc.37.38.csv")
incid.combined.38.39 <- read.csv("/home/david/benchmark_master_model/results/333_combined.phylo.100.calibration.parms.table.val.inc.38.39.csv")
incid.combined.39.40 <- read.csv("/home/david/benchmark_master_model/results/333_combined.phylo.100.calibration.parms.table.val.inc.39.40.csv")


incid.combineds <- rbind(incid.combined.35.36, incid.combined.36.37, 
                         incid.combined.37.38, incid.combined.38.39, 
                         incid.combined.39.40)

incid.combineds$incid_year <- c(rep("35-36", nrow(incid.combined.35.36)),
                                rep("36-37", nrow(incid.combined.36.37)),
                                rep("37-38", nrow(incid.combined.37.38)),
                                rep("38-39", nrow(incid.combined.38.39)),
                                rep("39-40", nrow(incid.combined.39.40)))

incid.combineds$Calibration <- rep("Combined", length(incid.combineds$incid_year))


# All together

incidence_all <- rbind(incid.classics, incid.phylos, incid.combineds, incid.master)

for(i in 1:nrow(incidence_all)){
  
  if(incidence_all$incid_year[i]=="35-36"){
    incidence_all$incid_year[i] <- 2013
  }else if(incidence_all$incid_year[i]=="36-37"){
    incidence_all$incid_year[i] <- 2014
  }else if(incidence_all$incid_year[i]=="37-38"){
    incidence_all$incid_year[i] <- 2015
  }else if(incidence_all$incid_year[i]=="38-39"){
    incidence_all$incid_year[i] <- 2016
  }else if(incidence_all$incid_year[i]=="39-40"){
    incidence_all$incid_year[i] <- 2017
  }
  
}

incidence_all <- rename(incidence_all, Year = incid_year, Incidence = mean)


write.csv(incidence_all, file = "/home/david/benchmark_master_model/results/PREDIC_INCIDENCE_PARAMS.csv")




# Plot

# plot.incidence <- ggplot(incidence_all, 
#                          aes(x=Year, y=Incidence, group=Gender, color=Gender))+
#   geom_line(size=.3) +
#   geom_point()+
#   facet_wrap(age_group ~ Calibration, ncol = 4, nrow = 6)+
#   theme(legend.position="bottom")



plot.incidence <- ggplot(incidence_all, 
                         aes(x=Year, y=Incidence, group=age_group, color=age_group))+
  geom_line(size=.3) +
  geom_point()+
  facet_wrap(Gender ~ Calibration, ncol = 4, nrow = 6)+
  xlab("Time (year)") + ylab("HIV incidence")+
  theme(legend.position="bottom")


ggsave("/home/david/benchmark_master_model/results/predic_incidence_params.pdf",
       plot.incidence,
       width = 25,
       height = 15,
       units = "cm")




# (2) With median of selected parameter values

# Epi-bio-behavioural calibration

incid.classic.bis.35.36 <- read.csv("/home/david/benchmark_master_model/results/111_2_classic.calibration.median.parms.table.val.inc.35.36.csv")
incid.classic.bis.36.37 <- read.csv("/home/david/benchmark_master_model/results/111_2_classic.calibration.median.parms.table.val.inc.36.37.csv")
incid.classic.bis.37.38 <- read.csv("/home/david/benchmark_master_model/results/111_2_classic.calibration.median.parms.table.val.inc.37.38.csv")
incid.classic.bis.38.39 <- read.csv("/home/david/benchmark_master_model/results/111_2_classic.calibration.median.parms.table.val.inc.38.39.csv")
incid.classic.bis.39.40 <- read.csv("/home/david/benchmark_master_model/results/111_2_classic.calibration.median.parms.table.val.inc.39.40.csv")

incid.classics.bis <- rbind(incid.classic.bis.35.36, incid.classic.bis.36.37, 
                            incid.classic.bis.37.38, incid.classic.bis.38.39, 
                            incid.classic.bis.39.40)

incid.classics.bis$incid_year <- c(rep("35-36", nrow(incid.classic.bis.35.36)),
                                   rep("36-37", nrow(incid.classic.bis.36.37)),
                                   rep("37-38", nrow(incid.classic.bis.37.38)),
                                   rep("38-39", nrow(incid.classic.bis.38.39)),
                                   rep("39-40", nrow(incid.classic.bis.39.40)))

incid.classics.bis$Calibration <- rep("HIV_Epi_Behavior", length(incid.classics.bis$incid_year))




# Phylo calibration

incid.phylo.bis.35.36 <- read.csv("/home/david/benchmark_master_model/results/222_2_phylo.100.calibration.median.params.table.val.inc.35.36.csv")
incid.phylo.bis.36.37 <- read.csv("/home/david/benchmark_master_model/results/222_2_phylo.100.calibration.median.params.table.val.inc.36.37.csv")
incid.phylo.bis.37.38 <- read.csv("/home/david/benchmark_master_model/results/222_2_phylo.100.calibration.median.params.table.val.inc.37.38.csv")
incid.phylo.bis.38.39 <- read.csv("/home/david/benchmark_master_model/results/222_2_phylo.100.calibration.median.params.table.val.inc.38.39.csv")
incid.phylo.bis.39.40 <- read.csv("/home/david/benchmark_master_model/results/222_2_phylo.100.calibration.median.params.table.val.inc.39.40.csv")

incid.phylos.bis <- rbind(incid.phylo.bis.35.36, incid.phylo.bis.36.37, 
                          incid.phylo.bis.37.38, incid.phylo.bis.38.39, 
                          incid.phylo.bis.39.40)

incid.phylos.bis$incid_year <- c(rep("35-36", nrow(incid.phylo.bis.35.36)),
                                 rep("36-37", nrow(incid.phylo.bis.36.37)),
                                 rep("37-38", nrow(incid.phylo.bis.37.38)),
                                 rep("38-39", nrow(incid.phylo.bis.38.39)),
                                 rep("39-40", nrow(incid.phylo.bis.39.40)))

incid.phylos.bis$Calibration <- rep("PhyloTree", length(incid.phylos.bis$incid_year))




# Combined calibration

incid.combined.bis.35.36 <- read.csv("/home/david/benchmark_master_model/results/333_2_combined.phylo.100.calibration.median.parms.table.val.inc.35.36.csv")
incid.combined.bis.36.37 <- read.csv("/home/david/benchmark_master_model/results/333_2_combined.phylo.100.calibration.median.parms.table.val.inc.36.37.csv")
incid.combined.bis.37.38 <- read.csv("/home/david/benchmark_master_model/results/333_2_combined.phylo.100.calibration.median.parms.table.val.inc.37.38.csv")
incid.combined.bis.38.39 <- read.csv("/home/david/benchmark_master_model/results/333_2_combined.phylo.100.calibration.median.parms.table.val.inc.38.39.csv")
incid.combined.bis.39.40 <- read.csv("/home/david/benchmark_master_model/results/333_2_combined.phylo.100.calibration.median.parms.table.val.inc.39.40.csv")


incid.combineds.bis <- rbind(incid.combined.bis.35.36, incid.combined.bis.36.37, 
                             incid.combined.bis.37.38, incid.combined.bis.38.39, 
                             incid.combined.bis.39.40)

incid.combineds.bis$incid_year <- c(rep("35-36", nrow(incid.combined.bis.35.36)),
                                    rep("36-37", nrow(incid.combined.bis.36.37)),
                                    rep("37-38", nrow(incid.combined.bis.37.38)),
                                    rep("38-39", nrow(incid.combined.bis.38.39)),
                                    rep("39-40", nrow(incid.combined.bis.39.40)))

incid.combineds.bis$Calibration <- rep("Combined", length(incid.combineds.bis$incid_year))





# All together

incidence_all.bis <- rbind(incid.classics.bis, incid.phylos.bis, incid.combineds.bis, incid.master)


for(i in 1:nrow(incidence_all.bis)){
  
  if(incidence_all.bis$incid_year[i]=="35-36"){
    incidence_all.bis$incid_year[i] <- 2013
  }else if(incidence_all.bis$incid_year[i]=="36-37"){
    incidence_all.bis$incid_year[i] <- 2014
  }else if(incidence_all.bis$incid_year[i]=="37-38"){
    incidence_all.bis$incid_year[i] <- 2015
  }else if(incidence_all.bis$incid_year[i]=="38-39"){
    incidence_all.bis$incid_year[i] <- 2016
  }else if(incidence_all.bis$incid_year[i]=="39-40"){
    incidence_all.bis$incid_year[i] <- 2017
  }
  
}

incidence_all.bis <- rename(incidence_all.bis, Year = incid_year, Incidence = mean)

write.csv(incidence_all.bis, file = "/home/david/benchmark_master_model/results/PREDIC_INCIDENCE_MED_PARAMS.csv")


# Plot

# plot.incidence.bis <- ggplot(incidence_all.bis, 
#                              aes(x=incid_year, y=median, group=Gender, color=Gender))+
#   geom_line(size=.3) +
#   geom_point()+
#   facet_wrap(age_group ~ Calibration, ncol = 4, nrow = 6)+
#   theme(legend.position="bottom")


plot.incidence.med <- ggplot(incidence_all.bis, 
                             aes(x=Year, y=Incidence, group=age_group, color=age_group))+
  geom_line(size=.3) +
  geom_point()+
  facet_wrap(Gender ~ Calibration, ncol = 4, nrow = 6)+
  xlab("Time (year)") + ylab("HIV incidence")+
  theme(legend.position="bottom")

ggsave("/home/david/benchmark_master_model/results/predic_incidence_med_params.pdf",
       plot.incidence.med,
       width = 25,
       height = 15,
       units = "cm")




# Table 5: Summary features/statistics -----------------


clas_vals <- round(as.numeric(mean.epi.behav.features), digits = 5)

clas_params <- names(mean.epi.behav.features)

phy_vals <- round(as.numeric(mean.phylo_MCAR_100.features), digits = 5)

phy_params <- names(mean.phylo_MCAR_100.features)

ALL.vals <- c(clas_vals, phy_vals)
ALL.params <- c(clas_params, phy_params)

ALL_Features <- data.frame(ALL.params, ALL.vals)
names(ALL_Features) <- c("Parameter", "Value")

ALL.clas.f <- data.frame(clas_params, clas_vals)
names(ALL.clas.f) <- c("Parameter", "Value")



write.csv(ALL.clas.f, file = "/home/david/benchmark_master_model/results/FEATURES_EPI_BEHAVIOUR.csv")



ALL.phy.f <- data.frame(phy_params, phy_vals)
names(ALL.phy.f) <- c("Parameter", "Value")

write.csv(ALL.phy.f, file = "/home/david/benchmark_master_model/results/FEATURES_PHYLO.csv")





# Comparing predictions in master model and after calibrations ----------------------


# Using mean of absolute relative error 


# 1//N * sum(x_i - x)/x
# x: parameter value from the master model 
# x_i: N parameter values from simulations


# onwardtransmi_df
# age_mixing_df
# incidence_all
# incidence_all.bis

# Calibration 1: 
# epi.metrics_after_calibration_run_params_hiv_epi_behavior # 406*69
# epi.metrics_after_calibration_run_median_params_epi_behavior # 848*69


# Calibration 2:
# epi.metrics_after_calibration_run_params_phylo_100 # 442*69
# epi.metrics_after_calibration_run_median_params_phylo # 909*69

# Calibration 3:
# epi.metrics_after_calibration_run_params_combined # 413*69
# epi.metrics_after_calibration_run_median_params_combined # 857*69


perform_df <- NULL
perform_df[[1]] <- onwardtransmi_df
perform_df[[2]] <- age_mixing_df
perform_df[[3]] <- incidence_all
perform_df[[4]] <- incidence_all.bis

perform_df[[5]] <- epi.metrics_after_calibration_run_params_hiv_epi_behavior
perform_df[[6]] <- epi.metrics_after_calibration_run_median_params_epi_behavior
perform_df[[7]] <- epi.metrics_after_calibration_run_params_phylo_100
perform_df[[8]] <- epi.metrics_after_calibration_run_median_params_phylo
perform_df[[9]] <- epi.metrics_after_calibration_run_params_combined
perform_df[[10]] <- epi.metrics_after_calibration_run_median_params_combined

saveRDS(perform_df, file="/home/david/benchmark_master_model/results/Performance.RDS")


performance_data_list <- readRDS("/home/david/benchmark_master_model/results/Performance.RDS")


onwardtransmi_df <- performance_data_list[[1]] 

age_mixing_df <- performance_data_list[[2]] 

incidence_all <- performance_data_list[[3]] 
incidence_all.bis <- performance_data_list[[4]] 

df.calib_1 <- performance_data_list[[5]]
df.calib_1_bis<- performance_data_list[[6]]
df.calib_2 <- performance_data_list[[7]]
df.calib_2_bis <- performance_data_list[[8]]
df.calib_3 <- performance_data_list[[9]]
df.calib_3_bis <- performance_data_list[[10]]



# Onward transmission --------


# Mean of onward transmission 

mean_onwardtransm_vec <- onwardtransmi_df[1,]

# MRE after calibration 1
mean_onwardtransm_error_calib_1_1 <- sum(abs(df.calib_1$metr.transm.av - as.numeric(mean_onwardtransm_vec[1]))/abs(as.numeric(mean_onwardtransm_vec[1])))/length(df.calib_1$metr.transm.av)

mean_onwardtransm_error_calib_1_2 <- sum(abs(df.calib_1_bis$metr.transm.av - as.numeric(mean_onwardtransm_vec[1]))/abs(as.numeric(mean_onwardtransm_vec[1])))/length(df.calib_1_bis$metr.transm.av)

# MRE after calibration 2
mean_onwardtransm_error_calib_2_1 <- sum(abs(df.calib_2$metr.transm.av - as.numeric(mean_onwardtransm_vec[1]))/abs(as.numeric(mean_onwardtransm_vec[1])))/length(df.calib_2$metr.transm.av)

mean_onwardtransm_error_calib_2_2 <- sum(abs(df.calib_2_bis$metr.transm.av - as.numeric(mean_onwardtransm_vec[1]))/abs(as.numeric(mean_onwardtransm_vec[1])))/length(df.calib_2_bis$metr.transm.av)

# MRE after calibration 3
mean_onwardtransm_error_calib_3_1 <- sum(abs(df.calib_3$metr.transm.av - as.numeric(mean_onwardtransm_vec[1]))/abs(as.numeric(mean_onwardtransm_vec[1])))/length(df.calib_3$metr.transm.av)

mean_onwardtransm_error_calib_3_2 <- sum(abs(df.calib_3_bis$metr.transm.av - as.numeric(mean_onwardtransm_vec[1]))/abs(as.numeric(mean_onwardtransm_vec[1])))/length(df.calib_3_bis$metr.transm.av)



# Median of onward transmission 

med_onwardtransm_vec <- onwardtransmi_df[2,]


# MRE after calibration 1
med_onwardtransm_error_calib_1_1 <- sum(abs(df.calib_1$metr.transm.med - as.numeric(med_onwardtransm_vec[1]))/abs(as.numeric(med_onwardtransm_vec[1])))/length(df.calib_1$metr.transm.med)

med_onwardtransm_error_calib_1_2 <- sum(abs(df.calib_1_bis$metr.transm.med - as.numeric(med_onwardtransm_vec[1]))/abs(as.numeric(med_onwardtransm_vec[1])))/length(df.calib_1_bis$metr.transm.med)

# MRE after calibration 2
med_onwardtransm_error_calib_2_1 <- sum(abs(df.calib_2$metr.transm.med - as.numeric(med_onwardtransm_vec[1]))/abs(as.numeric(med_onwardtransm_vec[1])))/length(df.calib_2$metr.transm.med)


med_onwardtransm_error_calib_2_2 <- sum(abs(df.calib_2_bis$metr.transm.med - as.numeric(med_onwardtransm_vec[1]))/abs(as.numeric(med_onwardtransm_vec[1])))/length(df.calib_2_bis$metr.transm.med)


# MRE after calibration 3
med_onwardtransm_error_calib_3_1 <- sum(abs(df.calib_3$metr.transm.med - as.numeric(med_onwardtransm_vec[1]))/abs(as.numeric(med_onwardtransm_vec[1])))/length(df.calib_3$metr.transm.med)

med_onwardtransm_error_calib_3_2 <- sum(abs(df.calib_3_bis$metr.transm.med - as.numeric(med_onwardtransm_vec[1]))/abs(as.numeric(med_onwardtransm_vec[1])))/length(df.calib_3_bis$metr.transm.med)


# Standard deviation of onward transmission 

sd_onwardtransm_vec <- onwardtransmi_df[3,]

# MRE after calibration 1
sd_onwardtransm_error_calib_1_1 <- sum(abs(df.calib_1$metr.transm.sd - as.numeric(sd_onwardtransm_vec[1]))/abs(as.numeric(sd_onwardtransm_vec[1])))/length(df.calib_1$metr.transm.sd)

sd_onwardtransm_error_calib_1_2 <- sum(abs(df.calib_1_bis$metr.transm.sd - as.numeric(sd_onwardtransm_vec[1]))/abs(as.numeric(sd_onwardtransm_vec[1])))/length(df.calib_1_bis$metr.transm.sd)

# MRE after calibration 2
sd_onwardtransm_error_calib_2_1 <- sum(abs(df.calib_2$metr.transm.sd - as.numeric(sd_onwardtransm_vec[1]))/abs(as.numeric(sd_onwardtransm_vec[1])))/length(df.calib_2$metr.transm.sd)

sd_onwardtransm_error_calib_2_2 <- sum(abs(df.calib_2_bis$metr.transm.sd - as.numeric(sd_onwardtransm_vec[1]))/abs(as.numeric(sd_onwardtransm_vec[1])))/length(df.calib_2_bis$metr.transm.sd)

# MRE after calibration 3
sd_onwardtransm_error_calib_3_1 <- sum(abs(df.calib_3$metr.transm.sd - as.numeric(sd_onwardtransm_vec[1]))/abs(as.numeric(sd_onwardtransm_vec[1])))/length(df.calib_3$metr.transm.sd)

sd_onwardtransm_error_calib_3_2 <- sum(abs(df.calib_3_bis$metr.transm.sd - as.numeric(sd_onwardtransm_vec[1]))/abs(as.numeric(sd_onwardtransm_vec[1])))/length(df.calib_3_bis$metr.transm.sd)




# MRE vector for Mean of onward transmission 

error_mean_onward_transmission <- c(mean_onwardtransm_error_calib_1_1, mean_onwardtransm_error_calib_1_2, 
                                    mean_onwardtransm_error_calib_2_1, mean_onwardtransm_error_calib_2_2, 
                                    mean_onwardtransm_error_calib_3_1, mean_onwardtransm_error_calib_3_2) 



amazina <- c("MRE_1_1", "MRE_1_2", "MRE_2_1", "MRE_2_2", "MRE_3_1", "MRE_3_2")


names(error_mean_onward_transmission) <- amazina


error_mean_onward_transmission <- round(error_mean_onward_transmission, digits = 5)



# MRE vector for Median of onward transmission 

error_med_onward_transmission <- c(med_onwardtransm_error_calib_1_1, med_onwardtransm_error_calib_1_2, 
                                   med_onwardtransm_error_calib_2_1, med_onwardtransm_error_calib_2_2, 
                                   med_onwardtransm_error_calib_3_1, med_onwardtransm_error_calib_3_2) 


names(error_med_onward_transmission) <- amazina


error_med_onward_transmission <- round(error_med_onward_transmission, digits = 5)



# MRE vector for Standard deviation of onward transmission 

error_sd_onward_transmission <- c(sd_onwardtransm_error_calib_1_1, sd_onwardtransm_error_calib_1_2, 
                                  sd_onwardtransm_error_calib_2_1, sd_onwardtransm_error_calib_2_2, 
                                  sd_onwardtransm_error_calib_3_1, sd_onwardtransm_error_calib_3_2) 

error_sd_onward_transmission <- as.numeric(error_sd_onward_transmission)


names(error_sd_onward_transmission) <- amazina


error_sd_onward_transmission <- round(error_sd_onward_transmission, digits = 5)


# MRE table

onward_transmission_error_table <- rbind(error_mean_onward_transmission, error_med_onward_transmission, error_sd_onward_transmission)


write.csv(onward_transmission_error_table, file = "/home/david/benchmark_master_model/results/ERROR_PREDIC_ONWARD_TRANSMISSIONS.csv")




# Age mixing patterns in partnership --------



# AAD 

AAD_vec <- age_mixing_df[1,]

# MRE after calibration 1
AAD_age.mix_error_calib_1_1 <- sum(abs(df.calib_1$metr.AAD.male - as.numeric(AAD_vec[1]))/abs(as.numeric(AAD_vec[1])))/length(df.calib_1$metr.AAD.male)

AAD_age.mix_error_calib_1_2 <- sum(abs(df.calib_1_bis$metr.AAD.male - as.numeric(AAD_vec[1]))/abs(as.numeric(AAD_vec[1])))/length(df.calib_1_bis$metr.AAD.male)

# MRE after calibration 2
AAD_age.mix_error_calib_2_1 <- sum(abs(df.calib_2$metr.AAD.male - as.numeric(AAD_vec[1]))/abs(as.numeric(AAD_vec[1])))/length(df.calib_2$metr.AAD.male)

AAD_age.mix_error_calib_2_2 <- sum(abs(df.calib_2_bis$metr.AAD.male - as.numeric(AAD_vec[1]))/abs(as.numeric(AAD_vec[1])))/length(df.calib_2_bis$metr.AAD.male)

# MRE after calibration 3
AAD_age.mix_error_calib_3_1 <- sum(abs(df.calib_3$metr.AAD.male - as.numeric(AAD_vec[1]))/abs(as.numeric(AAD_vec[1])))/length(df.calib_3$metr.AAD.male)

AAD_age.mix_error_calib_3_2 <- sum(abs(df.calib_3_bis$metr.AAD.male - as.numeric(AAD_vec[1]))/abs(as.numeric(AAD_vec[1])))/length(df.calib_3_bis$metr.AAD.male)


# SDAD 

SDAD_vec <- age_mixing_df[2,]

# MRE after calibration 1
SDAD_age.mix_error_calib_1_1 <- sum(abs(df.calib_1$metr.SDAD.male - as.numeric(SDAD_vec[1]))/abs(as.numeric(SDAD_vec[1])))/length(df.calib_1$metr.SDAD.male)

SDAD_age.mix_error_calib_1_2 <- sum(abs(df.calib_1_bis$metr.SDAD.male - as.numeric(SDAD_vec[1]))/abs(as.numeric(SDAD_vec[1])))/length(df.calib_1_bis$metr.SDAD.male)

# MRE after calibration 2
SDAD_age.mix_error_calib_2_1 <- sum(abs(df.calib_2$metr.SDAD.male - as.numeric(SDAD_vec[1]))/abs(as.numeric(SDAD_vec[1])))/length(df.calib_2$metr.SDAD.male)

SDAD_age.mix_error_calib_2_2 <- sum(abs(df.calib_2_bis$metr.SDAD.male - as.numeric(SDAD_vec[1]))/abs(as.numeric(SDAD_vec[1])))/length(df.calib_2_bis$metr.SDAD.male)

# MRE after calibration 3
SDAD_age.mix_error_calib_3_1 <- sum(abs(df.calib_3$metr.SDAD.male - as.numeric(SDAD_vec[1]))/abs(as.numeric(SDAD_vec[1])))/length(df.calib_3$metr.SDAD.male)

SDAD_age.mix_error_calib_3_2 <- sum(abs(df.calib_3_bis$metr.SDAD.male - as.numeric(SDAD_vec[1]))/abs(as.numeric(SDAD_vec[1])))/length(df.calib_3_bis$metr.SDAD.male)



# BSD 


BSD_vec <- age_mixing_df[3,]

# MRE after calibration 1
BSD_age.mix_error_calib_1_1 <- sum(abs(df.calib_1$metr.BSD.male - as.numeric(BSD_vec[1]))/abs(as.numeric(BSD_vec[1])))/length(df.calib_1$metr.BSD.male)

BSD_age.mix_error_calib_1_2 <- sum(abs(df.calib_1_bis$metr.BSD.male - as.numeric(BSD_vec[1]))/abs(as.numeric(BSD_vec[1])))/length(df.calib_1_bis$metr.BSD.male)

# MRE after calibration 2
BSD_age.mix_error_calib_2_1 <- sum(abs(df.calib_2$metr.BSD.male - as.numeric(BSD_vec[1]))/abs(as.numeric(BSD_vec[1])))/length(df.calib_2$metr.BSD.male)

BSD_age.mix_error_calib_2_2 <- sum(abs(df.calib_2_bis$metr.BSD.male - as.numeric(BSD_vec[1]))/abs(as.numeric(BSD_vec[1])))/length(df.calib_2_bis$metr.BSD.male)

# MRE after calibration 3
BSD_age.mix_error_calib_3_1 <- sum(abs(df.calib_3$metr.BSD.male - as.numeric(BSD_vec[1]))/abs(as.numeric(BSD_vec[1])))/length(df.calib_3$metr.BSD.male)

BSD_age.mix_error_calib_3_2 <- sum(abs(df.calib_3_bis$metr.BSD.male - as.numeric(BSD_vec[1]))/abs(as.numeric(BSD_vec[1])))/length(df.calib_3_bis$metr.BSD.male)


# WSD 

WSD_vec <- age_mixing_df[4,]

# MRE after calibration 1
WSD_age.mix_error_calib_1_1 <- sum(abs(df.calib_1$metr.WSD.male - as.numeric(WSD_vec[1]))/abs(as.numeric(WSD_vec[1])))/length(df.calib_1$metr.WSD.male)

WSD_age.mix_error_calib_1_2 <- sum(abs(df.calib_1_bis$metr.WSD.male - as.numeric(WSD_vec[1]))/abs(as.numeric(WSD_vec[1])))/length(df.calib_1_bis$metr.WSD.male)

# MRE after calibration 2
WSD_age.mix_error_calib_2_1 <- sum(abs(df.calib_2$metr.WSD.male - as.numeric(WSD_vec[1]))/abs(as.numeric(WSD_vec[1])))/length(df.calib_2$metr.WSD.male)

WSD_age.mix_error_calib_2_2 <- sum(abs(df.calib_2_bis$metr.WSD.male - as.numeric(WSD_vec[1]))/abs(as.numeric(WSD_vec[1])))/length(df.calib_2_bis$metr.WSD.male)

# MRE after calibration 3
WSD_age.mix_error_calib_3_1 <- sum(abs(df.calib_3$metr.WSD.male - as.numeric(WSD_vec[1]))/abs(as.numeric(WSD_vec[1])))/length(df.calib_3$metr.WSD.male)

WSD_age.mix_error_calib_3_2 <- sum(abs(df.calib_3_bis$metr.WSD.male - as.numeric(WSD_vec[1]))/abs(as.numeric(WSD_vec[1])))/length(df.calib_3_bis$metr.WSD.male)



# Slope 

slope_vec <- age_mixing_df[5,]

# MRE after calibration 1
slope_age.mix_error_calib_1_1 <- sum(abs(df.calib_1$metr.slope.male - as.numeric(slope_vec[1]))/abs(as.numeric(slope_vec[1])))/length(df.calib_1$metr.slope.male)

slope_age.mix_error_calib_1_2 <- sum(abs(df.calib_1_bis$metr.slope.male - as.numeric(slope_vec[1]))/abs(as.numeric(slope_vec[1])))/length(df.calib_1_bis$metr.slope.male)

# MRE after calibration 2
slope_age.mix_error_calib_2_1 <- sum(abs(df.calib_2$metr.slope.male - as.numeric(slope_vec[1]))/abs(as.numeric(slope_vec[1])))/length(df.calib_2$metr.slope.male)

slope_age.mix_error_calib_2_2 <- sum(abs(df.calib_2_bis$metr.slope.male - as.numeric(slope_vec[1]))/abs(as.numeric(slope_vec[1])))/length(df.calib_2_bis$metr.slope.male)

# MRE after calibration 3
slope_age.mix_error_calib_3_1 <- sum(abs(df.calib_3$metr.slope.male - as.numeric(slope_vec[1]))/abs(as.numeric(slope_vec[1])))/length(df.calib_3$metr.slope.male)

slope_age.mix_error_calib_3_2 <- sum(abs(df.calib_3_bis$metr.slope.male - as.numeric(slope_vec[1]))/abs(as.numeric(slope_vec[1])))/length(df.calib_3_bis$metr.slope.male)




# Intercept 

interc_vec <- age_mixing_df[6,]

# MRE after calibration 1
interc_age.mix_error_calib_1_1 <- sum(abs(df.calib_1$metr.intercept.male - as.numeric(interc_vec[1]))/abs(as.numeric(interc_vec[1])))/length(df.calib_1$metr.intercept.male)

interc_age.mix_error_calib_1_2 <- sum(abs(df.calib_1_bis$metr.intercept.male - as.numeric(interc_vec[1]))/abs(as.numeric(interc_vec[1])))/length(df.calib_1_bis$metr.intercept.male)

# MRE after calibration 2
interc_age.mix_error_calib_2_1 <- sum(abs(df.calib_2$metr.intercept.male - as.numeric(interc_vec[1]))/abs(as.numeric(interc_vec[1])))/length(df.calib_2$metr.intercept.male)

interc_age.mix_error_calib_2_2 <- sum(abs(df.calib_2_bis$metr.intercept.male - as.numeric(interc_vec[1]))/abs(as.numeric(interc_vec[1])))/length(df.calib_2_bis$metr.intercept.male)

# MRE after calibration 3
interc_age.mix_error_calib_3_1 <- sum(abs(df.calib_3$metr.intercept.male - as.numeric(interc_vec[1]))/abs(as.numeric(interc_vec[1])))/length(df.calib_3$metr.intercept.male)

interc_age.mix_error_calib_3_2 <- sum(abs(df.calib_3_bis$metr.intercept.male - as.numeric(interc_vec[1]))/abs(as.numeric(interc_vec[1])))/length(df.calib_3_bis$metr.intercept.male)




# MRE vector for AAD of ge mixing in partnership 

error_AAD_age.mix <- c(AAD_age.mix_error_calib_1_1, AAD_age.mix_error_calib_1_2,
                       AAD_age.mix_error_calib_2_1, AAD_age.mix_error_calib_2_2,
                       AAD_age.mix_error_calib_3_1, AAD_age.mix_error_calib_3_2)

error_AAD_age.mix <- round(error_AAD_age.mix, digits = 5)

names(error_AAD_age.mix) <- amazina



# MRE vector for SDAD of ge mixing in partnership

error_SDAD_age.mix <- c(SDAD_age.mix_error_calib_1_1, SDAD_age.mix_error_calib_1_2,
                        SDAD_age.mix_error_calib_2_1, SDAD_age.mix_error_calib_2_2,
                        SDAD_age.mix_error_calib_3_1, SDAD_age.mix_error_calib_3_2)


error_SDAD_age.mix <- round(error_SDAD_age.mix, digits = 5)

names(error_SDAD_age.mix) <- amazina



# MRE vector for BSD of ge mixing in partnership

error_BSD_age.mix <- c(BSD_age.mix_error_calib_1_1, BSD_age.mix_error_calib_1_2,
                       BSD_age.mix_error_calib_2_1, BSD_age.mix_error_calib_2_2,
                       BSD_age.mix_error_calib_3_1, BSD_age.mix_error_calib_3_2)


error_BSD_age.mix <- round(error_BSD_age.mix, digits = 5)

names(error_BSD_age.mix) <- amazina


# MRE vector for WSD of ge mixing in partnership

error_WSD_age.mix <- c(WSD_age.mix_error_calib_1_1, WSD_age.mix_error_calib_1_2,
                       WSD_age.mix_error_calib_2_1, WSD_age.mix_error_calib_2_2,
                       WSD_age.mix_error_calib_3_1, WSD_age.mix_error_calib_3_2)


error_WSD_age.mix <- round(error_WSD_age.mix, digits = 5)

names(error_WSD_age.mix) <- amazina


# MRE vector for Slope of ge mixing in partnership

error_slope_age.mix <- c(slope_age.mix_error_calib_1_1, slope_age.mix_error_calib_1_2,
                         slope_age.mix_error_calib_2_1, slope_age.mix_error_calib_2_2,
                         slope_age.mix_error_calib_3_1, slope_age.mix_error_calib_3_2)


error_slope_age.mix <- round(error_slope_age.mix, digits = 5)

names(error_slope_age.mix) <- amazina



# MRE vector for Intercept of ge mixing in partnership

error_interc_age.mix <- c(interc_age.mix_error_calib_1_1, interc_age.mix_error_calib_1_2,
                          interc_age.mix_error_calib_2_1, interc_age.mix_error_calib_2_2,
                          interc_age.mix_error_calib_3_1, interc_age.mix_error_calib_3_2)


error_interc_age.mix <- round(error_interc_age.mix, digits = 5)

names(error_interc_age.mix) <- amazina



# MRE table for age mixing

age_mixing_error_table <- rbind(error_AAD_age.mix, error_SDAD_age.mix,
                                error_BSD_age.mix, error_WSD_age.mix,
                                error_slope_age.mix, error_interc_age.mix)


write.csv(age_mixing_error_table, file = "/home/david/benchmark_master_model/results/ERROR_PREDIC_AGE_MIXING.csv")





# Temporal trend of incidence ---------


incidence_all_0 <- dplyr::filter(incidence_all, incidence_all$Calibration=="Benchmark")

incidence_all_0_m <- dplyr::filter(incidence_all_0, incidence_all_0$Gender=="men") # temporal incidence for men
incidence_all_0_w <- dplyr::filter(incidence_all_0, incidence_all_0$Gender=="women") # temporal incidece for women


# MRE after calibration 1

incidence_tab_calib_1_1 <- df.calib_1 %>%
  select(contains("incid.")) 

incidence_tab_calib_1_1_m <- incidence_tab_calib_1_1 %>%
  select(contains(".m."))

incidence_tab_calib_1_1_w <- incidence_tab_calib_1_1 %>%
  select(contains(".w."))


incidence_tab_calib_1_2 <- df.calib_1_bis %>%
  select(contains("incid.")) 

incidence_tab_calib_1_2_m <- incidence_tab_calib_1_2 %>%
  select(contains(".m."))

incidence_tab_calib_1_2_w <- incidence_tab_calib_1_2 %>%
  select(contains(".w."))


# For men

incidence_error_m_1_1 <- vector()

for(i in 1:ncol(incidence_tab_calib_1_1_m)){
  
  # MRE between mean value from master model of 
  # incidence at a given year for a given age group: incidence_all_0_m$Incidence[i] 
  # and vector of incidence values at athe same given year for a given age group: incidence_tab_calib_1_m[,i] 
  incidence_error_i <- 
    sum(abs(incidence_tab_calib_1_1_m[,i] - incidence_all_0_m$Incidence[i])/abs(incidence_all_0_m$Incidence[i]))/length(incidence_tab_calib_1_1_m[,i])
  incidence_error_m_1_1 <- c(incidence_error_m_1_1, incidence_error_i)
  
}


incidence_error_m_1_2 <- vector()

for(i in 1:ncol(incidence_tab_calib_1_2_m)){
  
  # MRE between mean value from master model of 
  # incidence at a given year for a given age group: incidence_all_0_m$Incidence[i] 
  # and vector of incidence values at athe same given year for a given age group: incidence_tab_calib_1_m[,i] 
  incidence_error_i <- sum(abs(incidence_tab_calib_1_2_m[,i] - incidence_all_0_m$Incidence[i])/abs(incidence_all_0_m$Incidence[i]))/length(incidence_tab_calib_1_2_m[,i])
  
  incidence_error_m_1_2 <- c(incidence_error_m_1_2, incidence_error_i)
  
}



# For women

incidence_error_w_1_1 <- vector()

for(i in 1:ncol(incidence_tab_calib_1_1_w)){
  
  # MRE between mean value from master model of 
  # incidence at a given year for a given age group: incidence_all_0_w$Incidence[i] 
  # and vector of incidence values at athe same given year for a given age group: incidence_tab_calib_1_w[,i] 
  incidence_error_i <- sum(abs(incidence_tab_calib_1_1_w[,i] - incidence_all_0_w$Incidence[i])/abs(incidence_all_0_w$Incidence[i]))/length(incidence_tab_calib_1_1_w[,i])
  
  incidence_error_w_1_1 <- c(incidence_error_w_1_1, incidence_error_i)
  
}


incidence_error_w_1_2 <- vector()

for(i in 1:ncol(incidence_tab_calib_1_2_w)){
  
  # MRE between mean value from master model of 
  # incidence at a given year for a given age group: incidence_all_0_w$Incidence[i] 
  # and vector of incidence values at athe same given year for a given age group: incidence_tab_calib_1_w[,i] 
  incidence_error_i <- sum(abs(incidence_tab_calib_1_2_w[,i] - incidence_all_0_w$Incidence[i])/abs(incidence_all_0_w$Incidence[i]))/length(incidence_tab_calib_1_2_w[,i])
  
  incidence_error_w_1_2 <- c(incidence_error_w_1_2, incidence_error_i)
  
}


# MRE after calibration 2

incidence_tab_calib_2_1 <- df.calib_2 %>%
  select(contains("incid.")) 

incidence_tab_calib_2_1_m <- incidence_tab_calib_2_1 %>%
  select(contains(".m."))

incidence_tab_calib_2_1_w <- incidence_tab_calib_2_1 %>%
  select(contains(".w."))


incidence_tab_calib_2_2 <- df.calib_2_bis %>%
  select(contains("incid.")) 

incidence_tab_calib_2_2_m <- incidence_tab_calib_2_2 %>%
  select(contains(".m."))

incidence_tab_calib_2_2_w <- incidence_tab_calib_2_2 %>%
  select(contains(".w."))


# For men

incidence_error_m_2_1 <- vector()

for(i in 1:ncol(incidence_tab_calib_2_1_m)){
  
  # MRE between mean value from master model of 
  # incidence at a given year for a given age group: incidence_all_0_m$Incidence[i] 
  # and vector of incidence values at athe same given year for a given age group: incidence_tab_calib_2_m[,i] 
  incidence_error_i <- sum(abs(incidence_tab_calib_2_1_m[,i] - incidence_all_0_m$Incidence[i])/abs(incidence_all_0_m$Incidence[i]))/length(incidence_tab_calib_2_1_m[,i])
  
  incidence_error_m_2_1 <- c(incidence_error_m_2_1, incidence_error_i)
  
}


incidence_error_m_2_2 <- vector()

for(i in 1:ncol(incidence_tab_calib_2_2_m)){
  
  # MRE between mean value from master model of 
  # incidence at a given year for a given age group: incidence_all_0_m$Incidence[i] 
  # and vector of incidence values at athe same given year for a given age group: incidence_tab_calib_2_m[,i] 
  incidence_error_i <- sum(abs(incidence_tab_calib_2_2_m[,i] - incidence_all_0_m$Incidence[i])/abs(incidence_all_0_m$Incidence[i]))/length(incidence_tab_calib_2_2_m[,i])
  
  incidence_error_m_2_2 <- c(incidence_error_m_2_2, incidence_error_i)
  
}



# For women

incidence_error_w_2_1 <- vector()

for(i in 1:ncol(incidence_tab_calib_2_1_w)){
  
  # MRE between mean value from master model of 
  # incidence at a given year for a given age group: incidence_all_0_w$Incidence[i] 
  # and vector of incidence values at athe same given year for a given age group: incidence_tab_calib_2_w[,i] 
  incidence_error_i <- sum(abs(incidence_tab_calib_2_1_w[,i] - incidence_all_0_w$Incidence[i])/abs(incidence_all_0_w$Incidence[i]))/length(incidence_tab_calib_2_1_w[,i])
  
  incidence_error_w_2_1 <- c(incidence_error_w_2_1, incidence_error_i)
  
}


incidence_error_w_2_2 <- vector()

for(i in 1:ncol(incidence_tab_calib_2_2_w)){
  
  # MRE between mean value from master model of 
  # incidence at a given year for a given age group: incidence_all_0_w$Incidence[i] 
  # and vector of incidence values at athe same given year for a given age group: incidence_tab_calib_2_w[,i] 
  incidence_error_i <- sum(abs(incidence_tab_calib_2_2_w[,i] - incidence_all_0_w$Incidence[i])/abs(incidence_all_0_w$Incidence[i]))/length(incidence_tab_calib_2_2_w[,i])
  
  incidence_error_w_2_2 <- c(incidence_error_w_2_2, incidence_error_i)
  
}



# MRE after calibration 3


incidence_tab_calib_3_1 <- df.calib_3 %>%
  select(contains("incid.")) 

incidence_tab_calib_3_1_m <- incidence_tab_calib_3_1 %>%
  select(contains(".m."))

incidence_tab_calib_3_1_w <- incidence_tab_calib_3_1 %>%
  select(contains(".w."))


incidence_tab_calib_3_2 <- df.calib_3_bis %>%
  select(contains("incid.")) 

incidence_tab_calib_3_2_m <- incidence_tab_calib_3_2 %>%
  select(contains(".m."))

incidence_tab_calib_3_2_w <- incidence_tab_calib_3_2 %>%
  select(contains(".w."))


# For men

incidence_error_m_3_1 <- vector()

for(i in 1:ncol(incidence_tab_calib_3_1_m)){
  
  # MRE between mean value from master model of 
  # incidence at a given year for a given age group: incidence_all_0_m$Incidence[i] 
  # and vector of incidence values at athe same given year for a given age group: incidence_tab_calib_3_m[,i] 
  incidence_error_i <- sum(abs(incidence_tab_calib_3_1_m[,i] - incidence_all_0_m$Incidence[i])/abs(incidence_all_0_m$Incidence[i]))/length(incidence_tab_calib_3_1_m[,i])
  
  incidence_error_m_3_1 <- c(incidence_error_m_3_1, incidence_error_i)
  
}


incidence_error_m_3_2 <- vector()

for(i in 1:ncol(incidence_tab_calib_3_2_m)){
  
  # MRE between mean value from master model of 
  # incidence at a given year for a given age group: incidence_all_0_m$Incidence[i] 
  # and vector of incidence values at athe same given year for a given age group: incidence_tab_calib_3_m[,i] 
  incidence_error_i <- sum(abs(incidence_tab_calib_3_2_m[,i] - incidence_all_0_m$Incidence[i])/abs(incidence_all_0_m$Incidence[i]))/length(incidence_tab_calib_3_2_m[,i])
  
  incidence_error_m_3_2 <- c(incidence_error_m_3_2, incidence_error_i)
  
}



# For women

incidence_error_w_3_1 <- vector()

for(i in 1:ncol(incidence_tab_calib_3_1_w)){
  
  # MRE between mean value from master model of 
  # incidence at a given year for a given age group: incidence_all_0_w$Incidence[i] 
  # and vector of incidence values at athe same given year for a given age group: incidence_tab_calib_3_w[,i] 
  incidence_error_i <- sum(abs(incidence_tab_calib_3_1_w[,i] - incidence_all_0_w$Incidence[i])/abs(incidence_all_0_w$Incidence[i]))/length(incidence_tab_calib_3_1_w[,i])
  
  incidence_error_w_3_1 <- c(incidence_error_w_3_1, incidence_error_i)
  
}


incidence_error_w_3_2 <- vector()

for(i in 1:ncol(incidence_tab_calib_3_2_w)){
  
  # MRE between mean value from master model of 
  # incidence at a given year for a given age group: incidence_all_0_w$Incidence[i] 
  # and vector of incidence values at athe same given year for a given age group: incidence_tab_calib_3_w[,i] 
  incidence_error_i <- sum(abs(incidence_tab_calib_3_2_w[,i] - incidence_all_0_w$Incidence[i])/abs(incidence_all_0_w$Incidence[i]))/length(incidence_tab_calib_3_2_w[,i])
  
  incidence_error_w_3_2 <- c(incidence_error_w_3_2, incidence_error_i)
  
}



men_incid_error <- data.frame(incidence_error_m_1_1, incidence_error_m_1_2, 
                              incidence_error_m_2_1, incidence_error_m_2_2,
                              incidence_error_m_3_1, incidence_error_m_3_2)

men_incid_error <- round(men_incid_error, digits = 5)

men_incid_error.df <- cbind(men_incid_error, incidence_all_0_m$age_group, incidence_all_0_m$Year, incidence_all_0_m$Gender)

names(men_incid_error.df) <- c("MRE_1_1", "MRE_1_2", "MRE_2_1", "MRE_2_2", "MRE_3_1", "MRE_3_2", "Age_group", "Year", "Gender")



women_incid_error <- data.frame(incidence_error_w_1_1, incidence_error_w_1_2, 
                                incidence_error_w_2_1, incidence_error_w_2_2,
                                incidence_error_w_3_1, incidence_error_w_3_2)

women_incid_error <- round(women_incid_error, digits = 5)

women_incid_error.df <- cbind(women_incid_error, incidence_all_0_m$age_group, incidence_all_0_m$Year, incidence_all_0_w$Gender)

names(women_incid_error.df) <- c("MRE_1_1", "MRE_1_2", "MRE_2_1", "MRE_2_2", "MRE_3_1", "MRE_3_2", "Age_group", "Year", "Gender")


incid_error <- rbind(men_incid_error.df, women_incid_error.df)


# MRE age mixing Table

write.csv(incid_error, file = "/home/david/benchmark_master_model/results/ERROR_PREDIC_INCIDENCE.csv")


write.csv(men_incid_error.df, file = "/home/david/benchmark_master_model/results/ERROR_PREDIC_INCIDENCE_MALES.csv")


write.csv(women_incid_error.df, file = "/home/david/benchmark_master_model/results/ERROR_PREDIC_INCIDENCE_FEMALES.csv")






# Comparing selected parameters' values in calibration scenarios with benchmark values ------------------


# Using absolute relative error

abc_res_calibration_hiv_epi_behavior <- readRDS("/home/david/benchmark_master_model/calibration_epi_behav/abc_res_calibration_hiv_epi_behavior.RDS")
abc_params_hiv_epi_behavior <- as.data.frame(abc_res_calibration_hiv_epi_behavior$unadj.values)
colnames(abc_params_hiv_epi_behavior) <- params.names

abc_res_calibration_phylo_100 <- readRDS("/home/david/benchmark_master_model/calibration_phylo/abc_res_calibration_phylo.RDS")
abc_params_phylo_100 <- as.data.frame(abc_res_calibration_phylo_100$unadj.values)
colnames(abc_params_phylo_100) <- params.names


abc_res_calibration_combined_phylo_100 <- readRDS("/home/david/benchmark_master_model/calibration_combined/abc_res_calibration_combined_phylo.RDS")
abc_params_combined <- as.data.frame(abc_res_calibration_combined_phylo_100$unadj.values)
colnames(abc_params_combined) <- params.names


# Compute relative error between parameter inputvector for benchmark and accepted parameter values in calibration



# Epi-behaviour


e.vec_epi_beha <- list()

for (i in 1:ncol(abc_params_hiv_epi_behavior)) {
  
  e.i <- (1/length(abc_params_hiv_epi_behavior[,i] )) * abs((abc_params_hiv_epi_behavior[,i] - inputvector[i])/inputvector[i])
  
  e.vec_epi_beha[[i]] <- e.i
  
}

e.vec_epi_beha <- as.data.frame(e.vec_epi_beha)

names(e.vec_epi_beha) <- params.names

# boxplot(e.vec_epi_beha)




# Phylo



e.vec_phylo <- list()

for (i in 1:ncol(abc_params_phylo_100)) {
  
  e.i <- (1/length(abc_params_phylo_100[,i] )) * abs((abc_params_phylo_100[,i] - inputvector[i])/inputvector[i])
  
  e.vec_phylo[[i]] <- e.i
  
}


e.vec_phylo <- as.data.frame(e.vec_phylo)

names(e.vec_phylo) <- params.names

# boxplot(e.vec_phylo)



# Combined


e.vec_combined <- list()

for (i in 1:ncol(abc_params_combined)) {
  
  e.i <- (1/length(abc_params_combined[,i] )) * abs((abc_params_combined[,i] - inputvector[i])/inputvector[i])
  
  e.vec_combined[[i]] <- e.i
  
}


e.vec_combined <- as.data.frame(e.vec_combined)

names(e.vec_combined) <- params.names

# boxplot(e.vec_combined)


e_a <- e.vec_epi_beha
e_b <- e.vec_phylo
e_c <- e.vec_combined


e_a$calibration <- rep("epi_behav", nrow(e_a))
e_b$calibration <- rep("phylo", nrow(e_b))
e_c$calibration <- rep("combined", nrow(e_c))


e_fun <- function(x=e_a) {
  
  calib_i <- unique(x$calibration)
  values_cols <- vector()
  params_cols <- vector()
  calibration <- vector()
  for(i in 1:ncol(x)-1){
    x_col <- x[,i]
    x_parm <- rep(names(x)[i], length(x_col))
    x_cal <- rep(calib_i, length(x_col))
    
    values_cols <- c(values_cols, x_col)
    params_cols <- c(params_cols, x_parm)
    calibration <- c(calibration, x_cal)
    
  }
  
  
  error_gplt <- cbind(values_cols, params_cols, calibration)
  
  error_gplt <- as.data.frame(error_gplt)
  
  return(error_gplt)
  
}


A <- as.data.frame(e_fun(x=e_a))
B <- as.data.frame(e_fun(x=e_b))
C <- as.data.frame(e_fun(x=e_c))


e_data <- rbind(A, B, C)
e_data$values_cols <- as.numeric(e_data$values_cols)
e_data$params_cols <- as.character(e_data$params_cols)
e_data$calibration <- as.character(e_data$calibration)

# grouped boxplot
p.params_error <- ggplot(e_data, aes(x=params_cols, y=values_cols,  fill=calibration)) + 
  geom_boxplot() + 
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5))+
  ylab("Error")+
  xlab("Parameters")


ggsave(filename = "000_error_selected_params.pdf",
       plot = p.params_error,
       path = "/home/david/benchmark_master_model/results",
       width = 25, height = 15, units = "cm")





# Comparing summary statistics associated to selected parameters in calibration scenarios ------------------

# Using absolute relative error 

abc_ss_hiv_epi_behavior <- as.data.frame(abc_res_calibration_hiv_epi_behavior$ss)
colnames(abc_ss_hiv_epi_behavior) <- names_ss_epi_behav

abc_ss_phylo_100 <- as.data.frame(abc_res_calibration_phylo_100$ss)
colnames(abc_ss_phylo_100) <- names_ss_phylo


abc_ss_combined <- as.data.frame(abc_res_calibration_combined_phylo_100$ss)
colnames(abc_ss_combined) <- c(names_ss_epi_behav, names_ss_phylo)
# 
# 
# 
# 
# ss_fun <- function(x=ss) {
#   
#   values_cols <- vector()
#   ss_cols <- vector()
#   
#   for(i in 1:ncol(x)){
#     x_col <- x[,i]
#     x_ss <- rep(names(x)[i], length(x_col))
#     
#     values_cols <- c(values_cols, x_col)
#     ss_cols <- c(ss_cols, x_ss)
#     
#   }
#   
#   
#   ss_gplt <- cbind(ss_cols, values_cols)
#   
#   ss_gplt <- as.data.frame(ss_gplt)
#   
#   return(ss_gplt)
#   
# }
# 


# Formatting for visualization



# Epi-Behav



ss.vec_epi_behav <- list()
for (i in 1:ncol(abc_ss_hiv_epi_behavior)) {
  
  e.i <- (1/length(abc_ss_hiv_epi_behavior[,i] )) * abs((abc_ss_hiv_epi_behavior[,i] - mean.epi.behav.features[i])/mean.epi.behav.features[i])
  
  ss.vec_epi_behav[[i]] <- e.i
  
}

ss.vec_epi_behav <- as.data.frame(ss.vec_epi_behav)

names(ss.vec_epi_behav) <- names_ss_epi_behav



dat.ss.vec_epi_behav <- melt(ss.vec_epi_behav, measure.vars=c(names(ss.vec_epi_behav)))


p.dat.ss.vec_epi_behav <- ggplot(dat.ss.vec_epi_behav) +
  geom_boxplot(aes(x=variable, y=value))+
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5))+
  ylab("Error")+
  xlab("Summary statistics")


ggsave(filename = "000_error_ss_epi_behav.pdf",
       plot = p.dat.ss.vec_epi_behav,
       path = "/home/david/benchmark_master_model/results",
       width = 16, height = 12, units = "cm")



# Phylo



ss.vec_phylo <- list()

for (i in 1:ncol(abc_ss_phylo_100)) {
  
  e.i <- (1/length(abc_ss_phylo_100[,i] )) * abs((abc_ss_phylo_100[,i] - mean.phylo_MCAR_100.features[i])/mean.phylo_MCAR_100.features[i])
  
  ss.vec_phylo[[i]] <- e.i
  
}

ss.vec_phylo <- as.data.frame(ss.vec_phylo)

names(ss.vec_phylo) <- names_ss_phylo



dat.ss.vec_phylo <- melt(ss.vec_phylo, measure.vars=c(names(ss.vec_phylo)))


p.dat.ss.vec_phylo <- ggplot(dat.ss.vec_phylo) +
  geom_boxplot(aes(x=variable, y=value))+
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5))+
  ylab("Error")+
  xlab("Summary statistics")


ggsave(filename = "000_error_ss_phylo.pdf",
       plot = p.dat.ss.vec_phylo,
       path = "/home/david/benchmark_master_model/results",
       width = 16, height = 12, units = "cm")



# Combined


mean_combined <- c(mean.epi.behav.features, mean.phylo_MCAR_100.features)


ss.vec_combined <- list()

for (i in 1:ncol(abc_ss_combined)) {
  
  e.i <- (1/length(abc_ss_combined[,i] )) * abs((abc_ss_combined[,i] - mean_combined[i])/mean_combined[i])
  
  ss.vec_combined[[i]] <- e.i
  
}

ss.vec_combined <- as.data.frame(ss.vec_combined)

names(ss.vec_combined) <- c(names_ss_epi_behav, names_ss_phylo)


dat.ss.vec_combined <- melt(ss.vec_combined, measure.vars=c(names(ss.vec_combined)))


p.dat.ss.vec_combined <- ggplot(dat.ss.vec_combined) +
  geom_boxplot(aes(x=variable, y=value))+
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5))+
  ylab("Error")+
  xlab("Summary statistics")


ggsave(filename = "000_error_ss_combined.pdf",
       plot = p.dat.ss.vec_combined,
       path = "/home/david/benchmark_master_model/results",
       width = 25, height = 15, units = "cm")


