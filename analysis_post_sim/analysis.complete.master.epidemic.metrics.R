
# To produce incidence trend


analysis.complete.master.epidemic.metrics <- function(keyword = "master.model.", 
                                                      epi.metric.df = df){
  
  
  
  # Incidence trend  ------------------
  
  
  incidence.df <- epi.metric.df %>%
    select(contains("incid.")) 
  
  
  # 35 - 36
  
  incid.35.36.df <- incidence.df %>%
    select(contains("metr.incid.35.36.")) 
  
  incid.35.36.m.15.24 <- quant.mean(incid.35.36.df$metr.incid.35.36.m.15.24)
  incid.35.36.m.25.29 <- quant.mean(incid.35.36.df$metr.incid.35.36.m.25.29)
  incid.35.36.m.30.34 <- quant.mean(incid.35.36.df$metr.incid.35.36.m.30.34)
  incid.35.36.m.35.39 <- quant.mean(incid.35.36.df$metr.incid.35.36.m.35.39)
  incid.35.36.m.40.44 <- quant.mean(incid.35.36.df$metr.incid.35.36.m.40.44)
  incid.35.36.m.45.49 <- quant.mean(incid.35.36.df$metr.incid.35.36.m.45.49)
  incid.35.36.w.15.24 <- quant.mean(incid.35.36.df$metr.incid.35.36.w.15.24)
  incid.35.36.w.25.29 <- quant.mean(incid.35.36.df$metr.incid.35.36.w.25.29)
  incid.35.36.w.30.34 <- quant.mean(incid.35.36.df$metr.incid.35.36.w.30.34)
  incid.35.36.w.35.39 <- quant.mean(incid.35.36.df$metr.incid.35.36.w.35.39)
  incid.35.36.w.40.44 <- quant.mean(incid.35.36.df$metr.incid.35.36.w.40.44)
  incid.35.36.w.45.49 <- quant.mean(incid.35.36.df$metr.incid.35.36.w.45.49)
  
  
  Gender <- c(rep("men", 6), rep("women", 6))
  age_group <- rep(c("15-24", "25-29", "30-34", "35-39", "40-44", "45-49"), 2)
  
  val.inc.35.36.F <- c(incid.35.36.m.15.24[2], incid.35.36.m.25.29[2], incid.35.36.m.30.34[2],
                       incid.35.36.m.35.39[2], incid.35.36.m.40.44[2], incid.35.36.m.45.49[2],
                       incid.35.36.w.15.24[2], incid.35.36.w.25.29[2], incid.35.36.w.30.34[2],
                       incid.35.36.w.35.39[2], incid.35.36.w.40.44[2], incid.35.36.w.45.49[2])
  
  val.inc.35.36.L <- c(incid.35.36.m.15.24[1], incid.35.36.m.25.29[1], incid.35.36.m.30.34[1],
                       incid.35.36.m.35.39[1], incid.35.36.m.40.44[1], incid.35.36.m.45.49[1],
                       incid.35.36.w.15.24[1], incid.35.36.w.25.29[1], incid.35.36.w.30.34[1],
                       incid.35.36.w.35.39[1], incid.35.36.w.40.44[1], incid.35.36.w.45.49[1])
  
  val.inc.35.36.U <- c(incid.35.36.m.15.24[3], incid.35.36.m.25.29[3], incid.35.36.m.30.34[3],
                       incid.35.36.m.35.39[3], incid.35.36.m.40.44[3], incid.35.36.m.45.49[3],
                       incid.35.36.w.15.24[3], incid.35.36.w.25.29[3], incid.35.36.w.30.34[3],
                       incid.35.36.w.35.39[3], incid.35.36.w.40.44[3], incid.35.36.w.45.49[3])
  
  
  val.inc.35.36 <- data.frame(val.inc.35.36.L, val.inc.35.36.F, val.inc.35.36.U, Gender, age_group)
  
  table.val.inc.35.36 <- data.frame(age_group, val.inc.35.36.L, val.inc.35.36.F, val.inc.35.36.U, Gender)
  
  
  names(table.val.inc.35.36) <- c("age_group", "lower.Q1", "mean", "upper.Q3", "Gender") # c("age_group", "lower.Q1", "median", "upper.Q3", "Gender")

  
  # table.val.inc.35.36 %>% 
  #   kable() %>% 
  #   kable_styling("striped")
  
  
  write.csv(table.val.inc.35.36, file = paste0("/home/david/benchmark_master_model/results/", keyword, "table.val.inc.35.36.csv"))
  
  
  
  
  # 36 - 37
  
  incid.36.37.df <- incidence.df %>%
    select(contains("metr.incid.36.37.")) 
  
  
  incid.36.37.m.15.24 <- quant.mean(incid.36.37.df$metr.incid.36.37.m.15.24)
  incid.36.37.m.25.29 <- quant.mean(incid.36.37.df$metr.incid.36.37.m.25.29)
  incid.36.37.m.30.34 <- quant.mean(incid.36.37.df$metr.incid.36.37.m.30.34)
  incid.36.37.m.35.39 <- quant.mean(incid.36.37.df$metr.incid.36.37.m.35.39)
  incid.36.37.m.40.44 <- quant.mean(incid.36.37.df$metr.incid.36.37.m.40.44)
  incid.36.37.m.45.49 <- quant.mean(incid.36.37.df$metr.incid.36.37.m.45.49)
  incid.36.37.w.15.24 <- quant.mean(incid.36.37.df$metr.incid.36.37.w.15.24)
  incid.36.37.w.25.29 <- quant.mean(incid.36.37.df$metr.incid.36.37.w.25.29)
  incid.36.37.w.30.34 <- quant.mean(incid.36.37.df$metr.incid.36.37.w.30.34)
  incid.36.37.w.35.39 <- quant.mean(incid.36.37.df$metr.incid.36.37.w.35.39)
  incid.36.37.w.40.44 <- quant.mean(incid.36.37.df$metr.incid.36.37.w.40.44)
  incid.36.37.w.45.49 <- quant.mean(incid.36.37.df$metr.incid.36.37.w.45.49)
  
  
  val.inc.36.37.F <- c(incid.36.37.m.15.24[2], incid.36.37.m.25.29[2], incid.36.37.m.30.34[2],
                       incid.36.37.m.35.39[2], incid.36.37.m.40.44[2], incid.36.37.m.45.49[2],
                       incid.36.37.w.15.24[2], incid.36.37.w.25.29[2], incid.36.37.w.30.34[2],
                       incid.36.37.w.35.39[2], incid.36.37.w.40.44[2], incid.36.37.w.45.49[2])
  
  val.inc.36.37.L <- c(incid.36.37.m.15.24[1], incid.36.37.m.25.29[1], incid.36.37.m.30.34[1],
                       incid.36.37.m.35.39[1], incid.36.37.m.40.44[1], incid.36.37.m.45.49[1],
                       incid.36.37.w.15.24[1], incid.36.37.w.25.29[1], incid.36.37.w.30.34[1],
                       incid.36.37.w.35.39[1], incid.36.37.w.40.44[1], incid.36.37.w.45.49[1])
  
  val.inc.36.37.U <- c(incid.36.37.m.15.24[3], incid.36.37.m.25.29[3], incid.36.37.m.30.34[3],
                       incid.36.37.m.35.39[3], incid.36.37.m.40.44[3], incid.36.37.m.45.49[3],
                       incid.36.37.w.15.24[3], incid.36.37.w.25.29[3], incid.36.37.w.30.34[3],
                       incid.36.37.w.35.39[3], incid.36.37.w.40.44[3], incid.36.37.w.45.49[3])
  
  
  val.inc.36.37 <- data.frame(val.inc.36.37.L, val.inc.36.37.F, val.inc.36.37.U, Gender, age_group)
  
  table.val.inc.36.37 <- data.frame(age_group, val.inc.36.37.L, val.inc.36.37.F, val.inc.36.37.U, Gender)
  
  
  names(table.val.inc.36.37) <- c("age_group", "lower.Q1", "mean", "upper.Q3", "Gender") # c("age_group", "lower.Q1", "median", "upper.Q3", "Gender")

  
  # table.val.inc.36.37 %>% 
  #   kable() %>% 
  #   kable_styling("striped")
  
  
  write.csv(table.val.inc.36.37, file = paste0("/home/david/benchmark_master_model/results/", keyword, "table.val.inc.36.37.csv"))
  
  
  
  
  
  # 37 - 38
  
  incid.37.38.df <- incidence.df %>%
    select(contains("metr.incid.37.38.")) 
  
  
  
  incid.37.38.m.15.24 <- quant.mean(incid.37.38.df$metr.incid.37.38.m.15.24)
  incid.37.38.m.25.29 <- quant.mean(incid.37.38.df$metr.incid.37.38.m.25.29)
  incid.37.38.m.30.34 <- quant.mean(incid.37.38.df$metr.incid.37.38.m.30.34)
  incid.37.38.m.35.39 <- quant.mean(incid.37.38.df$metr.incid.37.38.m.35.39)
  incid.37.38.m.40.44 <- quant.mean(incid.37.38.df$metr.incid.37.38.m.40.44)
  incid.37.38.m.45.49 <- quant.mean(incid.37.38.df$metr.incid.37.38.m.45.49)
  incid.37.38.w.15.24 <- quant.mean(incid.37.38.df$metr.incid.37.38.w.15.24)
  incid.37.38.w.25.29 <- quant.mean(incid.37.38.df$metr.incid.37.38.w.25.29)
  incid.37.38.w.30.34 <- quant.mean(incid.37.38.df$metr.incid.37.38.w.30.34)
  incid.37.38.w.35.39 <- quant.mean(incid.37.38.df$metr.incid.37.38.w.35.39)
  incid.37.38.w.40.44 <- quant.mean(incid.37.38.df$metr.incid.37.38.w.40.44)
  incid.37.38.w.45.49 <- quant.mean(incid.37.38.df$metr.incid.37.38.w.45.49)
  
  
  val.inc.37.38.F <- c(incid.37.38.m.15.24[2], incid.37.38.m.25.29[2], incid.37.38.m.30.34[2],
                       incid.37.38.m.35.39[2], incid.37.38.m.40.44[2], incid.37.38.m.45.49[2],
                       incid.37.38.w.15.24[2], incid.37.38.w.25.29[2], incid.37.38.w.30.34[2],
                       incid.37.38.w.35.39[2], incid.37.38.w.40.44[2], incid.37.38.w.45.49[2])
  
  val.inc.37.38.L <- c(incid.37.38.m.15.24[1], incid.37.38.m.25.29[1], incid.37.38.m.30.34[1],
                       incid.37.38.m.35.39[1], incid.37.38.m.40.44[1], incid.37.38.m.45.49[1],
                       incid.37.38.w.15.24[1], incid.37.38.w.25.29[1], incid.37.38.w.30.34[1],
                       incid.37.38.w.35.39[1], incid.37.38.w.40.44[1], incid.37.38.w.45.49[1])
  
  val.inc.37.38.U <- c(incid.37.38.m.15.24[3], incid.37.38.m.25.29[3], incid.37.38.m.30.34[3],
                       incid.37.38.m.35.39[3], incid.37.38.m.40.44[3], incid.37.38.m.45.49[3],
                       incid.37.38.w.15.24[3], incid.37.38.w.25.29[3], incid.37.38.w.30.34[3],
                       incid.37.38.w.35.39[3], incid.37.38.w.40.44[3], incid.37.38.w.45.49[3])
  
  
  val.inc.37.38 <- data.frame(val.inc.37.38.L, val.inc.37.38.F, val.inc.37.38.U, Gender, age_group)
  
  table.val.inc.37.38 <- data.frame(age_group, val.inc.37.38.L, val.inc.37.38.F, val.inc.37.38.U, Gender)
  
  
  names(table.val.inc.37.38) <- c("age_group", "lower.Q1", "mean", "upper.Q3", "Gender") # c("age_group", "lower.Q1", "median", "upper.Q3", "Gender")

  
  # table.val.inc.37.38 %>% 
  #   kable() %>% 
  #   kable_styling("striped")
  
  write.csv(table.val.inc.37.38, file = paste0("/home/david/benchmark_master_model/results/", keyword, "table.val.inc.37.38.csv"))
  
  
  
  
  # 38 - 39
  
  incid.38.39.df <- incidence.df %>%
    select(contains("metr.incid.38.39.")) 
  
  
  
  incid.38.39.m.15.24 <- quant.mean(incid.38.39.df$metr.incid.38.39.m.15.24)
  incid.38.39.m.25.29 <- quant.mean(incid.38.39.df$metr.incid.38.39.m.25.29)
  incid.38.39.m.30.34 <- quant.mean(incid.38.39.df$metr.incid.38.39.m.30.34)
  incid.38.39.m.35.39 <- quant.mean(incid.38.39.df$metr.incid.38.39.m.35.39)
  incid.38.39.m.40.44 <- quant.mean(incid.38.39.df$metr.incid.38.39.m.40.44)
  incid.38.39.m.45.49 <- quant.mean(incid.38.39.df$metr.incid.38.39.m.45.49)
  incid.38.39.w.15.24 <- quant.mean(incid.38.39.df$metr.incid.38.39.w.15.24)
  incid.38.39.w.25.29 <- quant.mean(incid.38.39.df$metr.incid.38.39.w.25.29)
  incid.38.39.w.30.34 <- quant.mean(incid.38.39.df$metr.incid.38.39.w.30.34)
  incid.38.39.w.35.39 <- quant.mean(incid.38.39.df$metr.incid.38.39.w.35.39)
  incid.38.39.w.40.44 <- quant.mean(incid.38.39.df$metr.incid.38.39.w.40.44)
  incid.38.39.w.45.49 <- quant.mean(incid.38.39.df$metr.incid.38.39.w.45.49)
  
  
  val.inc.38.39.F <- c(incid.38.39.m.15.24[2], incid.38.39.m.25.29[2], incid.38.39.m.30.34[2],
                       incid.38.39.m.35.39[2], incid.38.39.m.40.44[2], incid.38.39.m.45.49[2],
                       incid.38.39.w.15.24[2], incid.38.39.w.25.29[2], incid.38.39.w.30.34[2],
                       incid.38.39.w.35.39[2], incid.38.39.w.40.44[2], incid.38.39.w.45.49[2])
  
  val.inc.38.39.L <- c(incid.38.39.m.15.24[1], incid.38.39.m.25.29[1], incid.38.39.m.30.34[1],
                       incid.38.39.m.35.39[1], incid.38.39.m.40.44[1], incid.38.39.m.45.49[1],
                       incid.38.39.w.15.24[1], incid.38.39.w.25.29[1], incid.38.39.w.30.34[1],
                       incid.38.39.w.35.39[1], incid.38.39.w.40.44[1], incid.38.39.w.45.49[1])
  
  val.inc.38.39.U <- c(incid.38.39.m.15.24[3], incid.38.39.m.25.29[3], incid.38.39.m.30.34[3],
                       incid.38.39.m.35.39[3], incid.38.39.m.40.44[3], incid.38.39.m.45.49[3],
                       incid.38.39.w.15.24[3], incid.38.39.w.25.29[3], incid.38.39.w.30.34[3],
                       incid.38.39.w.35.39[3], incid.38.39.w.40.44[3], incid.38.39.w.45.49[3])
  
  
  val.inc.38.39 <- data.frame(val.inc.38.39.L, val.inc.38.39.F, val.inc.38.39.U, Gender, age_group)
  
  table.val.inc.38.39 <- data.frame(age_group, val.inc.38.39.L, val.inc.38.39.F, val.inc.38.39.U, Gender)
  
  
  names(table.val.inc.38.39) <- c("age_group", "lower.Q1", "mean", "upper.Q3", "Gender") # c("age_group", "lower.Q1", "median", "upper.Q3", "Gender")
  
  # table.val.inc.38.39 %>% 
  #   kable() %>% 
  #   kable_styling("striped")
  
  
  write.csv(table.val.inc.38.39, file = paste0("/home/david/benchmark_master_model/results/", keyword, "table.val.inc.38.39.csv"))
  
  
  
  
  # 39 - 40
  
  incid.39.40.df <- incidence.df %>%
    select(contains("metr.incid.39.40.")) 
  
  
  
  incid.39.40.m.15.24 <- quant.mean(incid.39.40.df$metr.incid.39.40.m.15.24)
  incid.39.40.m.25.29 <- quant.mean(incid.39.40.df$metr.incid.39.40.m.25.29)
  incid.39.40.m.30.34 <- quant.mean(incid.39.40.df$metr.incid.39.40.m.30.34)
  incid.39.40.m.35.39 <- quant.mean(incid.39.40.df$metr.incid.39.40.m.35.39)
  incid.39.40.m.40.44 <- quant.mean(incid.39.40.df$metr.incid.39.40.m.40.44)
  incid.39.40.m.45.49 <- quant.mean(incid.39.40.df$metr.incid.39.40.m.45.49)
  incid.39.40.w.15.24 <- quant.mean(incid.39.40.df$metr.incid.39.40.w.15.24)
  incid.39.40.w.25.29 <- quant.mean(incid.39.40.df$metr.incid.39.40.w.25.29)
  incid.39.40.w.30.34 <- quant.mean(incid.39.40.df$metr.incid.39.40.w.30.34)
  incid.39.40.w.35.39 <- quant.mean(incid.39.40.df$metr.incid.39.40.w.35.39)
  incid.39.40.w.40.44 <- quant.mean(incid.39.40.df$metr.incid.39.40.w.40.44)
  incid.39.40.w.45.49 <- quant.mean(incid.39.40.df$metr.incid.39.40.w.45.49)
  
  
  val.inc.39.40.F <- c(incid.39.40.m.15.24[2], incid.39.40.m.25.29[2], incid.39.40.m.30.34[2],
                       incid.39.40.m.35.39[2], incid.39.40.m.40.44[2], incid.39.40.m.45.49[2],
                       incid.39.40.w.15.24[2], incid.39.40.w.25.29[2], incid.39.40.w.30.34[2],
                       incid.39.40.w.35.39[2], incid.39.40.w.40.44[2], incid.39.40.w.45.49[2])
  
  val.inc.39.40.L <- c(incid.39.40.m.15.24[1], incid.39.40.m.25.29[1], incid.39.40.m.30.34[1],
                       incid.39.40.m.35.39[1], incid.39.40.m.40.44[1], incid.39.40.m.45.49[1],
                       incid.39.40.w.15.24[1], incid.39.40.w.25.29[1], incid.39.40.w.30.34[1],
                       incid.39.40.w.35.39[1], incid.39.40.w.40.44[1], incid.39.40.w.45.49[1])
  
  val.inc.39.40.U <- c(incid.39.40.m.15.24[3], incid.39.40.m.25.29[3], incid.39.40.m.30.34[3],
                       incid.39.40.m.35.39[3], incid.39.40.m.40.44[3], incid.39.40.m.45.49[3],
                       incid.39.40.w.15.24[3], incid.39.40.w.25.29[3], incid.39.40.w.30.34[3],
                       incid.39.40.w.35.39[3], incid.39.40.w.40.44[3], incid.39.40.w.45.49[3])
  
  
  val.inc.39.40 <- data.frame(val.inc.39.40.L, val.inc.39.40.F, val.inc.39.40.U, Gender, age_group)
  
  table.val.inc.39.40 <- data.frame(age_group, val.inc.39.40.L, val.inc.39.40.F, val.inc.39.40.U, Gender)
  
  
  names(table.val.inc.39.40) <- c("age_group", "lower.Q1", "mean", "upper.Q3", "Gender")
  
  # table.val.inc.39.40 %>% 
  #   kable() %>% 
  #   kable_styling("striped")
  
  write.csv(table.val.inc.39.40, file = paste0("/home/david/benchmark_master_model/results/", keyword, "table.val.inc.39.40.csv"))
  
  
  
  
  
  # Plots
  
  plot.15.24.df.males <- c(table.val.inc.35.36[1,3], table.val.inc.36.37[1,3], table.val.inc.37.38[1,3],
                           table.val.inc.38.39[1,3], table.val.inc.39.40[1,3])
  
  plot.15.24.df.females <- c(table.val.inc.35.36[7,3], table.val.inc.36.37[7,3], table.val.inc.37.38[7,3],
                             table.val.inc.38.39[7,3], table.val.inc.39.40[7,3])
  
  
  plot.25.29.df.males <- c(table.val.inc.35.36[2,3], table.val.inc.36.37[2,3], table.val.inc.37.38[2,3],
                           table.val.inc.38.39[2,3], table.val.inc.39.40[2,3])
  
  plot.25.29.df.females <- c(table.val.inc.35.36[8,3], table.val.inc.36.37[8,3], table.val.inc.37.38[8,3],
                             table.val.inc.38.39[8,3], table.val.inc.39.40[8,3])
  
  plot.30.34.df.males <- c(table.val.inc.35.36[3,3], table.val.inc.36.37[3,3], table.val.inc.37.38[3,3],
                           table.val.inc.38.39[3,3], table.val.inc.39.40[3,3])
  
  plot.30.34.df.females <- c(table.val.inc.35.36[9,3], table.val.inc.36.37[9,3], table.val.inc.37.38[9,3],
                             table.val.inc.38.39[9,3], table.val.inc.39.40[9,3])
  
  
  
  plot.35.39.df.males <- c(table.val.inc.35.36[4,3], table.val.inc.36.37[4,3], table.val.inc.37.38[4,3],
                           table.val.inc.38.39[4,3], table.val.inc.39.40[4,3])
  
  plot.35.39.df.females <- c(table.val.inc.35.36[10,3], table.val.inc.36.37[10,3], table.val.inc.37.38[10,3],
                             table.val.inc.38.39[10,3], table.val.inc.39.40[10,3])
  
  
  plot.40.44.df.males <- c(table.val.inc.35.36[5,3], table.val.inc.36.37[5,3], table.val.inc.37.38[5,3],
                           table.val.inc.38.39[5,3], table.val.inc.39.40[5,3])
  
  plot.40.44.df.females <- c(table.val.inc.35.36[11,3], table.val.inc.36.37[11,3], table.val.inc.37.38[11,3],
                             table.val.inc.38.39[11,3], table.val.inc.39.40[11,3])
  
  
  plot.45.49.df.males <- c(table.val.inc.35.36[6,3], table.val.inc.36.37[6,3], table.val.inc.37.38[6,3],
                           table.val.inc.38.39[6,3], table.val.inc.39.40[6,3])
  
  plot.45.49.df.females <- c(table.val.inc.35.36[12,3], table.val.inc.36.37[12,3], table.val.inc.37.38[12,3],
                             table.val.inc.38.39[12,3], table.val.inc.39.40[12,3])
  
  
  
  # Time <- c("35-36", "36-37", "37-38", "38-39", "39-40")
  
  Time <- c(2013, 2014, 2015, 2016, 2017)
  
  plot.15.24.df.m <- data.frame(plot.15.24.df.males, Time) 
  plot.15.24.df.f <- data.frame(plot.15.24.df.females, Time) 
  plot.25.29.df.m <- data.frame(plot.25.29.df.males, Time) 
  plot.25.29.df.f <- data.frame(plot.25.29.df.females, Time) 
  plot.30.34.df.m <- data.frame(plot.30.34.df.males, Time) 
  plot.30.34.df.f <- data.frame(plot.30.34.df.females, Time) 
  plot.35.39.df.m <- data.frame(plot.35.39.df.males, Time) 
  plot.35.39.df.f <- data.frame(plot.35.39.df.females, Time) 
  plot.40.44.df.m <- data.frame(plot.40.44.df.males, Time) 
  plot.40.44.df.f <- data.frame(plot.40.44.df.females, Time) 
  plot.45.49.df.m <- data.frame(plot.45.49.df.males, Time) 
  plot.45.49.df.f <- data.frame(plot.45.49.df.females, Time) 
  
  
  names(plot.15.24.df.m) <- c("Incidence", "Time")
  names(plot.15.24.df.f) <- c("Incidence", "Time") 
  names(plot.25.29.df.m) <- c("Incidence", "Time")
  names(plot.25.29.df.f) <- c("Incidence", "Time")
  names(plot.30.34.df.m) <- c("Incidence", "Time")
  names(plot.30.34.df.f) <- c("Incidence", "Time")
  names(plot.35.39.df.m) <- c("Incidence", "Time")
  names(plot.35.39.df.f) <- c("Incidence", "Time")
  names(plot.40.44.df.m) <- c("Incidence", "Time")
  names(plot.40.44.df.f) <- c("Incidence", "Time") 
  names(plot.45.49.df.m) <- c("Incidence", "Time") 
  names(plot.45.49.df.f) <- c("Incidence", "Time") 
  
  plot.15.24.df <- rbind(plot.15.24.df.m, plot.15.24.df.f)
  plot.25.29.df <- rbind(plot.25.29.df.m, plot.25.29.df.f)
  plot.30.34.df <- rbind(plot.30.34.df.m, plot.30.34.df.f)
  plot.35.39.df <- rbind(plot.35.39.df.m, plot.35.39.df.f)
  plot.40.44.df <- rbind(plot.40.44.df.m, plot.40.44.df.f)
  plot.45.49.df <- rbind(plot.45.49.df.m, plot.45.49.df.f)
  
  plot.15.24.df$Gender <- c(rep("men", 5), rep("women", 5))
  plot.25.29.df$Gender <- c(rep("men", 5), rep("women", 5))
  plot.30.34.df$Gender <- c(rep("men", 5), rep("women", 5))
  plot.35.39.df$Gender <- c(rep("men", 5), rep("women", 5))
  plot.40.44.df$Gender <- c(rep("men", 5), rep("women", 5))
  plot.45.49.df$Gender <- c(rep("men", 5), rep("women", 5))
  
  
  plot.df <- rbind(plot.15.24.df, plot.25.29.df, plot.30.34.df, 
                   plot.35.39.df, plot.40.44.df, plot.45.49.df)
  
  plot.df$age_group <- c(rep("15 - 24 years", 10), rep("25 - 29 years", 10), rep("30 - 34 years", 10),
                         rep("35 - 39 years", 10), rep("40 - 44 years", 10), rep("45 - 49 years", 10))
  
  plot.df <- rename(plot.df, Year = Time, gender = Gender)
  
  plot.incidence.trend <- ggplot(plot.df, aes(x=Year, y=Incidence, colour=gender, group = gender)) + 
    geom_line(size=1) +
    geom_point(size=2) + 
    facet_wrap(~age_group)+
    xlab("Time (year)") + ylab("HIV incidence")+
    theme(legend.position="bottom")
  
  ggsave(filename = paste0(keyword, "metric.plot.incidence.trend.pdf"),
         plot = plot.incidence.trend,
         path = "/home/david/benchmark_master_model/results",
         width = 22, height = 10, units = "cm")
  
  
  
  
  
  # Age mixing in partnerships ----------------
  
  agemixpartnership.df <- epi.metric.df %>% 
    select("metr.AAD.male", "metr.SDAD.male",  "metr.slope.male",  "metr.WSD.male", "metr.BSD.male" , "metr.intercept.male" )
  
  
  pop.AAD.male <- quant.mean(agemixpartnership.df$metr.AAD.male)
  pop.SDAD.male <- quant.mean(agemixpartnership.df$metr.SDAD.male)
  pop.slope.male <- quant.mean(agemixpartnership.df$metr.slope.male)
  pop.WSD.male <- quant.mean(agemixpartnership.df$metr.WSD.male)
  pop.BSD.male <- quant.mean(agemixpartnership.df$metr.BSD.male)
  pop.intercept.male <- quant.mean(agemixpartnership.df$metr.intercept.male)
  
  age.mix.F <- c(pop.AAD.male[2], pop.SDAD.male[2], pop.BSD.male[2],
                 pop.WSD.male[2], pop.slope.male[2], pop.intercept.male[2])
  
  age.mix.F <- round(age.mix.F, digits = 3)
  
  age.mix.L <- c(pop.AAD.male[1], pop.SDAD.male[1], pop.BSD.male[1],
                 pop.WSD.male[1], pop.slope.male[1], pop.intercept.male[1])
  
  age.mix.L <- round(age.mix.L, digits = 3)
  
  age.mix.U <- c(pop.AAD.male[3], pop.SDAD.male[3], pop.BSD.male[3],
                 pop.WSD.male[3], pop.slope.male[3], pop.intercept.male[3])
  
  age.mix.U <- round(age.mix.U, digits = 3)
  
  param.name <- c("AAD.male", "SDAD.male", "BSD.male" , "WSD.male", "slope.male",   "intercept.male" ) 
  
  age.mixing.pop <- data.frame(param.name, age.mix.L, age.mix.F, age.mix.U)
  
  colnames(age.mixing.pop) <- c("param", "lower.Q1", "mean", "upper.Q3") # c("param", "lower.Q1", "med", "upper.Q3")
  
  # age.mixing.pop %>% 
  #   kable() %>% 
  #   kable_styling("striped") 
  
  
  write.csv(age.mixing.pop, file = paste0("/home/david/benchmark_master_model/results/", keyword, "age.mixing.pop.csv"))
  
  
  
  
  # Onwards transmissions ---------------
  
  onwardstransmission.df <- epi.metric.df %>%
    select(contains(".transm.")) 
  
  
  onwardstransmission.df <- na.omit(onwardstransmission.df)
  
  
  transm.mean <- quant.mean(onwardstransmission.df$metr.transm.av)
  transm.med <- quant.mean(onwardstransmission.df$metr.transm.med)
  transm.sd <- quant.mean(onwardstransmission.df$metr.transm.sd)
  
  transm.F <- c(transm.mean[2], transm.med[2], transm.sd[2])
  
  transm.F <- round(transm.F, digits = 3)
  
  transm.L <- c(transm.mean[1], transm.med[1], transm.sd[1])
  
  transm.L <- round(transm.L, digits = 3)
  
  transm.U <- c(transm.mean[3], transm.med[3], transm.sd[3])
  
  transm.U <- round(transm.U, digits = 3)
  
  params.name <- c("Mean", "Median", "SD") 
  
  onwards.transmission <- data.frame(params.name, transm.L, transm.F, transm.U)
  
  colnames(onwards.transmission) <- c("param", "lower.Q1", "mean", "upper.Q3")
  
  # onwards.transmission %>% 
  #   kable() %>% 
  #   kable_styling("striped") 
  
  write.csv(onwards.transmission, file = paste0("/home/david/benchmark_master_model/results/", keyword, "onwards.transmission.csv"))
  
  
  
}



