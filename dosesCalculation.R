rm(list = ls(all = TRUE))

# declare the project root
here::i_am("dosesCalculation.R")

library(tidyverse)
source(here::here("utils.R"))


# Input parameters --------------------------------------------------------

numOfDosesTab <- tibble("rate_pla" = numeric(), "n_target_cases" = numeric(), "AbSampleSize" = numeric(),
                        "nullHR" = numeric(), "altHR" = numeric(), "tau" = numeric(), 
                        "minNumOfDoses" = numeric(),
                        "Q1NumOfDoses" = numeric(),
                        "meanNumOfDoses" = numeric(),
                        "Q2NumOfDoses" = numeric(),
                        "Q3NumOfDoses" = numeric(),
                        "maxNumOfDoses" = numeric())

# 2-arm designs and low-dose vs pla comparisons in 3-arm designs
n_target_cases <- c(21, 22, 40, 135)
rate_pla <- 0.026
nullHR <- 0.75
altHR <- c(0.15, 0.15, 0.25, 0.35)
rate_cens <- 0.075
tau <- c(1, rep(1.5, 3))
iter <- 2000

for (i in 1:4){
  df <- oper_chars_eff_phase_dosesCalculation_twoArms(n_target_cases = n_target_cases[i],
                                                      rate_pla = rate_pla,
                                                      nullHR = nullHR,
                                                      altHR = altHR[i],
                                                      rate_cens = rate_cens,
                                                      tau = tau[i],
                                                      iter = iter)
  
  numOfDosesTab <- add_row(.data = numOfDosesTab, "rate_pla" = rate_pla,
                           "n_target_cases" = n_target_cases[i], 
                           "AbSampleSize" = df$n_enrolled[1] / 2,
                           "nullHR" = nullHR, "altHR" = altHR[i], "tau" = tau[i], 
                           "minNumOfDoses" = min(df$numOfDoses),
                           "Q1NumOfDoses" = quantile(df$numOfDoses, probs = 0.25),
                           "meanNumOfDoses" = mean(df$numOfDoses),
                           "Q2NumOfDoses" = quantile(df$numOfDoses, probs = 0.5),
                           "Q3NumOfDoses" = quantile(df$numOfDoses, probs = 0.75),
                           "maxNumOfDoses" = max(df$numOfDoses))
}

numOfDosesTab