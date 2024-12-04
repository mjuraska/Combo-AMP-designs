# declare the project root
here::i_am("comboAMP_exponentialModel.R")

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
n_target_cases <- 17
rate_pla <- 0.032
nullHR <- 0.7
altHR <- 0.1
rate_cens <- 0.075
tau <- 1
iter <- 2000


df <- oper_chars_eff_phase_dosesCalculation_twoArms(n_target_cases = n_target_cases,
                           rate_pla = rate_pla,
                           nullHR = nullHR,
                           altHR = altHR,
                           rate_cens = rate_cens,
                           tau = tau,
                           iter = iter)

numOfDosesTab <- add_row(.data = numOfDosesTab, "rate_pla" = rate_pla,
                         "n_target_cases" = n_target_cases, 
                         "AbSampleSize" = df$n_enrolled[1]/2,
                         "nullHR" = nullHR, "altHR" = altHR, "tau" = tau, 
                         "minNumOfDoses" = min(df$numOfDoses),
                         "Q1NumOfDoses" = quantile(df$numOfDoses, probs = 0.25),
                         "meanNumOfDoses" = mean(df$numOfDoses),
                         "Q2NumOfDoses" = quantile(df$numOfDoses, probs = 0.5),
                         "Q3NumOfDoses" = quantile(df$numOfDoses, probs = 0.75),
                         "maxNumOfDoses" = max(df$numOfDoses))