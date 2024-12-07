rm(list = ls(all = TRUE))

# declare the project root
here::i_am("dosesCalculation.R")

library(tidyverse)
source(here::here("utils.R"))


# Input parameters --------------------------------------------------------

numOfDosesTab <- tibble("rate_pla" = numeric(), "n_target_cases" = numeric(), 
                        "n_ab" = numeric(), "n_pla" = numeric(),
                        "nullHR" = numeric(), "altHR" = numeric(), "tau" = numeric(), 
                        "minNumOfDoses" = numeric(),
                        "Q1NumOfDoses" = numeric(),
                        "meanNumOfDoses" = numeric(),
                        "Q2NumOfDoses" = numeric(),
                        "Q3NumOfDoses" = numeric(),
                        "maxNumOfDoses" = numeric())

# 2-arm designs and low-dose vs pla comparisons in 3-arm designs
# n_target_cases <- c(23, 23, 43, 135)
# rate_pla <- 0.026
# nullHR <- 0.75
# altHR <- c(0.15, 0.15, 0.25, 0.35)
# rate_cens <- 0.075
# p_ab <- 0.5
# p_pla <- 0.5
# n_enroll_4m <- c(1000, 1000, 800, 800)
# tau <- c(1, rep(1.5, 3))
# iter <- 2000

# high-dose vs pla comparisons in 3-arm designs
n_target_cases <- c(37, 108)
rate_pla <- 0.026
nullHR <- 0.75
altHR <- c(0.15, 0.15)
rate_cens <- 0.075
p_ab <- 1/3
p_pla <- 2/3
n_enroll_4m <- c(600, 600)
tau <- c(1.5, 1.5)
iter <- 2000

for (i in 1:length(n_target_cases)){
  df <- oper_chars_eff_phase_dosesCalculation_twoArms(n_target_cases = n_target_cases[i],
                                                      rate_pla = rate_pla,
                                                      nullHR = nullHR,
                                                      altHR = altHR[i],
                                                      rate_cens = rate_cens,
                                                      p_ab = p_ab,
                                                      p_pla = p_pla,
                                                      n_enroll_4m = n_enroll_4m[i],
                                                      tau = tau[i],
                                                      iter = iter)
  
  numOfDosesTab <- add_row(.data = numOfDosesTab, "rate_pla" = rate_pla,
                           "n_target_cases" = n_target_cases[i], 
                           "n_ab" = df$n_ab[1],
                           "n_pla" = df$n_pla[1],
                           "nullHR" = nullHR, "altHR" = altHR[i], "tau" = tau[i], 
                           "minNumOfDoses" = min(df$numOfDoses),
                           "Q1NumOfDoses" = quantile(df$numOfDoses, probs = 0.25),
                           "meanNumOfDoses" = mean(df$numOfDoses),
                           "Q2NumOfDoses" = quantile(df$numOfDoses, probs = 0.5),
                           "Q3NumOfDoses" = quantile(df$numOfDoses, probs = 0.75),
                           "maxNumOfDoses" = max(df$numOfDoses))
}

numOfDosesTab