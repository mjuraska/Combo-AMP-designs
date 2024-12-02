# declare the project root
here::i_am("comboAMP_exponentialModel.R")

library(tidyverse)
source(here::here("utils.R"))


# Input parameters --------------------------------------------------------

n_target_cases <- 40
rate_pla <- 0.032
nullHR <- 0.7
altHR <- 0.25
rate_cens <- 0.075
tau <- 78 / 52
iter <- 2000
# correlates expansion phase
n_to_enroll <- 4000
n_target_cases_ab <- 35
# 2-arm design
# altHR_high <- altHR
# 3-arm design
altHR_high <- 0.15


# Run the simulation ------------------------------------------------------

df <- oper_chars_eff_phase(n_target_cases = n_target_cases,
                           rate_pla = rate_pla,
                           nullHR = nullHR,
                           altHR = altHR,
                           rate_cens = rate_cens,
                           tau = tau,
                           iter = iter)

# power
mean(df$wald_pval <= 0.025)
mean(df$cuminc_pval <= 0.025)

# check the time when the target number of events is accrued
summary(df$analysisTime)

# are there iterations where the target event count was reached before the
# enrollment was complete?
summary(df$not_enrolled)

# 2-arm design
# n_on_study <- df$n_enrolled[1] - 
#   n_target_cases - df$n_enrolled[1] * tau * rate_cens
# 3-arm design
n_on_study <- (5 / 4) * df$n_enrolled[1] - 
  n_target_cases - altHR_high * (n_target_cases - df$n1[1]) / 2 - 
  (5 / 4) * df$n_enrolled[1] * tau * rate_cens

# 2-arm design
# n_obs_cases_ab <- mean(df$n_cases_ab)
# 3-arm design
n_obs_cases_ab <- mean(df$n_cases_ab) + 
  altHR_high * (n_target_cases - mean(df$n_cases_ab)) / 2

df2 <- duration_corr_exp_phase(n_on_study = n_on_study,
                               n_to_enroll = n_to_enroll,
                               n_obs_cases_ab = n_obs_cases_ab,
                               n_target_cases_ab = n_target_cases_ab,
                               rate_pla = rate_pla,
                               altHR = altHR_high,
                               rate_cens = rate_cens,
                               iter = iter)

# check the time when the target number of events is accrued in the Ab arm
summary(df2$analysisTime)
#Correlates expansion sample size
round(df$n_enrolled[1] - n_target_cases - df$n_enrolled[1] * tau * rate_cens + n_to_enroll, 0)
# are there iterations where the target event count in the Ab arm was reached
# before the enrollment was complete?
summary(df2$not_enrolled)

# event splits
# df <- df %>%
#   group_by(n_cases_pla, n_cases_ab) %>%
#   summarise(prob = n() / max(iter))
# df
