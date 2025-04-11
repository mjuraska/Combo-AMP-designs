# declare the project root
here::i_am("comboAMP_exponentialModel.R")

library(tidyverse)
source(here::here("utils.R"))


# Input parameters --------------------------------------------------------

# 2-arm design
design <- "3arm"
n_target_cases <- 43
rate_pla <- 0.026
nullHR <- 0.75
altHR <- 0.25
rate_cens <- 0.075
# p_ab <- 0.5
# p_pla <- 0.5
n_enroll_4m <- 1000
tau <- 1.5
size <- 0.025
iter <- 1000
path <- "output"

# 3-arm design
p_ab_h <- 0.2
p_ab_l <- 0.4
p_pla <- 0.4
altHR_high <- 0.15

# correlates expansion phase
n_to_enroll <- 5000
n_target_cases_ab <- 35


# Run the simulation ------------------------------------------------------

df <- oper_chars_eff_phase(n_target_cases = n_target_cases,
                           rate_pla = rate_pla,
                           nullHR = nullHR,
                           altHR = altHR,
                           rate_cens = rate_cens,
                           p_ab = p_ab,
                           p_pla = p_pla,
                           n_enroll_4m = n_enroll_4m,
                           tau = tau,
                           iter = iter)

plot_time_to_analysis(df, path = path)
plot_fu_time_eff_phase(df, path = path)
plot_case_split_eff_phase(df, path = path)

df <- oper_chars_eff_phase_3arm(n_target_cases = n_target_cases,
                                rate_pla = rate_pla,
                                nullHR = nullHR,
                                altHR_l = altHR,
                                altHR_h = altHR_high,
                                rate_cens = rate_cens,
                                p_ab_h = p_ab_h,
                                p_ab_l = p_ab_l,
                                p_pla = p_pla,
                                n_enroll_4m = n_enroll_4m,
                                tau = tau,
                                iter = iter)

plot_case_ab_h_eff_phase(df, path = path)
plot_n_doses(df, var_name = "n_doses_ab_h", path = path,
             x_lab = "Number of Administered High Doses of Ab\nby Primary Analysis",
             title = "High-Dose Ab Arm", file_name = "n_doses_ab_h.pdf")
plot_n_doses(df, var_name = "n_doses_ab_l", path = path,
             x_lab = "Number of Administered Low Doses of Ab\nby Primary Analysis",
             title = "Low-Dose Ab Arm", file_name = "n_doses_ab_l.pdf")

# power
mean(df$wald_pval <= size)
mean(df$meanEventTime)
# check the time when the target number of events is accrued
summary(df$analysisTime)

# are there iterations where the target event count was reached before the
# enrollment was complete?
summary(df$not_enrolled)

if (design == "2arm"){
  # 2-arm design
  n_on_study <- df$n_enrolled[1] - n_target_cases - df$n_enrolled[1] * tau * rate_cens
  n_obs_cases_ab <- mean(df$n_cases_ab)
} else {
  # 3-arm design
  n_on_study <- (5 / 4) * df$n_enrolled[1] - n_target_cases - altHR_high * (n_target_cases - df$n1[1]) / 2 - (5 / 4) * df$n_enrolled[1] * tau * rate_cens  
  n_obs_cases_ab <- mean(df$n_cases_ab) + altHR_high * (n_target_cases - mean(df$n_cases_ab)) / 2
}

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
n_on_study  + n_to_enroll
# are there iterations where the target event count in the Ab arm was reached
# before the enrollment was complete?
summary(df2$not_enrolled)


# Get summaries from simulation results -----------------------------------

#distribution of follow-up time in the primary efficacy assessment
ans_primary <- c(mean(df$minEventTime), mean(df$Q1EventTime), mean(df$meanEventTime), mean(df$medEventTime), mean(df$Q3EventTime), mean(df$maxEventTime))
names(ans_primary) <- c("min", "Q1", "mean", "med", "Q3", "max")
round(ans_primary, 2)
# distribution of follow-up time in the correlate expansion
ans_expansion <- c(mean(df2$minEventTime), mean(df2$Q1EventTime), mean(df2$meanEventTime), mean(df2$Q3EventTime), mean(df2$maxEventTime))
names(ans_expansion) <- c("min", "Q1", "mean", "Q3", "max")
round(ans_expansion, 2)

