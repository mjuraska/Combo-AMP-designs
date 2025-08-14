# declare the project root
here::i_am("runComboAMP.R")

library(tidyverse)
source(here::here("code/param.R"))
source(here::here("param/utils.R"))


# Find the target endpoint count ------------------------------------------

registerDoParallel(cores = n_cores)

for (n in n_target_cases:250){
  df <- oper_chars_eff_phase(compare = "h", nullHR = nullHR, altHR_h = altHR_h,
                             altHR_l = altHR_l, info_fractions = info_fractions,
                             alpha_1sided = alpha_1sided,
                             n_target_cases = n,
                             rate_pla = rate_pla,  rate_cens = rate_cens,
                             p_ab_h = p_ab_h, p_ab_l = p_ab_l, p_pla = p_pla,
                             n_enroll_m = n_enroll_m, tau = tau, iter = iter,
                             n_cores = n_cores, verbose = FALSE)
  power <- mean(df$reject_H0)
  cat("n =", n, "; power =", power, "\n")

  if (power >= pwr){
    break
  }
}


# Get sample size and endpoint splits -------------------------------------

df <- oper_chars_eff_phase(compare = "h", nullHR = nullHR, altHR_h = altHR_h,
                           altHR_l = altHR_l, info_fractions = info_fractions, 
                           alpha_1sided = alpha_1sided, 
                           n_target_cases = n_target_cases,
                           rate_pla = rate_pla,  rate_cens = rate_cens,
                           p_ab_h = p_ab_h, p_ab_l = p_ab_l, p_pla = p_pla,
                           n_enroll_m = n_enroll_m, tau = tau, iter = iter,
                           n_cores = n_cores, verbose = TRUE)
mean(df$reject_H0)

plot_time_to_end_stage1(df, path = path)
# plot_fu_time_eff_phase(df, path = path)
plot_h_vs_p_case_split_end_stage1(df, path = path)
plot_l_case_count_end_stage1(df, path = path)
# plot_n_doses(df, var_name = "n_doses_ab_h", path = path,
#              x_lab = "Number of Administered High Doses of Ab\nby Primary Analysis",
#              title = "High-Dose Ab Arm", file_name = "n_doses_ab_h.pdf")
# plot_n_doses(df, var_name = "n_doses_ab_l", path = path,
#              x_lab = "Number of Administered Low Doses of Ab\nby Primary Analysis",
#              title = "Low-Dose Ab Arm", file_name = "n_doses_ab_l.pdf")

# power
mean(df$pval <= size)
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

