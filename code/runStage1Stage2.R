source(here::here("code/param.R"))
source(here::here("code/utils.R"))


# Run both Stage 1 and Stage 2 --------------------------------------------

registerDoParallel(cores = n_cores)

df <- run_stage1_stage2(compare = "h", nullHR = nullHR, altHR_h = altHR_h,
                        altHR_l = altHR_l, nullHR_h_l = nullHR_h_l, 
                        info_fractions_stage1 = info_fractions, 
                        alpha_1sided = alpha_1sided, 
                        n_target_cases = n_target_cases,
                        n_stage2_target_cases_h_l = n_stage2_target_cases_h_l,
                        rate_pla = rate_pla,  rate_cens = rate_cens,
                        p_ab_h = p_ab_h, p_ab_l = p_ab_l, p_pla = p_pla,
                        n_enroll_m = n_enroll_m, tau = tau, iter = iter,
                        n_cores = n_cores, verbose = TRUE)

saveRDS(df, here::here(path, "runStage1Stage2.rds"))


# Plot results ------------------------------------------------------------

# df <- readRDS(here::here(path, "runStage1Stage2.rds")) %>%
df <- df %>%
  mutate(n_cases_ab = n_cases_ab_h + n_cases_ab_l,
         ind_trial_complete_end_stage1 = as.numeric(init_stage2 == 0 | (init_stage2 == 1 & stage2_stop_IA == 1)))

mean(df$n_cases_ab)
mean(df$n_cases_ab >= 36)

# probability of stopping at the end of Stage 1
mean(df$ind_trial_complete_end_stage1)

# probability of rejecting Stage 2 H0 if Stage 2 is initiated
mean(df$stage2_reject_H0, na.rm=TRUE)

plot_case_count_end_stage2(df, varname = "n_cases_ab", path = path)
plot_case_count_end_stage2(df, varname = "n_cases_ab_l", path = path)
plot_case_count_end_stage2(df, varname = "n_cases_ab_h", path = path)
plot_time_to_complete_trial(df, path)
plot_stage2_IA_count(df, path)
