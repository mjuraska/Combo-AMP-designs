rm(list=ls(all=TRUE))
# declare the project root
here::i_am("comboAMP_exponentialModel.R")

library(tidyverse)
source(here::here("utils3.R"))

# Input parameters --------------------------------------------------------
#0.032 or 0.026 for placebo incidence rate
rate_pla <- 0.026
nullHR <- 0.7
rate_cens <- 0.075

iter <- 2000
# correlates expansion phase
n_target_cases_ab <- 35


# null hypothesis
nullHR <- 0.7

# alternative hypothesis
altHRs <- c(0.1, 0.15, 0.1, 0.15, 0.25, 0.35)
taus <- c(1, 1, 1.5, 1.5, 1.5, 1.5)
ans <- list()
for(i in 1:length(altHRs)){
  altHR <- altHRs[i]
  tau <- taus[i]
  power_wald_v <- NULL
  power_cuminc_v <- NULL
  totalSampleSize <- NULL
  nAbSampleSize <- NULL
  meanDuration <- NULL
  meanDurationCorrelates <- NULL
  expansionSampleSize <- NULL
  correlatesNotEnrolled <- NULL
  if(altHR %in% c(0.1, 0.15)){
    n_target_cases_v <- seq(15, 35, 1)
    n_to_enroll <- 5000
  }else if (altHR == 0.25){
    n_target_cases_v <- seq(40, 70, 1)
    n_to_enroll <- 0
  }else if (altHR == 0.35){
    n_target_cases_v <- seq(80, 250, 1)
    n_to_enroll <- 0
  }
  
  for(j in 1:length(n_target_cases_v)){
    n_target_cases <- n_target_cases_v[j]
    df <- oper_chars_eff_phase(n_target_cases = n_target_cases_v[j],
                               rate_pla = rate_pla,
                               nullHR = nullHR,
                               altHR = altHR,
                               rate_cens = rate_cens,
                               tau = tau,
                               iter = iter)
    power_wald_v[j] <- mean(df$wald_pval <= 0.025)
    power_cuminc_v[j]  <- mean(df$cuminc_pval <= 0.025, na.rm = TRUE) 
    totalSampleSize[j] <- mean(df$n_enrolled)
    nAbSampleSize[j] <- round(mean(df$n1), 0)
    meanDuration[j] <- mean(df$analysisTime)
    expansionSampleSize[j] <-  round(df$n_enrolled[1] - n_target_cases - df$n_enrolled[1] * tau * rate_cens + n_to_enroll, 0)
    
    df2 <- duration_corr_exp_phase(n_on_study = df$n_enrolled[1] - n_target_cases - df$n_enrolled[1] * tau * rate_cens,
                                   n_to_enroll = n_to_enroll,
                                   n_obs_cases_ab = mean(df$n_cases_ab),
                                   n_target_cases_ab = n_target_cases_ab,
                                   rate_pla = rate_pla,
                                   altHR = altHR,
                                   rate_cens = rate_cens,
                                   iter = iter)
    meanDurationCorrelates[j] <- mean(df2$analysisTime)
    correlatesNotEnrolled[j] <- mean(df2$not_enrolled)
  }
  ans[[i]] <- data.frame(altHR = altHR, tau = tau, rate_pla = rate_pla, n_target_cases_v = n_target_cases_v, power_wald = power_wald_v, power_cuminc = power_cuminc_v,
                       totalSampleSize = totalSampleSize, nAbSampleSize = nAbSampleSize,
                       meanDuration = meanDuration, expansionSampleSize = expansionSampleSize,
                       meanDurationCorrelates = meanDurationCorrelates, correlatesNotEnrolled = correlatesNotEnrolled)
}

saveRDS(ans, here::here(paste0("power_expoentialModel_iterated_", ifelse(rate_pla == 0.032, 1, 2), "utils3.rds")))
tmp <- ans[[1]]




cbind(seq(90, 110, 1), ans[5])
