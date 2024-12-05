rm(list=ls(all=TRUE))
# declare the project root
here::i_am("comboAMP_exponentialModel.R")

library(tidyverse)
source(here::here("utils.R"))

# Input parameters --------------------------------------------------------
#to change from slide from slide
####################################
nullHR <- 0.7
size = 0.025
rate_pla <- 0.026
####################################


rate_cens <- 0.075
iter <- 2000
# correlates expansion phase
n_target_cases_ab <- 35

# alternative hypothesis
altHRs <- c(0.1, 0.15, 0.1, 0.15, 0.25, 0.25)
altHighHRs <- c(0.1, 0.15, 0.1, 0.15, 0.1, 0.15)
taus <- c(1, 1, 1.5, 1.5, 1.5, 1.5)
designs <- c("2arm", "2arm", "2arm", "2arm", "3arm", "3arm")
ans <- list()
for(i in 1:length(altHRs)){
  altHR <- altHRs[i]
  altHR_high <- altHighHRs[i]
  tau <- taus[i]
  p_ab <- 0.5
  p_pla <- 0.5
  design <- designs[i]
  n_enroll_4m <- ifelse(design == "2arm", 1000, 800)
  power_wald_v <- NULL
  totalSampleSize <- NULL
  nAbSampleSize <- NULL
  meanDuration <- NULL
  meanDurationCorrelates <- NULL
  expansionSampleSize <- NULL
  if(design == "2arm"){
    n_target_cases_v <- seq(15, 40, 1)
    n_to_enroll <- 5000
  }else if (design == "3arm"){
    n_target_cases_v <- seq(40, 60, 1)
    n_to_enroll <- 4000
  }
  
  for(j in 1:length(n_target_cases_v)){
    n_target_cases <- n_target_cases_v[j]
    df <- oper_chars_eff_phase(n_target_cases = n_target_cases_v[j],
                               rate_pla = rate_pla,
                               nullHR = nullHR,
                               altHR = altHR,
                               rate_cens = rate_cens,
                               p_ab = p_ab,
                               p_pla = p_pla,
                               n_enroll_4m = n_enroll_4m,
                               tau = tau,
                               iter = iter,
                               size = size*2)
    
    power_wald_v[j] <- mean(df$wald_pval <= size)
    totalSampleSize[j] <- mean(df$n_enrolled)
    nAbSampleSize[j] <- round(mean(df$n1), 0)
    meanDuration[j] <- mean(df$analysisTime)
    expansionSampleSize[j] <-  round(df$n_enrolled[1] - n_target_cases - df$n_enrolled[1] * tau * rate_cens + n_to_enroll, 0)
    
    if(design == "3arm"){
      n_on_study <- (5 / 4) * df$n_enrolled[1] - n_target_cases - altHR_high * (n_target_cases - df$n1[1]) / 2 - (5 / 4) * df$n_enrolled[1] * tau * rate_cens  
      n_obs_cases_ab <- mean(df$n_cases_ab) + altHR_high * (n_target_cases - mean(df$n_cases_ab)) / 2

    }else{
      n_on_study <- df$n_enrolled[1] - n_target_cases - df$n_enrolled[1] * tau * rate_cens
      n_obs_cases_ab <- mean(df$n_cases_ab)
    }
   
   
    df2 <- duration_corr_exp_phase(n_on_study = n_on_study,
                                   n_to_enroll = n_to_enroll,
                                   n_obs_cases_ab = n_obs_cases_ab,
                                   n_target_cases_ab = n_target_cases_ab,
                                   rate_pla = rate_pla,
                                   altHR = altHR_high,
                                   rate_cens = rate_cens,
                                   iter = iter)
    meanDurationCorrelates[j] <- mean(df2$analysisTime)
    expansionSampleSize <- n_on_study + n_to_enroll
    
    #stop the loop if power is reached for both methods
    if(sum(power_wald_v>=0.9) > 0){
      break
    }
  }
  ans[[i]] <- data.frame(altHR = altHR, tau = tau, rate_pla = rate_pla, n_target_cases_v = n_target_cases_v[1:j], power_wald = power_wald_v, 
                       totalSampleSize = totalSampleSize, nAbSampleSize = nAbSampleSize,
                       meanDuration = meanDuration, expansionSampleSize = expansionSampleSize,
                       meanDurationCorrelates = meanDurationCorrelates)
}

saveRDS(ans, here::here(paste0("power_expoentialModel_iterated_slide2.rds")))
