# declare the project root
here::i_am("comboAMP_exponentialModel.R")

library(tidyverse)
source(here::here("utils.R"))

# Input parameters --------------------------------------------------------

rate_pla <- 0.032
nullHR <- 0.7
rate_cens <- 0.075

iter <- 2000
# correlates expansion phase
n_to_enroll <- 5000
n_target_cases_ab <- 35


# null hypothesis
nullHR <- 0.7

# alternative hypothesis
altHRs <- c(0.1, 0.15, 0.25, 0.35)

ans <- list()
for(altHR in altHRs){
  power_wald_v <- NULL
  power_cuminc_v <- NULL

  if(altHR %in% c(0.1, 0.15)){
    n_target_cases_v <- seq(15, 35, 1)
    tau <- 1
  }else if (altHR == 0.25){
    n_target_cases_v <- seq(45, 60, 1)
    tau <- 1.5
  }else if (altHR == 0.35){
    n_target_cases_v <- seq(90, 110, 1)
    tau <- 1.5
  }
  for(j in 1:length(n_target_cases_v)){
    df <- oper_chars_eff_phase(n_target_cases = n_target_cases_v[j],
                               rate_pla = rate_pla,
                               nullHR = nullHR,
                               altHR = altHR,
                               rate_cens = rate_cens,
                               tau = tau,
                               iter = iter)
    power_wald_v[j] <- mean(df$wald_pval <= 0.025)
    power_cuminc_v[j]  <- mean(df$cuminc_pval <= 0.025, na.rm = TRUE) 
  }
  ans[[paste(altHR)]] <- data.frame(power_wald = power_wald_v, power_cuminc = power_cuminc_v)
}

saveRDS(ans, here::here("power_expoentialModel_iterated.rds"))

