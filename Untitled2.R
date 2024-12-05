# declare the project root
here::i_am("comboAMP_exponentialModel.R")

library(tidyverse)
source(here::here("utils.R"))


# Input parameters --------------------------------------------------------

n_target_cases_v <- c(25)
ans <- list()
for(i in 1:length(n_target_cases_v)){
  n_target_cases <- n_target_cases_v[i]
  rate_pla <- 0.026
  nullHR <- 0.70
  altHR <- 0.15
  rate_cens <- 0.075
  tau <- 1
  iter <- 2000
  # correlates expansion phase
  n_to_enroll <- 5000
  n_target_cases_ab <- 35
  # 2-arm design
  altHR_high <- altHR
  # 3-arm design
  #altHR_high <- 0.15
  
  ans[[i]]<- simul(n_target_cases = n_target_cases,
                   rate_pla = rate_pla,
                   nullHR = nullHR,
                   altHR = altHR,
                   altHRhigh = altHRhigh,
                   rate_cens = rate_cens,
                   tau = tau,
                   iter = iter,
                   n_target_cases_ab = n_target_cases_ab)
  
}

ans[[1]]$power
ans[[2]]$power
ans[[3]]$power