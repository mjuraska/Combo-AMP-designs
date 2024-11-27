# declare the project root
here::i_am("comboAMP_exponentialModel.R")

library(tidyverse)
source(here::here("utils.R"))

df <- oper_chars_eff_phase(n_target_cases = 12,
                           rate_pla = 0.032,
                           nullHR = 0.7,
                           altHR = 0.1,
                           rate_cens = 0.075,
                           tau = 1,
                           iter = 1000)


# power
mean(sapply(out, "[[", "wald_pval") <= 0.025)
# mean(sapply(out, "[[", "logrank_pval") <= 0.025)

mean(sapply(out, "[[", "cuminc_pval") <= 0.025, na.rm = TRUE) #some tests have NA since the treatment group has <=1 case
cuminc_pval <- sapply(out, "[[", "cuminc_pval")
# check the time when the target number of events is accrued
mean(sapply(out, "[[", "analysisTime"))

summary(sapply(out, "[[", "notEnrolled"))


# splits
# df <- data.frame(pla = sapply(out, "[[", "pla_events"),
#                  vax = sapply(out, "[[", "vax_events")) %>%
#   group_by(pla, vax) %>%
#   summarise(prob = n() / iter)
# df
