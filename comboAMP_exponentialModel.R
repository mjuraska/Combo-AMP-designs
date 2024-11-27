# declare the project root
here::i_am("comboAMP_exponentialModel.R")

library(survival)
library(tidyverse)
library(DescTools)
source(here::here("functions_eventSplits_sampleSize.R"))
source(here::here("cuminc_functions.R"))

# event target
eventTarget <- 20

# null hypothesis
nullHR <- 0.7

# alternative hypothesis
altHR <- 0.10

# total sample size
n <- N(eventTarget, p1 = 0.5, p0 = 0.5, rate1 = 0.032 * altHR, rate0 = 0.032, 
       rateC = 0.075, tau = 1)
n
# sample size per arm under 1:1 allocation
n <- ceiling(n / 2)

n1(eventTarget, hr = altHR, p1 = 0.5, p0 = 0.5)

# number of iterations
iter <- 2000

# enrollment rate: 1000 participants/4 months
enrollPeriod <- 2 * n * (4 / 12) / 1000

# upper bound of 95% CI for the HR
out <- lapply(1:iter, function(i){
  set.seed(i)
  
  enrollTime <- runif(2 * n, max = enrollPeriod)
  tx <- rep(0:1, each = n)
  tm <- c(rexp(n, 0.032), rexp(n, 0.032 * altHR))
  cens <- rexp(2 * n, 0.075)
  eventTime <- pmin(tm, cens)
  eventInd <- as.numeric(tm <= cens)
  calTime <- enrollTime + eventTime
  analysisTime <- sort(calTime[eventInd == 1])[eventTarget]
  eventInd <- ifelse(calTime > analysisTime, 0, eventInd)
  calTime <- pmin(calTime, analysisTime)
  eventTime <- pmax(calTime - enrollTime, 0)
  
  notEnrolled <- sum(eventTime == 0)
  tx <- tx[eventTime > 0]
  eventInd <- eventInd[eventTime > 0]
  eventTime <- eventTime[eventTime > 0]
  
  split <- as.numeric(tapply(eventInd, tx, sum))
  
  # check power
  #Binomial score Test when the number of cases is very low for the treatment group
  if(split[2] <=1){
    df <- table(tx, eventInd)
    rownames(df) <- c("placebo","vaccine")
    colnames(df) <- c("nonEvent", "Event")
    df2 <- df[c(2,1),c(2,1)]
    CIscore <- RelRisk(df2, method = "score", conf.level = 0.95)
    score_pval <- ifelse(CIscore["upr.ci"] < 0.3, 0.001, 1)
  }
  
  # 1-sided Wald test
  sfit <- summary(coxph(Surv(eventTime, eventInd) ~ tx)) 
  #warning is given when the number of cases is zero for the treatment group
  stat <- (sfit$coef[1, 1] - log(nullHR)) / sfit$coef[1, 3]
  wald_pval <- pnorm(stat)
  if(split[2] == 0){wald_pval <- score_pval}
  # 1-sided log-rank test: doesn't have correct size for testing
  # H0: HR >= nullHR, where nullHR < 1; it does have correct size for nullHR = 1
  # LR <- survdiff(Surv(eventTime, eventInd) ~ tx)
  # HR <- (LR$obs[2] / LR$exp[2]) / (LR$obs[1] / LR$exp[1])
  # pval <- pchisq(LR$chisq, length(LR$n)-1, lower.tail = FALSE) / 2
  # logrank_pval <- ifelse(HR < nullHR, pval, 1 - pval)
  
  # 1-sided wald test based on nelson-aalen estimator of the cumulative incidence
  data <- data.frame("eventTime" = eventTime, "eventInd" = eventInd, "tx" = tx)
  cuminc_est <- naCumInc( data = data, futimeVar = "eventTime", eventVar = "eventInd", 
                          groupVar = "tx", censor = list( minAtRisk=150))
  
  cuminc_PE_est <- EffCIR(cuminc_est, refLvl = "0", cmpLvl="1", nullHypEff=0.3, test = "oneSided")
  cuminc_pval <- cuminc_PE_est$tests$pvalue
  #NA is given when the number of cases is less than or equal to 1 for the treatment group sometimes
  if(split[2] <= 1 & is.na(cuminc_pval)){cuminc_pval <- score_pval}
  
  return(list(wald_pval = wald_pval, cuminc_pval = cuminc_pval, 
              #CIR = cuminc_PE_est$CIR, df = data,
              pla_events = split[1], vax_events = split[2],
              analysisTime = analysisTime, notEnrolled = notEnrolled))
})

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
