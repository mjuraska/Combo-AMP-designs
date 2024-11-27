# compute total events for adequate power under a given design alternative
# (Schoenfeld, 1983, Biometrics)
((qnorm(0.975) + qnorm(0.9))^2) / ((1/4) * (log(0.7) - log(0.25))^2)

# get number of events in the numerator given
# total events 'n', hazard ratio 'hr', and probabilities 'p1' and 'p0' of
# assignment to Z=1 and Z=0, respectively;
# this is useful for obtaining event splits by arm
# for a given total number of events
n1 <- function(n, hr, p1=0.5, p0=0.5){
  p <- p1 / p0
  return(hr * p * n / (1 + p * hr))
}

# compute total target sample size 'N' for a 2-arm trial as a function of
# total events 'n', probability 'p1' of being assigned treatment Z=1,
# probability 'p0' of being assigned treatment Z=0, 
# and exponential rates (per PYR) for T|Z=1 ('rate1'),
# T|Z=0 ('rate0') and censoring time C ('rateC');
# 'tau' is the duration of follow-up (in years)
N <- function(n, p1=0.5, p0=0.5, rate1, rate0, rateC, tau=72/52){
  r1 <- rate1 / (rate1 + rateC)
  pEvent1 <- r1 - r1 * exp(-(rate1 + rateC) * tau)
  r0 <- rate0 / (rate0 + rateC)
  pEvent0 <- r0 - r0 * exp(-(rate0 + rateC) * tau)
  return(n / (p1 * pEvent1 + p0 * pEvent0))
}

#' @param n_on_study number of originally enrolled participants at risk at the
#'   time of cross-over
#' @param n_to_enroll number of participants to be additionally enrolled to
#'   speed up time to accrual of \code{n_target_cases_ab} cases
#' @param n_obs_cases_ab number of cases already observed in the Ab arm at the
#'   start of cross-over
#' @param n_target_cases_ab target number of cases in the Ab arm
#' @param rate_pla exponential rate (per PYR) for T|Z=placebo
#' @param altHR design alternative HR
#' @param rate_cens exponential rate (per PYR) for C
duration_corr_exp_phase <- function(n_on_study, n_to_enroll, 
                                    n_obs_cases_ab, n_target_cases_ab, 
                                    rate_pla, altHR, rate_cens, iter){
  # enrollment rate: 1000 participants/4 months
  enrollPeriod <- n_to_enroll * (4 / 12) / 1000
  
  rate_ab <- rate_pla * altHR
  
  df <- plyr::ldply(1:iter, function(i){
    set.seed(i)
    
    enrollTime <- c(rep(0, n_on_study), runif(n_to_enroll, max = enrollPeriod))
    tm <- rexp(n_on_study + n_to_enroll, rate_ab)
    cens <- rexp(n_on_study + n_to_enroll, rate_cens)
    eventTime <- pmin(tm, cens)
    eventInd <- as.numeric(tm <= cens)
    calTime <- enrollTime + eventTime
    analysisTime <- sort(calTime[eventInd == 1])[n_target_cases_ab]
    eventInd <- ifelse(calTime > analysisTime, 0, eventInd)
    calTime <- pmin(calTime, analysisTime)
    eventTime <- pmax(calTime - enrollTime, 0)
    notEnrolled <- sum(eventTime == 0)
    
    return(data.frame(iter = i, analysisTime = analysisTime, 
                      notEnrolled = notEnrolled))
  })
  
  return(df)
}