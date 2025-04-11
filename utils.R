library(survival)
library(DescTools)
source(here::here("cuminc_functions.R"))

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

# 2-arm designs
#' @param n_enroll_4m an integer value specifying the enrollment rate as the
#'   number of participants enrolled per 4 months
oper_chars_eff_phase <- function(n_target_cases, rate_pla, nullHR, altHR, 
                                 rate_cens, p_ab = 0.5, p_pla = 0.5, 
                                 n_enroll_4m = 1000, tau, iter, 
                                 minAtRisk = 150, size = 0.10){
  
  # total sample size
  n <- N(n_target_cases, p1 = p_ab, p0 = p_pla, rate1 = rate_pla * altHR, 
         rate0 = rate_pla, rateC = rate_cens, tau = tau)
  
  # sample size in each arm
  n_pla <- ceiling(n * p_pla)
  n_ab <- ceiling(n * p_ab)
  cat("Total sample size for the 2 arms included in the comparison:", n_pla + n_ab, "\n")
  
  # enrollment rate: 'n_enroll_4m' participants / 4 months
  enrollPeriod <- (n_pla + n_ab) * (4 / 12) / n_enroll_4m
  
  n1 <- n1(n_target_cases, hr = altHR, p1 = p_ab, p0 = p_pla)
  #cat("Expected number of events in the Ab arm (version 1):", n1, "\n")
  
  rate1 <- rate_pla * altHR
  r1 <- rate1 / (rate1 + rate_cens)
  pEvent1 <- r1 - r1 * exp(-(rate1 + rate_cens) * tau)
  cat("Expected number of events in the Ab arm (version 2):", n_ab * pEvent1, "\n")

  df <- plyr::ldply(1:iter, function(i){
    set.seed(i)
    
    enrollTime <- runif(n_pla + n_ab, max = enrollPeriod)
    FPFI <- min(enrollTime)
    tx <- rep(0:1, c(n_pla, n_ab))
    tm <- c(rexp(n_pla, rate_pla), rexp(n_ab, rate_pla * altHR))
    cens <- rexp(n_pla + n_ab, rate_cens)
    eventTime <- pmin(tm, cens)
    eventInd <- as.numeric(tm <= cens)
    calTime <- enrollTime + eventTime
    analysisTime <- sort(calTime[eventInd == 1])[n_target_cases]
    eventInd <- ifelse(calTime > analysisTime, 0, eventInd)
    calTime <- pmin(calTime, analysisTime)
    eventTime <- pmax(calTime - enrollTime, 0)
    
    not_enrolled <- sum(eventTime == 0)
    tx <- tx[eventTime > 0]
    eventInd <- eventInd[eventTime > 0]
    eventTime <- eventTime[eventTime > 0]
    
    split <- as.numeric(tapply(eventInd, tx, sum))
    
    # binomial score test 
    df <- table(tx, eventInd)
    rownames(df) <- c("placebo", "Ab")
    colnames(df) <- c("nonEvent", "event")
    df2 <- df[c(2, 1), c(2, 1)]
    CIscore <- RelRisk(df2, method = "score", conf.level = 1 - size)
    score_pval <- ifelse(CIscore["upr.ci"] < nullHR, 0.001, 1)
    
    if(split[2] == 0){
      wald_pval <- score_pval
    } else {
      # 1-sided Wald test
      sfit <- summary(suppressWarnings(coxph(Surv(eventTime, eventInd) ~ tx)))
      #warning is given when the number of cases is zero for the treatment group
      stat <- (sfit$coef[1, 1] - log(nullHR)) / sfit$coef[1, 3]
      wald_pval <- pnorm(stat)
    }
    
    # 1-sided log-rank test: doesn't have correct size for testing
    # H0: HR >= nullHR, where nullHR < 1; it does have correct size for nullHR = 1
    # LR <- survdiff(Surv(eventTime, eventInd) ~ tx)
    # HR <- (LR$obs[2] / LR$exp[2]) / (LR$obs[1] / LR$exp[1])
    # pval <- pchisq(LR$chisq, length(LR$n)-1, lower.tail = FALSE) / 2
    # logrank_pval <- ifelse(HR < nullHR, pval, 1 - pval)
    
    # # 1-sided Wald test based on the Nelson-Aalen estimator for the cumulative
    # # incidence
    # data <- data.frame("eventTime" = eventTime, "eventInd" = eventInd, "tx" = tx)
    # cuminc_est <- naCumInc( data = data, futimeVar = "eventTime", eventVar = "eventInd", 
    #                         groupVar = "tx", censor = list( minAtRisk = minAtRisk, time = max(eventTime[eventInd == 1])))
    # #    If *both* 'minAtRisk' and 'time'=t are specified, the time that satisfies
    # #    the 'minAtRisk' criteria will be determined, and then the *smaller* of that time and
    # #    time 't' will be used for censoring (and *that* value returned via 'censorTime'
    # cuminc_PE_est <- EffCIR(cuminc_est, refLvl = "0", cmpLvl="1", nullHypEff=1 - nullHR, test = "oneSided")
    # cuminc_pval <- cuminc_PE_est$tests$pvalue
    # # NA is given when the number of cases is less than or equal to 1 for the treatment group sometimes
    # if(is.na(cuminc_pval)){cuminc_pval <- score_pval}
    
    analysisTime <- analysisTime - FPFI
    
    return(data.frame(iter = i, n_enrolled = n_pla + n_ab,  n1 = n1, 
                      pEvent1 = pEvent1, analysisTime = analysisTime, 
                      n_cases_pla = split[1], n_cases_ab = split[2],
                      wald_pval = wald_pval, 
                      # cuminc_pval = cuminc_pval, 
                      not_enrolled = not_enrolled, 
                      meanEventTime = mean(eventTime),
                      Q1EventTime = quantile(eventTime, probs = 0.25),
                      medEventTime = quantile(eventTime, probs = 0.5),
                      Q3EventTime = quantile(eventTime, probs = 0.75),
                      maxEventTime = max(eventTime),
                      minEventTime = min(eventTime)))
  })
  
  return(df)
}

plot_time_to_analysis <- function(df, path){
  m <- mean(df$analysisTime)
  breaks <- c(pretty(df$analysisTime), m)
  labels <- c(pretty(df$analysisTime), round(m, 1))
  
  p <- ggplot(df, aes(x = analysisTime, y = ..density..)) +
    geom_histogram(fill = "cornsilk", color = "gray60") +
    geom_density() +
    geom_vline(xintercept = m, linetype = "dashed") +
    annotate("text", x = m, y = Inf, hjust = -0.07, vjust = 1.5, 
             label = "Mean", size = 2.8) +
    scale_x_continuous(breaks = breaks, labels = labels) +
    xlab("Time (Years) from FPFI to Primary Analysis") +
    ylab("Density") +
    theme_bw() +
    theme(panel.border = element_blank())
  
  ggsave(here::here(path, "time_to_analysis.pdf"), plot = p, 
         width = 4.5, height = 4.5)
  
  return(invisible(NULL))
}

plot_fu_time_eff_phase <- function(df, path){
  m <- mean(df$meanEventTime)
  breaks <- c(pretty(df$meanEventTime), m)
  labels <- c(pretty(df$meanEventTime), "")
  
  p <- ggplot(df, aes(x = meanEventTime, y = ..density..)) +
    geom_histogram(fill = "cornsilk", color = "gray60") +
    geom_density() +
    geom_vline(xintercept = m, linetype = "dashed") +
    annotate("text", x = m, y = Inf, hjust = -0.07, vjust = 1.5, 
             label = "Mean", size = 2.8) +
    scale_x_continuous(breaks = breaks, labels = labels) +
    xlab("Average Participant Follow-up Duration (Years)\nin Efficacy Phase") +
    ylab("Density") +
    theme_bw() +
    theme(panel.border = element_blank())
  
  ggsave(here::here(path, "fu_time_eff_phase.pdf"), plot = p, 
         width = 4.5, height = 4.5)
  
  return(invisible(NULL))
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
  
  #cat("The number on study at the start of cross-over:", n_on_study, "\n")
  
  # enrollment rate: 1000 participants / 4 months
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
    not_enrolled <- sum(eventTime == 0)
  
    newlyenrolledInd <- c(rep(0, n_on_study), rep(1, n_to_enroll))
    return(data.frame(iter = i, analysisTime = analysisTime, 
                      meanEventTime = mean(eventTime[newlyenrolledInd == 1]),
                      Q1EventTime = quantile(eventTime[newlyenrolledInd == 1], probs = 0.25),
                      Q3EventTime = quantile(eventTime[newlyenrolledInd == 1], probs = 0.75),
                      maxEventTime = max(eventTime[newlyenrolledInd == 1]),
                      minEventTime = min(eventTime[newlyenrolledInd == 1]),
                      not_enrolled = not_enrolled))
  })
  
  return(df)
}


oper_chars_eff_phase_dosesCalculation_twoArms <- function(n_target_cases, rate_pla, nullHR, altHR, 
                                 rate_cens, p_ab = 0.5, p_pla = 0.5, n_enroll_4m = 1000, tau, iter){
  
  # total sample size
  n <- N(n_target_cases, p1 = p_ab, p0 = p_pla, rate1 = rate_pla * altHR, 
         rate0 = rate_pla, rateC = rate_cens, tau = tau)
  
  # sample size in each arm
  n_pla <- ceiling(n * p_pla)
  n_ab <- ceiling(n * p_ab)
  
  # enrollment rate: 'n_enroll_4m' participants / 4 months
  enrollPeriod <- (n_ab + n_pla) * (4 / 12) / n_enroll_4m
  
  n1 <- n1(n_target_cases, hr = altHR, p1 = p_ab, p0 = p_pla)
  
  rate1 <- rate_pla * altHR
  r1 <- rate1 / (rate1 + rate_cens)
  pEvent1 <- r1 - r1 * exp(-(rate1 + rate_cens) * tau)
  
  df <- plyr::ldply(1:iter, function(i){
    set.seed(i)
    
    enrollTime <- runif(n_pla + n_ab, max = enrollPeriod)
    tx <- rep(0:1, c(n_pla, n_ab))
    tm <- c(rexp(n_pla, rate_pla), rexp(n_ab, rate_pla * altHR))
    cens <- rexp(n_pla + n_ab, rate_cens)
    eventTime <- pmin(tm, cens)
    eventInd <- as.numeric(tm <= cens)
    calTime <- enrollTime + eventTime
    analysisTime <- sort(calTime[eventInd == 1])[n_target_cases]
    eventInd <- ifelse(calTime > analysisTime, 0, eventInd)
    calTime <- pmin(calTime, analysisTime)
    eventTime <- pmax(calTime - enrollTime, 0)
    
    not_enrolled <- sum(eventTime == 0)
    tx <- tx[eventTime > 0]
    eventInd <- eventInd[eventTime > 0]
    eventTime <- eventTime[eventTime > 0]
    
    split <- as.numeric(tapply(eventInd, tx, sum))
    
    # calculating the total number of doses given; focus on Ab arm; doses are given every 6 months, including enrollment
    m <- ceiling(eventTime[tx == 1] / 0.5)

    return(data.frame(iter = i, n_pla = n_pla, n_ab = n_ab, n1 = n1, numOfDoses = sum(m)))
  })
  
  return(df)
}

