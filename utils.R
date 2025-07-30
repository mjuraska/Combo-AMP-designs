library(tidyverse)
library(survival)
library(foreach)
library(iterators)
library(parallel)
library(doParallel)
library(doRNG)
# source(here::here("cuminc_functions.R"))

# compute total events for adequate power under a given design alternative
# (Schoenfeld, 1983, Biometrics)


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

primary_test <- function(df, nullHR, alpha_1sided){
  tab <- with(df, table(tx, eventInd))
  
  if(tab[2, 2] == 0){
    # convert 'tab' to numeric before passing it to RelRisk() to avoid integer
    # overflow
    tab <- matrix(as.numeric(tab[c(2, 1), c(2, 1)]), nrow = 2)
    CIscore <- DescTools::RelRisk(tab, method = "score", 
                                  conf.level = 1 - 2 * alpha_1sided)
    pval <- ifelse(CIscore["upr.ci"] < nullHR, 1e-9, 1)
  } else {
    sfit <- summary(coxph(Surv(eventTime, eventInd) ~ tx, 
                          data = mutate(df, tx = as.numeric(tx > 0))))
    stat <- (sfit$coef[1, 1] - log(nullHR)) / sfit$coef[1, 3]
    pval <- pnorm(stat)
  }
  
  return(pval)
}

perform_stage1_analysis <- function(df, at_case_count, nullHR, alpha_1sided){
  anal_time <- sort(df %>% 
                      filter(eventInd == 1, prim_comp == 1) %>% 
                      pull(calTime))[at_case_count]
  df <- df %>%
    mutate(eventInd = if_else(calTime > anal_time, 0, eventInd),
           calTime = pmin(calTime, anal_time),
           eventTime = pmax(calTime - enrollTime, 0)) %>%
    filter(eventTime > 0)
  
  pval <- primary_test(df %>% filter(prim_comp == 1), nullHR = nullHR, 
                       alpha_1sided = alpha_1sided)
  
  return(list(df = df, anal_time = anal_time, pval = pval))
}

#' @param compare Is the primary test for the comparison of high-dose mAb vs.
#'   placebo (\code{"h"}) or low-dose mAb vs. placebo (\code{"l"})?
#' @param n_target_cases number of cases, pooling over the mAb arm specified by
#'   \code{compare} and the placebo arm, triggering the final primary test
#' @param nullHR the null hypothesis HR for the comparison defined by \code{compare}
#' @param info_fractions function assumes a single interim analysis
oper_chars_eff_phase <- function(compare = c("h", "l"), nullHR, altHR_h, 
                                 altHR_l, method = "obf", 
                                 info_fractions = c(0.5, 1), 
                                 alpha_1sided = 0.025, 
                                 n_target_cases, rate_pla, 
                                 rate_cens, p_ab_h = 0.2, 
                                 p_ab_l = 0.4, p_pla = 0.4, 
                                 n_enroll_m = NULL, tau, iter, n_cores,
                                 verbose = TRUE){
  compare <- match.arg(compare)
  
  out <- foreach(i = 1:iter) %dorng% {
    source(here::here("utils.R"))
    
    p_ab_prim <- ifelse(compare == "h", p_ab_h, p_ab_l)
    altHR_prim <- ifelse(compare == "h", altHR_h, altHR_l)
    
    alpha_seq <- ldbounds::ldBounds(alpha = alpha_1sided, t = info_fractions,
                                    iuse = 1, sides = 1)$nom.alpha
    
    # 1 Ab arm + placebo sample size 
    n_2arm <- N(n_target_cases, p1 = p_ab_prim / (p_ab_prim + p_pla),
                p0 = p_pla / (p_ab_prim + p_pla), rate1 = rate_pla * altHR_prim,
                rate0 = rate_pla, rateC = rate_cens, tau = tau)
    n_total <- n_2arm / (p_ab_prim + p_pla)
    
    # sample size in each arm
    n_ab_h <- ceiling(n_total * p_ab_h)
    n_ab_l <- ceiling(n_total * p_ab_l)
    n_pla <- ceiling(n_total * p_pla)
    
    n_total <- n_ab_h + n_ab_l + n_pla
    
    # enrollment rate: 'n_enroll_m' participants / month
    if (is.null(n_enroll_m)){
      enrollPeriod <- 1.5
    } else {
      enrollPeriod <- (n_total / n_enroll_m) / 12
    }
    
    # get expected number of cases in each of high and low-dose mAb arm
    p_case_ab <- sapply(c(altHR_h, altHR_l), function(altHR){
      rate <- rate_pla * altHR
      r <- rate / (rate + rate_cens)
      return(r - r * exp(-(rate + rate_cens) * tau))
    })
    e_n_cases_ab <- round(c(n_ab_h, n_ab_l) * p_case_ab, digits = 0)
    e_n_cases_pla <- n_target_cases - ifelse(compare == "h", e_n_cases_ab[1], e_n_cases_ab[2])
    
    if (i == 1){
      msg <- paste0("Total sample size for the 2 arms included in the primary comparison:", 
                    ifelse(compare == "h", n_ab_h, n_ab_l) + n_pla, "\n",
                    "Total sample size for all 3 arms:", n_total, "\n",
                    "Expected number of events in the high-dose Ab arm:", e_n_cases_ab[1], "\n",
                    "Expected number of events in the low-dose Ab arm:", e_n_cases_ab[2], "\n",
                    "Expected number of events in the placebo arm:", e_n_cases_pla, "\n")  
    } else {
      msg <- NA
    }
    
    set.seed(i)
    
    enrollTime <- runif(n_total, max = enrollPeriod)
    FPFI <- min(enrollTime)
    enrollTime <- enrollTime - FPFI
    tx <- rep(0:2, c(n_pla, n_ab_l, n_ab_h))
    prim_comp <- as.numeric(tx != ifelse(compare == "h", 1, 2))
    tm <- c(rexp(n_pla, rate = rate_pla), 
            rexp(n_ab_l, rate = rate_pla * altHR_l),
            rexp(n_ab_h, rate = rate_pla * altHR_h))
    cens <- rexp(n_total, rate = rate_cens)
    eventTime <- pmin(tm, cens)
    eventInd <- as.numeric(tm <= cens)
    calTime <- enrollTime + eventTime
    df <- data.frame(enrollTime, tx, prim_comp, eventTime, eventInd, calTime)
    
    # interim analysis
    a <- perform_stage1_analysis(
      df, 
      at_case_count = ceiling(info_fractions[1] * n_target_cases),
      nullHR = nullHR,
      alpha_1sided = alpha_seq[1]
      )
    
    if (a$pval <= alpha_seq[1]){
      split <- as.numeric(with(a$df, tapply(eventInd, tx, sum)))
      dat <- data.frame(iter = i, stop_at_IA = 1, anal_time = a$anal_time, 
                        pval = a$pval, reject_H0 = 1, n_cases_pla = split[1], 
                        n_cases_ab_l = split[2], n_cases_ab_h = split[3], 
                        mean_eventTime = mean(a$df$eventTime))
      if (verbose){
        return(list(dat = dat, msg = msg))
      } else {
        return(dat)  
      }
      
    } else {
      # proceed to final analysis
      a <- perform_stage1_analysis(
        df, 
        at_case_count = n_target_cases,
        nullHR = nullHR,
        alpha_1sided = alpha_seq[2]
      )
      
      split <- as.numeric(with(a$df, tapply(eventInd, tx, sum)))
      
      dat <- data.frame(iter = i, stop_at_IA = 0, anal_time = a$anal_time, 
                        pval = a$pval, 
                        reject_H0 = as.numeric(a$pval <= alpha_seq[2]), 
                        n_cases_pla = split[1], 
                        n_cases_ab_l = split[2], 
                        n_cases_ab_h = split[3], 
                        mean_eventTime = mean(a$df$eventTime))
      if (verbose){
        return(list(dat = dat, msg = msg))
      } else {
        return(dat)  
      }
    }
  }
  
  if (verbose){
    cat(na.omit(sapply(out, "[[", "msg")))
    out <- do.call(rbind, lapply(out, "[[", "dat"))
  } else {
    out <- do.call(rbind, out)  
  }
  
  return(out)
}

plot_time_to_end_stage1 <- function(df, path){
  m <- mean(df$anal_time)
  breaks <- c(1.5, 2, 2.5, 3, 3.5, m)
  labels <- c(1.5, 2, 2.5, 3, 3.5, round(m, 1))
  
  p <- ggplot(df, aes(x = anal_time, y = ..density..)) +
    geom_histogram(fill = "cornsilk", color = "gray60") +
    geom_density() +
    geom_vline(xintercept = m, linetype = "dashed") +
    annotate("text", x = m, y = Inf, hjust = -0.07, vjust = 1.5, 
             label = "Mean", size = 2.8) +
    scale_x_continuous(breaks = breaks, labels = labels, minor_breaks = NULL) +
    xlab("Time (Years) from FPFI to End of Stage 1") +
    ylab("Density") +
    theme_bw() +
    theme(panel.border = element_blank())
  
  ggsave(here::here(path, "time_to_end_stage1.pdf"), plot = p, 
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

plot_h_vs_p_case_split_end_stage1 <- function(df, path){
  df <- df %>%
    group_by(n_cases_ab_h, n_cases_pla) %>%
    summarise(p = n() / nrow(df)) %>%
    ungroup() %>%
    mutate(split = paste0(n_cases_ab_h, ":", n_cases_pla),
           split = factor(split, levels = split))
  
  p <- ggplot(df, aes(x = split, y = p)) +
    geom_col() +
    xlab("High-Dose Ab:Placebo Endpoint Split at End of Stage 1") +
    ylab("Probability") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
          panel.border = element_blank())
  
  ggsave(here::here(path, paste0("h_vs_p_case_split_end_stage1.pdf")), 
         plot = p, width = 6, height = 4.5)
  
  return(invisible(NULL))
}

plot_l_case_count_end_stage1 <- function(df, path){
  df <- df %>%
    group_by(n_cases_ab_l) %>%
    summarise(p = n() / nrow(df)) %>%
    ungroup()
  
  p <- ggplot(df, aes(x = factor(n_cases_ab_l), y = p)) +
    geom_col() +
    xlab("Low-Dose Ab Arm Endpoint Count at End of Stage 1") +
    ylab("Probability") +
    theme_bw() +
    theme(panel.border = element_blank())
  
  ggsave(here::here(path, paste0("l_case_count_end_stage1.pdf")), 
         plot = p, width = 6, height = 4.5)
  
  return(invisible(NULL))
}

plot_n_doses <- function(df, var_name, path,
                         x_lab = "Number of Received Low Doses of Ab by Primary Analysis",
                         title = "Low-Dose Arm",
                         file_name){
  m <- mean(df[, var_name])
  # breaks <- c(pretty(df[, var_name]), m)
  # labels <- c(pretty(df[, var_name]), round(m, 1))
  
  p <- ggplot(df, aes(x = .data[[var_name]], y = ..density..)) +
    geom_histogram(fill = "cornsilk", color = "gray60") +
    geom_density() +
    geom_vline(xintercept = m, linetype = "dashed") +
    annotate("text", x = m, y = Inf, hjust = -0.07, vjust = 1.5, 
             label = paste0("Mean = ", round(m, 0)), size = 2.8) +
    scale_x_continuous(breaks = pretty(df[, var_name])) +
    xlab(x_lab) +
    ylab("Density") +
    ggtitle(title) +
    theme_bw() +
    theme(panel.border = element_blank())
  
  ggsave(here::here(path, file_name), plot = p, 
         width = 4.4, height = 4.5)
  
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

