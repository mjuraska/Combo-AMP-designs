library(tidyverse)
library(survival)
library(foreach)
library(iterators)
library(parallel)
library(doParallel)
library(doRNG)
library(sievePH)
library(patchwork)
library(data.table)

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

perform_test <- function(df, nullHR, alpha_1sided){
  tab <- with(df, table(tx, eventInd))
  
  if(tab[2, 2] == 0){
    # convert 'tab' to numeric before passing it to RelRisk() to avoid integer
    # overflow
    tab <- matrix(as.numeric(tab[c(2, 1), c(2, 1)]), nrow = 2)
    CIscore <- DescTools::RelRisk(tab, method = "score", 
                                  conf.level = 1 - 2 * alpha_1sided)
    pval <- ifelse(CIscore["upr.ci"] < nullHR, 1e-9, 1)
  } else {
    sfit <- summary(coxph(Surv(eventTime, eventInd) ~ tx, data = df))
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
  
  pval <- perform_test(df %>% 
                         filter(prim_comp == 1) %>% 
                         mutate(tx = as.numeric(tx > 0)), 
                       nullHR = nullHR, 
                       alpha_1sided = alpha_1sided)
  
  reject_H0 <- as.numeric(pval <= alpha_1sided)
  
  return(list(df = df, at_case_count = at_case_count, anal_time = anal_time, 
              pval = pval, reject_H0 = reject_H0))
}

perform_stage2_analysis <- function(df, at_case_count, nullHR, alpha_1sided){
  anal_time <- sort(df %>% 
                      filter(eventInd == 1) %>% 
                      pull(calTime))[at_case_count]
  df <- df %>%
    mutate(eventInd = if_else(calTime > anal_time, 0, eventInd),
           calTime = pmin(calTime, anal_time),
           eventTime = pmax(calTime - enrollTime, 0)) %>%
    filter(eventTime > 0)
  
  pval <- perform_test(df %>% mutate(tx = tx_stage2 - 1), 
                       nullHR = nullHR, 
                       alpha_1sided = alpha_1sided)
  
  reject_H0 <- as.numeric(pval <= alpha_1sided)
  
  return(list(df = df, anal_time = anal_time, 
              pval = pval, reject_H0 = reject_H0))
}

get_cum_events <- function(tm, df){
  idx <- which(df$calTime <= tm)
  if (length(idx) == 0){
    return(0)
  } else {
    return(df$cumEvents[max(idx)])
  }
}

get_alpha_stage2 <- function(df, stage1_anal_time, stage2_target_h_l,
                             alpha_1sided){
  df <- df %>%
    arrange(calTime) %>%
    mutate(cumEvents = cumsum(eventInd))
  
  # Stage 2 analysis times in 6-month steps
  a_times <- seq(stage1_anal_time, max(df$calTime), by = 0.5)
  
  cum_events <- sapply(a_times, get_cum_events, df = df)
  
  cutoff <- which(cum_events >= stage2_target_h_l)[1]
  if (!is.na(cutff)){
    cum_events <- cum_events[1:cutoff]  
  }
  
  nom_alpha <- ldbounds::ldBounds(alpha = alpha_1sided, 
                                  t = cum_events / dplyr::last(cum_events),
                                  iuse = 1, sides = 1)$nom.alpha
  return(list(nom_alpha = nom_alpha, cum_events = cum_events))
}

#' @param compare Is the primary test for the comparison of high-dose mAb vs.
#'   placebo (\code{"h"}) or low-dose mAb vs. placebo (\code{"l"})?
#' @param n_target_cases number of cases, pooling over the mAb arm specified by
#'   \code{compare} and the placebo arm, triggering the final primary test
#' @param nullHR the null hypothesis HR for the comparison defined by \code{compare}
#' @param info_fractions function assumes a single interim analysis
run_stage1 <- function(compare = c("h", "l"), nullHR, altHR_h, 
                       altHR_l, method = "obf", 
                       info_fractions = c(0.5, 1), 
                       alpha_1sided = 0.025, 
                       n_target_cases, rate_pla, 
                       rate_cens, p_ab_h = 0.2, 
                       p_ab_l = 0.4, p_pla = 0.4, 
                       n_enroll_m = NULL, tau, iter, n_cores,
                       verbose = TRUE, seed = 3291){
  compare <- match.arg(compare)
  
  set.seed(seed)
  
  out <- foreach(i = 1:iter) %dorng% {
    source(here::here("code/utils.R"))
    
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
                        reject_H0 = a$reject_H0, 
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

#' @param compare Is the primary test for the comparison of high-dose mAb vs.
#'   placebo (\code{"h"}) or low-dose mAb vs. placebo (\code{"l"})?
#' @param n_target_cases number of cases, pooling over the mAb arm specified by
#'   \code{compare} and the placebo arm, triggering the final primary test
#' @param nullHR the null hypothesis HR for the comparison defined by \code{compare}
#' @param info_fractions function assumes a single interim analysis
run_stage1_stage2 <- function(compare = c("h", "l"), nullHR, altHR_h, 
                              altHR_l, nullHR_h_l, method = "obf", 
                              info_fractions = c(0.5, 1), 
                              alpha_1sided = 0.025, 
                              n_target_cases, n_stage2_target_cases_h_l, rate_pla, 
                              rate_cens, p_ab_h = 0.2, 
                              p_ab_l = 0.4, p_pla = 0.4, 
                              n_enroll_m = NULL, tau, iter, n_cores,
                              verbose = TRUE, seed = 3291){
  compare <- match.arg(compare)
  
  set.seed(seed)
  
  out <- foreach(i = 1:iter) %dorng% {
    source(here::here("code/utils.R"))
    
    p_ab_prim <- ifelse(compare == "h", p_ab_h, p_ab_l)
    altHR_prim <- ifelse(compare == "h", altHR_h, altHR_l)
    
    alpha1_seq <- ldbounds::ldBounds(alpha = alpha_1sided, t = info_fractions,
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
    
    # stage 1 interim analysis
    a1 <- perform_stage1_analysis(
      df, 
      at_case_count = ceiling(info_fractions[1] * n_target_cases),
      nullHR = nullHR,
      alpha_1sided = alpha1_seq[1]
    )
    
    if (a1$reject_H0 == 0){
      # proceed to stage 1 final analysis
      a1 <- perform_stage1_analysis(
        df, 
        at_case_count = n_target_cases,
        nullHR = nullHR,
        alpha_1sided = alpha1_seq[2]
      )
    }
    
    if (a1$reject_H0 == 1){
      # stage 2 is initiated; perform crossover
      df <- df %>%
        mutate(enrollTime = if_else(tx == 0, a1$anal_time, enrollTime),
               tx_stage2 = if_else(tx == 0, rbinom(n_total, 1, 0.5) + 1, tx),
               eventTime = if_else(tx == 0 & tx_stage2 == 1, rexp(n_total, rate = rate_pla * altHR_l), eventTime),
               eventTime = if_else(tx == 0 & tx_stage2 == 2, rexp(n_total, rate = rate_pla * altHR_h), eventTime),
               cens = rexp(n_total, rate = rate_cens),
               eventTime = if_else(tx == 0, pmin(eventTime, cens), eventTime),
               eventInd = if_else(tx == 0, as.numeric(eventTime <= cens), eventInd),
               calTime = if_else(tx == 0, enrollTime + eventTime, calTime))
      
      stage2_IA <- get_alpha_stage2(df, stage1_anal_time = a1$anal_time,
                                    stage2_target_h_l = n_stage2_target_cases_h_l,
                                    alpha_1sided = alpha_1sided)
      alpha2_seq <- stage2_IA$nom_alpha
      cum_events <- stage2_IA$cum_events
      
      # the first stage 2 analysis
      a2 <- perform_stage2_analysis(
        df, 
        at_case_count = cum_events[1],
        nullHR = nullHR_h_l,
        alpha_1sided = alpha2_seq[1]
      )
      
      k <- 2
      while (a2$reject_H0 == 0 & k <= length(cum_events)){
        a2 <- perform_stage2_analysis(
          df, 
          at_case_count = cum_events[k],
          nullHR = nullHR_h_l,
          alpha_1sided = alpha2_seq[k]
        )
        k <- k + 1
      }
      
      split <- as.numeric(with(a2$df, tapply(eventInd, tx_stage2, sum)))
      
      dat <- data.frame(iter = i, init_stage2 = 1, stage2_stop_time = a2$anal_time, 
                        stage2_pval = a2$pval, 
                        stage2_reject_H0 = a2$reject_H0, 
                        n_cases_ab_pla = sum(a1$df %>% filter(tx == 0) %>% pull(eventInd)),
                        n_cases_ab_l = split[1], 
                        n_cases_ab_h = split[2])
      
      if (verbose){
        return(list(dat = dat, msg = msg))
      } else {
        return(dat)  
      }
      
    } else {
      # stage 2 is not initiated and the trial is over
      # TO DO: finish the output
      split <- as.numeric(with(a1$df, tapply(eventInd, tx, sum)))
      
      dat <- data.frame(iter = i, init_stage2 = 0, stage2_stop_time = NA, 
                        stage2_pval = NA, 
                        stage2_reject_H0 = NA, 
                        n_cases_pla = split[1], 
                        n_cases_ab_l = split[2], 
                        n_cases_ab_h = split[3])
      
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

est_pe_by_log10ic80 <- function(n_total, n_enroll_m = NULL, p_pla, p_ab_l, p_ab_h, 
                                rate_pla, altHR_l, altHR_h, rate_cens,
                                n_target_cases_h_l, dens, beta, conf_level = 0.95, 
                                iter, seed = 3291){
  set.seed(seed)
  
  # dens$y <- dens$y[dens$x < 1]
  # dens$x <- dens$x[dens$x < 1]
  
  out <- foreach(i = 1:iter) %dorng% {
    source(here::here("code/utils.R"))
    
    # sample size in each arm
    n_ab_h <- round(n_total * p_ab_h, 0)
    n_ab_l <- round(n_total * p_ab_l, 0)
    n_pla <- round(n_total * p_pla, 0)
    n_total <- n_ab_h + n_ab_l + n_pla
    
    altHR <- (altHR_h * p_ab_h + altHR_l * p_ab_l) / (p_ab_h + p_ab_l)
    
    # enrollment rate: 'n_enroll_m' participants / month
    if (is.null(n_enroll_m)){
      enrollPeriod <- 1.5
    } else {
      enrollPeriod <- (n_total / n_enroll_m) / 12
    }
    
    # simulate trial data
    enrollTime <- runif(n_total, max = enrollPeriod)
    FPFI <- min(enrollTime)
    enrollTime <- enrollTime - FPFI
    tx <- rep(0:2, c(n_pla, n_ab_l, n_ab_h))
    tm <- c(rexp(n_pla, rate = rate_pla), 
            rexp(n_ab_l, rate = rate_pla * altHR_l),
            rexp(n_ab_h, rate = rate_pla * altHR_h))
    cens <- rexp(n_total, rate = rate_cens)
    eventTime <- pmin(tm, cens)
    eventInd <- as.numeric(tm <= cens)
    # sample from the support using the estimated density as weights
    log10_ic80_comb <- 
      c(sample(dens$x, size = n_pla, replace = TRUE, prob = dens$y),
        sample(dens$x, size = n_ab_l, replace = TRUE, 
               prob = dens$y * exp(beta * (dens$x - 1) - log(altHR_l))),
        sample(dens$x, size = n_ab_h, replace = TRUE, 
               prob = dens$y * exp(beta * (dens$x - 1) - log(altHR_h))))
    # log10_ic80_comb <- pmin(log10_ic80_comb, 1)
    log10_ic80_comb <- ifelse(eventInd == 1, log10_ic80_comb, NA)
    calTime <- enrollTime + eventTime
    df <- data.frame(enrollTime, tx, eventTime, eventInd, log10_ic80_comb, calTime)
    
    anal_times <- sort(df %>% 
                         filter(tx != 0, eventInd == 1) %>% 
                         pull(calTime))
    
    if (length(anal_times) < n_target_cases_h_l){
      # the target event count has not been reached
      return(NA)
    } else {
      anal_time <- anal_times[n_target_cases_h_l]
      
      # apply data cut at anal_time
      df <- df %>%
        mutate(eventInd = if_else(calTime > anal_time, 0, eventInd),
               calTime = pmin(calTime, anal_time),
               eventTime = pmax(calTime - enrollTime, 0),
               log10_ic80_comb = if_else(eventInd == 1, log10_ic80_comb, NA),
               tx = as.numeric(tx != 0)) %>%
        filter(eventTime > 0)
      
      fit <- sievePH(eventTime = df$eventTime, eventInd = df$eventInd,
                     mark = df$log10_ic80_comb, tx = df$tx)
      
      markRng <- range(df$log10_ic80_comb, na.rm = TRUE)
      markGrid <- seq(markRng[1], markRng[2], length.out = 500)
      # store a data frame with columns mark, TE, LB, UB
      sfit <- try(summary(fit, markGrid = markGrid, contrast = "te", 
                          sieveAlternative = "oneSided", confLevel = conf_level)$te)
      
      if (inherits(sfit, "try-error")){
        return("summary.sievePH-error")
      } else {
        sfit$iter <- i
        return(sfit)
      }
    }
  }
  
  return(out)
}

#' @param df_ic80 data frame of IC80s against 704 placebo viruses
#' @param h a numeric vector of bnAb-specific multiplicative constants by which
#'   IC80s are multiplied to define a new regimen
#' @param df_c data frame of individual-level concentrations
#' 
#' @return a data frame with a combination PT80 for each ptid, day, virus
get_comb_pt80 <- function(df_ic80, h = c(1, 1, 1), df_c){
  ab <- c("PGDM", "PGT", "VRC")
  
  # IC80s of a new regimen
  df_ic80 <- df_ic80 %>% 
    mutate(!!ab[1] := h[1] * .data[[ab[1]]],
           !!ab[2] := h[2] * .data[[ab[2]]],
           !!ab[3] := h[3] * .data[[ab[3]]])
  
  mat_ic80 <- as.matrix(df_ic80 %>% 
                          arrange(isolate) %>% 
                          select(all_of(ab)))
  
  mat_c <- as.matrix(df_c %>% 
                       arrange(id, time) %>% 
                       select(all_of(paste0("conc_", ab))))
  
  k <- 1
  out <- vector("list", length(unique(df_c$id)))
  for (i in unique(df_c$id)){
    idx <- which(df_c$id == i)
    comb_pt80 <- (1 / mat_ic80) %*% t(mat_c[idx, ])
    out[[k]] <- data.frame(id = rep(df_c$id[idx], each = nrow(df_ic80)),
                           time = rep(df_c$time[idx], each = nrow(df_ic80)),
                           virus = rep(df_ic80$isolate, times = length(idx)),
                           comb_pt80 = as.vector(comb_pt80))
    k <- k + 1
  }
  
  return(bind_rows(out))
}

#' @param df_pe data frame outputted by est_pe_by_log10ic80()
#' @param scale_c average effective concentration-at-exposure scaling constant
#' @param df_comb_pt80 data frame outputted by get_comb_pt80()
#' 
#' @return a numeric vector with the predicted PE point estimate and CI
predict_pe <- function(df_pe, scale_c, df_comb_pt80, approxfun_rule = 1){
  df_pe <- df_pe %>% 
    mutate(pt80 = scale_c / 10^mark) %>%
    # approxfun() below requires sorted pt80
    arrange(pt80)
  
  # create interpolation functions once
  f_pe <- approxfun(x = df_pe$pt80, y = df_pe$TE, rule = approxfun_rule)
  f_lb <- approxfun(x = df_pe$pt80, y = df_pe$LB, rule = approxfun_rule)
  
  df_pred_pe <- df_comb_pt80 %>%
    mutate(pe = f_pe(comb_pt80),
           lb = f_lb(comb_pt80)) %>%
    group_by(id, time) %>%
    summarise(mean_pe = mean(pe),
              mean_lb = mean(lb),
              .groups = "drop") %>%
    group_by(time) %>%
    summarise(med_pe = as.numeric(quantile(mean_pe, prob = 0.5)),
              lb_A = as.numeric(quantile(mean_pe, prob = 0.025)),
              lb_B = as.numeric(quantile(mean_lb, prob = 0.5)),
              .groups = "drop")
  
  return(c(ptEst_pe = mean(df_pred_pe$med_pe), 
           lb_pe_A = mean(df_pred_pe$lb_A),
           lb_pe_B = mean(df_pred_pe$lb_B)))
}

# the fastest
predict_pe_dt <- function(df_pe, scale_c, df_comb_pt80, approxfun_rule = 1){
  setDT(df_pe)
  setDT(df_comb_pt80)
  
  # Compute pt80 and sort (approxfun needs sorted x)
  df_pe[, pt80 := scale_c / (10^mark)]
  setorder(df_pe, pt80)
  
  # create interpolation functions once
  f_pe <- approxfun(x = df_pe$pt80, y = df_pe$TE, rule = approxfun_rule)
  f_lb <- approxfun(x = df_pe$pt80, y = df_pe$LB, rule = approxfun_rule)
  f_ub <- approxfun(x = df_pe$pt80, y = df_pe$UB, rule = approxfun_rule)
  
  df_comb_pt80[, `:=`(logrr = log(1 - f_pe(comb_pt80)),
                      lb = log(1 - f_ub(comb_pt80)),
                      ub = log(1 - f_lb(comb_pt80)))]
  
  df_summary <- df_comb_pt80[, .(mean_pe = mean(1 - exp(logrr), na.rm = TRUE),
                                 mean_logrr = mean(logrr, na.rm = TRUE),
                                 mean_lb = mean(lb, na.rm = TRUE),
                                 mean_ub = mean(ub, na.rm = TRUE)), 
                             by = .(id, time)]
  
  df_pred <- df_summary[, .(med_pe = as.numeric(quantile(mean_pe, prob = 0.5, na.rm = TRUE)),
                            med_logrr = as.numeric(quantile(mean_logrr, prob = 0.5, na.rm = TRUE)),
                            lb_pe_A = as.numeric(quantile(mean_pe, prob = 0.025, na.rm = TRUE)),
                            ub_pe_A = as.numeric(quantile(mean_pe, prob = 0.975, na.rm = TRUE)),
                            lb_logrr_cA = as.numeric(quantile(mean_logrr, prob = 0.025, na.rm = TRUE)),
                            ub_logrr_cA = as.numeric(quantile(mean_logrr, prob = 0.975, na.rm = TRUE)),
                            lb_logrr_B = as.numeric(quantile(mean_lb, prob = 0.5, na.rm = TRUE)),
                            ub_logrr_B = as.numeric(quantile(mean_ub, prob = 0.5, na.rm = TRUE))), 
                        by = time]
  
  return(c(ptEst_pe_A = mean(df_pred$med_pe, na.rm = TRUE),
           ptEst_pe_cAB = 1 - exp(mean(df_pred$med_logrr, na.rm = TRUE)),
           lb_pe_A  = mean(df_pred$lb_pe_A, na.rm = TRUE),
           ub_pe_A  = mean(df_pred$ub_pe_A, na.rm = TRUE),
           lb_pe_cA  = 1 - exp(mean(df_pred$ub_logrr_cA, na.rm = TRUE)),
           ub_pe_cA  = 1 - exp(mean(df_pred$lb_logrr_cA, na.rm = TRUE)),
           lb_pe_B  = 1 - exp(mean(df_pred$ub_logrr_B, na.rm = TRUE)),
           ub_pe_B  = 1 - exp(mean(df_pred$lb_logrr_B, na.rm = TRUE))))
}

predict_pe_parallel <- function(df_pe, scale_c, df_comb_pt80, 
                                approxfun_rule = 1, n_cores){
  df_pe <- df_pe %>% 
    mutate(pt80 = scale_c / 10^mark) %>%
    # approxfun() below requires sorted pt80
    arrange(pt80)
  
  # create interpolation functions once
  f_pe <- approxfun(x = df_pe$pt80, y = df_pe$TE, rule = approxfun_rule)
  f_lb <- approxfun(x = df_pe$pt80, y = df_pe$LB, rule = approxfun_rule)
  
  # split df_comb_pt80 by time for parallel processing
  time_groups <- split(df_comb_pt80, df_comb_pt80$time)
  
  df_pred_pe <- foreach(time_val = names(time_groups), .combine = bind_rows) %dorng% {
    library(tidyverse)
    
    df_group <- time_groups[[time_val]]
    
    df_group <- df_group %>%
      mutate(pe = f_pe(comb_pt80),
             lb = f_lb(comb_pt80)) %>%
      group_by(id, time) %>%
      summarise(mean_pe = mean(pe, na.rm = TRUE),
                mean_lb = mean(lb, na.rm = TRUE),
                .groups = "drop")
    
    return(summarise(df_group,
                     time = unique(time),
                     med_pe = as.numeric(quantile(mean_pe, prob = 0.5, na.rm = TRUE)),
                     lb_A = as.numeric(quantile(mean_pe, prob = 0.025, na.rm = TRUE)),
                     lb_B = as.numeric(quantile(mean_lb, prob = 0.5, na.rm = TRUE))))
  }
  
  return(c(ptEst_pe = mean(df_pred_pe$med_pe, na.rm = TRUE),
           lb_pe_A = mean(df_pred_pe$lb_A, na.rm = TRUE),
           lb_pe_B = mean(df_pred_pe$lb_B, na.rm = TRUE)))
}

# slower than predict_pe_parallel()
predict_pe_parallel_dt <- function(df_pe, scale_c, df_comb_pt80, 
                                   approxfun_rule = 1, n_cores){
  # convert to data.table
  setDT(df_pe)
  setDT(df_comb_pt80)
  
  df_pe[, pt80 := scale_c / (10^mark)]
  setorder(df_pe, pt80)
  
  # create interpolation functions once
  f_pe <- approxfun(x = df_pe$pt80, y = df_pe$TE, rule = approxfun_rule)
  f_lb <- approxfun(x = df_pe$pt80, y = df_pe$LB, rule = approxfun_rule)
  
  # split df_comb_pt80 by time for parallel processing
  time_groups <- split(df_comb_pt80, by = "time", keep.by = TRUE)
  
  df_pred_pe <- foreach(time_val = names(time_groups), .combine = rbind, 
                        .packages = "data.table") %dorng% {
    df_group <- time_groups[[time_val]]
    
    df_group[, pe := f_pe(comb_pt80)]
    df_group[, lb := f_lb(comb_pt80)]
    
    df_group <- df_group[, .(mean_pe = mean(pe, na.rm = TRUE), 
                             mean_lb = mean(lb, na.rm = TRUE)), 
                         by = .(id, time)]
    
    return(
      data.table(time = unique(df_group$time),
                 med_pe = as.numeric(quantile(df_group$mean_pe, prob = 0.5, na.rm = TRUE)),
                 lb_A   = as.numeric(quantile(df_group$mean_pe, prob = 0.025, na.rm = TRUE)),
                 lb_B   = as.numeric(quantile(df_group$mean_lb, prob = 0.5, na.rm = TRUE)))
      )
  }
  
  return(c(ptEst_pe = mean(df_pred_pe[, med_pe], na.rm = TRUE),
           lb_pe_A = mean(df_pred_pe[, lb_A], na.rm = TRUE),
           lb_pe_B = mean(df_pred_pe[, lb_B], na.rm = TRUE)))
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

