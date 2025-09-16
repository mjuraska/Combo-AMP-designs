rm(list = ls(all = TRUE))

library(tidyverse)
library(grid)

source(here::here("code/param.R"))
source(here::here("code/utils.R"))


# Get potency weights of Combo-AMP bnAbs for 704 placebo viruses ----------

d_ic80 <- read.csv(here::here("data", ic80_file)) %>%
  select(-VRC01) %>%
  rename_with(~ case_when(grepl("PGDM", .x) ~ "PGDM",
                          grepl("PGT", .x) ~ "PGT",
                          grepl("VRC", .x) ~ "VRC",
                          TRUE ~ .x)) %>%
  mutate(comb_ic80 = 1 / ((1 / PGDM) + (1 / PGT) + (1 / VRC)),
         log10_comb_ic80 = log10(comb_ic80),
         w_PGDM = comb_ic80 / PGDM,
         w_PGT = comb_ic80 / PGT,
         w_VRC = comb_ic80 / VRC,
         w_sum = w_PGDM + w_PGT + w_VRC)


# Get scaling constant using pop-level concentrations ---------------------

d_pop_conc <- plyr::ldply(pop_conc_files, function(f){
  read.csv(here::here("data", f))
}) %>% 
  group_by(drug, dose) %>%
  summarise(med = quantile(Cc_ss, prob = 0.5),
            gm = exp(mean(log(Cc_ss))))

gm_pop_conc_l <- exp(mean(log(drop(as.matrix(d_ic80 %>% select(w_PGDM, w_PGT, w_VRC)) %*% 
                              (d_pop_conc %>% filter(dose == "IV 0.4g") %>% pull(gm))))))
gm_pop_conc_h <- exp(mean(log(drop(as.matrix(d_ic80 %>% select(w_PGDM, w_PGT, w_VRC)) %*% 
                              (d_pop_conc %>% filter(dose != "IV 0.4g") %>% pull(gm))))))
gm_pop_conc <- exp(mean(log(c(gm_pop_conc_h, gm_pop_conc_l))))


# Get scaling constant using ind-level concentrations ---------------------

d_ind_conc <- plyr::ldply(1:length(ind_conc_files), function(i){
  read.csv(here::here("data", ind_conc_files[i])) %>%
    rename_with(~ case_when(grepl("PGDM", .x) ~ "conc_PGDM",
                            grepl("PGT", .x) ~ "conc_PGT",
                            grepl("VRC", .x) ~ "conc_VRC",
                            TRUE ~ .x)) %>%
    mutate(dose = if_else(i == 1, "low", "high"))
}) %>%
  mutate(id = if_else(dose == "low", id, id + max(id[dose == "low"])))

gm_ind_conc <- d_ind_conc %>%
  group_by(id) %>%
  summarise(gm_PGDM = exp(mean(log(conc_PGDM))),
            gm_PGT = exp(mean(log(conc_PGT))),
            gm_VRC = exp(mean(log(conc_VRC))),
            .groups = "drop")

m_w <- as.matrix(d_ic80 %>% select(w_PGDM, w_PGT, w_VRC))

gm_ind_conc <- exp(mean(log(
  sapply(1:NROW(gm_ind_conc), function(r){
  conc <- as.numeric(gm_ind_conc[r, paste0("gm_", c("PGDM", "PGT", "VRC"))])
  return(mean(drop(m_w %*% conc)))
})
)))


# Plot potency-based weights of the 3 bnAbs -------------------------------

p <- ggtern::ggtern(d_ic80, aes(x = w_PGDM, y = w_PGT, z = w_VRC)) +
  geom_point(size = 3, alpha = 0.5) +
  scale_T_continuous(breaks = seq(0, 1, 0.1), labels = seq(0, 1, 0.1)) +
  scale_L_continuous(breaks = seq(0, 1, 0.1), labels = seq(0, 1, 0.1)) +
  scale_R_continuous(breaks = seq(0, 1, 0.1), labels = seq(0, 1, 0.1)) +
  theme(
    clip = "off"   # allow labels outside the panel
  ) +
  labs(T = expression(w[1]), L = expression(w[2]), R = expression(w[3])) +
  theme_bw()

ggsave(here::here(file.path(path, "potencyBasedWeights.pdf")), plot = p, height = 5, width = 5)


# Get density ratio model coefficients under PE(log10(IC80)=1) = 0 --------

# density estimate of log10(IC80_comb) against 704 placebo viruses
dens <- density(d_ic80$log10_comb_ic80, n = 1000)
plot(dens)

f_beta <- function(beta, dens, altHR){
  x <- dens$x 
  y <- dens$y * exp(beta * (dens$x - 1) - log(altHR))
  dx <- diff(x)
  area <- sum(dx * (y[-1] + y[-length(y)]) / 2)
  return(area - 1)
}

# overall altHR (high + low vs. placebo) under 1:1 randomization
altHR <- 0.5 * (altHR_h + altHR_l)

# unique root
beta <- uniroot(f_beta, interval = c(0.5, 2), dens = dens, altHR = altHR, tol = 1e-9)$root

v <- seq(-3, 1, length.out = 200)
ve <- 1 - exp(beta * (v - 1))
df_ve <- data.frame(v = v, ve = ve)

dens <- density(d_ic80$log10_comb_ic80, n = 300)
df_dens <- data.frame(tx = rep(0:1, each = length(dens$x)),
                      v = rep(dens$x, 2), 
                      dens = c(dens$y, 
                               dens$y * exp(beta * (dens$x - 1) - log(altHR)))) %>%
  filter(v >= -3, v <= 1)


# Plot PE by log10(IC80) --------------------------------------------------

p <- ggplot() +
  geom_hline(yintercept = 1 - altHR, linetype = "dashed") +
  annotate("text", x = min(df_ve$v), y = 0.8, hjust = 0.05, vjust = -0.6,
           label = "Overall PE (Ab high + Ab low vs. placebo)", size = 2.3) +
  geom_line(aes(x = v, y = dens, group = factor(tx), color = factor(tx)), 
            data = df_dens) +
  scale_color_manual(values = c("blue1", "magenta"), 
                     labels = c("Placebo", "Ab high + Ab low"),
                     name = "PDF") +
  geom_line(aes(x = v, y = ve), data = df_ve, linewidth = 1.2) +
  scale_x_continuous(breaks = -3:1, labels = 10^(-3:1)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.25), 
                     labels = seq(0, 1, by = 0.25) * 100) +
  labs(x = "Combination IC80 of Combo-AMP Regimen\nagainst Autologous Virus", 
       y = "Prevention Efficacy (%)") +
  theme_bw() + 
  theme(legend.position = "bottom")

ggsave(here::here(file.path(path, "truePEbyLog10combIC80_dataFromLily.pdf")), plot = p, height = 5, width = 4.9)


# Plot PE by PT80 ---------------------------------------------------------

# emailed by Lily on 8/13/25
# df_conc <- data.frame(ab = rep(c("PGDM1400LS", "PGT121LS", "VRC07.523LS"), each = 2),
#                       dose = rep(c("IV 0.4g", "IV 1.6g"), 3),
#                       med_conc = c(33.30190776, 133.207631, 48.21225924, 192.849037, 19.17685883, 153.4148706))
# d <- df_conc %>% filter(dose == "IV 0.4g") %>% select(ab, med_conc) %>%
#   mutate(gm_ic80 = c(exp(mean(log(dat$ic80_ab3))), exp(mean(log(dat$ic80_ab2))), exp(mean(log(dat$ic80_ab1)))),
#          pt80 = med_conc / gm_ic80)
# gm_ic80_comb <- 10^(mean(dat$log10_ic80_comb))
# conc <- gm_ic80_comb * sum(d$pt80)

x_values <- c(5, 10, 20, 50, 100, 200, 500, 2000, 5000, 20000)
x_breaks <- -log10(gm_ind_conc / x_values)

p <- ggplot() +
  geom_hline(yintercept = 1 - altHR, linetype = "dashed") +
  annotate("text", x = max(-df_ve$v), y = 0.8, hjust = 1.05, vjust = 1.6,
           label = "Overall PE (Ab high + Ab low vs. placebo)", size = 2.3) +
  geom_line(aes(x = -v, y = dens, group = factor(tx), color = factor(tx)),
            data = df_dens) +
  scale_color_manual(values = c("blue1", "magenta"),
                     labels = c("Placebo", "Ab high + Ab low"),
                     name = "PDF") +
  geom_line(aes(x = -v, y = ve), data = df_ve, linewidth = 1.2) +
  # scale_x_continuous(breaks = -1:3, labels = conc / 10^(1:-3)) +
  scale_x_continuous(breaks = x_breaks, labels = x_values, minor_breaks = NULL) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.25), 
                     labels = seq(0, 1, by = 0.25) * 100) +
  labs(x = "Combination PT80 of Combo-AMP Regimen\nagainst Autologous Virus", 
       y = "Prevention Efficacy (%)") +
  theme_bw() + 
  theme(legend.position = "bottom")

ggsave(here::here(file.path(path, "truePEbycombPT80_dataFromLily.pdf")), plot = p, height = 5, width = 4.9)


# Run the simulation ------------------------------------------------------

true_pe <- NA
thres <- NA
plot_est_pe <- FALSE
plot_pred_pe <- TRUE
plot_cp_pe <- FALSE
p <- list()
pred_pe <- list()
n_h_l <- 36

dens <- density(d_ic80$log10_comb_ic80, n = 1000)

registerDoParallel(cores = n_cores)

h <- list(c(1, 1, 1), 
          c(0.7, 1, 1), 
          c(1, 0.7, 1), 
          c(1, 1, 0.7),
          c(0.7, 0.7, 1),
          c(0.7, 1, 0.7),
          c(1, 0.7, 0.7),
          c(0.7, 0.7, 0.7))

start_time <- Sys.time()

for (j in 1:length(h)){
  # get PT80s of a new regimen
  d_pt80 <- get_comb_pt80(d_ic80, 
                          h = h[[j]], 
                          df_c = d_ind_conc %>% 
                            filter(id %in% 1001:1400))
  
  # a list of data frames
  l_pe <- est_pe_by_log10ic80(n_total = n_total, p_pla = p_pla, p_ab_l = p_ab_l, 
                              p_ab_h = p_ab_h, rate_pla = rate_pla, 
                              altHR_l = altHR_l, altHR_h = altHR_h, 
                              rate_cens = rate_cens, n_target_cases_h_l = n_h_l,
                              dens = dens, beta = beta, iter = iter)
  
  idx <- which(sapply(l_pe, function(x){ is.character(x) | is.na(x) }))
  if (length(idx) > 0){ l_pe <- l_pe[-idx] }
  
  if (plot_est_pe){
    df <- do.call(rbind, l_pe)
    
    p[[j]] <- ggplot() +
      geom_line(aes(x = -mark, y = TE, group = iter), data = df, alpha = 0.05) +
      geom_line(aes(x = -v, y = ve), data = df_ve, linewidth = 1.2, color = "red") +
      coord_cartesian(xlim = c(-1, 3), ylim = c(0, 1)) +
      # scale_x_continuous(breaks = -1:3, labels = 20 / 10^(1:-3)) +
      scale_x_continuous(breaks = x_breaks, labels = x_values, minor_breaks = NULL) +
      scale_y_continuous(breaks = seq(0, 1, by = 0.25), 
                         labels = seq(0, 1, by = 0.25) * 100) +
      labs(x = "Combination PT80 of Combo-AMP Regimen\nagainst Autologous Virus", 
           y = "Est. Prevention Efficacy (%)",
           title = paste0("High-Ab + Low-Ab Endpoint Count = ", n_h_l)) +
      theme_bw()
    
    ggsave(here::here(file.path(path, paste0("estPEbycombPT80_target_h_l=", n_h_l, ".pdf"))), 
           plot = p[[j]], height = 5, width = 5)
  }
  
  # predict_pe_dt() is faster than predict_pe(), predict_pe_parallel(), or
  # predict_pe_parallel_dt()
  pred_pe[[j]] <- lapply(l_pe, predict_pe_dt, scale_c = gm_ind_conc, 
                         df_comb_pt80 = d_pt80, approxfun_rule = 1)
  pred_pe[[j]] <- as.data.frame(do.call(rbind, pred_pe[[j]]))
  pred_pe[[j]]$n_h_l <- n_h_l
  pred_pe[[j]]$j <- j
  pred_pe[[j]]$h_PGDM <- h[[j]][1]
  pred_pe[[j]]$h_PGT <- h[[j]][2]
  pred_pe[[j]]$h_VRC <- h[[j]][3]
}

pred_pe <- bind_rows(pred_pe)

if (!is.na(true_pe)){
  pred_pe <- pred_pe %>%
    mutate(indCover_pe_A = as.numeric(lb_pe_A < true_pe & ub_pe_A > true_pe),
           indCover_pe_cA = as.numeric(lb_pe_cA < true_pe & ub_pe_cA > true_pe),
           indCover_pe_B = as.numeric(lb_pe_B < true_pe & ub_pe_B > true_pe))  
}

pred_pe_fname <- paste0("predPE_highDose_n_h_l=", n_h_l, "_h=vary_",
                    format(Sys.time(), "%d%b%Y_%H%M"), ".rds")
saveRDS(pred_pe, here::here(path, pred_pe_fname))
diff_time <- Sys.time() - start_time
diff_time

if (!is.na(true_pe)){
  pred_pe %>%
    group_by(n_h_l) %>%
    summarise(cp_A = mean(indCover_pe_A),
              cp_cA = mean(indCover_pe_cA),
              cp_B = mean(indCover_pe_B),
              .groups = "drop")
  
  if (!is.na(thres)){
    pred_pe %>%
      group_by(n_h_l) %>%
      summarise(p_ll95pe_above_thres = mean(lb_pe_B > thres), 
                .groups = "drop")  
  }
}


# Plot MC median estimates and CIs across different 'h' vectors -----------

df <- readRDS(here::here(path, pred_pe_fname)) %>%
  group_by(j) %>%
  summarise(med_ptEst_pe_B = as.numeric(quantile(ptEst_pe_cAB, prob = 0.5)),
            med_lb_pe_B = as.numeric(quantile(lb_pe_B, prob = 0.5)),
            med_ub_pe_B = as.numeric(quantile(ub_pe_B, prob = 0.5)),
            .groups = "drop")

h_labs <- sapply(h, function(vec) paste(vec, collapse = "\n"))

p <- ggplot(df, aes(x = factor(j), y = med_ptEst_pe_B)) +
  geom_errorbar(aes(ymin = med_lb_pe_B, ymax = med_ub_pe_B), width = 0.4) +
  geom_point(size = 2) +
  scale_x_discrete(labels = h_labs) +
  coord_cartesian(ylim = c(0.75, 1), clip = "off") +
  # scale_y_continuous(breaks = seq(0.7, 1, by = 0.1),
  #                    labels = paste0(seq(0.7, 1, by = 0.1) * 100, "%")) +
  labs(x = "Multiplicative Constants for IC80s of the 3 New bnAbs\nagainst All HVTN 704 Placebo Viruses",
       y = "Monte-Carlo Median Estimate and 95% CI for PE") +
  theme_bw() +
  theme(axis.text.x = element_text(lineheight = 0.9),
        plot.margin = margin(20, 5, 5, 5)) +
  annotation_custom(
    textGrob(paste(c("h1", "h2", "h3"), collapse = "\n"),
             x = unit(0.02, "npc"),
             y = unit(-0.085, "npc"),
             just = "right",
             gp = gpar(fontsize = 7.5, lineheight = 1.1))
  )
p
ggsave(here::here(path, paste0(unlist(strsplit(pred_pe_fname, "\\."))[1], ".pdf")), 
       plot = p, height = 4, width = 5.5)


# Plot multiple panels of estimated PE by comb PT80 -----------------------

if (plot_est_pe){
  combined_p <- Reduce(`+`, p) + plot_layout(ncol = 3)
  ggsave(here::here(path, "estPEbycombPT80_dataFromLily.pdf"), plot = combined_p, 
         height = 5, width = 15)  
}


# Plot MC med estimates, CIs, coverage prob if true PE known --------------

if (plot_cp_pe){
  df <- readRDS(here::here(path, pred_pe_fname)) %>%
    group_by(n_h_l) %>%
    summarise(med_ptEst_pe_A = as.numeric(quantile(ptEst_pe_A, prob = 0.5)),
              med_ptEst_pe_cA = as.numeric(quantile(ptEst_pe_cAB, prob = 0.5)),
              med_ptEst_pe_B = as.numeric(quantile(ptEst_pe_cAB, prob = 0.5)),
              med_lb_pe_A = as.numeric(quantile(lb_pe_A, prob = 0.5)),
              med_ub_pe_A = as.numeric(quantile(ub_pe_A, prob = 0.5)),
              med_lb_pe_cA = as.numeric(quantile(lb_pe_cA, prob = 0.5)),
              med_ub_pe_cA = as.numeric(quantile(ub_pe_cA, prob = 0.5)),
              med_lb_pe_B = as.numeric(quantile(lb_pe_B, prob = 0.5)),
              med_ub_pe_B = as.numeric(quantile(ub_pe_B, prob = 0.5)),
              cp_A = mean(indCover_pe_A),
              cp_cA = mean(indCover_pe_cA),
              cp_B = mean(indCover_pe_B),
              .groups = "drop") %>%
    pivot_longer(cols = -n_h_l,
                 names_to = c(".value", "method"),
                 names_pattern = "(.*)_(A|cA|B)$") %>%
    mutate(n_h_l = factor(n_h_l),
           method = forcats::fct_recode(factor(method, 
                                               levels = c("A", "cA", "B")),
                                        "Method A" = "A",
                                        "Bias-Corrected Method A" = "cA",
                                        "Method B" = "B"),
           cp_fmt = sapply(cp, function(x) format(round(x, 2), nsmall = 2)))
  
  p <- ggplot(df, aes(x = n_h_l, y = med_ptEst_pe)) +
    geom_hline(yintercept = true_pe, linetype = "dashed", size = 0.8, color = "red") +
    geom_errorbar(aes(ymin = med_lb_pe, ymax = med_ub_pe), width = 0.4) +
    geom_point(size = 2) +
    geom_text(aes(label = cp_fmt, y = med_lb_pe),
              vjust = 1.7, size = 3) +
    coord_cartesian(ylim = c(0.7, 1)) +
    scale_y_continuous(breaks = seq(0.7, 1, by = 0.1),
                       labels = paste0(seq(0.7, 1, by = 0.1) * 100, "%")) +
    facet_grid(. ~ method, scales = "fixed") +
    labs(x = "Combo-AMP High-Ab + Low-Ab Endpoint Count",
         y = "Monte-Carlo Median Estimate and 95% CI for PE") +
    theme_bw()
  p
  ggsave(here::here(path, paste0("bias_cp_", unlist(strsplit(pred_pe_fname, "\\."))[1], ".pdf")), 
         plot = p, height = 4, width = 6)
}


# Plot MC distributions of est and LL95 for PE ----------------------------

if (plot_pred_pe){
  df <- readRDS(here::here(path, pred_pe_fname)) %>%
    filter(n_h_l == 36)
  p <- list()
  p[[1]] <- ggplot(df, aes(x = factor(n_h_l), y = ptEst_pe_B)) +
    geom_boxplot(color = "black", width = 0.5, lwd = 0.6, outlier.shape = 1, 
                 outlier.alpha = 0.5) +
    # geom_hline(yintercept = 1 - altHR, linetype = "dashed") +
    coord_cartesian(ylim = c(0, 1)) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.2),
                       labels = paste0(seq(0, 1, by = 0.2) * 100, "%")) +
    labs(x = "Combo-AMP High-Ab + Low-Ab\nEndpoint Count",
         y = "Monte-Carlo Sampling Distribution of\nPoint Estimates of PE",
         title = "Algorithm B:\nPoint Estimate of PE") +
    theme_bw()
  
  p[[2]] <- ggplot(df, aes(x = factor(n_h_l), y = lb_pe_B)) +
    geom_hline(yintercept = true_pe, linetype = "dashed", size = 0.8, color = "red") +
    geom_hline(yintercept = 0.7, linetype = "dotted", size = 1.2, color = "blue") +
    geom_boxplot(color = "black", width = 0.5, lwd = 0.6, outlier.shape = 1, 
                 outlier.alpha = 0.5) +
    coord_cartesian(ylim = c(0.4, 1)) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.1),
                       labels = paste0(seq(0, 1, by = 0.1) * 100, "%")) +
    labs(x = "Combo-AMP\nHigh-Ab + Low-Ab\nEndpoint Count",
         y = "Monte-Carlo Sampling Distribution of\nLower 95% Confidence Limit for PE") +
    theme_bw()
  p[[2]]
  ggsave(here::here(path, paste0("LL95PE_B_", unlist(strsplit(pred_pe_fname, "\\."))[1], ".pdf")), 
         plot = p[[2]], height = 5 * 0.9, width = 3 * 0.9)
  
  p[[3]] <- ggplot(df, aes(x = factor(n_h_l), y = lb_pe_cA)) +
    geom_boxplot(color = "black", width = 0.5, lwd = 0.6, outlier.shape = 1, 
                 outlier.alpha = 0.5) +
    coord_cartesian(ylim = c(0, 1)) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.2),
                       labels = paste0(seq(0, 1, by = 0.2) * 100, "%")) +
    labs(x = "Combo-AMP High-Ab + Low-Ab\nEndpoint Count",
         y = "Monte-Carlo Sampling Distribution of\nLower 95% Uncertainty Limit for PE",
         title = "Bias-Corrected Algorithm A:\nLower 95% Uncertainty Limit") +
    theme_bw()
  
  # combined_p <- Reduce(`+`, p) + plot_layout(ncol = 3)
  # ggsave(here::here(path, "MCdistribPredPE_highDose_h=1_1_1_v4.pdf"), plot = combined_p, 
  #        height = 5, width = 10)
}


