library(tidyverse)
library(np)

source(here::here("code/param.R"))
source(here::here("code/utils.R"))


# Get IC80 of Combo-AMP regimen against 704 placebo viruses ---------------

# pubid_pla <- read.csv("/trials/vaccine/p704/analysis/efficacy/adata/amp_survival_wk104_neut.csv") %>%
#   filter(protocol == "HVTN 704", tx == "C3", hiv1event == 1) %>%
#   pull(pub_id)
# 
# dat <- read.csv("/trials/vaccine/p704/analysis/manuscripts/NeutTiterBiomarker/data/VTN704_Non_Par_breakthrough_NAb_20210331.csv") %>%
#   select(protnum, isolate, poscrit, titer, mab_name) %>%
#   mutate(pub_id = sapply(strsplit(isolate, split = "_"), function(x){ paste0(substring(x[1], first = 2), "-", x[2]) })) %>%
#   filter(pub_id %in% pubid_pla, poscrit == 80, mab_name %in% c("VRC07-523LS", "PGT121.414LS", "PGDM1400")) %>%
#   mutate(titer = as.numeric(gsub(">>", "", titer)),
#          mab_name = recode(mab_name, "VRC07-523LS" = "1", "PGT121.414LS" = "2", 
#                            "PGDM1400" = "3")) %>%
#   pivot_wider(names_from = "mab_name", values_from = "titer", names_prefix = "ic80_ab") %>%
#   mutate(log10_ic80_comb = log10(1 / ((1 / ic80_ab1) + (1 / ic80_ab2) + (1 / ic80_ab3))))

# files emailed by Lily on 8/13/25 as an alternative to the data approach above
d_ic80 <- read.csv(here::here("data", ic80_file)) %>%
  rename_with(~ gsub("\\.", "", .x)) %>%
  mutate(comb_ic80 = 1 / ((1 / PGDM1400) + (1 / PGT121414LS) + (1 / VRC07523LS)),
         log10_comb_ic80 = log10(comb_ic80),
         w_PGDM = comb_ic80 / PGDM1400,
         w_PGT = comb_ic80 / PGT121414LS,
         w_VRC = comb_ic80 / VRC07523LS,
         w_sum = w_PGDM + w_PGT + w_VRC)

d_pop_conc <- plyr::ldply(pop_conc_files, function(f){
  read.csv(here::here("data", f))
}) %>% 
  group_by(drug, dose) %>%
  summarise(med = quantile(Cc_ss, prob = 0.5),
            gm = exp(mean(log(Cc_ss))))

gm_conc_l <- exp(mean(log(drop(as.matrix(d_ic80 %>% select(w_PGDM, w_PGT, w_VRC)) %*% 
                              (d_pop_conc %>% filter(dose == "IV 0.4g") %>% pull(gm))))))
gm_conc_h <- exp(mean(log(drop(as.matrix(d_ic80 %>% select(w_PGDM, w_PGT, w_VRC)) %*% 
                              (d_pop_conc %>% filter(dose != "IV 0.4g") %>% pull(gm))))))
# this is the average effective concentration at exposure scaling constant
gm_conc <- exp(mean(log(c(gm_conc_h, gm_conc_l))))

d_ind_conc <- read.csv(here::here("data", ind_conc_file)) %>%
  # for consistency with IC80 file
  rename(conc_PGDM1400 = conc_PGDM1400LS)


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

x_values <- c(10, 20, 50, 100, 200, 500, 2000, 5000, 20000)
x_breaks <- -log10(gm_conc / x_values)

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

# get PT80s of a new regimen
d_pt80 <- get_comb_pt80(d_ic80, c(1, 1, 1), d_ind_conc)

dens <- density(d_ic80$log10_comb_ic80, n = 1000)

registerDoParallel(cores = n_cores)

plot_est_pe <- FALSE
plot_pred_pe <- TRUE
p <- list()
n_h_l <- c(25, 50, 75)

df_pred_pe <- plyr::ldply(1:length(n_h_l), function(j){
  # a list of data frames
  l_pe <- est_pe_by_log10ic80(n_total = n_total, p_pla = p_pla, p_ab_l = p_ab_l, 
                              p_ab_h = p_ab_h, rate_pla = rate_pla, 
                              altHR_l = altHR_l, altHR_h = altHR_h, 
                              rate_cens = rate_cens, n_target_cases_h_l = n_h_l[j],
                              dens = dens, beta = beta, iter = iter)
  
  idx <- which(sapply(l_pe, function(x){ is.character(x) }))
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
           title = paste0("High-Ab + Low-Ab Endpoint Count = ", n_h_l[j])) +
      theme_bw()
    
    ggsave(here::here(file.path(path, paste0("estPEbycombPT80_target_h_l=", n_h_l[j], ".pdf"))), 
           plot = p[[k]], height = 5, width = 5)
  }
  
  out <- lapply(l_pe, predict_pe, scale_c = gm_conc, df_comb_pt80 = d_pt80)
  out <- as.data.frame(do.call(rbind, out))
  colnames(out) <- c("pe", "lb")
  out$n_h_l <- n_h_l[j]
  
  return(out)
})

if (plot_est_pe){
  combined_p <- Reduce(`+`, p) + plot_layout(ncol = 3)
  ggsave(here::here(path, "estPEbycombPT80_dataFromLily.pdf"), plot = combined_p, 
         height = 5, width = 15)  
}

if (plot_pred_pe){
  p <- ggplot(df_pred_pe, aes(x = factor(n_h_l), y = pe)) +
    geom_boxplot(color = "black", width = 0.5, lwd = 0.6) +
    labs(x = "High-Ab + Low-Ab Endpoint Count",
         y = "Monte-Carlo Sampling Distribution of Point Estimates of PE")
  p
  
  p <- ggplot(df_pred_pe, aes(x = factor(n_h_l), y = lb)) +
    geom_boxplot(color = "black", width = 0.5, lwd = 0.6) +
    labs(x = "High-Ab + Low-Ab Endpoint Count",
         y = "Monte-Carlo Sampling Distribution of Lower 95\% Uncertainty Limit for PE")
  p
}
