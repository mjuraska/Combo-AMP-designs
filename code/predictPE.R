library(tidyverse)
library(np)

source(here::here("code/param.R"))
source(here::here("code/utils.R"))


# Get IC80 of Combo-AMP regimen against 704 placebo viruses ---------------

pubid_pla <- read.csv("/trials/vaccine/p704/analysis/efficacy/adata/amp_survival_wk104_neut.csv") %>%
  filter(protocol == "HVTN 704", tx == "C3", hiv1event == 1) %>%
  pull(pub_id)

dat <- read.csv("/trials/vaccine/p704/analysis/manuscripts/NeutTiterBiomarker/data/VTN704_Non_Par_breakthrough_NAb_20210331.csv") %>%
  select(protnum, isolate, poscrit, titer, mab_name) %>%
  mutate(pub_id = sapply(strsplit(isolate, split = "_"), function(x){ paste0(substring(x[1], first = 2), "-", x[2]) })) %>%
  filter(pub_id %in% pubid_pla, poscrit == 80, mab_name %in% c("VRC07-523LS", "PGT121.414LS", "PGDM1400")) %>%
  mutate(titer = as.numeric(gsub(">>", "", titer)),
         mab_name = recode(mab_name, "VRC07-523LS" = "1", "PGT121.414LS" = "2", 
                           "PGDM1400" = "3")) %>%
  pivot_wider(names_from = "mab_name", values_from = "titer", names_prefix = "ic80_ab") %>%
  mutate(log10_ic80_comb = log10(1 / ((1 / ic80_ab1) + (1 / ic80_ab2) + (1 / ic80_ab3))))

# file emailed by Lily on 8/13/25 as an alternative to the data approach above
# dat2 <- read.csv(here::here("data/704_placebo_IC80_vrc01_vrc07_pgt121_pgdm1400.csv")) %>%
#   mutate(log10_ic80_comb = log10(1 / ((1 / PGDM1400) + (1 / PGT121.414LS) + (1 / VRC07.523LS)))) %>%
#   select(log10_ic80_comb)

# Get density ratio model coefficients under PE(log10(IC80)=1) = 0 --------

# density estimate of log10(IC80_comb) against 704 placebo viruses
dens <- density(dat$log10_ic80_comb, n = 1000)
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

dens <- density(dat$log10_ic80_comb, n = 300)
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

ggsave(here::here(file.path(path, "truePEbyLog10combIC80.pdf")), plot = p, height = 5, width = 4.9)


# Plot PE by PT80 ---------------------------------------------------------

# emailed by Lily on 8/13/25
df_conc <- data.frame(ab = rep(c("PGDM1400LS", "PGT121LS", "VRC07.523LS"), each = 2),
                      dose = rep(c("IV 0.4g", "IV 1.6g"), 3),
                      med_conc = c(33.30190776, 133.207631, 48.21225924, 192.849037, 19.17685883, 153.4148706))
d <- df_conc %>% filter(dose == "IV 0.4g") %>% select(ab, med_conc) %>%
  mutate(gm_ic80 = c(exp(mean(log(dat$ic80_ab3))), exp(mean(log(dat$ic80_ab2))), exp(mean(log(dat$ic80_ab1)))),
         pt80 = med_conc / gm_ic80)
gm_ic80_comb <- 10^(mean(dat$log10_ic80_comb))
conc <- gm_ic80_comb * sum(d$pt80)

x_values <- c(1, 2, 5, 10, 20, 50, 100, 200, 500, 2000, 10000)
x_breaks <- -log10(conc / x_values)

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

ggsave(here::here(file.path(path, "truePEbycombPT80.pdf")), plot = p, height = 5, width = 4.9)


# Run the simulation ------------------------------------------------------

registerDoParallel(cores = n_cores)

dens <- density(dat$log10_ic80_comb, n = 1000)

p <- list()
k <- 1
for (n_h_l in c(25, 50, 75)){
  l <- est_PE_by_PT80(n_total = n_total, p_pla = p_pla, p_ab_l = p_ab_l, 
                      p_ab_h = p_ab_h, rate_pla = rate_pla, altHR_l = altHR_l,
                      altHR_h = altHR_h, rate_cens = rate_cens, n_target_cases_h_l = n_h_l,
                      dens = dens, beta = beta, iter = iter)
  
  idx <- which(sapply(l, function(x){ is.character(x) }))
  if (length(idx) > 0){ l <- l[-idx] }
  df <- do.call(rbind, l)
  
  p[[k]] <- ggplot() +
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
         plot = p[[k]], height = 5, width = 5)
  k <- k + 1
}

combined_p <- Reduce(`+`, p) + plot_layout(ncol = 3)
ggsave(here::here(path, "estPEbycombPT80.pdf"), plot = combined_p, 
       height = 5, width = 15)

