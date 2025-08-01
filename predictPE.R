library(tidyverse)
library(np)

source(here::here("param.R"))
source(here::here("utils.R"))


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
  mutate(log10_ic80_comb = log10(1 / ((1 / ic80_ab1) + (1 / ic80_ab2) + (1 / ic80_ab3)))) %>%
  select(log10_ic80_comb)


# Get density ratio model coefficients under PE(0) = 0 --------------------

# density estimate of log10(IC80_comb) in placebo arm
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

p <- ggplot() +
  geom_hline(yintercept = 1 - altHR, linetype = "dashed") +
  annotate("text", x = min(df$v), y = 0.8, hjust = 0.05, vjust = -0.6,
           label = "Overall PE (Ab high + Ab low vs. placebo)", size = 2.3) +
  geom_line(aes(x = v, y = dens, group = factor(tx), color = factor(tx)), 
            data = df_dens) +
  scale_color_manual(values = c("blue1", "magenta"), 
                     labels = c("Placebo", "Ab high + Ab low"),
                     name = "PDF") +
  geom_line(aes(x = v, y = ve), data = df_ve, size = 1) +
  scale_x_continuous(breaks = -3:1, labels = 10^(-3:1)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.25), 
                     labels = seq(0, 1, by = 0.25) * 100) +
  labs(x = "Combination IC80 of Combo-AMP Regimen", 
       y = "Prevention Efficacy (%)") +
  theme_bw() + 
  theme(legend.position = "bottom")

ggsave(here::here(file.path(path, "truePEbyLog10combIC80.pdf")), plot = p, height = 5, width = 5)


# Sample from the support using the estimated density as weights
samples <- sample(dens$x, size = 1000, replace = TRUE, prob = dens$y)
