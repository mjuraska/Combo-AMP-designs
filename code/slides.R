rm(list = ls(all = TRUE))

library(tidyverse)
library(grid)

source(here::here("code/param.R"))
source(here::here("code/utils.R"))

d_ind_conc <- read.csv(here::here("data", ind_conc_files[2])) %>%
  rename_with(~ case_when(grepl("PGDM", .x) ~ "bnAb1",
                          grepl("PGT", .x) ~ "bnAb2",
                          grepl("VRC", .x) ~ "bnAb3",
                          TRUE ~ .x)) %>%
  pivot_longer(cols = starts_with("bnAb"),
               names_to = "ab",
               values_to = "conc") %>%
  mutate(ab = factor(ab, levels = paste0("bnAb", 1:3))) %>%
  filter(id < 201)

p <- ggplot(d_ind_conc, aes(x = time, y = log10(conc), group = id)) +
  geom_line(alpha = 0.2) +
  facet_grid(. ~ ab, scales = "fixed") +
  scale_x_continuous(breaks = c(seq(0, 7 * 23, by = 28), 7 * 26), 
                     labels = c(seq(0, 7 * 23, by = 28) / 7, 26)) +
  scale_y_continuous(breaks = 0:3, labels = 10^(0:3)) +
  labs(x = "Weeks since Administration",
       y = expression("Serum Concentration (" * mu * "g/ml)")) +
  theme_bw()
p
ggsave(here::here(path, "ind_conc.pdf"), plot = p, height = 4, width = 5)


d_ic80 <- read.csv(here::here("data", ic80_file)) %>%
  select(-VRC01) %>%
  rename_with(~ case_when(grepl("PGDM", .x) ~ "bnAb1",
                          grepl("PGT", .x) ~ "bnAb2",
                          grepl("VRC", .x) ~ "bnAb3",
                          TRUE ~ .x)) %>%
  pivot_longer(cols = starts_with("bnAb"),
               names_to = "ab",
               values_to = "ic80") %>%
  mutate(ab = factor(ab, levels = paste0("bnAb", 1:3)),
         x_jit = as.numeric(ab) + runif(n(), -0.10, 0.10))

set.seed(123)
p <- ggplot(d_ic80, aes(y = log10(ic80))) +
  geom_boxplot(aes(x = ab), width = 0.5, lwd = 0.6, outlier.shape = NA, alpha = 0) +
  geom_line(aes(x = x_jit, group = isolate), alpha = 0.3, color = "gray50") + 
  geom_point(aes(x = x_jit)) +
  scale_y_continuous(breaks = c(-2:1, log10(50)), labels = c(10^(-2:1), 50)) +
  labs(x = NULL, 
       y = expression(IC[80] * " (" * mu * "g/ml) against Virus Panel")) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 11))
p
ggsave(here::here(path, "ic80.pdf"), plot = p, height = 4, width = 5)

