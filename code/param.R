# Input parameters --------------------------------------------------------

# compute total events for adequate power under a given design alternative
# (Schoenfeld, 1983, Biometrics) as a starting value
ceiling(((qnorm(0.975) + qnorm(0.9))^2) / ((1/4) * (log(0.4) - log(0.1))^2))

n_target_cases <- 43
rate_pla <- 0.0075
nullHR <- 0.4
altHR_h <- 0.1
altHR_l <- 0.3
rate_cens <- 0.075
p_ab_h <- 1 / 3
p_ab_l <- 1 / 3
p_pla <- 1 / 3
n_enroll_m <- NULL
tau <- 1.5
alpha_1sided <- 0.025
info_fractions <- c(0.75, 1)
iter <- 1000
n_cores <- 16
pwr <- 0.9
path <- "output/2025-07-29"

n_total <- 10821

pop_conc_files <- c("sim_PK_vrc07_fixDose0.4g_wk26.csv",
                    "sim_PK_vrc07_fixDose3.2g_wk26.csv",
                    "sim_PK_pgt121ls_fixDose0.4g_wk26.csv",
                    "sim_PK_pgt121ls_fixDose1.6g_wk26.csv",
                    "sim_PK_pgdm1400ls_fixDose0.4g_wk26.csv",
                    "sim_PK_pgdm1400ls_fixDose1.6g_wk26.csv")
ind_conc_files <- c("PK_1000individal_fixDose_0.4vrc07_0.4pgt121ls_0.4pgdm1400ls_over26Wk.csv",
                    "PK_1000individal_fixDose_3.2vrc07_1.6pgt121ls_1.6pgdm1400ls_over26Wk.csv")
ic80_file <- "704_placebo_IC80_vrc01_vrc07_pgt121_pgdm1400.csv"