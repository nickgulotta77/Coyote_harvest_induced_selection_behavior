################################################################################
# Manuscript title: Social status mediates harvest-induced selection on 
#                   behavioral type and predictability in Coyotes (Canis latrans)
#
# Script title:     Double Hierarchical Generalized Linear Model (DHGLM) Analysis
#                   Behavioral Type, Predictability (rIIV), and Survival in Coyotes
# Authors: Nick A Gulotta, Joey H. Hinton, Michael J. Chamberlain
# Date:    06/08/2026
# Journal: Ecology Letters
#
# Description:
#   This script fits Bayesian DHGLMs to estimate among-individual differences
#   in behavioral type (mean habitat use) and residual intra-individual
#   variation (rIIV; predictability) in distance to shrub and hardwood
#   landcover for resident and transient coyotes. Individual-level
#   random effects (BLUPs) from each model are then linked to
#   survival using a Bayesian model-averaging approach that propagates
#   uncertainty from the DHGLM through to the survival analysis.
################################################################################


# ------------------------------------------------------------------------------
# 1. Load Required Packages
# ------------------------------------------------------------------------------

require(parallel)
require(remotes)
require(ggplot2)
require(knitr)
require(readr)
require(ggeffects)
require(dplyr)
require(tidybayes)
require(coda)
require(lubridate)
require(posterior)
require(brms)
require(tidyr)
require(psych)


# ------------------------------------------------------------------------------
# 2. Custom ggplot2 Theme
# ------------------------------------------------------------------------------

theme_turkey <- function(){
  theme_bw() +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          axis.text.x = element_text(size=12),
          axis.line.x = element_line(color = "black"),
          axis.line.y = element_line(color = "black"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank(),
          plot.title = element_text(size = 15),
          legend.text = element_text(size = 10),
          legend.title = element_blank(),
          legend.key = element_rect(colour = NA, fill = NA),
          legend.background = element_rect(color = "black",
                                           fill = "transparent",
                                           size = 4, linetype = "blank"))
}


# ------------------------------------------------------------------------------
# 3. Configure Stan Backend and Set Random Seed
# ------------------------------------------------------------------------------

library(cmdstanr)

options(mc.cores = 12,
        brms.backend = "cmdstanr")

bayes_seed <- 1234


# ------------------------------------------------------------------------------
# 4. Load and Prepare Data
# ------------------------------------------------------------------------------

data <- read_csv("TriState_daily_Landcover_FINAL.csv")

# Retain only daily locations with >= 4 GPS fixes to ensure reliable
# daily landcover estimates
data <- data %>%
  group_by(n_points) %>%
  filter(n_points >= 4)


# ------------------------------------------------------------------------------
# 5. Transform and Standardize Response Variables
# ------------------------------------------------------------------------------

# --- Shrub landcover ---
# Square-root transform to improve normality, then standardize (mean = 0, SD = 1)
hist(data$Shrub)
describe(data$Shrub)

data$Shrub_s <- sqrt(data$Shrub)
hist(data$Shrub_s)
describe(data$Shrub_s)

data$shrub_c <- scale(data$Shrub_s, center = T, scale = T)
hist(data$shrub_c)
describe(data$shrub_c)

# --- Hardwood landcover ---
# Square-root transform to improve normality, then standardize (mean = 0, SD = 1)
hist(data$Hardwoods)
describe(data$Hardwoods)

data$Hardwoods_s <- sqrt(data$Hardwoods)
hist(data$Hardwoods_s)
describe(data$Hardwoods_s)

data$Hardwoods_c <- scale(data$Hardwoods_s, center = T, scale = T)
hist(data$Hardwoods_c)
describe(data$Hardwoods_c)

# ------------------------------------------------------------------------------
# 5a.Factorize and reorder fixed effects
# ------------------------------------------------------------------------------

data$Year<-as.factor(data$Year)
data$Sex<-as.factor(data$Sex)
data$StudySite<-as.factor(data$StudySite)
data$Status_clean<-as.factor(data$Status_clean)
data$Season = factor(data$Season, levels = c( "Spring", "Summer", "Fall", "Winter"))

# ------------------------------------------------------------------------------
# 6. Summarize Monitoring Data
# ------------------------------------------------------------------------------

library(dplyr)
library(lubridate)

df <- data %>% mutate(Date = as.Date(Date))

# Number of unique individuals per year and study site
ids_per_year <- df %>%
  group_by(Year, StudySite) %>%
  summarise(
    n_individuals = n_distinct(ID),
    .groups = "drop"
  )
print(ids_per_year)

ids_per_year_site <- df %>%
  group_by(Year, StudySite) %>%
  summarise(n_individuals = n_distinct(ID), .groups = "drop")
print(ids_per_year_site)

# Identify individuals monitored across more than one year
multi_year_ids <- df %>%
  group_by(ID) %>%
  summarise(
    n_years = n_distinct(Year),
    years   = paste(sort(unique(Year)), collapse = ", "),
    .groups = "drop"
  ) %>%
  filter(n_years > 1)
print(multi_year_ids)

# Identify individuals with > 365 days between first and last GPS fix
over_1yr_span <- df %>%
  group_by(ID) %>%
  summarise(
    first_date = min(Date, na.rm = TRUE),
    last_date  = max(Date,  na.rm = TRUE),
    span_days  = as.integer(last_date - first_date) + 1,
    .groups = "drop"
  ) %>%
  filter(span_days > 365)
print(over_1yr_span)

# Summarize number of daily observations per individual
daily_counts <- data %>%
  group_by(ID) %>%
  summarise(n_days = n())
daily_counts

daily_counts %>%
  summarise(
    min_days    = min(n_days),
    max_days    = max(n_days),
    mean_days   = mean(n_days),
    median_days = median(n_days)
  )


# ------------------------------------------------------------------------------
# 7. DHGLM Model Specifications and MCMC Settings
# ------------------------------------------------------------------------------

library(brms)
options(brms.backend = "cmdstanr")
options(mc.cores = parallel::detectCores())

# MCMC settings
chains  <- 2
iter    <- 4000
warmup  <- 2000
cores   <- 4

# Count of individuals per dispersal status
ct <- data %>%
  ungroup() %>%
  distinct(ID, Status_clean) %>%
  count(Status_clean)

ct


# ------------------------------------------------------------------------------
# 8A. DHGLM - Hardwood Landcover
# ------------------------------------------------------------------------------
# The DHGLM simultaneously models:
#   (1) mean behavioral expression (behavioral type) via the location submodel
#   (2) residual intra-individual variation (rIIV; predictability) via the
#       dispersion submodel (sigma)
# Random effects of individual ID in both submodels allow estimation of
# among-individual differences in both behavioral type and predictability

double_model1 = bf(Hardwoods_c ~ Sex + Status_clean*Season + StudySite + Year + (1 | ID),
                   sigma ~ Sex + Status_clean*Season + (1|ID))

hardwood_brm <- brm(double_model1,
                    data = data, family = "gaussian",
                    iter = 4000,
                    warmup = 1000,
                    chains = 4,
                    threads = threading(2),
                    init = "random",
                    cores = 15,
                    seed = bayes_seed,
                    backend = "cmdstanr",
                    control = list(adapt_delta = 0.99, max_treedepth = 15))

 #saveRDS(hardwood_brm, "hardwood_riiv_FINAL_check.rds")
hardwood_brm <- read_rds("hardwood_riiv_FINAL_check.rds")

summary(hardwood_brm)
plot(conditional_effects(hardwood_brm, dpar = "sigma"))
plot(conditional_effects(hardwood_brm))
describe_posterior(hardwood_brm)

# Extract and summarize posterior distributions for dispersion submodel parameters
library(brms)
library(bayestestR)
library(posterior)
library(dplyr)

draws <- as_draws_df(hardwood_brm)
sigma_pars <- grep("^b_sigma_", names(draws), value = TRUE)
sigma_results <- describe_posterior(draws[, sigma_pars])
sigma_results

# --- Variance partitioning for hardwood behavioral type ---

# Among-individual variance (behavioral type)
var.ID <- posterior_samples(hardwood_brm)$"sd_ID__Intercept"^2
mean(var.ID); HPDinterval(as.mcmc(var.ID), 0.95)

# Residual (within-individual) variance
var.res <- exp(posterior_samples(hardwood_brm)$"b_sigma_Intercept")^2
mean(var.res); HPDinterval(as.mcmc(var.res), 0.95)

# Repeatability (proportion of total variance explained by among-individual differences)
rep <- var.ID / (var.ID + var.res)
mean(rep); HPDinterval(as.mcmc(rep), 0.95)

# Among-individual variation in predictability (rIIV)
var.rIIV <- (posterior_samples(hardwood_brm)$"sd_ID__sigma_Intercept")
mean(var.rIIV); HPDinterval(as.mcmc(var.rIIV), 0.95)

# Coefficient of variation in predictability (CVP)
log.norm.res <- exp(posterior_samples(hardwood_brm)$"sd_ID__sigma_Intercept"^2)
CVP <- sqrt(log.norm.res - 1)
mean(CVP); HPDinterval(as.mcmc(CVP), 0.95)


# ------------------------------------------------------------------------------
# 8B. DHGLM - Shrub Landcover
# ------------------------------------------------------------------------------

double_model2 = bf(shrub_c ~ Sex + Status_clean*Season + StudySite + Year + (1 | ID),
                   sigma ~ Sex + Status_clean*Season + (1|ID))

shrub_brm <- brm(double_model2,
                 data = data, family = "gaussian",
                 iter = 4000,
                 warmup = 1000,
                 chains = 4,
                 threads = threading(2),
                 init = "random",
                 cores = 15,
                 seed = bayes_seed,
                 backend = "cmdstanr",
                 control = list(adapt_delta = 0.99, max_treedepth = 15))

 saveRDS(shrub_brm, "shrub_riiv_FINAL_check.rds")
#shrub_brm <- read_rds("shrub_riiv_FINAL_check.rds")

summary(shrub_brm)
plot(conditional_effects(shrub_brm, dpar = "sigma"))
plot(conditional_effects(shrub_brm))
describe_posterior(shrub_brm)

# Extract and summarize posterior distributions for dispersion submodel parameters
draws <- as_draws_df(shrub_brm)
sigma_pars <- grep("^b_sigma_", names(draws), value = TRUE)
sigma_results <- describe_posterior(draws[, sigma_pars])
sigma_results

# --- Variance partitioning for shrub behavioral type ---

# Among-individual variance (behavioral type)
var.ID <- posterior_samples(shrub_brm)$"sd_ID__Intercept"^2
mean(var.ID); HPDinterval(as.mcmc(var.ID), 0.95)

# Residual (within-individual) variance
var.res <- exp(posterior_samples(shrub_brm)$"b_sigma_Intercept")^2
mean(var.res); HPDinterval(as.mcmc(var.res), 0.95)

# Repeatability
rep <- var.ID / (var.ID + var.res)
mean(rep); HPDinterval(as.mcmc(rep), 0.95)

# Among-individual variation in predictability (rIIV)
var.rIIV <- (posterior_samples(shrub_brm)$"sd_ID__sigma_Intercept")
mean(var.rIIV); HPDinterval(as.mcmc(var.rIIV), 0.95)

# Coefficient of variation in predictability (CVP)
log.norm.res <- exp(posterior_samples(shrub_brm)$"sd_ID__sigma_Intercept"^2)
CVP <- sqrt(log.norm.res - 1)
mean(CVP); HPDinterval(as.mcmc(CVP), 0.95)


# ==============================================================================
# 9. Fitness Analyses: Linking Behavioral BLUPs to Annual Survival
#
#    Approach: For each posterior draw of individual-level random effects
#    (BLUPs), we fit a binomial GLM regressing annual survival on BLUP values.
#    Repeating this across 2000 posterior draws propagates uncertainty from
#    the DHGLM through to the survival analysis, producing a posterior
#    distribution of log-odds ratios that accounts for measurement error in
#    individual behavioral estimates.
# ==============================================================================


# ------------------------------------------------------------------------------
# 9A. Shrub - Transients - Behavioral Type
# ------------------------------------------------------------------------------

require(dplyr)
variables(shrub_brm)
draws_df <- brms::as_draws_df(shrub_brm)

# Extract posterior draws of individual-level random effects (behavioral type BLUPs)
blup_cols <- grep("^r_ID\\[", names(draws_df), value = TRUE)
meta_cols  <- intersect(c(".chain", ".draw", ".iteration"), names(draws_df))

o1 <- draws_df %>%
  dplyr::select(dplyr::all_of(c(blup_cols, meta_cols))) %>%
  dplyr::slice_sample(n = 2000) %>%
  dplyr::select(-dplyr::any_of(c(".chain", ".draw", ".iteration")))

# Load survival data and filter to transients, excluding unknown/non-natural causes of death
sur <- read_csv("survival_meta_02_17_2026.csv")
survival_final <- sur %>%
  filter(Status_clean %in% c("resident", "transient")) %>%
  group_by(COD) %>%
  filter(!COD %in% c("Unknown", "Vehicle", "Car", "Distemper"))

require(bayestestR)

s_ga <- survival_final %>%
  group_by(Status_clean) %>%
  filter(Status_clean != "resident")

# Bayesian model-averaging loop: fit one GLM per posterior draw of BLUPs
FR2_lm_list <- list()
require(data.table)
require(tidyr)
require(stringr)

for (iter in seq_len(nrow(o1))) {
  blups_FR2 <- select(as_tibble(o1), contains("r_ID"))[iter, ]
  blups_FR  <- tibble(
    ID       = str_remove(colnames(blups_FR2), "r_ID"),
    blups_FR = as.numeric(blups_FR2))
  data_FR2_lm          <- merge(blups_FR, s_ga, by = "ID")
  FR2_lm_list[[iter]]  <- glm(Survival ~ blups_FR, data = data_FR2_lm, family = binomial)
}

# Extract log-odds ratios across all GLMs to form posterior distribution
harv_ga_da <- as.mcmc(
  sapply(FR2_lm_list, function(x) { summary(x)$coefficients[2, 1] })
)

plot(harv_ga_da)
harv_ga_da_sur <- data.frame(harv_ga_da)
describe_posterior(harv_ga_da_sur)

# Extract mean BLUP per individual for plotting
library(data.table)
FR <- select(as_tibble(o1), contains("r_ID"))

FR_long <- FR %>%
  pivot_longer(everything(), names_to = c("ID")) %>%
  tibble()

FR_long$ID <- str_remove_all(FR_long$ID, "r_ID")
FR_long     <- tibble(FR_long)

FR <- FR_long %>%
  group_by(ID) %>%
  summarise(BLUP = mean(value))

FR <- left_join(s_ga, FR)

# Generate predicted survival curves from each GLM and summarize
predictions_df <- data.frame()

for (i in seq_along(FR2_lm_list)) {
  model                   <- FR2_lm_list[[i]]
  data_FR2_lm$prediction  <- predict(model, newdata = data_FR2_lm, type = "response")
  predictions_df          <- rbind(predictions_df, data.frame(model = rep(i, nrow(data_FR2_lm)), data_FR2_lm))
}

predictions_df$BLUP_raw <- predictions_df$blups_FR

# Summarize predictions: mean and 95% credible intervals per individual
df_summary <- predictions_df %>%
  group_by(ID) %>%
  summarise(mean     = mean(prediction, na.rm = TRUE),
            ci_lower = quantile(prediction, 0.025),
            ci_upper = quantile(prediction, 0.975),
            BLUP_raw = mean(BLUP_raw))

dr          <- FR %>% select(ID, BLUP, Survival)
df_summary1 <- left_join(df_summary, dr)
df_summary1$BLUP_r <- df_summary1$BLUP

# Back-transform standardized BLUPs to original scale for plotting
df_summary1$BLUP_raw1 <- df_summary1$BLUP_raw * 3.37 + 10.3   # revert square-root transformation
df_summary1$BLUP_raw2 <- df_summary1$BLUP_raw1^2               # revert to original units
df_summary1$BLUP      <- df_summary1$BLUP * 3.37 + 10.3
df_summary1$BLUP      <- df_summary1$BLUP^2

# Plot: behavioral type (shrub) vs. survival probability - transients
require(ggplot2)
p33 <- ggplot(df_summary1, aes(x = BLUP_raw2, y = mean)) +
  geom_line(color="#3E4A4F", linewidth=5) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper),fill="#3E4A4F",linetype=0, alpha = 0.25, linewidth=0) +
  geom_point(aes(x=BLUP_raw2, y= Survival), color= "#3E4A4F",alpha=0.6,size=8)+ theme_turkey()+
  xlab("Distance to shrub landcover (m)")+
  ylab("Survival Probaility")+
  theme_turkey()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14),plot.margin = margin(r=20, l=15))+
  scale_y_continuous(n.breaks = 5)+
  ggtitle(label ="", subtitle = "") +
  theme(legend.title = element_text(size=74, face="bold"), legend.position = "bottom") +
  theme(axis.title = element_text(size = 70, face = "bold"),
        plot.title = element_text(size=72),
        legend.text = element_blank(),
        legend.position = "none",
        axis.ticks.length = unit(.5, "cm"),
        axis.ticks = element_line(linewidth = 3, color = "black"),
        panel.border = element_rect(colour ="black", fill=NA, size=2.05),
        axis.text.y   = element_text(size=68,colour = "black"),
        axis.text.x   = element_text(size=68, colour = "black"),
        panel.grid.major = element_blank(),
        plot.margin = margin(r=55, l=20, b=35, t=35),
        panel.background = element_blank(),
        panel.grid.minor = element_blank())
p33


# ------------------------------------------------------------------------------
# 9B. Shrub - Transients - Predictability (rIIV)
# ------------------------------------------------------------------------------

require(dplyr)
variables(shrub_brm)
draws_df <- brms::as_draws_df(shrub_brm)

# Extract posterior draws of individual-level random effects from dispersion submodel (rIIV BLUPs)
blup_cols <- grep("^r_ID__sigma\\[", names(draws_df), value = TRUE)
meta_cols  <- intersect(c(".chain", ".draw", ".iteration"), names(draws_df))

o1 <- draws_df %>%
  dplyr::select(dplyr::all_of(c(blup_cols, meta_cols))) %>%
  dplyr::slice_sample(n = 2000) %>%
  dplyr::select(-dplyr::any_of(c(".chain", ".draw", ".iteration")))

# Load and filter survival data
sur            <- read_csv("survival_meta_02_17_2026.csv")
survival_final <- sur %>%
  filter(Status_clean %in% c("resident", "transient")) %>%
  group_by(COD) %>%
  filter(!COD %in% c("Unknown", "Vehicle", "Car", "Distemper"))

require(bayestestR)
s_ga <- survival_final %>%
  group_by(Status_clean) %>%
  filter(Status_clean != "resident")

# Bayesian model-averaging loop
FR2_lm_list <- list()
require(data.table)
require(tidyr)
require(stringr)

for (iter in seq_len(nrow(o1))) {
  blups_FR2 <- select(as_tibble(o1), contains("r_ID__sigma"))[iter, ]
  blups_FR  <- tibble(
    ID       = str_remove(colnames(blups_FR2), "r_ID__sigma"),
    blups_FR = as.numeric(blups_FR2))
  data_FR2_lm         <- merge(blups_FR, s_ga, by = "ID")
  FR2_lm_list[[iter]] <- glm(Survival ~ blups_FR, data = data_FR2_lm, family = binomial)
}

# Extract log-odds ratios
harv_ga_da <- as.mcmc(
  sapply(FR2_lm_list, function(x) { summary(x)$coefficients[2, 1] })
)

plot(harv_ga_da)
harv_ga_da_sur <- data.frame(harv_ga_da)
describe_posterior(harv_ga_da_sur)

# Extract mean BLUP per individual for plotting
library(data.table)
FR <- select(as_tibble(o1), contains("r_ID__sigma"))

FR_long     <- FR %>% pivot_longer(everything(), names_to = c("ID")) %>% tibble()
FR_long$ID  <- str_remove_all(FR_long$ID, "r_ID__sigma")
FR_long     <- tibble(FR_long)

FR <- FR_long %>%
  group_by(ID) %>%
  summarise(BLUP = mean(value))

FR <- left_join(s_ga, FR)

# Generate predicted survival curves and summarize
predictions_df <- data.frame()

for (i in seq_along(FR2_lm_list)) {
  model                   <- FR2_lm_list[[i]]
  data_FR2_lm$prediction  <- predict(model, newdata = data_FR2_lm, type = "response")
  predictions_df          <- rbind(predictions_df, data.frame(model = rep(i, nrow(data_FR2_lm)), data_FR2_lm))
}

predictions_df$BLUP_raw <- predictions_df$blups_FR

df_summary <- predictions_df %>%
  group_by(ID) %>%
  summarise(mean     = mean(prediction, na.rm = TRUE),
            ci_lower = quantile(prediction, 0.025),
            ci_upper = quantile(prediction, 0.975),
            BLUP_raw = mean(BLUP_raw))

dr          <- FR %>% select(ID, BLUP, Survival)
df_summary1 <- left_join(df_summary, dr)
df_summary1$BLUP_r <- df_summary1$BLUP

# Plot: predictability (shrub rIIV) vs. survival probability - transients
require(ggplot2)
p33 <- ggplot(df_summary1, aes(x = BLUP_raw, y = mean)) +
  geom_line(color="#3E4A4F", linewidth=5) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper),fill="#3E4A4F",linetype=0, alpha = 0.25, linewidth=0) +
  geom_point(aes(x=BLUP_raw, y= Survival), color="#3E4A4F",alpha=0.6,size=8)+ theme_turkey()+
  xlab("Predictability (rIIV)\nDistance to shrub landcover")+
  ylab("Survival Probaility")+
  theme_turkey()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14),plot.margin = margin(r=20, l=15))+
  scale_y_continuous(n.breaks = 5)+
  ggtitle(label ="", subtitle = "") +
  theme(legend.title = element_text(size=74, face="bold"), legend.position = "bottom") +
  theme(axis.title = element_text(size = 70, face = "bold"),
        plot.title = element_text(size=72),
        legend.text = element_blank(),
        legend.position = "none",
        axis.ticks.length = unit(.5, "cm"),
        axis.ticks = element_line(linewidth = 3, color = "black"),
        panel.border = element_rect(colour ="black", fill=NA, size=2.05),
        axis.text.y   = element_text(size=68,colour = "black"),
        axis.text.x   = element_text(size=68, colour = "black"),
        panel.grid.major = element_blank(),
        plot.margin = margin(r=55, l=20, b=35, t=35),
        panel.background = element_blank(),
        panel.grid.minor = element_blank())
p33


# ------------------------------------------------------------------------------
# 9C. Shrub - Residents - Behavioral Type
# ------------------------------------------------------------------------------

require(dplyr)
variables(shrub_brm)
draws_df <- brms::as_draws_df(shrub_brm)

# Extract behavioral type BLUPs
blup_cols <- grep("^r_ID\\[", names(draws_df), value = TRUE)
meta_cols  <- intersect(c(".chain", ".draw", ".iteration"), names(draws_df))

o1 <- draws_df %>%
  dplyr::select(dplyr::all_of(c(blup_cols, meta_cols))) %>%
  dplyr::slice_sample(n = 2000) %>%
  dplyr::select(-dplyr::any_of(c(".chain", ".draw", ".iteration")))

# Load and filter survival data to residents
sur            <- read_csv("survival_meta_02_17_2026.csv")
survival_final <- sur %>%
  filter(Status_clean %in% c("resident", "transient")) %>%
  group_by(COD) %>%
  filter(!COD %in% c("Unknown", "Vehicle", "Car", "Distemper"))

require(bayestestR)
s_ga <- survival_final %>%
  group_by(Status_clean) %>%
  filter(Status_clean != "transient")

# Remove individual with incomplete data
s_ga <- s_ga %>% filter(ID2 != "AL45F")
table(s_ga$Survival)

# Bayesian model-averaging loop
FR2_lm_list <- list()
require(data.table)
require(tidyr)
require(stringr)

for (iter in seq_len(nrow(o1))) {
  blups_FR2 <- select(as_tibble(o1), contains("r_ID"))[iter, ]
  blups_FR  <- tibble(
    ID       = str_remove(colnames(blups_FR2), "r_ID"),
    blups_FR = as.numeric(blups_FR2))
  data_FR2_lm         <- merge(blups_FR, s_ga, by = "ID")
  FR2_lm_list[[iter]] <- glm(Survival ~ blups_FR, data = data_FR2_lm, family = binomial)
}

# Extract log-odds ratios
harv_ga_da <- as.mcmc(
  sapply(FR2_lm_list, function(x) { summary(x)$coefficients[2, 1] })
)

plot(harv_ga_da)
harv_ga_da_sur <- data.frame(harv_ga_da)
describe_posterior(harv_ga_da_sur)

# Extract mean BLUP per individual for plotting
library(data.table)
FR <- select(as_tibble(o1), contains("r_ID"))

FR_long     <- FR %>% pivot_longer(everything(), names_to = c("ID")) %>% tibble()
FR_long$ID  <- str_remove_all(FR_long$ID, "r_ID")
FR_long     <- tibble(FR_long)

FR <- FR_long %>%
  group_by(ID) %>%
  summarise(BLUP = mean(value))

FR <- left_join(s_ga, FR)

# Generate predicted survival curves and summarize
predictions_df <- data.frame()

for (i in seq_along(FR2_lm_list)) {
  model                   <- FR2_lm_list[[i]]
  data_FR2_lm$prediction  <- predict(model, newdata = data_FR2_lm, type = "response")
  predictions_df          <- rbind(predictions_df, data.frame(model = rep(i, nrow(data_FR2_lm)), data_FR2_lm))
}

predictions_df$BLUP_raw <- predictions_df$blups_FR

df_summary <- predictions_df %>%
  group_by(ID) %>%
  summarise(mean     = mean(prediction, na.rm = TRUE),
            ci_lower = quantile(prediction, 0.025),
            ci_upper = quantile(prediction, 0.975),
            BLUP_raw = mean(BLUP_raw))

dr          <- FR %>% select(ID, BLUP, Survival)
df_summary1 <- left_join(df_summary, dr)
df_summary1$BLUP_r <- df_summary1$BLUP

# Back-transform to original scale
df_summary1$BLUP_raw1 <- df_summary1$BLUP_raw * 3.37 + 10.3
df_summary1$BLUP_raw2 <- df_summary1$BLUP_raw1^2
df_summary1$BLUP      <- df_summary1$BLUP * 3.37 + 10.3
df_summary1$BLUP      <- df_summary1$BLUP^2

# Plot: behavioral type (shrub) vs. survival probability - residents
require(ggplot2)
p33 <- ggplot(df_summary1, aes(x = BLUP_raw2, y = mean)) +
  geom_line(color="#A98F64", linewidth=5) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper),fill="#A98F64",linetype=0, alpha = 0.25, linewidth=0) +
  geom_point(aes(x=BLUP_raw2, y= Survival), color="#A98F64",alpha=0.6,size=8)+ theme_turkey()+
  xlab("Distance to shrub landcover (m)")+
  ylab("Survival Probaility")+
  theme_turkey()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14),plot.margin = margin(r=20, l=15))+
  scale_y_continuous(n.breaks = 5)+
  ggtitle(label ="", subtitle = "") +
  theme(legend.title = element_text(size=74, face="bold"), legend.position = "bottom") +
  theme(axis.title = element_text(size = 70, face = "bold"),
        plot.title = element_text(size=72),
        legend.text = element_blank(),
        legend.position = "none",
        axis.ticks.length = unit(.5, "cm"),
        axis.ticks = element_line(linewidth = 3, color = "black"),
        panel.border = element_rect(colour ="black", fill=NA, size=2.05),
        axis.text.y   = element_text(size=68,colour = "black"),
        axis.text.x   = element_text(size=68, colour = "black"),
        panel.grid.major = element_blank(),
        plot.margin = margin(r=55, l=20, b=35, t=35),
        panel.background = element_blank(),
        panel.grid.minor = element_blank())
p33


# ------------------------------------------------------------------------------
# 9D. Shrub - Residents - Predictability (rIIV)
#     Note: quadratic term included to account for potential nonlinear
#     relationship between predictability and survival
# ------------------------------------------------------------------------------

require(dplyr)
variables(shrub_brm)
draws_df <- brms::as_draws_df(shrub_brm)

# Extract rIIV BLUPs from dispersion submodel
blup_cols <- grep("^r_ID__sigma\\[", names(draws_df), value = TRUE)
meta_cols  <- intersect(c(".chain", ".draw", ".iteration"), names(draws_df))

o1 <- draws_df %>%
  dplyr::select(dplyr::all_of(c(blup_cols, meta_cols))) %>%
  dplyr::slice_sample(n = 2000) %>%
  dplyr::select(-dplyr::any_of(c(".chain", ".draw", ".iteration")))

# Load and filter survival data to residents
sur            <- read_csv("survival_meta_02_17_2026.csv")
survival_final <- sur %>%
  filter(Status_clean %in% c("resident", "transient")) %>%
  group_by(COD) %>%
  filter(!COD %in% c("Unknown", "Vehicle", "Car", "Distemper"))

require(bayestestR)
s_ga <- survival_final %>%
  group_by(Status_clean) %>%
  filter(Status_clean != "transient")
s_ga <- s_ga %>% filter(ID2 != "AL45F")

# Bayesian model-averaging loop with quadratic term
FR2_lm_list <- list()
require(data.table)
require(tidyr)
require(stringr)

for (iter in seq_len(nrow(o1))) {
  blups_FR2 <- select(as_tibble(o1), contains("r_ID__sigma"))[iter, ]
  blups_FR  <- tibble(
    ID       = str_remove(colnames(blups_FR2), "r_ID__sigma"),
    blups_FR = as.numeric(blups_FR2))
  data_FR2_lm         <- merge(blups_FR, s_ga, by = "ID")
  FR2_lm_list[[iter]] <- glm(Survival ~ blups_FR + I(blups_FR^2), data = data_FR2_lm, family = binomial)
}

# Extract quadratic coefficient (index 3) across all GLMs
harv_ga_da <- as.mcmc(
  sapply(FR2_lm_list, function(x) { summary(x)$coefficients[3, 1] })
)

plot(harv_ga_da)
harv_ga_da_sur <- data.frame(harv_ga_da)
describe_posterior(harv_ga_da_sur)

# Build mean BLUP per individual and predict nonlinear survival curve
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(ggplot2)

FR <- as_tibble(o1) %>%
  select(contains("r_ID__sigma")) %>%
  pivot_longer(cols = everything(), names_to = "ID", values_to = "value") %>%
  mutate(ID = str_remove_all(ID, "r_ID__sigma")) %>%
  group_by(ID) %>%
  summarise(BLUP = mean(value, na.rm = TRUE), .groups = "drop")

FR <- left_join(s_ga, FR, by = "ID")

# Predict over BLUP grid and summarize across posterior draws
blup_grid <- data.frame(
  blups_FR = seq(min(FR$BLUP, na.rm = TRUE),
                 max(FR$BLUP, na.rm = TRUE),
                 length.out = 110)
)

pred_mat <- sapply(FR2_lm_list, function(m) {
  predict(m, newdata = blup_grid, type = "response")
})

pred_summ <- data.frame(
  BLUP_raw = blup_grid$blups_FR,
  mean     = apply(pred_mat, 1, mean, na.rm = TRUE),
  ci_lower = apply(pred_mat, 1, quantile, probs = 0.025, na.rm = TRUE),
  ci_upper = apply(pred_mat, 1, quantile, probs = 0.975, na.rm = TRUE)
)

# Plot: predictability (shrub rIIV) vs. survival probability - residents (nonlinear)
p33 <- ggplot(pred_summ, aes(x = BLUP_raw, y = mean)) +
  geom_line(color = "#A98F64", linewidth = 5) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper),
              fill = "#A98F64", linetype = 0, alpha = 0.25, linewidth = 0) +
  geom_point(data = FR, aes(x = BLUP, y = Survival),
             color = "#A98F64", alpha = 0.6, size = 8) +
  xlab("Predictability (rIIV)\nDistance to shrub landcover") +
  ylab("Survival Probability") +
  theme_turkey()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14),plot.margin = margin(r=20, l=15))+
  scale_y_continuous(n.breaks = 5)+
  ggtitle(label ="", subtitle = "") +
  theme(legend.title = element_text(size=74, face="bold"), legend.position = "bottom") +
  theme(axis.title = element_text(size = 70, face = "bold"),
        plot.title = element_text(size=72),
        legend.text = element_blank(),
        legend.position = "none",
        axis.ticks.length = unit(.5, "cm"),
        axis.ticks = element_line(linewidth = 3, color = "black"),
        panel.border = element_rect(colour ="black", fill=NA, size=2.05),
        axis.text.y   = element_text(size=68,colour = "black"),
        axis.text.x   = element_text(size=68, colour = "black"),
        panel.grid.major = element_blank(),
        plot.margin = margin(r=55, l=20, b=35, t=35),
        panel.background = element_blank(),
        panel.grid.minor = element_blank())
p33


# ------------------------------------------------------------------------------
# 9E. Hardwood - Transients - Behavioral Type
# ------------------------------------------------------------------------------

require(dplyr)
variables(hardwood_brm)
draws_df <- brms::as_draws_df(hardwood_brm)

# Extract behavioral type BLUPs
blup_cols <- grep("^r_ID\\[", names(draws_df), value = TRUE)
meta_cols  <- intersect(c(".chain", ".draw", ".iteration"), names(draws_df))

o1 <- draws_df %>%
  dplyr::select(dplyr::all_of(c(blup_cols, meta_cols))) %>%
  dplyr::slice_sample(n = 2000) %>%
  dplyr::select(-dplyr::any_of(c(".chain", ".draw", ".iteration")))

# Load and filter survival data to transients
sur            <- read_csv("survival_meta_02_17_2026.csv")
survival_final <- sur %>%
  filter(Status_clean %in% c("resident", "transient")) %>%
  group_by(COD) %>%
  filter(!COD %in% c("Unknown", "Vehicle", "Car", "Distemper"))

require(bayestestR)
s_ga <- survival_final %>%
  group_by(Status_clean) %>%
  filter(Status_clean != "resident")
s_ga <- s_ga %>% filter(ID2 != "AL45F")

# Bayesian model-averaging loop
FR2_lm_list <- list()
require(data.table)
require(tidyr)
require(stringr)

for (iter in seq_len(nrow(o1))) {
  blups_FR2 <- select(as_tibble(o1), contains("r_ID"))[iter, ]
  blups_FR  <- tibble(
    ID       = str_remove(colnames(blups_FR2), "r_ID"),
    blups_FR = as.numeric(blups_FR2))
  data_FR2_lm         <- merge(blups_FR, s_ga, by = "ID")
  FR2_lm_list[[iter]] <- glm(Survival ~ blups_FR, data = data_FR2_lm, family = binomial)
}

# Extract log-odds ratios
harv_ga_da <- as.mcmc(
  sapply(FR2_lm_list, function(x) { summary(x)$coefficients[2, 1] })
)

plot(harv_ga_da)
harv_ga_da_sur <- data.frame(harv_ga_da)
describe_posterior(harv_ga_da_sur)

# Extract mean BLUP per individual for plotting
library(data.table)
FR <- select(as_tibble(o1), contains("r_ID"))

FR_long     <- FR %>% pivot_longer(everything(), names_to = c("ID")) %>% tibble()
FR_long$ID  <- str_remove_all(FR_long$ID, "r_ID")
FR_long     <- tibble(FR_long)

FR <- FR_long %>%
  group_by(ID) %>%
  summarise(BLUP = mean(value))

FR <- left_join(s_ga, FR)

# Generate predicted survival curves and summarize
predictions_df <- data.frame()

for (i in seq_along(FR2_lm_list)) {
  model                   <- FR2_lm_list[[i]]
  data_FR2_lm$prediction  <- predict(model, newdata = data_FR2_lm, type = "response")
  predictions_df          <- rbind(predictions_df, data.frame(model = rep(i, nrow(data_FR2_lm)), data_FR2_lm))
}

predictions_df$BLUP_raw <- predictions_df$blups_FR

df_summary <- predictions_df %>%
  group_by(ID) %>%
  summarise(mean     = mean(prediction, na.rm = TRUE),
            ci_lower = quantile(prediction, 0.025),
            ci_upper = quantile(prediction, 0.975),
            BLUP_raw = mean(BLUP_raw))

dr          <- FR %>% select(ID, BLUP, Survival)
df_summary1 <- left_join(df_summary, dr)
df_summary1$BLUP_r <- df_summary1$BLUP

# Back-transform to original scale
df_summary1$BLUP_raw1 <- df_summary1$BLUP_raw * 3.24 + 7.93
df_summary1$BLUP_raw2 <- df_summary1$BLUP_raw1^2
df_summary1$BLUP      <- df_summary1$BLUP * 3.24 + 7.93
df_summary1$BLUP      <- df_summary1$BLUP^2

# Plot: behavioral type (hardwood) vs. survival probability - transients
require(ggplot2)
p33 <- ggplot(df_summary1, aes(x = BLUP_raw2, y = mean)) +
  geom_line(color="#3E4A4F", linewidth=5) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper),fill="#3E4A4F",linetype=0, alpha = 0.25, linewidth=0) +
  geom_point(aes(x=BLUP_raw2, y= Survival), color="#3E4A4F",alpha=0.6,size=8)+ theme_turkey()+
  xlab("Distance to hardwood landcover (m)")+
  ylab("Survival Probaility")+
  theme_turkey()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14),plot.margin = margin(r=20, l=15))+
  scale_y_continuous(n.breaks = 5)+scale_x_continuous(breaks =c(30,60,90, 120))+
  ggtitle(label ="", subtitle = "") +
  theme(legend.title = element_text(size=74, face="bold"), legend.position = "bottom") +
  theme(axis.title = element_text(size = 70, face = "bold"),
        plot.title = element_text(size=72),
        legend.text = element_blank(),
        legend.position = "none",
        axis.ticks.length = unit(.5, "cm"),
        axis.ticks = element_line(linewidth = 3, color = "black"),
        panel.border = element_rect(colour ="black", fill=NA, size=2.05),
        axis.text.y   = element_text(size=68,colour = "black"),
        axis.text.x   = element_text(size=68, colour = "black"),
        panel.grid.major = element_blank(),
        plot.margin = margin(r=55, l=20, b=35, t=35),
        panel.background = element_blank(),
        panel.grid.minor = element_blank())
p33


# ------------------------------------------------------------------------------
# 9F. Hardwood - Transients - Predictability (rIIV)
# ------------------------------------------------------------------------------

require(dplyr)
variables(hardwood_brm)
draws_df <- brms::as_draws_df(hardwood_brm)

# Extract rIIV BLUPs from dispersion submodel
blup_cols <- grep("^r_ID__sigma\\[", names(draws_df), value = TRUE)
meta_cols  <- intersect(c(".chain", ".draw", ".iteration"), names(draws_df))

o1 <- draws_df %>%
  dplyr::select(dplyr::all_of(c(blup_cols, meta_cols))) %>%
  dplyr::slice_sample(n = 2000) %>%
  dplyr::select(-dplyr::any_of(c(".chain", ".draw", ".iteration")))

# Load and filter survival data to transients
sur            <- read_csv("survival_meta_02_17_2026.csv")
survival_final <- sur %>%
  filter(Status_clean %in% c("resident", "transient")) %>%
  group_by(COD) %>%
  filter(!COD %in% c("Unknown", "Vehicle", "Car", "Distemper"))

require(bayestestR)
s_ga <- survival_final %>%
  group_by(Status_clean) %>%
  filter(Status_clean != "resident")

# Bayesian model loop
FR2_lm_list <- list()
require(data.table)
require(tidyr)
require(stringr)

for (iter in seq_len(nrow(o1))) {
  blups_FR2 <- select(as_tibble(o1), contains("r_ID__sigma"))[iter, ]
  blups_FR  <- tibble(
    ID       = str_remove(colnames(blups_FR2), "r_ID__sigma"),
    blups_FR = as.numeric(blups_FR2))
  data_FR2_lm         <- merge(blups_FR, s_ga, by = "ID")
  FR2_lm_list[[iter]] <- glm(Survival ~ blups_FR, data = data_FR2_lm, family = binomial)
}

# Extract log-odds ratios
harv_ga_da <- as.mcmc(
  sapply(FR2_lm_list, function(x) { summary(x)$coefficients[2, 1] })
)

plot(harv_ga_da)
harv_ga_da_sur <- data.frame(harv_ga_da)
describe_posterior(harv_ga_da_sur)

# Extract mean BLUP per individual for plotting
library(data.table)
FR <- select(as_tibble(o1), contains("r_ID__sigma"))

FR_long     <- FR %>% pivot_longer(everything(), names_to = c("ID")) %>% tibble()
FR_long$ID  <- str_remove_all(FR_long$ID, "r_ID__sigma")
FR_long     <- tibble(FR_long)

FR <- FR_long %>%
  group_by(ID) %>%
  summarise(BLUP = mean(value))

FR <- left_join(s_ga, FR)

# Generate predicted survival curves and summarize
predictions_df <- data.frame()

for (i in seq_along(FR2_lm_list)) {
  model                   <- FR2_lm_list[[i]]
  data_FR2_lm$prediction  <- predict(model, newdata = data_FR2_lm, type = "response")
  predictions_df          <- rbind(predictions_df, data.frame(model = rep(i, nrow(data_FR2_lm)), data_FR2_lm))
}

predictions_df$BLUP_raw <- predictions_df$blups_FR

df_summary <- predictions_df %>%
  group_by(ID) %>%
  summarise(mean     = mean(prediction, na.rm = TRUE),
            ci_lower = quantile(prediction, 0.025),
            ci_upper = quantile(prediction, 0.975),
            BLUP_raw = mean(BLUP_raw))

dr          <- FR %>% select(ID, BLUP, Survival)
df_summary1 <- left_join(df_summary, dr)
df_summary1$BLUP_r <- df_summary1$BLUP

# Plot: predictability (hardwood rIIV) vs. survival probability - transients
require(ggplot2)
p33 <- ggplot(df_summary1, aes(x = BLUP_raw, y = mean)) +
  geom_line(color="#3E4A4F", linewidth=5) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper),fill="#3E4A4F",linetype=0, alpha = 0.25, linewidth=0) +
  geom_point(aes(x=BLUP_raw, y= Survival), color="#3E4A4F",alpha=0.6,size=8)+ theme_turkey()+
  xlab("Predictability (rIIV)\nDistance to hardwood landcover")+
  ylab("Survival Probaility")+
  theme_turkey()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14),plot.margin = margin(r=20, l=15))+
  scale_y_continuous(n.breaks = 5)+
  ggtitle(label ="", subtitle = "") +
  theme(legend.title = element_text(size=74, face="bold"), legend.position = "bottom") +
  theme(axis.title = element_text(size = 70, face = "bold"),
        plot.title = element_text(size=72),
        legend.text = element_blank(),
        legend.position = "none",
        axis.ticks.length = unit(.5, "cm"),
        axis.ticks = element_line(linewidth = 3, color = "black"),
        panel.border = element_rect(colour ="black", fill=NA, size=2.05),
        axis.text.y   = element_text(size=68,colour = "black"),
        axis.text.x   = element_text(size=68, colour = "black"),
        panel.grid.major = element_blank(),
        plot.margin = margin(r=55, l=20, b=35, t=35),
        panel.background = element_blank(),
        panel.grid.minor = element_blank())
p33


# ------------------------------------------------------------------------------
# 9G. Hardwood - Residents - Behavioral Type
# ------------------------------------------------------------------------------

require(dplyr)
variables(hardwood_brm)
draws_df <- brms::as_draws_df(hardwood_brm)

# Extract behavioral type BLUPs
blup_cols <- grep("^r_ID\\[", names(draws_df), value = TRUE)
meta_cols  <- intersect(c(".chain", ".draw", ".iteration"), names(draws_df))

o1 <- draws_df %>%
  dplyr::select(dplyr::all_of(c(blup_cols, meta_cols))) %>%
  dplyr::slice_sample(n = 2000) %>%
  dplyr::select(-dplyr::any_of(c(".chain", ".draw", ".iteration")))

# Load and filter survival data to residents
sur            <- read_csv("survival_meta_02_17_2026.csv")
survival_final <- sur %>%
  filter(Status_clean %in% c("resident", "transient")) %>%
  group_by(COD) %>%
  filter(!COD %in% c("Unknown", "Vehicle", "Car", "Distemper"))

require(bayestestR)
s_ga <- survival_final %>%
  group_by(Status_clean) %>%
  filter(Status_clean != "transient")
s_ga <- s_ga %>% filter(ID2 != "AL45F")

# Bayesian model-averaging loop
FR2_lm_list <- list()
require(data.table)
require(tidyr)
require(stringr)

for (iter in seq_len(nrow(o1))) {
  blups_FR2 <- select(as_tibble(o1), contains("r_ID"))[iter, ]
  blups_FR  <- tibble(
    ID       = str_remove(colnames(blups_FR2), "r_ID"),
    blups_FR = as.numeric(blups_FR2))
  data_FR2_lm         <- merge(blups_FR, s_ga, by = "ID")
  FR2_lm_list[[iter]] <- glm(Survival ~ blups_FR, data = data_FR2_lm, family = binomial)
}

# Extract log-odds ratios
harv_ga_da <- as.mcmc(
  sapply(FR2_lm_list, function(x) { summary(x)$coefficients[2, 1] })
)

plot(harv_ga_da)
harv_ga_da_sur <- data.frame(harv_ga_da)
describe_posterior(harv_ga_da_sur)

# Extract mean BLUP per individual for plotting
library(data.table)
FR <- select(as_tibble(o1), contains("r_ID"))

FR_long     <- FR %>% pivot_longer(everything(), names_to = c("ID")) %>% tibble()
FR_long$ID  <- str_remove_all(FR_long$ID, "r_ID")
FR_long     <- tibble(FR_long)

FR <- FR_long %>%
  group_by(ID) %>%
  summarise(BLUP = mean(value))

FR <- left_join(s_ga, FR)

# Generate predicted survival curves and summarize
predictions_df <- data.frame()

for (i in seq_along(FR2_lm_list)) {
  model                   <- FR2_lm_list[[i]]
  data_FR2_lm$prediction  <- predict(model, newdata = data_FR2_lm, type = "response")
  predictions_df          <- rbind(predictions_df, data.frame(model = rep(i, nrow(data_FR2_lm)), data_FR2_lm))
}

predictions_df$BLUP_raw <- predictions_df$blups_FR

df_summary <- predictions_df %>%
  group_by(ID) %>%
  summarise(mean     = mean(prediction, na.rm = TRUE),
            ci_lower = quantile(prediction, 0.025),
            ci_upper = quantile(prediction, 0.975),
            BLUP_raw = mean(BLUP_raw))

dr          <- FR %>% select(ID, BLUP, Survival)
df_summary1 <- left_join(df_summary, dr)
df_summary1$BLUP_r <- df_summary1$BLUP

# Back-transform to original scale
df_summary1$BLUP_raw1 <- df_summary1$BLUP_raw * 3.24 + 7.93
df_summary1$BLUP_raw2 <- df_summary1$BLUP_raw1^2
df_summary1$BLUP      <- df_summary1$BLUP * 3.24 + 7.93
df_summary1$BLUP      <- df_summary1$BLUP^2

# Plot: behavioral type (hardwood) vs. survival probability - residents
require(ggplot2)
p33 <- ggplot(df_summary1, aes(x = BLUP_raw2, y = mean)) +
  geom_line(color="#A98F64", linewidth=5) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper),fill="#A98F64",linetype=0, alpha = 0.25, linewidth=0) +
  geom_point(aes(x=BLUP_raw2, y= Survival), color="#A98F64",alpha=0.6,size=8)+ theme_turkey()+
  xlab("Distance to hardwood landcover (m)")+
  ylab("Survival Probaility")+
  theme_turkey()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14),plot.margin = margin(r=20, l=15))+
  scale_y_continuous(n.breaks = 5)+
  ggtitle(label ="", subtitle = "") +
  theme(legend.title = element_text(size=74, face="bold"), legend.position = "bottom") +
  theme(axis.title = element_text(size = 70, face = "bold"),
        plot.title = element_text(size=72),
        legend.text = element_blank(),
        legend.position = "none",
        axis.ticks.length = unit(.5, "cm"),
        axis.ticks = element_line(linewidth = 3, color = "black"),
        panel.border = element_rect(colour ="black", fill=NA, size=2.05),
        axis.text.y   = element_text(size=68,colour = "black"),
        axis.text.x   = element_text(size=68, colour = "black"),
        panel.grid.major = element_blank(),
        plot.margin = margin(r=55, l=20, b=35, t=35),
        panel.background = element_blank(),
        panel.grid.minor = element_blank())
p33


# ------------------------------------------------------------------------------
# 9H. Hardwood - Residents - Predictability (rIIV)
#     Note: quadratic term included to account for potential nonlinear
#     relationship between predictability and survival
# ------------------------------------------------------------------------------

require(dplyr)
variables(hardwood_brm)
draws_df <- brms::as_draws_df(hardwood_brm)

# Extract rIIV BLUPs from dispersion submodel
blup_cols <- grep("^r_ID__sigma\\[", names(draws_df), value = TRUE)
meta_cols  <- intersect(c(".chain", ".draw", ".iteration"), names(draws_df))

o1 <- draws_df %>%
  dplyr::select(dplyr::all_of(c(blup_cols, meta_cols))) %>%
  dplyr::slice_sample(n = 2000) %>%
  dplyr::select(-dplyr::any_of(c(".chain", ".draw", ".iteration")))

# Load and filter survival data to residents
sur            <- read_csv("survival_meta_02_17_2026.csv")
survival_final <- sur %>%
  filter(Status_clean %in% c("resident", "transient")) %>%
  group_by(COD) %>%
  filter(!COD %in% c("Unknown", "Vehicle", "Car", "Distemper"))

require(bayestestR)
s_ga <- survival_final %>%
  group_by(Status_clean) %>%
  filter(Status_clean != "transient")

# Bayesian model-averaging loop with quadratic term
FR2_lm_list <- list()
require(data.table)
require(tidyr)
require(stringr)

for (iter in seq_len(nrow(o1))) {
  blups_FR2 <- select(as_tibble(o1), contains("r_ID__sigma"))[iter, ]
  blups_FR  <- tibble(
    ID       = str_remove(colnames(blups_FR2), "r_ID__sigma"),
    blups_FR = as.numeric(blups_FR2))
  data_FR2_lm         <- merge(blups_FR, s_ga, by = "ID")
  FR2_lm_list[[iter]] <- glm(Survival ~ blups_FR + I(blups_FR^2), data = data_FR2_lm, family = binomial)
}

# Extract quadratic coefficient (index 3) across all GLMs
harv_ga_da <- as.mcmc(
  sapply(FR2_lm_list, function(x) { summary(x)$coefficients[3, 1] })
)

plot(harv_ga_da)
harv_ga_da_sur <- data.frame(harv_ga_da)
describe_posterior(harv_ga_da_sur)

# Build mean BLUP per individual and predict nonlinear survival curve
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(ggplot2)

FR <- as_tibble(o1) %>%
  select(contains("r_ID__sigma")) %>%
  pivot_longer(cols = everything(), names_to = "ID", values_to = "value") %>%
  mutate(ID = str_remove_all(ID, "r_ID__sigma")) %>%
  group_by(ID) %>%
  summarise(BLUP = mean(value, na.rm = TRUE), .groups = "drop")

FR <- left_join(s_ga, FR, by = "ID")

# Predict over BLUP grid and summarize across posterior draws
blup_grid <- data.frame(
  blups_FR = seq(min(FR$BLUP, na.rm = TRUE),
                 max(FR$BLUP, na.rm = TRUE),
                 length.out = 110)
)

pred_mat <- sapply(FR2_lm_list, function(m) {
  predict(m, newdata = blup_grid, type = "response")
})

pred_summ <- data.frame(
  BLUP_raw = blup_grid$blups_FR,
  mean     = apply(pred_mat, 1, mean, na.rm = TRUE),
  ci_lower = apply(pred_mat, 1, quantile, probs = 0.025, na.rm = TRUE),
  ci_upper = apply(pred_mat, 1, quantile, probs = 0.975, na.rm = TRUE)
)

# Plot: predictability (hardwood rIIV) vs. survival probability - residents (nonlinear)
p33 <- ggplot(pred_summ, aes(x = BLUP_raw, y = mean)) +
  geom_line(color = "#A98F64", linewidth = 5) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper),
              fill = "#A98F64", linetype = 0, alpha = 0.25, linewidth = 0) +
  geom_point(data = FR, aes(x = BLUP, y = Survival),
             color = "#A98F64", alpha = 0.6, size = 8) +
  xlab("Predictability (rIIV)\nDistance to hardwood landcover") +
  ylab("Survival Probability") +
  theme_turkey()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14),plot.margin = margin(r=20, l=15))+
  scale_y_continuous(n.breaks = 5)+
  ggtitle(label ="", subtitle = "") +
  theme(legend.title = element_text(size=74, face="bold"), legend.position = "bottom") +
  theme(axis.title = element_text(size = 70, face = "bold"),
        plot.title = element_text(size=72),
        legend.text = element_blank(),
        legend.position = "none",
        axis.ticks.length = unit(.5, "cm"),
        axis.ticks = element_line(linewidth = 3, color = "black"),
        panel.border = element_rect(colour ="black", fill=NA, size=2.05),
        axis.text.y   = element_text(size=68,colour = "black"),
        axis.text.x   = element_text(size=68, colour = "black"),
        panel.grid.major = element_blank(),
        plot.margin = margin(r=55, l=20, b=35, t=35),
        panel.background = element_blank(),
        panel.grid.minor = element_blank())
p33

################################################################################
# End of Script
################################################################################
