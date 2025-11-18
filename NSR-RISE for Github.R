############################################################
# Biomarker–Epilepsy Analysis 
# Code for data availability / reproducibility
# Assumes a de-identified, analysis-ready dataset
############################################################

## 0. Packages ----
library(dplyr)
library(tidyr)
library(ggplot2)
library(tableone)
library(sandwich)
library(lmtest)

## 1. Import analysis dataset ----
# One row per participant; wide-format biomarkers (e.g., bioX_T1, bioX_T2, bioX_T3)
# Replace with actual file name / path in your repository.
dat <- read.csv("analytic_biomarker_dataset.csv")

## Expected minimal variables:
# id                    : participant ID
# epilepsy_24m          : 0/1 outcome (epilepsy by ~24 months)
# eeg_days_ge3          : 0/1 covariate (≥3 days of EEG seizures)
# worst_eeg_severe      : 0/1 covariate (severely abnormal EEG background)
# abnl_neuro_exam       : 0/1 covariate (abnormal neurologic exam at discharge)
# biomarker_*           : numeric biomarker concentrations (T1–T3, or single timepoint)

## 2. Define biomarker variables ----
# Here we assume all biomarker columns start with "bio_".
# Modify the pattern to match your variable naming.
biomarker_vars <- grep("^bio_", names(dat), value = TRUE)

# Optional: if you have separate timepoints, you can keep them all
# and use a "sample" indicator in long format (see below).

## 3. Basic descriptive table for clinical variables ----
vars_table <- c("eeg_days_ge3", "worst_eeg_severe", "abnl_neuro_exam")
dat$epilepsy_24m_factor <- factor(dat$epilepsy_24m, labels = c("No Epilepsy", "Epilepsy"))

tbl1 <- CreateTableOne(
  data  = dat,
  vars  = vars_table,
  strata = "epilepsy_24m_factor",
  test  = TRUE
)
print(tbl1)

## 4. Optional: reshape biomarkers to long format ----
# (useful for plotting trajectories over time)
# Here we assume variables like bio_GFAP_T1, bio_GFAP_T2, bio_GFAP_T3.
# Adjust the separator / patterns to match your dataset.
dat_long <- dat |>
  pivot_longer(
    cols = starts_with("bio_"),
    names_to = c("biomarker", "timepoint"),
    names_sep = "_T",
    values_to = "value"
  )

# Example plot for one biomarker:
ggplot(
  dat_long |> filter(biomarker == "bio_GFAP"),
  aes(x = timepoint, y = value, color = epilepsy_24m_factor, group = interaction(epilepsy_24m_factor, id))
) +
  geom_line(alpha = 0.3) +
  stat_summary(fun = mean, geom = "line", size = 1.2, aes(group = epilepsy_24m_factor)) +
  labs(x = "Timepoint", y = "Biomarker concentration", color = "Epilepsy at 24 mo")

## 5. Helper: robust Poisson regression per biomarker ----
# We use a log-link Poisson model with robust standard errors to estimate
# relative risk of epilepsy per log2 unit increase in biomarker.

fit_biomarker <- function(bio_var, data) {
  df <- data
  
  # log2-transform; add small constant if needed to avoid log(0)
  df$bio_log2 <- log2(df[[bio_var]] + 1e-6)
  
  # unadjusted model
  m_unadj <- glm(epilepsy_24m ~ bio_log2, family = poisson, data = df)
  vc_unadj <- vcovHC(m_unadj, type = "HC0")
  ct_unadj <- coeftest(m_unadj, vcov = vc_unadj)["bio_log2", ]
  
  rr_unadj  <- exp(ct_unadj["Estimate"])
  se_unadj  <- ct_unadj["Std. Error"]
  rr_l_unadj <- exp(ct_unadj["Estimate"] - qnorm(0.975) * se_unadj)
  rr_u_unadj <- exp(ct_unadj["Estimate"] + qnorm(0.975) * se_unadj)
  p_unadj    <- ct_unadj["Pr(>|z|)"]
  
  # adjusted model (example covariates)
  m_adj <- glm(
    epilepsy_24m ~ eeg_days_ge3 + worst_eeg_severe + abnl_neuro_exam + bio_log2,
    family = poisson,
    data = df
  )
  vc_adj <- vcovHC(m_adj, type = "HC0")
  ct_adj <- coeftest(m_adj, vcov = vc_adj)["bio_log2", ]
  
  rr_adj  <- exp(ct_adj["Estimate"])
  se_adj  <- ct_adj["Std. Error"]
  rr_l_adj <- exp(ct_adj["Estimate"] - qnorm(0.975) * se_adj)
  rr_u_adj <- exp(ct_adj["Estimate"] + qnorm(0.975) * se_adj)
  p_adj    <- ct_adj["Pr(>|z|)"]
  
  tibble::tibble(
    biomarker = bio_var,
    rr_unadj  = rr_unadj,
    ci_l_unadj = rr_l_unadj,
    ci_u_unadj = rr_u_unadj,
    p_unadj   = p_unadj,
    rr_adj    = rr_adj,
    ci_l_adj  = rr_l_adj,
    ci_u_adj  = rr_u_adj,
    p_adj     = p_adj
  )
}

## 6. Run models for all biomarkers and apply FDR correction ----
library(purrr)
library(tibble)

results <- map_dfr(biomarker_vars, fit_biomarker, data = dat) |>
  mutate(
    p_unadj_fdr = p.adjust(p_unadj, method = "fdr"),
    p_adj_fdr   = p.adjust(p_adj,   method = "fdr")
  )

print(results)

## 7. Save results table for the manuscript ----
write.csv(results, "biomarker_epilepsy_results.csv", row.names = FALSE)

## 8. (Optional) simple forest plot for a subset of biomarkers ----
top_biomarkers <- results |>
  arrange(p_adj_fdr) |>
  slice(1:10)

ggplot(top_biomarkers,
       aes(x = biomarker, y = rr_adj, ymin = ci_l_adj, ymax = ci_u_adj)) +
  geom_pointrange() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  coord_flip() +
  labs(
    x = "Biomarker",
    y = "Adjusted RR for epilepsy (95% CI)",
    title = "Top biomarkers associated with epilepsy at ~24 months"
  )