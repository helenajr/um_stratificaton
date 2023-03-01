# Load libraries --------------------------------
library(tidyverse)
library(lubridate)
library(plotROC)
library(patchwork)
library(ROCit)
library(scales)
library(here)
library(PropCIs)
source(here("scripts", "Functions.R"))

# At this point LUMPO scores have been calculated, added and the dataset has been filtered to only
# include those wth 5 years of follow up or outcome within 5 years

# Load datasets -------------------------------------
roc_data <- readRDS(here("data", "roc_data.RDS"))

str(roc_data)

roc_data <- roc_data %>%
  mutate(
    area = lbd * uh,
    score = score / 100,
    with3 = as_factor(if_else(is.na(chr3), "No result", "Chr3 result")),
    m3_str = if_else(chr3 == "L" | is.na(chr3), "Surveillance", "No surveillance"),
    stage_2a = if_else(stage %in% c("I"), "No surveillance", "Surveillance"),
    stage_2b = if_else(stage %in% c("I", "IIA"), "No surveillance", "Surveillance"),
    stage_3a = if_else(stage %in% c("I", "IIA", "IIB"), "No surveillance", "Surveillance")
  )

roc_data2 <- roc_data %>%
  filter(with3 == "Chr3 result")

roc_data3 <- roc_data %>%
  filter(with3 != "Chr3 result")

roc_data_ngs <- readRDS(here("data", "roc_data_ngs.RDS")) %>% # ngs = no genetics service
  mutate(score = score / 100)

# Prevalence of outcome (outcome == "Mets") -----------------------------------
prevalence <- calc_prevalence(roc_data)
ci_prev <- calc_ci_prev(roc_data)

prevalence2 <- calc_prevalence(roc_data2)
ci_prev2 <- calc_ci_prev(roc_data2)

prevalence3 <- calc_prevalence(roc_data3)
ci_prev3 <- calc_ci_prev(roc_data3)

# Table 2: Descriptive statistics (full dataset)----------------------------------------------------
roc_data %>%
  count(outcome)

median(roc_data$ageatpm)
range(roc_data$ageatpm)

roc_data %>%
  count(gender)

median(roc_data$lbd)
range(roc_data$lbd)

median(roc_data$uh)
range(roc_data$uh)

roc_data %>%
  count(cbi)

roc_data %>%
  count(eospread)

roc_data %>%
  count(epithelioid)

roc_data %>%
  count(loops)

roc_data %>%
  mutate(mitcount = case_when(
    mitcount <= 1 ~ "0-1",
    mitcount >= 2 & mitcount <= 3 ~ "2-3",
    mitcount >= 4 & mitcount <= 7 ~ "4-7",
    mitcount > 7 ~ "7+"
  )) %>%
  count(mitcount)

roc_data %>%
  count(chr3)

roc_data %>%
  count(chr8q)

roc_data %>%
  count(primarytreatment)

# Figure 1: Distribution plots (full dataset)----------------------------------------------------------
dist_lumpo <- roc_data %>%
  ggplot(aes(MM)) +
  stat_bin(binwidth = 0.05, boundary = 0, color = "black", fill = "black", alpha = 0.4) +
  geom_vline(xintercept = 0.1, color = "red", linetype = "dashed") +
  labs(
    x = "LUMPO 5-year MAM (binwidth = 0.05)",
    title = "A: LUMPO III"
  ) +
  theme_classic()

dist_ajcc <- roc_data %>%
  mutate(stage = factor(stage, levels = c("I", "IIA", "IIB", "IIIA", "IIIB", "IIIC"))) %>%
  ggplot(aes(stage)) +
  geom_bar(color = "purple", fill = "purple", alpha = 0.5) +
  geom_vline(xintercept = 1.5, color = "red", linetype = "dashed") +
  labs(title = "B: AJCC stage") +
  theme_classic()

dist_lpm <- roc_data %>%
  ggplot(aes(score)) +
  stat_bin(binwidth = 0.1, boundary = 0, color = "darkgoldenrod", fill = "#E69F00") +
  labs(
    x = "LPM 5-year MAM (binwidth = 0.1)",
    title = "C: LPM"
  ) +
  theme_classic()

bar_gen <- roc_data %>%
  mutate(
    chr3 = if_else(chr3 == "L", "Monosomy", "Disomy"),
    chr3 = if_else(is.na(chr3), "Unknown", chr3)
  ) %>%
  ggplot(aes(chr3)) +
  geom_bar(color = "pink3", fill = "pink") +
  labs(x = "Chromosome 3 status", title = "D: Genetics") +
  theme_classic()

dist_lumpo + dist_ajcc + dist_lpm + bar_gen

# Sensitivity and specificity point estimates (Figures 2, 3 and 5, Table 3, 5 and 6)-------------------------------------------------

# Monosomy 3 tables (full dataset and each subpopulation)
m3_tibble <- tibble(
  sensitivity = calc_sens2(roc_data, m3_str),
  specificity = calc_spec2(roc_data, m3_str),
  TP = prevalence * sensitivity,
  TN = (1 - prevalence) * specificity,
  FN = prevalence - TP,
  FP = (1 - prevalence) - TN,
  NPV = TN / (TN + FN),
  PPV = TP / (FP + TP),
  surveillance = TP + FP,
  sens_lower = calc_sens_ci2(roc_data, m3_str)[[1]][1],
  sens_upper = calc_sens_ci2(roc_data, m3_str)[[1]][2],
  spec_lower = calc_spec_ci2(roc_data, m3_str)[[1]][1],
  spec_upper = calc_spec_ci2(roc_data, m3_str)[[1]][2],
  inv_spec = 1 - specificity,
  inv_spec_upper = 1 - spec_upper,
  inv_spec_lower = 1 - spec_lower
)

m3_tibble2 <- tibble(
  sensitivity = calc_sens2(roc_data2, m3_str),
  specificity = calc_spec2(roc_data2, m3_str),
  TP = prevalence2 * sensitivity,
  TN = (1 - prevalence2) * specificity,
  FN = prevalence2 - TP,
  FP = (1 - prevalence2) - TN,
  NPV = TN / (TN + FN),
  PPV = TP / (FP + TP),
  surveillance = TP + FP,
  sens_lower = calc_sens_ci2(roc_data2, m3_str)[[1]][1],
  sens_upper = calc_sens_ci2(roc_data2, m3_str)[[1]][2],
  spec_lower = calc_spec_ci2(roc_data2, m3_str)[[1]][1],
  spec_upper = calc_spec_ci2(roc_data2, m3_str)[[1]][2]
) %>%
  mutate(
    inv_spec = 1 - specificity,
    inv_spec_upper = 1 - spec_upper,
    inv_spec_lower = 1 - spec_lower
  )

m3_tibble3 <- tibble(
  sensitivity = calc_sens2(roc_data3, m3_str),
  specificity = calc_spec2(roc_data3, m3_str),
  TP = prevalence3 * sensitivity,
  TN = (1 - prevalence3) * specificity,
  FN = prevalence3 - TP,
  FP = (1 - prevalence3) - TN,
  NPV = TN / (TN + FN),
  PPV = TP / (FP + TP),
  surveillance = TP + FP,
  sens_lower = calc_sens_ci2(roc_data3, m3_str)[[1]][1],
  sens_upper = calc_sens_ci2(roc_data3, m3_str)[[1]][2],
  spec_lower = calc_sens_ci2(roc_data3, m3_str)[[1]][1],
  spec_upper = calc_sens_ci2(roc_data3, m3_str)[[1]][2],
  inv_spec = 1 - specificity,
  inv_spec_upper = 1 - spec_upper,
  inv_spec_lower = 1 - spec_lower
)

# AJCC tables (full dataset and each subpopulation, ngs dataset same as full)
stage2a_tibble <- tibble(
  sensitivity = calc_sens2(roc_data, stage_2a),
  specificity = calc_spec2(roc_data, stage_2a),
  TP = prevalence * sensitivity,
  TN = (1 - prevalence) * specificity,
  FN = prevalence - TP,
  FP = (1 - prevalence) - TN,
  NPV = TN / (TN + FN),
  PPV = TP / (FP + TP),
  surveillance = TP + FP,
  sens_lower = calc_sens_ci2(roc_data, stage_2a)[[1]][1],
  sens_upper = calc_sens_ci2(roc_data, stage_2a)[[1]][2],
  spec_lower = calc_spec_ci2(roc_data, stage_2a)[[1]][1],
  spec_upper = calc_spec_ci2(roc_data, stage_2a)[[1]][2],
  inv_spec = 1 - specificity,
  inv_spec_upper = 1 - spec_upper,
  inv_spec_lower = 1 - spec_lower
)

stage2a_tibble2 <- tibble(
  sensitivity = calc_sens2(roc_data2, stage_2a),
  specificity = calc_spec2(roc_data2, stage_2a),
  TP = prevalence2 * sensitivity,
  TN = (1 - prevalence2) * specificity,
  FN = prevalence2 - TP,
  FP = (1 - prevalence2) - TN,
  NPV = TN / (TN + FN),
  PPV = TP / (FP + TP),
  surveillance = TP + FP,
  sens_lower = calc_sens_ci2(roc_data2, stage_2a)[[1]][1],
  sens_upper = calc_sens_ci2(roc_data2, stage_2a)[[1]][2],
  spec_lower = calc_spec_ci2(roc_data2, stage_2a)[[1]][1],
  spec_upper = calc_spec_ci2(roc_data2, stage_2a)[[1]][2],
  inv_spec = 1 - specificity,
  inv_spec_upper = 1 - spec_upper,
  inv_spec_lower = 1 - spec_lower
)

stage2a_tibble3 <- tibble(
  sensitivity = calc_sens2(roc_data3, stage_2a),
  specificity = calc_spec2(roc_data3, stage_2a),
  TP = prevalence3 * sensitivity,
  TN = (1 - prevalence3) * specificity,
  FN = prevalence3 - TP,
  FP = (1 - prevalence3) - TN,
  NPV = TN / (TN + FN),
  PPV = TP / (FP + TP),
  surveillance = TP + FP,
  sens_lower = calc_sens_ci2(roc_data3, stage_2a)[[1]][1],
  sens_upper = calc_sens_ci2(roc_data3, stage_2a)[[1]][2],
  spec_lower = calc_sens_ci2(roc_data3, stage_2a)[[1]][1],
  spec_upper = calc_sens_ci2(roc_data3, stage_2a)[[1]][2],
  inv_spec = 1 - specificity,
  inv_spec_upper = 1 - spec_upper,
  inv_spec_lower = 1 - spec_lower
)

stage2b_tibble <- tibble(
  sensitivity = calc_sens2(roc_data, stage_2b),
  specificity = calc_spec2(roc_data, stage_2b),
  TP = prevalence * sensitivity,
  TN = (1 - prevalence) * specificity,
  FN = prevalence - TP,
  FP = (1 - prevalence) - TN,
  NPV = TN / (TN + FN),
  PPV = TP / (FP + TP),
  surveillance = TP + FP,
  sens_lower = calc_sens_ci2(roc_data, stage_2b)[[1]][1],
  sens_upper = calc_sens_ci2(roc_data, stage_2b)[[1]][2],
  spec_lower = calc_spec_ci2(roc_data, stage_2b)[[1]][1],
  spec_upper = calc_spec_ci2(roc_data, stage_2b)[[1]][2],
  inv_spec = 1 - specificity,
  inv_spec_upper = 1 - spec_upper,
  inv_spec_lower = 1 - spec_lower
)

stage2b_tibble2 <- tibble(
  sensitivity = calc_sens2(roc_data2, stage_2b),
  specificity = calc_spec2(roc_data2, stage_2b),
  TP = prevalence2 * sensitivity,
  TN = (1 - prevalence2) * specificity,
  FN = prevalence2 - TP,
  FP = (1 - prevalence2) - TN,
  NPV = TN / (TN + FN),
  PPV = TP / (FP + TP),
  surveillance = TP + FP,
  sens_lower = calc_sens_ci2(roc_data2, stage_2b)[[1]][1],
  sens_upper = calc_sens_ci2(roc_data2, stage_2b)[[1]][2],
  spec_lower = calc_spec_ci2(roc_data2, stage_2b)[[1]][1],
  spec_upper = calc_spec_ci2(roc_data2, stage_2b)[[1]][2],
  inv_spec = 1 - specificity,
  inv_spec_upper = 1 - spec_upper,
  inv_spec_lower = 1 - spec_lower
  )

stage2b_tibble3 <- tibble(
  sensitivity = calc_sens2(roc_data3, stage_2b),
  specificity = calc_spec2(roc_data3, stage_2b),
  TP = prevalence3 * sensitivity,
  TN = (1 - prevalence3) * specificity,
  FN = prevalence3 - TP,
  FP = (1 - prevalence3) - TN,
  NPV = TN / (TN + FN),
  PPV = TP / (FP + TP),
  surveillance = TP + FP,
  sens_lower = calc_sens_ci2(roc_data3, stage_2b)[[1]][1],
  sens_upper = calc_sens_ci2(roc_data3, stage_2b)[[1]][2],
  spec_lower = calc_sens_ci2(roc_data3, stage_2b)[[1]][1],
  spec_upper = calc_sens_ci2(roc_data3, stage_2b)[[1]][2],
  inv_spec = 1 - specificity,
  inv_spec_upper = 1 - spec_upper,
  inv_spec_lower = 1 - spec_lower
  )

stage3a_tibble <- tibble(
  sensitivity = calc_sens2(roc_data, stage_3a),
  specificity = calc_spec2(roc_data, stage_3a),
  TP = prevalence * sensitivity,
  TN = (1 - prevalence) * specificity,
  FN = prevalence - TP,
  FP = (1 - prevalence) - TN,
  NPV = TN / (TN + FN),
  PPV = TP / (FP + TP),
  surveillance = TP + FP,
  sens_lower = calc_sens_ci2(roc_data, stage_3a)[[1]][1],
  sens_upper = calc_sens_ci2(roc_data, stage_3a)[[1]][2],
  spec_lower = calc_spec_ci2(roc_data, stage_3a)[[1]][1],
  spec_upper = calc_spec_ci2(roc_data, stage_3a)[[1]][2],
  inv_spec = 1 - specificity,
  inv_spec_upper = 1 - spec_upper,
  inv_spec_lower = 1 - spec_lower
)

stage3a_tibble2 <- tibble(
  sensitivity = calc_sens2(roc_data2, stage_3a),
  specificity = calc_spec2(roc_data2, stage_3a),
  TP = prevalence2 * sensitivity,
  TN = (1 - prevalence2) * specificity,
  FN = prevalence2 - TP,
  FP = (1 - prevalence2) - TN,
  NPV = TN / (TN + FN),
  PPV = TP / (FP + TP),
  surveillance = TP + FP,
  sens_lower = calc_sens_ci2(roc_data2, stage_3a)[[1]][1],
  sens_upper = calc_sens_ci2(roc_data2, stage_3a)[[1]][2],
  spec_lower = calc_spec_ci2(roc_data2, stage_3a)[[1]][1],
  spec_upper = calc_spec_ci2(roc_data2, stage_3a)[[1]][2],
  inv_spec = 1 - specificity,
  inv_spec_upper = 1 - spec_upper,
  inv_spec_lower = 1 - spec_lower
  )

stage3a_tibble3 <- tibble(
  sensitivity = calc_sens2(roc_data3, stage_3a),
  specificity = calc_spec2(roc_data3, stage_3a),
  TP = prevalence3 * sensitivity,
  TN = (1 - prevalence3) * specificity,
  FN = prevalence3 - TP,
  FP = (1 - prevalence3) - TN,
  NPV = TN / (TN + FN),
  PPV = TP / (FP + TP),
  surveillance = TP + FP,
  sens_lower = calc_sens_ci2(roc_data3, stage_3a)[[1]][1],
  sens_upper = calc_sens_ci2(roc_data3, stage_3a)[[1]][2],
  spec_lower = calc_sens_ci2(roc_data3, stage_3a)[[1]][1],
  spec_upper = calc_sens_ci2(roc_data3, stage_3a)[[1]][2],
  inv_spec = 1 - specificity,
  inv_spec_upper = 1 - spec_upper,
  inv_spec_lower = 1 - spec_lower
  )

# LUMPO tables
cuts_tibble_mm <- tibble(cutoff = c(0.5, 0.3, 0.2, 0.1, 0.05)) %>%
  rowwise() %>%
  mutate(
    sensitivity = calc_sens1(roc_data, MM, cutoff),
    specificity = calc_spec1(roc_data, MM, cutoff),
    TP = prevalence * sensitivity,
    TN = (1 - prevalence) * specificity,
    FN = prevalence - TP,
    FP = (1 - prevalence) - TN,
    NPV = TN / (TN + FN),
    PPV = TP / (FP + TP),
    surveillance = TP + FP,
    sens_lower = calc_sens_ci1(roc_data, MM, cutoff)[[1]][1],
    sens_upper = calc_sens_ci1(roc_data, MM, cutoff)[[1]][2],
    spec_lower = calc_spec_ci1(roc_data, MM, cutoff)[[1]][1],
    spec_upper = calc_spec_ci1(roc_data, MM, cutoff)[[1]][2],
    inv_spec = 1 - specificity,
    inv_spec_upper = 1 - spec_upper,
    inv_spec_lower = 1 - spec_lower
  )

cuts_tibble_mm2 <- tibble(cutoff = c(0.5, 0.3, 0.2, 0.15, 0.07)) %>%
  rowwise() %>%
  mutate(
    sensitivity = calc_sens1(roc_data2, MM, cutoff),
    specificity = calc_spec1(roc_data2, MM, cutoff),
    TP = prevalence2 * sensitivity,
    TN = (1 - prevalence2) * specificity,
    FN = prevalence2 - TP,
    FP = (1 - prevalence2) - TN,
    NPV = TN / (TN + FN),
    PPV = TP / (FP + TP),
    surveillance = TP + FP,
    sens_lower = calc_sens_ci1(roc_data2, MM, cutoff)[[1]][1],
    sens_upper = calc_sens_ci1(roc_data2, MM, cutoff)[[1]][2],
    spec_lower = calc_spec_ci1(roc_data2, MM, cutoff)[[1]][1],
    spec_upper = calc_spec_ci1(roc_data2, MM, cutoff)[[1]][2],
    inv_spec = 1 - specificity,
    inv_spec_upper = 1 - spec_upper,
    inv_spec_lower = 1 - spec_lower
  )

cuts_tibble_mm3 <- tibble(cutoff = c(0.2, 0.1, 0.045)) %>%
  rowwise() %>%
  mutate(
    sensitivity = calc_sens1(roc_data3, MM, cutoff),
    specificity = calc_spec1(roc_data3, MM, cutoff),
    TP = prevalence3 * sensitivity,
    TN = (1 - prevalence3) * specificity,
    FN = prevalence3 - TP,
    FP = (1 - prevalence3) - TN,
    NPV = TN / (TN + FN),
    PPV = TP / (FP + TP),
    surveillance = TP + FP,
    sens_lower = calc_sens_ci1(roc_data3, MM, cutoff)[[1]][1],
    sens_upper = calc_sens_ci1(roc_data3, MM, cutoff)[[1]][2],
    spec_lower = calc_spec_ci1(roc_data3, MM, cutoff)[[1]][1],
    spec_upper = calc_spec_ci1(roc_data3, MM, cutoff)[[1]][2],
    inv_spec = 1 - specificity,
    inv_spec_upper = 1 - spec_upper,
    inv_spec_lower = 1 - spec_lower
  )

cuts_tibble_mm_ngs <- tibble(cutoff = c(0.3, 0.2, 0.1, 0.07, 0.05)) %>% # LUMPO summary table
  rowwise() %>%
  mutate(
    sensitivity = calc_sens1(roc_data_ngs, MM, cutoff),
    specificity = calc_spec1(roc_data_ngs, MM, cutoff),
    TP = prevalence * sensitivity,
    TN = (1 - prevalence) * specificity,
    FN = prevalence - TP,
    FP = (1 - prevalence) - TN,
    NPV = TN / (TN + FN),
    PPV = TP / (FP + TP),
    surveillance = TP + FP,
    sens_lower = calc_sens_ci1(roc_data_ngs, MM, cutoff)[[1]][1],
    sens_upper = calc_sens_ci1(roc_data_ngs, MM, cutoff)[[1]][2],
    spec_lower = calc_spec_ci1(roc_data_ngs, MM, cutoff)[[1]][1],
    spec_upper = calc_spec_ci1(roc_data_ngs, MM, cutoff)[[1]][2],
    inv_spec = 1 - specificity,
    inv_spec_upper = 1 - spec_upper,
    inv_spec_lower = 1 - spec_lower
  )

# LPM tables
cuts_tibble_par <- tibble(cutoff = c(0.3, 0.2, 0.1, 0.05, 0.04, 0.03)) %>% # LPM
  rowwise() %>%
  mutate(
    sensitivity = calc_sens1(roc_data, score, cutoff),
    specificity = calc_spec1(roc_data, score, cutoff),
    TP = prevalence * sensitivity,
    TN = (1 - prevalence) * specificity,
    FN = prevalence - TP,
    FP = (1 - prevalence) - TN,
    NPV = TN / (TN + FN),
    PPV = TP / (FP + TP),
    surveillance = TP + FP,
    sens_lower = calc_sens_ci1(roc_data, score, cutoff)[[1]][1],
    sens_upper = calc_sens_ci1(roc_data, score, cutoff)[[1]][2],
    spec_lower = calc_spec_ci1(roc_data, score, cutoff)[[1]][1],
    spec_upper = calc_spec_ci1(roc_data, score, cutoff)[[1]][2],
    inv_spec = 1 - specificity,
    inv_spec_upper = 1 - spec_upper,
    inv_spec_lower = 1 - spec_lower
  )

cuts_tibble_par2 <- tibble(cutoff = c(0.2, 0.15, 0.14, 0.13, 0.12, 0.11, 0.05)) %>%
  rowwise() %>%
  mutate(
    sensitivity = calc_sens1(roc_data2, score, cutoff),
    specificity = calc_spec1(roc_data2, score, cutoff),
    TP = prevalence2 * sensitivity,
    TN = (1 - prevalence2) * specificity,
    FN = prevalence2 - TP,
    FP = (1 - prevalence2) - TN,
    NPV = TN / (TN + FN),
    PPV = TP / (FP + TP),
    surveillance = TP + FP,
    sens_lower = calc_sens_ci1(roc_data2, score, cutoff)[[1]][1],
    sens_upper = calc_sens_ci1(roc_data2, score, cutoff)[[1]][2],
    spec_lower = calc_spec_ci1(roc_data2, score, cutoff)[[1]][1],
    spec_upper = calc_spec_ci1(roc_data2, score, cutoff)[[1]][2],
    inv_spec = 1 - specificity,
    inv_spec_upper = 1 - spec_upper,
    inv_spec_lower = 1 - spec_lower
  )

cuts_tibble_par3 <- tibble(cutoff = c(0.1, 0.07, 0.06, 0.05, 0.04, 0.03)) %>%
  rowwise() %>%
  mutate(
    sensitivity = calc_sens1(roc_data3, score, cutoff),
    specificity = calc_spec1(roc_data3, score, cutoff),
    TP = prevalence3 * sensitivity,
    TN = (1 - prevalence3) * specificity,
    FN = prevalence3 - TP,
    FP = (1 - prevalence3) - TN,
    NPV = TN / (TN + FN),
    PPV = TP / (FP + TP),
    surveillance = TP + FP,
    sens_lower = calc_sens_ci1(roc_data3, score, cutoff)[[1]][1],
    sens_upper = calc_sens_ci1(roc_data3, score, cutoff)[[1]][2],
    spec_lower = calc_spec_ci1(roc_data3, score, cutoff)[[1]][1],
    spec_upper = calc_spec_ci1(roc_data3, score, cutoff)[[1]][2],
    inv_spec = 1 - specificity,
    inv_spec_upper = 1 - spec_upper,
    inv_spec_lower = 1 - spec_lower
  )

cuts_tibble_par_ngs <- tibble(cutoff = c(0.3, 0.2, 0.1, 0.07, 0.05)) %>%
  rowwise() %>%
  mutate(
    sensitivity = calc_sens1(roc_data_ngs, score, cutoff),
    specificity = calc_spec1(roc_data_ngs, score, cutoff),
    TP = prevalence * sensitivity,
    TN = (1 - prevalence) * specificity,
    FN = prevalence - TP,
    FP = (1 - prevalence) - TN,
    NPV = TN / (TN + FN),
    PPV = TP / (FP + TP),
    surveillance = TP + FP,
    sens_lower = calc_sens_ci1(roc_data_ngs, score, cutoff)[[1]][1],
    sens_upper = calc_sens_ci1(roc_data_ngs, score, cutoff)[[1]][2],
    spec_lower = calc_spec_ci1(roc_data_ngs, score, cutoff)[[1]][1],
    spec_upper = calc_spec_ci1(roc_data_ngs, score, cutoff)[[1]][2],
  )

# AUC with CIs for each dataset-------------------------------------------------------
roc_lumpo <- rocit(
  score = roc_data$MM,
  class = roc_data$outcome, negref = "No mets"
)
paste("roc_data lumpo", auc_label(roc_lumpo))

roc_lumpo2 <- rocit(
  score = roc_data2$MM,
  class = roc_data2$outcome, negref = "No mets"
)
paste("roc_data2 lumpo", auc_label(roc_lumpo2))

roc_lumpo3 <- rocit(
  score = roc_data3$MM,
  class = roc_data3$outcome, negref = "No mets"
)
paste("roc_data3 lumpo", auc_label(roc_lumpo3))

roc_lumpo_ngs <- rocit(
  score = roc_data_ngs$MM,
  class = roc_data_ngs$outcome, negref = "No mets"
)
paste("roc_data_ngs lumpo", auc_label(roc_lumpo_ngs))

roc_lpm <- rocit(
  score = roc_data$score, class = roc_data$outcome,
  negref = "No mets"
)
paste("roc_data lpm", auc_label(roc_lpm))

roc_lpm2 <- rocit(
  score = roc_data2$score, class = roc_data2$outcome,
  negref = "No mets"
)
paste("roc_data2 lpm", auc_label(roc_lpm2))

roc_lpm3 <- rocit(
  score = roc_data3$score, class = roc_data3$outcome,
  negref = "No mets"
)
paste("roc_data3 lpm", auc_label(roc_lpm3))

roc_lpm_ngs <- rocit(
  score = roc_data_ngs$score, class = roc_data_ngs$outcome,
  negref = "No mets"
)
paste("roc_data_ngs lpm", auc_label(roc_lpm_ngs))

# Figure 2-----------------------------------------------------------------------
longdata <- melt_roc(roc_data, "outcome_binary", c("MM", "score")) %>%
  mutate(
    name = if_else(name == "MM", "LUMPO III", "LPM"),
    name = factor(name, levels = c("LUMPO III", "LPM"))
  )

fig2 <- ggplot(longdata, aes(d = D, m = M, color = name)) +
  scale_color_manual(values = c("black", "#E69F00")) +
  annotate(
    geom = "rect", xmin = m3_tibble$inv_spec_upper, xmax = m3_tibble$inv_spec_lower, ymin = m3_tibble$sens_lower, ymax = m3_tibble$sens_upper,
    fill = "pink", alpha = 0.5
  ) +
  annotate(
    geom = "rect", xmin = stage2a_tibble$inv_spec_upper, xmax = stage2a_tibble$inv_spec_lower, ymin = stage2a_tibble$sens_lower, ymax = stage2a_tibble$sens_upper,
    fill = "purple", alpha = 0.5
  ) +
  annotate(
    geom = "rect", xmin = stage2b_tibble$inv_spec_upper, xmax = stage2b_tibble$inv_spec_lower, ymin = stage2b_tibble$sens_lower, ymax = stage2b_tibble$sens_upper,
    fill = "purple", alpha = 0.5
  ) +
  annotate(
    geom = "rect", xmin = stage3a_tibble$inv_spec_upper, xmax = stage3a_tibble$inv_spec_lower, ymin = stage3a_tibble$sens_lower, ymax = stage3a_tibble$sens_upper,
    fill = "purple", alpha = 0.5
  ) +
  annotate(
    geom = "rect", xmin = cuts_tibble_mm$inv_spec_upper[1], xmax = cuts_tibble_mm$inv_spec_lower[1],
    ymin = cuts_tibble_mm$sens_lower[1], ymax = cuts_tibble_mm$sens_upper[1],
    fill = "black", alpha = 0.4
  ) +
  annotate(
    geom = "rect", xmin = cuts_tibble_mm$inv_spec_upper[3], xmax = cuts_tibble_mm$inv_spec_lower[3],
    ymin = cuts_tibble_mm$sens_lower[3], ymax = cuts_tibble_mm$sens_upper[3],
    fill = "black", alpha = 0.4
  ) +
  annotate(
    geom = "rect", xmin = cuts_tibble_mm$inv_spec_upper[4], xmax = cuts_tibble_mm$inv_spec_lower[4],
    ymin = cuts_tibble_mm$sens_lower[4], ymax = cuts_tibble_mm$sens_upper[4],
    fill = "black", alpha = 0.4
  ) +
  annotate(
    geom = "rect", xmin = cuts_tibble_mm$inv_spec_upper[5], xmax = cuts_tibble_mm$inv_spec_lower[5],
    ymin = cuts_tibble_mm$sens_lower[5], ymax = cuts_tibble_mm$sens_upper[5],
    fill = "black", alpha = 0.4
  ) +
  geom_roc(n.cuts = 0, labels = FALSE) +
  labs(x = "1-specificity", y = "sensitivity") +
  geom_point(x = 1 - m3_tibble$specificity, y = m3_tibble$sensitivity, color = "#F8766D") +
  geom_point(x = stage2a_tibble$inv_spec, y = stage2a_tibble$sensitivity, color = "purple") +
  geom_point(x = stage2b_tibble$inv_spec, y = stage2b_tibble$sensitivity, color = "purple") +
  geom_point(x = stage3a_tibble$inv_spec, y = stage3a_tibble$sensitivity, color = "purple") +
  geom_point(x = cuts_tibble_mm$inv_spec[1], y = cuts_tibble_mm$sensitivity[1], color = "black") +
  annotate("text",
    x = cuts_tibble_mm$inv_spec[1] - 0.05, y = cuts_tibble_mm$sensitivity[1],
    label = cuts_tibble_mm$cutoff[1], size = 4
  ) +
  geom_point(x = cuts_tibble_mm$inv_spec[3], y = cuts_tibble_mm$sensitivity[3], color = "black") +
  annotate("text",
    x = cuts_tibble_mm$inv_spec[3] - 0.06, y = cuts_tibble_mm$sensitivity[3],
    label = cuts_tibble_mm$cutoff[3], size = 4
  ) +
  geom_point(x = cuts_tibble_mm$inv_spec[4], y = cuts_tibble_mm$sensitivity[4], color = "black") +
  annotate("text",
    x = cuts_tibble_mm$inv_spec[4], y = cuts_tibble_mm$sensitivity[4] + 0.06,
    label = cuts_tibble_mm$cutoff[4], size = 4
  ) +
  geom_point(x = cuts_tibble_mm$inv_spec[5], y = cuts_tibble_mm$sensitivity[5], color = "black") +
  annotate("text",
    x = cuts_tibble_mm$inv_spec[5], y = cuts_tibble_mm$sensitivity[5] + 0.05,
    label = cuts_tibble_mm$cutoff[5], size = 4
  ) +
  annotate("text",
    x = stage2a_tibble$inv_spec + 0.13, y = stage2a_tibble$sensitivity, label = ">= IIA", size = 4,
    color = "mediumorchid4"
  ) +
  annotate("text",
    x = stage2b_tibble$inv_spec + 0.13, y = stage2b_tibble$sensitivity, label = ">= IIB", size = 4,
    color = "mediumorchid4"
  ) +
  annotate("text",
    x = stage3a_tibble$inv_spec + 0.13, y = stage3a_tibble$sensitivity, label = ">= IIIA", size = 4,
    color = "mediumorchid4"
  ) +
  geom_point(x = 0.68, y = 0.08, color = "#F8766D") +
  annotate("text", x = 0.9, y = 0.08, label = "Monosomy 3", size = 3) +
  geom_point(x = 0.68, y = 0.03, color = "purple") +
  annotate("text", x = 0.9, y = 0.03, label = "AJCC system", size = 3) +
  theme_classic() +
  theme(legend.position = c(0.8, 0.25)) +
  theme(legend.title = element_blank()) +
  theme(legend.background = element_blank())

fig2

# Figure 3
longdata2 <- melt_roc(roc_data2, "outcome_binary", c("MM", "score")) %>%
  mutate(
    name = if_else(name == "MM", "LUMPO III", "LPM"),
    name = factor(name, levels = c("LUMPO III", "LPM"))
  )

fig3a <- ggplot(longdata2, aes(d = D, m = M, color = name)) +
  scale_color_manual(values = c("black", "#E69F00")) +
  annotate(
    geom = "rect", xmin = m3_tibble2$inv_spec_upper, xmax = m3_tibble2$inv_spec_lower,
    ymin = m3_tibble2$sens_lower, ymax = m3_tibble2$sens_upper,
    fill = "pink", alpha = 0.5
  ) +
  annotate(
    geom = "rect", xmin = stage2a_tibble2$inv_spec_upper, xmax = stage2a_tibble2$inv_spec_lower,
    ymin = stage2a_tibble2$sens_lower, ymax = stage2a_tibble2$sens_upper,
    fill = "purple", alpha = 0.5
  ) +
  annotate(
    geom = "rect", xmin = stage2b_tibble2$inv_spec_upper, xmax = stage2b_tibble2$inv_spec_lower,
    ymin = stage2b_tibble2$sens_lower, ymax = stage2b_tibble2$sens_upper,
    fill = "purple", alpha = 0.5
  ) +
  annotate(
    geom = "rect", xmin = stage3a_tibble2$inv_spec_upper, xmax = stage3a_tibble2$inv_spec_lower,
    ymin = stage3a_tibble2$sens_lower, ymax = stage3a_tibble2$sens_upper,
    fill = "purple", alpha = 0.5
  ) +
  annotate(
    geom = "rect", xmin = cuts_tibble_mm2$inv_spec_upper[1], xmax = cuts_tibble_mm2$inv_spec_lower[1],
    ymin = cuts_tibble_mm2$sens_lower[1], ymax = cuts_tibble_mm2$sens_upper[1],
    fill = "black", alpha = 0.4
  ) +
  annotate(
    geom = "rect", xmin = cuts_tibble_mm2$inv_spec_upper[4], xmax = cuts_tibble_mm2$inv_spec_lower[4],
    ymin = cuts_tibble_mm2$sens_lower[4], ymax = cuts_tibble_mm2$sens_upper[4],
    fill = "black", alpha = 0.4
  ) +
  annotate(
    geom = "rect", xmin = cuts_tibble_mm2$inv_spec_upper[5], xmax = cuts_tibble_mm2$inv_spec_lower[5],
    ymin = cuts_tibble_mm2$sens_lower[5], ymax = cuts_tibble_mm2$sens_upper[5],
    fill = "black", alpha = 0.4
  ) +
  geom_roc(n.cuts = 0, labels = FALSE) +
  labs(x = "1-specificity", y = "sensitivity", title = "A: Subpopulation with chr3 result") +
  geom_point(x = 1 - m3_tibble2$specificity, y = m3_tibble2$sensitivity, color = "#F8766D") +
  geom_point(x = stage2a_tibble2$inv_spec, y = stage2a_tibble2$sensitivity, color = "purple") +
  geom_point(x = stage2b_tibble2$inv_spec, y = stage2b_tibble2$sensitivity, color = "purple") +
  geom_point(x = stage3a_tibble2$inv_spec, y = stage3a_tibble2$sensitivity, color = "purple") +
  geom_point(x = cuts_tibble_mm2$inv_spec[1], y = cuts_tibble_mm2$sensitivity[1], color = "black") +
  annotate("text",
           x = cuts_tibble_mm2$inv_spec[1] - 0.07, y = cuts_tibble_mm2$sensitivity[1],
           label = cuts_tibble_mm2$cutoff[1], size = 4
  ) +
  geom_point(x = cuts_tibble_mm2$inv_spec[4], y = cuts_tibble_mm2$sensitivity[4], color = "black") +
  annotate("text",
           x = cuts_tibble_mm2$inv_spec[4], y = cuts_tibble_mm2$sensitivity[4] + 0.06,
           label = cuts_tibble_mm2$cutoff[4], size = 4
  ) +
  geom_point(x = cuts_tibble_mm2$inv_spec[5], y = cuts_tibble_mm2$sensitivity[5], color = "black") +
  annotate("text",
           x = cuts_tibble_mm2$inv_spec[5], y = cuts_tibble_mm2$sensitivity[5] + 0.05,
           label = cuts_tibble_mm2$cutoff[5], size = 4
  ) +
  annotate("text",
           x = stage2a_tibble2$inv_spec + 0.15, y = stage2a_tibble2$sensitivity, label = ">= IIA", size = 4,
           color = "mediumorchid4"
  ) +
  annotate("text",
           x = stage2b_tibble2$inv_spec + 0.15, y = stage2b_tibble2$sensitivity, label = ">= IIB", size = 4,
           color = "mediumorchid4"
  ) +
  annotate("text",
           x = stage3a_tibble2$inv_spec + 0.15, y = stage3a_tibble2$sensitivity, label = ">= IIIA", size = 4,
           color = "mediumorchid4"
  ) +
  geom_point(x = 0.68, y = 0.08, color = "#F8766D") +
  annotate("text", x = 0.9, y = 0.08, label = "Monosomy 3", size = 3) +
  geom_point(x = 0.68, y = 0.03, color = "purple") +
  annotate("text", x = 0.9, y = 0.03, label = "AJCC system", size = 3) +
  theme_classic() +
  theme(legend.position = c(0.8, 0.25)) +
  theme(legend.title = element_blank()) +
  theme(legend.background = element_blank())

fig3a

longdata3 <- melt_roc(roc_data3, "outcome_binary", c("MM", "score")) %>%
  mutate(
    name = if_else(name == "MM", "LUMPO III", "LPM"),
    name = factor(name, levels = c("LUMPO III", "LPM"))
  )

fig3b <- ggplot(longdata3, aes(d = D, m = M, color = name)) +
  scale_color_manual(values = c("black", "#E69F00")) +
  annotate(
    geom = "rect", xmin = stage2a_tibble3$inv_spec_upper, xmax = stage2a_tibble3$inv_spec_lower,
    ymin = stage2a_tibble3$sens_lower, ymax = stage2a_tibble3$sens_upper,
    fill = "purple", alpha = 0.5
  ) +
  annotate(
    geom = "rect", xmin = stage2b_tibble3$inv_spec_upper, xmax = stage2b_tibble3$inv_spec_lower,
    ymin = stage2b_tibble3$sens_lower, ymax = stage2b_tibble3$sens_upper,
    fill = "purple", alpha = 0.5
  ) +
  annotate(
    geom = "rect", xmin = stage3a_tibble3$inv_spec_upper, xmax = stage3a_tibble3$inv_spec_lower,
    ymin = stage3a_tibble3$sens_lower, ymax = stage3a_tibble3$sens_upper,
    fill = "purple", alpha = 0.5
  ) +
  annotate(
    geom = "rect", xmin = cuts_tibble_mm3$inv_spec_upper[1], xmax = cuts_tibble_mm3$inv_spec_lower[1],
    ymin = cuts_tibble_mm3$sens_lower[1], ymax = cuts_tibble_mm3$sens_upper[1],
    fill = "black", alpha = 0.4
  ) +
  annotate(
    geom = "rect", xmin = cuts_tibble_mm3$inv_spec_upper[2], xmax = cuts_tibble_mm3$inv_spec_lower[2],
    ymin = cuts_tibble_mm3$sens_lower[2], ymax = cuts_tibble_mm3$sens_upper[2],
    fill = "black", alpha = 0.4
  ) +
  annotate(
    geom = "rect", xmin = cuts_tibble_mm3$inv_spec_upper[3], xmax = cuts_tibble_mm3$inv_spec_lower[3],
    ymin = cuts_tibble_mm3$sens_lower[3], ymax = cuts_tibble_mm3$sens_upper[3],
    fill = "black", alpha = 0.4
  ) +
  geom_roc(n.cuts = 0, labels = FALSE) +
  labs(x = "1-specificity", y = "sensitivity", title = "B: Subpopulation with unknown chr3") +
  geom_point(x = stage2a_tibble3$inv_spec, y = stage2a_tibble3$sensitivity, color = "purple") +
  geom_point(x = stage2b_tibble3$inv_spec, y = stage2b_tibble3$sensitivity, color = "purple") +
  geom_point(x = stage3a_tibble3$inv_spec, y = stage3a_tibble3$sensitivity, color = "purple") +
  geom_point(x = cuts_tibble_mm3$inv_spec[1], y = cuts_tibble_mm3$sensitivity[1], color = "black") +
  annotate("text",
           x = cuts_tibble_mm3$inv_spec[1] - 0.07, y = cuts_tibble_mm3$sensitivity[1],
           label = cuts_tibble_mm3$cutoff[1], size = 4
  ) +
  geom_point(x = cuts_tibble_mm3$inv_spec[2], y = cuts_tibble_mm3$sensitivity[2], color = "black") +
  annotate("text",
           x = cuts_tibble_mm3$inv_spec[2] - 0.07, y = cuts_tibble_mm3$sensitivity[2],
           label = cuts_tibble_mm3$cutoff[2], size = 4
  ) +
  geom_point(x = cuts_tibble_mm3$inv_spec[3], y = cuts_tibble_mm3$sensitivity[3], color = "black") +
  annotate("text",
           x = cuts_tibble_mm3$inv_spec[3], y = cuts_tibble_mm3$sensitivity[3] + 0.11,
           label = cuts_tibble_mm3$cutoff[3], size = 4
  ) +
  annotate("text",
           x = stage2a_tibble3$inv_spec - 0.14, y = stage2a_tibble3$sensitivity, label = ">= IIA", size = 4,
           color = "mediumorchid4"
  ) +
  annotate("text",
           x = stage2b_tibble3$inv_spec - 0.13, y = stage2b_tibble3$sensitivity, label = ">= IIB", size = 4,
           color = "mediumorchid4"
  ) +
  annotate("text",
           x = stage3a_tibble3$inv_spec + 0.13, y = stage3a_tibble3$sensitivity, label = ">= IIIA", size = 4,
           color = "mediumorchid4"
  ) +
  theme_classic() +
  theme(legend.position = "none") +
  theme(legend.title = element_blank()) +
  theme(legend.background = element_blank())

fig3b

# Figure 5
longdata <- melt_roc(roc_data_ngs, "outcome_binary", c("MM", "score")) %>%
  mutate(
    name = if_else(name == "MM", "LUMPO III", "LPM"),
    name = factor(name, levels = c("LUMPO III", "LPM"))
  )

fig5 <- ggplot(longdata, aes(d = D, m = M, color = name)) +
  scale_color_manual(values = c("black", "#E69F00")) +
  annotate(
    geom = "rect", xmin = stage2a_tibble$inv_spec_upper, xmax = stage2a_tibble$inv_spec_lower, ymin = stage2a_tibble$sens_lower, ymax = stage2a_tibble$sens_upper,
    fill = "purple", alpha = 0.5
  ) +
  annotate(
    geom = "rect", xmin = stage2b_tibble$inv_spec_upper, xmax = stage2b_tibble$inv_spec_lower, ymin = stage2b_tibble$sens_lower, ymax = stage2b_tibble$sens_upper,
    fill = "purple", alpha = 0.5
  ) +
  annotate(
    geom = "rect", xmin = stage3a_tibble$inv_spec_upper, xmax = stage3a_tibble$inv_spec_lower, ymin = stage3a_tibble$sens_lower, ymax = stage3a_tibble$sens_upper,
    fill = "purple", alpha = 0.5
  ) +
  annotate(
    geom = "rect", xmin = cuts_tibble_mm_ngs$inv_spec_upper[1], xmax = cuts_tibble_mm_ngs$inv_spec_lower[1],
    ymin = cuts_tibble_mm_ngs$sens_lower[1], ymax = cuts_tibble_mm_ngs$sens_upper[1],
    fill = "black", alpha = 0.4
  ) +
  annotate(
    geom = "rect", xmin = cuts_tibble_mm_ngs$inv_spec_upper[3], xmax = cuts_tibble_mm_ngs$inv_spec_lower[3],
    ymin = cuts_tibble_mm_ngs$sens_lower[3], ymax = cuts_tibble_mm_ngs$sens_upper[3],
    fill = "black", alpha = 0.4
  ) +
  annotate(
    geom = "rect", xmin = cuts_tibble_mm_ngs$inv_spec_upper[4], xmax = cuts_tibble_mm_ngs$inv_spec_lower[4],
    ymin = cuts_tibble_mm_ngs$sens_lower[4], ymax = cuts_tibble_mm_ngs$sens_upper[4],
    fill = "black", alpha = 0.4
  ) +
  geom_roc(n.cuts = 0, labels = FALSE) +
  labs(x = "1-specificity", y = "sensitivity") +
  geom_point(x = stage2a_tibble$inv_spec, y = stage2a_tibble$sensitivity, color = "purple") +
  geom_point(x = stage2b_tibble$inv_spec, y = stage2b_tibble$sensitivity, color = "purple") +
  geom_point(x = stage3a_tibble$inv_spec, y = stage3a_tibble$sensitivity, color = "purple") +
  geom_point(x = cuts_tibble_mm_ngs$inv_spec[1], y = cuts_tibble_mm_ngs$sensitivity[1], color = "black") +
  annotate("text",
           x = cuts_tibble_mm_ngs$inv_spec[1] - 0.07, y = cuts_tibble_mm_ngs$sensitivity[1],
           label = cuts_tibble_mm_ngs$cutoff[1], size = 4
  ) +
  geom_point(x = cuts_tibble_mm_ngs$inv_spec[3], y = cuts_tibble_mm_ngs$sensitivity[3], color = "black") +
  annotate("text",
           x = cuts_tibble_mm_ngs$inv_spec[3] - 0.08, y = cuts_tibble_mm_ngs$sensitivity[3],
           label = cuts_tibble_mm_ngs$cutoff[3], size = 4
  ) +
  geom_point(x = cuts_tibble_mm_ngs$inv_spec[4], y = cuts_tibble_mm_ngs$sensitivity[4], color = "black") +
  annotate("text",
           x = cuts_tibble_mm_ngs$inv_spec[4], y = cuts_tibble_mm_ngs$sensitivity[4] + 0.06,
           label = cuts_tibble_mm_ngs$cutoff[4], size = 4
  ) +
  annotate("text",
           x = stage2a_tibble$inv_spec + 0.13, y = stage2a_tibble$sensitivity, label = ">= IIA", size = 4,
           color = "mediumorchid4"
  ) +
  annotate("text",
           x = stage2b_tibble$inv_spec + 0.13, y = stage2b_tibble$sensitivity, label = ">= IIB", size = 4,
           color = "mediumorchid4"
  ) +
  annotate("text",
           x = stage3a_tibble$inv_spec + 0.13, y = stage3a_tibble$sensitivity, label = ">= IIIA", size = 4,
           color = "mediumorchid4"
  ) +
  geom_point(x = 0.72, y = 0.04, color = "purple") +
  annotate("text", x = 0.9, y = 0.03, label = "AJCC system", size = 3) +
  theme_classic() +
  theme(legend.position = c(0.8, 0.25)) +
  theme(legend.title = element_blank()) +
  theme(legend.background = element_blank())

fig5

#### Supplementary figure 3 ------------------------

prev_tibble <- tibble(
  subpopulation = c("Chr3 result", "No result"),
  prevalence = c(prevalence2, prevalence3),
  prev_low = c(ci_prev2[[1]][[1]], ci_prev3[[1]][[1]]),
  prev_high = c(ci_prev2[[1]][[2]], ci_prev3[[1]][[2]])
)

figs3a <- prev_tibble %>% ggplot(aes(x = subpopulation, y = prevalence, ymax = prev_high, ymin = prev_low)) +
  geom_bar(stat = "identity") +
  geom_errorbar(width = 0.1) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_classic() +
  labs(title = "A: Prevalence of outcome")

figs3b <- roc_data %>%
  ggplot(aes(with3, area)) +
  geom_boxplot() +
  labs(x = NULL, y = "tumour size (mm^2)", title = "B: Tumour size in subpopulations") +
  theme_classic()

figs3 <- figs3a + figs3b
figs3

#### Summary of Figure 4 lumpo strategy
roc_data <- roc_data %>% mutate(
  fig4_str = case_when(
    is.na(chr3) & MM >= 0.045 ~ "Surveillance",
    is.na(chr3) & MM < 0.045 ~ "No surveillance",
    !is.na(chr3) & MM >= 0.07 ~ "Surveillance",
    !is.na(chr3) & MM < 0.07 ~ "No surveillance"
  )
)

fig4_lumpo <- tibble(
  sensitivity = calc_sens2(roc_data, fig4_str),
  specificity = calc_spec2(roc_data, fig4_str),
  TP = prevalence * sensitivity,
  TN = (1 - prevalence) * specificity,
  FN = prevalence - TP,
  FP = (1 - prevalence) - TN,
  NPV = TN / (TN + FN),
  PPV = TP / (FP + TP),
  surveillance = TP + FP,
  sens_lower = calc_sens_ci2(roc_data, fig4_str)[[1]][1],
  sens_upper = calc_sens_ci2(roc_data, fig4_str)[[1]][2],
  spec_lower = calc_spec_ci2(roc_data, fig4_str)[[1]][1],
  spec_upper = calc_spec_ci2(roc_data, fig4_str)[[1]][2]
) %>%
  mutate(
    inv_spec = 1 - specificity,
    inv_spec_upper = 1 - spec_upper,
    inv_spec_lower = 1 - spec_lower
  )

# Plot 3 strategies
strategies <- tibble(
  strategy = factor(c("LUMPO III", "Monosomy 3", "AJCC"), levels = c("LUMPO III", "Monosomy 3", "AJCC")),
  sensitivity = c(fig4_lumpo$sensitivity, m3_tibble$sensitivity, stage2a_tibble$sensitivity),
  specificity = c(fig4_lumpo$specificity, m3_tibble$specificity, stage2a_tibble$specificity),
  sens_upper = c(fig4_lumpo$sens_upper, m3_tibble$sens_upper, stage2a_tibble$sens_upper),
  sens_lower = c(fig4_lumpo$sens_lower, m3_tibble$sens_lower, stage2a_tibble$sens_lower),
  spec_upper = c(fig4_lumpo$spec_upper, m3_tibble$spec_upper, stage2a_tibble$spec_upper),
  spec_lower = c(fig4_lumpo$spec_lower, m3_tibble$spec_lower, stage2a_tibble$spec_lower)
) %>%
  mutate(
    inv_spec = 1 - specificity,
    inv_spec_upper = 1 - spec_upper,
    inv_spec_lower = 1 - spec_lower
  )

fig4 <- strategies %>% ggplot(aes(x = inv_spec, y = sensitivity)) +
  geom_hline(aes(yintercept = 0.95), linetype = "dashed", color = "grey", alpha = 0.5) +
  geom_rect(aes(
    xmin = inv_spec_upper, xmax = inv_spec_lower, ymin = sens_lower,
    ymax = sens_upper, fill = strategy
  ), alpha = 0.5) +
  geom_point(aes(color = strategy)) +
  geom_abline(aes(intercept = 0, slope = 1), color = "white") +
  scale_x_continuous("1-specificity", limits = c(0, 1)) +
  scale_y_continuous("sensitivity", limits = c(0, 1)) +
  scale_color_manual(values = c("black", "#F8766D", "purple")) +
  scale_fill_manual(values = c("black", "pink", "purple")) +
  theme(legend.position = (c(0.85, 0.25))) +
  labs(
    title = "B: Strategy summary"
  ) +
  theme_classic() +
  theme(legend.position = c(0.8, 0.25)) +
  theme(legend.title = element_blank()) +
  theme(legend.background = element_blank())

fig4

### Calculation of cost savings --------------------------------------
calc_saving(211.24, 200, 0.51, 0.38)

calc_saving(211.24, 200, 0.65, 0.44)

calc_saving(211.24, 200, 0.55, 0.37)

### Save all outputs
ggsave(here("outputs", "fig1.png"), plot = (dist_lumpo + dist_ajcc), width = 8.27, height = 4, units = "in")
ggsave(here("outputs", "fig2.png"), plot = fig2, width = 5, height = 5, units = "in")
ggsave(here("outputs", "fig3.png"), plot = (fig3a + fig3b), width = 8.27, height = 4, units = "in")
ggsave(here("outputs", "fig5.png"), plot = fig5, width = 5, height = 5, units = "in")
ggsave(here("outputs", "figs3.png"), plot = (figs3), width = 8.27, height = 4, units = "in")
ggsave(here("outputs", "fig4.png"), plot = fig4, width = 5, height = 5, units = "in")



