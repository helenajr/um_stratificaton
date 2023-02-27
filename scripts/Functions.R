library(tidyverse)
library(lubridate)
library(ROCit)

calc_prevalence <- function(df) {
  total <- nrow(df)

  df1 <- df %>%
   filter(outcome == "Mets")

  positives <- nrow(df1)

  prevalence <- positives / total

  return(prevalence)
}

calc_ci_prev <- function(df) {
  total <- nrow(df)

  df1 <- df %>%
    filter(outcome == "Mets")

  positives <- nrow(df1)

  ci <- exactci(positives, total, conf.level = 0.95)

  return(ci)
}

# Sensitivity functions ---------------------------------------------------------------------------
# Sensitivity calculation where classification is based on a numerical column with cutoff
calc_sens1 <- function(df, col, cut) {
  mycol <- enquo(col)
  df1 <- df %>%
    filter(outcome == "Mets") %>%
    mutate(classi = if_else(UQ(mycol) >= cut, "positive", "negative"))

  df_pos <- df1 %>%
    filter(classi == "positive")

  df_neg <- df1 %>%
    filter(classi == "negative")

  tp <- nrow(df_pos)
  fn <- nrow(df_neg)

  sensitivity <- tp / (tp + fn)

  return(sensitivity)
}

# Calculate sensitivity where classification is contained in a character column
calc_sens2 <- function(df, col) {
  mycol <- enquo(col)
  df1 <- df %>%
    filter(outcome == "Mets") %>%
    mutate(classi = if_else(UQ(mycol) == "Surveillance", "positive", "negative"))

  df_pos <- df1 %>%
    filter(classi == "positive")

  df_neg <- df1 %>%
    filter(classi == "negative")

  tp <- nrow(df_pos)
  fn <- nrow(df_neg)

  sensitivity <- tp / (tp + fn)

  return(sensitivity)
}

# Specificity functions -----------------------------------------------------------------
calc_spec1 <- function(df, col, cut) {
  mycol <- enquo(col)
  df1 <- df %>%
    filter(outcome == "No mets") %>%
    mutate(classi = if_else(UQ(mycol) >= cut, "positive", "negative"))

  df_pos <- df1 %>%
    filter(classi == "positive")

  df_neg <- df1 %>%
    filter(classi == "negative")

  tn <- nrow(df_neg)
  fp <- nrow(df_pos)

  specificity <- tn / (tn + fp)

  return(specificity)
}

calc_spec2 <- function(df, col) {
  mycol <- enquo(col)
  df1 <- df %>%
    filter(outcome == "No mets") %>%
    mutate(classi = if_else(UQ(mycol) == "Surveillance", "positive", "negative"))

  df_pos <- df1 %>%
    filter(classi == "positive")

  df_neg <- df1 %>%
    filter(classi == "negative")

  tn <- nrow(df_neg)
  fp <- nrow(df_pos)

  specificity <- tn / (tn + fp)

  return(specificity)
}

# Confidence interval functions -------------------------------------------------------
calc_sens_ci1 <- function(df, col, cut) {
  mycol <- enquo(col)
  df1 <- df %>%
    filter(outcome == "Mets") %>%
    mutate(classi = if_else(UQ(mycol) >= cut, "positive", "negative"))

  df_pos <- df1 %>%
    filter(classi == "positive")

  df_neg <- df1 %>%
    filter(classi == "negative")

  tp <- nrow(df_pos)
  fn <- nrow(df_neg)

  ci <- exactci(tp, (tp + fn), conf.level = 0.95)

  return(ci)
}

calc_spec_ci1 <- function(df, col, cut) {
  mycol <- enquo(col)
  df1 <- df %>%
    filter(outcome == "No mets") %>%
    mutate(classi = if_else(UQ(mycol) >= cut, "positive", "negative"))

  df_pos <- df1 %>%
    filter(classi == "positive")

  df_neg <- df1 %>%
    filter(classi == "negative")

  tn <- nrow(df_neg)
  fp <- nrow(df_pos)

  ci <- exactci(tn, (tn + fp), conf.level = 0.95)

  return(ci)
}

calc_sens_ci2 <- function(df, col) {
  mycol <- enquo(col)
  df1 <- df %>%
    filter(outcome == "Mets") %>%
    mutate(classi = if_else(UQ(mycol) == "Surveillance", "positive", "negative"))

  df_pos <- df1 %>%
    filter(classi == "positive")

  df_neg <- df1 %>%
    filter(classi == "negative")

  tp <- nrow(df_pos)
  fn <- nrow(df_neg)

  ci <- exactci(tp, (tp + fn), conf.level = 0.95)

  return(ci)
}

calc_spec_ci2 <- function(df, col) {
  mycol <- enquo(col)
  df1 <- df %>%
    filter(outcome == "No mets") %>%
    mutate(classi = if_else(UQ(mycol) == "Surveillance", "positive", "negative"))

  df_pos <- df1 %>%
    filter(classi == "positive")

  df_neg <- df1 %>%
    filter(classi == "negative")

  tn <- nrow(df_neg)
  fp <- nrow(df_pos)

  ci <- exactci(tn, (tn + fp), conf.level = 0.95)

  return(ci)
}

calc_chr3_sens_ci95 <- function(df) {
  df1 <- df %>%
    filter(outcome == "Mets") %>%
    mutate(classi = if_else(chr3 == "L" | is.na(chr3), "positive", "negative"))

  df_pos <- df1 %>%
    filter(classi == "positive")

  df_neg <- df1 %>%
    filter(classi == "negative")

  tp <- nrow(df_pos)
  fn <- nrow(df_neg)

  ci <- exactci(tp, (tp + fn), conf.level = 0.95)

  return(ci)
}

calc_chr3_spec_ci95 <- function(df) {
  df1 <- df %>%
    filter(outcome == "No mets") %>%
    mutate(classi = if_else(chr3 == "L" | is.na(chr3), "positive", "negative"))

  df_pos <- df1 %>%
    filter(classi == "positive")

  df_neg <- df1 %>%
    filter(classi == "negative")

  tn <- nrow(df_neg)
  fp <- nrow(df_pos)

  ci <- exactci(tn, (tn + fp), conf.level = 0.95)

  return(ci)
}

calc_nlumpo_sens_ci95 <- function(df) {
  df1 <- df %>%
    filter(outcome == "Mets") %>%
    mutate(classi = if_else(new_str_lumpo == "Surveillance", "positive", "negative"))

  df_pos <- df1 %>%
    filter(classi == "positive")

  df_neg <- df1 %>%
    filter(classi == "negative")

  tp <- nrow(df_pos)
  fn <- nrow(df_neg)

  ci <- exactci(tp, (tp + fn), conf.level = 0.95)

  return(ci)
}

calc_nlumpo_spec_ci95 <- function(df) {
  df1 <- df %>%
    filter(outcome == "No mets") %>%
    mutate(classi = if_else(new_str_lumpo == "Surveillance", "positive", "negative"))

  df_pos <- df1 %>%
    filter(classi == "positive")

  df_neg <- df1 %>%
    filter(classi == "negative")

  tn <- nrow(df_neg)
  fp <- nrow(df_pos)

  ci <- exactci(tn, (tn + fp), conf.level = 0.95)

  return(ci)
}

calc_nlumpo2_sens_ci95 <- function(df) {
  df1 <- df %>%
    filter(outcome == "Mets") %>%
    mutate(classi = if_else(new_str_lumpo2 == "Surveillance", "positive", "negative"))

  df_pos <- df1 %>%
    filter(classi == "positive")

  df_neg <- df1 %>%
    filter(classi == "negative")

  tp <- nrow(df_pos)
  fn <- nrow(df_neg)

  ci <- exactci(tp, (tp + fn), conf.level = 0.95)

  return(ci)
}

calc_nlumpo2_spec_ci95 <- function(df) {
  df1 <- df %>%
    filter(outcome == "No mets") %>%
    mutate(classi = if_else(new_str_lumpo2 == "Surveillance", "positive", "negative"))

  df_pos <- df1 %>%
    filter(classi == "positive")

  df_neg <- df1 %>%
    filter(classi == "negative")

  tn <- nrow(df_neg)
  fp <- nrow(df_pos)

  ci <- exactci(tn, (tn + fp), conf.level = 0.95)

  return(ci)
}

calc_nsimple_sens_ci95 <- function(df) {
  df1 <- df %>%
    filter(outcome == "Mets") %>%
    mutate(classi = if_else(new_str_simple == "Surveillance", "positive", "negative"))

  df_pos <- df1 %>%
    filter(classi == "positive")

  df_neg <- df1 %>%
    filter(classi == "negative")

  tp <- nrow(df_pos)
  fn <- nrow(df_neg)

  ci <- exactci(tp, (tp + fn), conf.level = 0.95)

  return(ci)
}

calc_nsimple_spec_ci95 <- function(df) {
  df1 <- df %>%
    filter(outcome == "No mets") %>%
    mutate(classi = if_else(new_str_simple == "Surveillance", "positive", "negative"))

  df_pos <- df1 %>%
    filter(classi == "positive")

  df_neg <- df1 %>%
    filter(classi == "negative")

  tn <- nrow(df_neg)
  fp <- nrow(df_pos)

  ci <- exactci(tn, (tn + fp), conf.level = 0.95)

  return(ci)
}

calc_stage2a_sens_ci95 <- function(df) {
  df1 <- df %>%
    filter(outcome == "Mets") %>%
    mutate(classi = if_else(stage_2a == "Surveillance", "positive", "negative"))

  df_pos <- df1 %>%
    filter(classi == "positive")

  df_neg <- df1 %>%
    filter(classi == "negative")

  tp <- nrow(df_pos)
  fn <- nrow(df_neg)

  ci <- exactci(tp, (tp + fn), conf.level = 0.95)

  return(ci)
}

calc_stage2a_spec_ci95 <- function(df) {
  df1 <- df %>%
    filter(outcome == "No mets") %>%
    mutate(classi = if_else(stage_2a == "Surveillance", "positive", "negative"))

  df_pos <- df1 %>%
    filter(classi == "positive")

  df_neg <- df1 %>%
    filter(classi == "negative")

  tn <- nrow(df_neg)
  fp <- nrow(df_pos)

  ci <- exactci(tn, (tn + fp), conf.level = 0.95)

  return(ci)
}

calc_stage2b_sens_ci95 <- function(df) {
  df1 <- df %>%
    filter(outcome == "Mets") %>%
    mutate(classi = if_else(stage_2b == "Surveillance", "positive", "negative"))

  df_pos <- df1 %>%
    filter(classi == "positive")

  df_neg <- df1 %>%
    filter(classi == "negative")

  tp <- nrow(df_pos)
  fn <- nrow(df_neg)

  ci <- exactci(tp, (tp + fn), conf.level = 0.95)

  return(ci)
}

calc_stage2b_spec_ci95 <- function(df) {
  df1 <- df %>%
    filter(outcome == "No mets") %>%
    mutate(classi = if_else(stage_2b == "Surveillance", "positive", "negative"))

  df_pos <- df1 %>%
    filter(classi == "positive")

  df_neg <- df1 %>%
    filter(classi == "negative")

  tn <- nrow(df_neg)
  fp <- nrow(df_pos)

  ci <- exactci(tn, (tn + fp), conf.level = 0.95)

  return(ci)
}

calc_stage3a_sens_ci95 <- function(df) {
  df1 <- df %>%
    filter(outcome == "Mets") %>%
    mutate(classi = if_else(stage_3a == "Surveillance", "positive", "negative"))

  df_pos <- df1 %>%
    filter(classi == "positive")

  df_neg <- df1 %>%
    filter(classi == "negative")

  tp <- nrow(df_pos)
  fn <- nrow(df_neg)

  ci <- exactci(tp, (tp + fn), conf.level = 0.95)

  return(ci)
}

calc_stage3a_spec_ci95 <- function(df) {
  df1 <- df %>%
    filter(outcome == "No mets") %>%
    mutate(classi = if_else(stage_3a == "Surveillance", "positive", "negative"))

  df_pos <- df1 %>%
    filter(classi == "positive")

  df_neg <- df1 %>%
    filter(classi == "negative")

  tn <- nrow(df_neg)
  fp <- nrow(df_pos)

  ci <- exactci(tn, (tn + fp), conf.level = 0.95)

  return(ci)
}

calc_stage3b_sens_ci95 <- function(df) {
  df1 <- df %>%
    filter(outcome == "Mets") %>%
    mutate(classi = if_else(stage_3b == "Surveillance", "positive", "negative"))

  df_pos <- df1 %>%
    filter(classi == "positive")

  df_neg <- df1 %>%
    filter(classi == "negative")

  tp <- nrow(df_pos)
  fn <- nrow(df_neg)

  ci <- exactci(tp, (tp + fn), conf.level = 0.95)

  return(ci)
}

calc_stage3b_spec_ci95 <- function(df) {
  df1 <- df %>%
    filter(outcome == "No mets") %>%
    mutate(classi = if_else(stage_3b == "Surveillance", "positive", "negative"))

  df_pos <- df1 %>%
    filter(classi == "positive")

  df_neg <- df1 %>%
    filter(classi == "negative")

  tn <- nrow(df_neg)
  fp <- nrow(df_pos)

  ci <- exactci(tn, (tn + fp), conf.level = 0.95)

  return(ci)
}

calc_score_sens_ci95 <- function(df, cut) {
  df1 <- df %>%
    filter(outcome == "Mets") %>%
    mutate(classi = if_else( score >= cut, "positive", "negative"))

  df_pos <- df1 %>%
    filter(classi == "positive")

  df_neg <- df1 %>%
    filter(classi == "negative")

  tp <- nrow(df_pos)
  fn <- nrow(df_neg)

  ci <- exactci(tp, (tp + fn), conf.level = 0.95)

  return(ci)
}

calc_score_spec_ci95 <- function(df, cut) {
  df1 <- df %>%
    filter(outcome == "No mets") %>%
    mutate(classi = if_else(score >= cut, "positive", "negative"))

  df_pos <- df1 %>%
    filter(classi == "positive")

  df_neg <- df1 %>%
    filter(classi == "negative")

  tn <- nrow(df_neg)
  fp <- nrow(df_pos)

  ci <- exactci(tn, (tn + fp), conf.level = 0.95)

  return(ci)
}

calc_m3tnm_sens_ci95 <- function(df) {
  df1 <- df %>%
    filter(outcome == "Mets") %>%
    mutate(classi = if_else(m3tnm_str == "Surveillance", "positive", "negative"))

  df_pos <- df1 %>%
    filter(classi == "positive")

  df_neg <- df1 %>%
    filter(classi == "negative")

  tp <- nrow(df_pos)
  fn <- nrow(df_neg)

  ci <- exactci(tp, (tp + fn), conf.level = 0.95)

  return(ci)
}

calc_m3tnm_spec_ci95 <- function(df) {
  df1 <- df %>%
    filter(outcome == "No mets") %>%
    mutate(classi = if_else(m3tnm_str == "Surveillance", "positive", "negative"))

  df_pos <- df1 %>%
    filter(classi == "positive")

  df_neg <- df1 %>%
    filter(classi == "negative")

  tn <- nrow(df_neg)
  fp <- nrow(df_pos)

  ci <- exactci(tn, (tn + fp), conf.level = 0.95)

  return(ci)
}
# 95% CI for medians from https://www.statology.org/confidence-interval-for-median/ ----------------------
calc_med_ci <- function(df) {
  n <- nrow(df)
  q <- 0.5
  z <- 1.96

  med <- median(df$t_to_mets)
  j <- ceiling(n*q - 1.96*sqrt((n*q)*(1-q)))
  k <- ceiling(n*q + 1.96*sqrt((n*q)*(1-q)))

  df <- df %>% arrange(t_to_mets)

  lower <- df$t_to_mets[j]
  upper <- df$t_to_mets[k]
  return(c(lower, upper))
}

# General sensitivity, specificity and confidence interval functions
calc_sens <- function(df, col) {
  mycol <- enquo(col)
  df1 <- df %>%
    filter(outcome == "Mets")

  df_pos <- df1 %>%
    filter(UQ(mycol) == "Surveillance")

  df_neg <- df1 %>%
    filter(UQ(mycol) == "No surveillance")

  tp <- nrow(df_pos)
  fn <- nrow(df_neg)

  sensitivity <- tp / (tp + fn)

  return(sensitivity)
}


calc_spec <- function(df, col) {
  mycol <- enquo(col)
  df1 <- df %>%
    filter(outcome == "No mets")

  df_pos <- df1 %>%
    filter(UQ(mycol) == "Surveillance")

  df_neg <- df1 %>%
    filter(UQ(mycol) == "No surveillance")

  tn <- nrow(df_neg)
  fp <- nrow(df_pos)

  specificity <- tn / (tn + fp)

  return(specificity)
}

calc_spec_ci95 <- function(df, col) {
  mycol <- enquo(col)
  df1 <- df %>%
    filter(outcome == "No mets")

  df_pos <- df1 %>%
    filter(UQ(mycol) == "Surveillance")

  df_neg <- df1 %>%
    filter(UQ(mycol) == "No surveillance")

  tn <- nrow(df_neg)
  fp <- nrow(df_pos)

  ci <- exactci(tn, (tn + fp), conf.level = 0.95)

  return(ci)
}

calc_sens_ci95 <- function(df, col) {
  mycol <- enquo(col)
  df1 <- df %>%
    filter(outcome == "Mets")

  df_pos <- df1 %>%
    filter(UQ(mycol) == "Surveillance")

  df_neg <- df1 %>%
    filter(UQ(mycol) == "No surveillance")

  tp <- nrow(df_pos)
  fn <- nrow(df_neg)

  ci <- exactci(tp, (tp + fn), conf.level = 0.95)

  return(ci)
}

# Calculate cost savings ----------------------------
calc_scansave <- function(n, prevalence, spec1, spec2) {

}

# Creates a plot label from the roc analysis output
auc_label <- function(roc_results) {

  auc <- round(ciAUC(roc_results)$AUC, 2)
  cilow <- round(ciAUC(roc_results)$lower, 2)
  ciup <- round(ciAUC(roc_results)$upper, 2)

  label <- paste("AUC =", auc, "(", cilow, "-" ,ciup, ")")
  return(label)
}
