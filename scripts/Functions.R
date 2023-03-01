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

# 95% CI for medians from https://www.statology.org/confidence-interval-for-median/ ----------------------
calc_med_ci <- function(df) {
  n <- nrow(df)
  q <- 0.5
  z <- 1.96

  med <- median(df$t_to_mets)
  j <- ceiling(n * q - 1.96 * sqrt((n * q) * (1 - q)))
  k <- ceiling(n * q + 1.96 * sqrt((n * q) * (1 - q)))

  df <- df %>% arrange(t_to_mets)

  lower <- df$t_to_mets[j]
  upper <- df$t_to_mets[k]
  return(c(lower, upper))
}

# Calculate cost savings ----------------------------
calc_saving <- function(cost_scan, pop_size, spec1, spec2) {
  prev <- prevalence * pop_size
  tn1 <- (pop_size-prev)*spec1
  fp1 <- (pop_size-prev)-tn1
  tn2 <- (pop_size-prev)*spec2
  fp2 <- (pop_size-prev)-tn2
  
  cost_5yscan <- cost_scan*10
  
  people_saved <- round(fp1 - fp2)
  money_saved <- people_saved * cost_5yscan
  results <- list(people_saved, money_saved)
  return(results)
}

# Creates a plot label from the roc analysis output
auc_label <- function(roc_results) {
  auc <- round(ciAUC(roc_results)$AUC, 2)
  cilow <- round(ciAUC(roc_results)$lower, 2)
  ciup <- round(ciAUC(roc_results)$upper, 2)

  label <- paste("AUC =", auc, "(", cilow, "-", ciup, ")")
  return(label)
}
