library(tidyverse)
library(here)

# Read in data
data <- readRDS(here("data", "data_lp.RDS"))
data_ngs <- readRDS(here("data", "data_lp_ngs.RDS"))

#Reshape the data for roc analysis
roc_data <- data %>%
  mutate(outcome = case_when(((status == "D" & causeofdeath == "MM" & lfuyears <= 5) |
                                (!is.na(onsetofmm) & lfumets <= 5)) ~ "Mets",
                             (lfuyears >= 5 & (is.na(onsetofmm) | lfumets >= 5)) ~ "No mets"),
         outcome_binary = if_else(outcome == "Mets", 1, 0)) %>%
  filter(!is.na(outcome)) %>%
  select(-oob, -diagnosis, -dateofpm, -secondarytreatment, -followup, -genetictest, -chr1,
         -ch6p, -ch6q, -chr8p, -nbap1, -comments, -onsetofmm, -lfumonths, -lfuyears,
         -lfumets, -dateofpmyear, -secondary, -lbd_cat, -age_cat, -t, -t2, -with3)

saveRDS(roc_data, here("data","roc_data.RDS"))

roc_data_ngs <- data_ngs %>%
  mutate(outcome = case_when(((status == "D" & causeofdeath == "MM" & lfuyears <= 5) |
                                (!is.na(onsetofmm) & lfumets <= 5)) ~ "Mets",
                             (lfuyears >= 5 & (is.na(onsetofmm) | lfumets >= 5)) ~ "No mets"),
         outcome_binary = if_else(outcome == "Mets", 1, 0)) %>%
  filter(!is.na(outcome)) %>%
  select(-oob, -diagnosis, -dateofpm, -secondarytreatment, -followup, -genetictest, -chr1,
         -ch6p, -ch6q, -chr8p, -nbap1, -comments, -onsetofmm, -lfumonths, -lfuyears,
         -lfumets, -dateofpmyear, -secondary, -lbd_cat, -age_cat, -t, -t2, -with3)

saveRDS(roc_data_ngs, here("data","roc_data_ngs.RDS"))
