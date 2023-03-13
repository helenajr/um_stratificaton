# This script adds the 5 year metastatic mortality scores from the lumpo model,
# the parsimonious model
library(tidyverse)
library(here)
source(here("models", "lumpo3CumInc.R"))

# Read in data
data <- readRDS(here("data", "processed_data.RDS"))

# Add LUMPO scores ---------------------------------------------------------------------
#Select only the necessary columns and put in the correct order to input into the model
data_temp <- data %>% select(ageatpm, gender, lbd, cbi, eospread, uh, epithelioid, loops, mitcount, chr3, chr8q)

#Convert vectors to correct data types
data_temp <- data_temp %>% mutate(gender = case_when(gender == "M" ~ 1,
                                                     gender == "F" ~ 0),
                                  cbi = if_else(cbi == "Y", 1, 0),
                                  eospread = if_else(eospread == "Y", 1, 0),
                                  epithelioid = if_else(epithelioid == "Y", 1, 0),
                                  loops = case_when(loops == "Y" ~ 1,
                                                    loops == "N" ~ 0),
                                  mitcount = if_else(mitcount < 2, 1, if_else(mitcount < 4, 2, if_else(mitcount < 8, 3, 4))),
                                  chr3 = case_when(chr3 == "L" ~ 1,
                                                   chr3 %in% c("G", "N") ~ 0),
                                  chr8q = case_when(chr8q == "G" ~ 1,
                                                    chr8q %in% c("L", "N") ~ 0))

MM <- vector("double", length = nrow(data))

for (i in seq_along(MM)) {
  #Extract values from a row into a list
  test <- slice(data_temp, i)
  val <- c(test[[1]], test[[2]], test[[3]], test[[4]], test[[5]], test[[6]], test[[7]],
           test[[8]], test[[9]], test[[10]], test[[11]])

  #Use values as input to LUMPO function
  result <- lumpo3CumInc(val)

  #Select correct output of LUMPO function
  MM[i] <- result[[4]][753]
}

#Round to 3 sf
MM <- signif(MM, 3)

#Add vectors to dataframe
data_l <- cbind(data, MM)

# Add parsimonious (simple) model ------------------------------------------
parsim <- read_csv(here("models", "parsim.csv")) %>%
  mutate(chr3 = case_when(chr3 == "disomy" ~ "N",
                          chr3 == "monosomy" ~ "L")) %>%
  rename(lbd_cat = lbd)

data_lp <- data_l %>%
  mutate(lbd_cat = case_when(lbd < 10.1 ~ 1,
                             lbd >= 10.1 & lbd < 12.1 ~ 2,
                             lbd >= 12.1 & lbd < 14.1 ~ 3,
                             lbd >= 14.1 & lbd < 16.1 ~ 4,
                             lbd >= 16.1 & lbd < 18.1 ~ 5,
                             lbd >= 18.1 & lbd < 28.1 ~ 6),
         age_cat = if_else(ageatpm < 81, "young", "old")) %>%
  left_join(parsim)

### Add classification under the AJCC system
data_lp <- data_lp %>%
  mutate(t = case_when((lbd <= 9  & uh <= 6) | (lbd > 9 & lbd <= 12 & uh <=3) ~ "T1",
                       (lbd <= 9  & uh > 6 & uh <= 9) | (lbd > 9 & lbd <= 12 & uh > 3 & uh <=9) |
                         (lbd > 12 & lbd <= 15 & uh <= 6) | (lbd > 15 & lbd <= 18 & uh <= 3) ~ "T2",
                       (lbd > 3 & lbd <= 12 & uh > 9 & uh <= 15) | (lbd > 12 & lbd <= 15 & uh > 6 & uh <= 15) |
                         (lbd > 15 & lbd <= 18 & uh > 3 & uh <= 12) ~ "T3",
                       (lbd > 12 & lbd <= 15 & uh > 15) | (lbd > 15 & lbd <=18 & uh > 12) | lbd > 18 ~ "T4"),
         t2 = case_when(t == "T1" & cbi == "N" & eospread == "N" ~ "T1a",
                        t == "T1" & cbi == "Y" & eospread == "N" ~ "T1b",
                        t == "T1" & cbi == "N" & eospread == "Y" ~ "T1c",
                        t == "T1" & cbi == "Y" & eospread == "Y" ~ "T1d",
                        t == "T2" & cbi == "N" & eospread == "N" ~ "T2a",
                        t == "T2" & cbi == "Y" & eospread == "N" ~ "T2b",
                        t == "T2" & cbi == "N" & eospread == "Y" ~ "T2c",
                        t == "T2" & cbi == "Y" & eospread == "Y" ~ "T2d",
                        t == "T3" & cbi == "N" & eospread == "N" ~ "T3a",
                        t == "T3" & cbi == "Y" & eospread == "N" ~ "T3b",
                        t == "T3" & cbi == "N" & eospread == "Y" ~ "T3c",
                        t == "T3" & cbi == "Y" & eospread == "Y" ~ "T3d",
                        t == "T4" & cbi == "N" & eospread == "N" ~ "T4a",
                        t == "T4" & cbi == "Y" & eospread == "N" ~ "T4b",
                        t == "T4" & cbi == "N" & eospread == "Y" ~ "T4c",
                        t == "T4" & cbi == "Y" & eospread == "Y" ~ "T4d"),
         stage = case_when(t2 == "T1a" ~ "I",
                           t2 %in% c("T1b", "T1c", "T1d", "T2a") ~ "IIA",
                           t2 %in% c("T2b", "T3a") ~ "IIB",
                           t2 %in% c("T2c", "T2d", "T3b", "T3c", "T4a") ~ "IIIA",
                           t2 %in% c("T3d", "T4b", "T4c") ~ "IIIB",
                           t2 == "T4d" ~ "IIIC"))

saveRDS(data_lp, here("data", "data_lp.RDS"))

# Selection of scores should be manually checked against website after running this script -----
test_l <- data_lp %>%
  select(ageatpm, gender, lbd, cbi, eospread, uh, epithelioid, loops, mitcount, chr3, chr8q, MM, score) %>%
  slice(100:106)

# Now run test_parsim_score to check join

# Do the same for the dataset without genetic data----------------------------------------

# Read in data
data_ngs <- readRDS(here("data", "processed_data.RDS")) %>%
  mutate(chr3 = NA,
         chr8q = NA) #This removes all the genetic data

# Add LUMPO scores ---------------------------------------------------------------------
#Select only the necessary columns and put in the correct order to input into the model
data_temp <- data_ngs %>% select(ageatpm, gender, lbd, cbi, eospread, uh, epithelioid, loops, mitcount, chr3, chr8q)

#Convert vectors to correct data types
data_temp <- data_temp %>% mutate(gender = case_when(gender == "M" ~ 1,
                                                     gender == "F" ~ 0),
                                  cbi = if_else(cbi == "Y", 1, 0),
                                  eospread = if_else(eospread == "Y", 1, 0),
                                  epithelioid = if_else(epithelioid == "Y", 1, 0),
                                  loops = case_when(loops == "Y" ~ 1,
                                                    loops == "N" ~ 0),
                                  mitcount = if_else(mitcount < 2, 1, if_else(mitcount < 4, 2, if_else(mitcount < 8, 3, 4))),
                                  chr3 = case_when(chr3 == "L" ~ 1,
                                                   chr3 %in% c("G", "N") ~ 0),
                                  chr8q = case_when(chr8q == "G" ~ 1,
                                                    chr8q %in% c("L", "N") ~ 0))

MM <- vector("double", length = nrow(data_ngs))

for (i in seq_along(MM)) {
  #Extract values from a row into a list
  test <- slice(data_temp, i)
  val <- c(test[[1]], test[[2]], test[[3]], test[[4]], test[[5]], test[[6]], test[[7]],
           test[[8]], test[[9]], test[[10]], test[[11]])
  
  #Use values as input to LUMPO function
  result <- lumpo3CumInc(val)
  
  #Select correct output of LUMPO function
  MM[i] <- result[[4]][753]
}

#Round to 3 sf
MM <- signif(MM, 3)

#Add vectors to dataframe
data_l <- cbind(data_ngs, MM)

# Add parsimonious (simple) model ------------------------------------------
parsim <- read_csv(here("models", "parsim.csv")) %>%
  mutate(chr3 = case_when(chr3 == "disomy" ~ "N",
                          chr3 == "monosomy" ~ "L")) %>%
  rename(lbd_cat = lbd)

data_lp <- data_l %>%
  mutate(lbd_cat = case_when(lbd < 10.1 ~ 1,
                             lbd >= 10.1 & lbd < 12.1 ~ 2,
                             lbd >= 12.1 & lbd < 14.1 ~ 3,
                             lbd >= 14.1 & lbd < 16.1 ~ 4,
                             lbd >= 16.1 & lbd < 18.1 ~ 5,
                             lbd >= 18.1 & lbd < 28.1 ~ 6),
         age_cat = if_else(ageatpm < 81, "young", "old")) %>%
  left_join(parsim)

### Add classification under the AJCC system
data_lp <- data_lp %>%
  mutate(t = case_when((lbd <= 9  & uh <= 6) | (lbd > 9 & lbd <= 12 & uh <=3) ~ "T1",
                       (lbd <= 9  & uh > 6 & uh <= 9) | (lbd > 9 & lbd <= 12 & uh > 3 & uh <=9) |
                         (lbd > 12 & lbd <= 15 & uh <= 6) | (lbd > 15 & lbd <= 18 & uh <= 3) ~ "T2",
                       (lbd > 3 & lbd <= 12 & uh > 9 & uh <= 15) | (lbd > 12 & lbd <= 15 & uh > 6 & uh <= 15) |
                         (lbd > 15 & lbd <= 18 & uh > 3 & uh <= 12) ~ "T3",
                       (lbd > 12 & lbd <= 15 & uh > 15) | (lbd > 15 & lbd <=18 & uh > 12) | lbd > 18 ~ "T4"),
         t2 = case_when(t == "T1" & cbi == "N" & eospread == "N" ~ "T1a",
                        t == "T1" & cbi == "Y" & eospread == "N" ~ "T1b",
                        t == "T1" & cbi == "N" & eospread == "Y" ~ "T1c",
                        t == "T1" & cbi == "Y" & eospread == "Y" ~ "T1d",
                        t == "T2" & cbi == "N" & eospread == "N" ~ "T2a",
                        t == "T2" & cbi == "Y" & eospread == "N" ~ "T2b",
                        t == "T2" & cbi == "N" & eospread == "Y" ~ "T2c",
                        t == "T2" & cbi == "Y" & eospread == "Y" ~ "T2d",
                        t == "T3" & cbi == "N" & eospread == "N" ~ "T3a",
                        t == "T3" & cbi == "Y" & eospread == "N" ~ "T3b",
                        t == "T3" & cbi == "N" & eospread == "Y" ~ "T3c",
                        t == "T3" & cbi == "Y" & eospread == "Y" ~ "T3d",
                        t == "T4" & cbi == "N" & eospread == "N" ~ "T4a",
                        t == "T4" & cbi == "Y" & eospread == "N" ~ "T4b",
                        t == "T4" & cbi == "N" & eospread == "Y" ~ "T4c",
                        t == "T4" & cbi == "Y" & eospread == "Y" ~ "T4d"),
         stage = case_when(t2 == "T1a" ~ "I",
                           t2 %in% c("T1b", "T1c", "T1d", "T2a") ~ "IIA",
                           t2 %in% c("T2b", "T3a") ~ "IIB",
                           t2 %in% c("T2c", "T2d", "T3b", "T3c", "T4a") ~ "IIIA",
                           t2 %in% c("T3d", "T4b", "T4c") ~ "IIIB",
                           t2 == "T4d" ~ "IIIC"))

saveRDS(data_lp, here("data", "data_lp_ngs.RDS"))
