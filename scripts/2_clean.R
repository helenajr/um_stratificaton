library(tidyverse)
library(lubridate)

# Excludes those treated pre-2007, diagnoses other than choroidal melanoma, those with mets at primary treatment,
# those with insufficient information to use lumpo through the website and corrects typos

# Read in data
data <- readRDS(here("data", "assembled_data.RDS"))

# Corrects all identified typos ----------------------------------------------------------
data <- data %>%
  mutate(causeofdeath = if_else(rowid == 1958, "Other", causeofdeath),
         dateofpm = if_else(rowid == 2187, "04/04/2016", dateofpm),
         dateofpm = if_else(rowid == 1344, "11/09/2012", dateofpm),
         ageatpm = ifelse(rowid == 1066, 22, ageatpm),
         gender = if_else(rowid == 1066, "M", gender),
         ageatpm = ifelse(rowid == 1369, 57, ageatpm),
         gender = if_else(rowid == 1369, "F", gender),
         ageatpm = ifelse(rowid == 1606, 66, ageatpm),
         uh = ifelse(rowid == 2142, 8.8, uh)) %>%
  filter(rowid != 235)

# Puts columns in correct format, correct follow up date --------------------------------
data <- data %>%
  rename(followup = dateoflastfollowupdeath,
         mitcount = mitoticcount40hpf) %>%
  mutate(dateofpm = dmy(dateofpm),
         followup = dmy(followup),
         followup = if_else(followup == as_date("2020-05-01") & status != "D", as_date("2018-09-01"), followup),
         followup = if_else(followup < dateofpm, dateofpm, followup),
         ageatpm = as.numeric(ageatpm),
         lbd = as.numeric(lbd),
         uh = as.numeric(uh),
         mitcount = as.numeric(mitcount),
         onsetofmm = dmy(onsetofmm),
         lfumonths = dateofpm %--% followup / dmonths(1),
         lfuyears = dateofpm %--% followup /dyears(1),
         lfumets = dateofpm %--% onsetofmm / dyears(1),
         dateofpmyear = year(dateofpm)) %>%
  mutate(primarytreatment = str_to_lower(primarytreatment),
         primarytreatment = str_replace_all(primarytreatment, " ", ""),
         primarytreatment = fct_collapse(primarytreatment, `endoresection+prxt` = c("endoresectionandprxt",
                                                                                    "endoresectionandplaque",
                                                                                    "endoresection+prxt",
                                                                                    "endoresection+plaque",
                                                                                    "endoresection")),
         primarytreatment = fct_collapse(primarytreatment, `localresection+prxt` = c("localresection",
                                                                                     "localresection+prxt")),
         secondary = if_else(is.na(secondarytreatment), "None", "Secondary Treatment")) %>%
  mutate(causeofdeath = if_else(status == "D" & is.na(causeofdeath), "Unknown", causeofdeath),
         causeofdeath = as.factor(causeofdeath),
         causeofdeath = fct_collapse(causeofdeath, MM = c("MM", "Brain mets"),
                                     Other = c("Other", "Myocardial infarction")),
         causeofdeath = factor(causeofdeath, levels = c("MM", "Other", "Unknown"))) %>%
  mutate(chr3 = na_if(chr3, "U"),
         chr3 = if_else(chr3 == "G", "N", chr3),
         chr8q = na_if(chr8q, "U"),
         chr8q = if_else(chr8q == "L", "N", chr8q),
         with3 = if_else(is.na(chr3), "No result", "Chr3 result"))

# Applies exclusion criteria --------------------------------------------------------------------
data <- data %>%
  filter(!is.na(lbd) & !is.na(uh) & !is.na(cbi) &
           !is.na(eospread) & diagnosis == "Choroidal melanoma" & !is.na(ageatpm) & !is.na(gender)) %>%
  filter(dateofpmyear > 2006,
         (is.na(lfumets) | lfumets > 0)) %>%
  filter(!primarytreatment %in% c("ttt", "avastin", "vitrectomy", "observation", "excisionbiopsy",
                                  "treatmentelsewhere")) #These are not standard treatments for CM


# Save the dataset
saveRDS(data, here("data", "processed_data.RDS"))
