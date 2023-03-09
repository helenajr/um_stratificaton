library(tidyverse)
library(here)

#This function removes column 27 if present (otherise causes problems) -----
cond_remove <- function(df) {
  colnames <- colnames(df)
  if (any(str_detect(colnames, "27"))) {
    df <- df %>% select(-"27")
  } else {
    df <- df
  }
  return(df)
}

# Reads in source data files ---------------------------------------------
files <- list.files(here("data", "csvs"))

# Corrects column naming and puts data from all csvs into a single data frame
df <- tibble()
for (i in seq_along(files)) {
  new_df <- read_csv(str_c("data/csvs/", files[i])) %>%
    select(1:27) %>%
    rename_with(~str_replace_all(., "\n|/|\\(mm\\)|cells|ID|No|\\.| ", ""),
                everything()) %>%
    rename_with(~str_to_lower(.)) %>%
    cond_remove(.) %>%
    mutate(across(where(is.numeric), ~as.character(.)))

  df <- bind_rows(df, new_df)
}

df %>%
  rowid_to_column(var = "rowid") %>%
  saveRDS(here("data", "assembled_data.RDS"))


