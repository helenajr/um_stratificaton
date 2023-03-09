1_Assemble
================
Helena
2023-03-09

## Assembly

This markdown file explains the code in 1_assemble.R. The raw dataset is
a series of csv files, one csv file for each year (2007-2016). The first
step is to combine these sets into one big table. Like most data which
has been manually collated in excel there is inconsistency among the
column names which makes this tricky… The following script works
reproducibly with the csv files used for the paper but may need tweaking
to be effective on other files. I could have just manually fixed the
problems with the csvs but I thought this would be more fun (and a
better way of doing it if you had loads of csv files).

``` r
#Load required packages
library(tidyverse)
library(here)
```

A fun quirk of the csvs is some of them have an extra column, which when
reading the files in in the loop below, causes the datasets without it
to have a blank column called ‘27’, which causes problems when you try
and combine the datasets, so I wrote a function to deal with this if it
is present.

``` r
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
```

This makes a list of the data file names:

``` r
files <- list.files(here("data", "csvs"))
```

Make a blank tibble imaginatively called df.

``` r
df <- tibble()
```

This loop reads in the data, one file at a time, into a table called
new_df. It then selects the correct columns and alters the names of
columns to make sure they are consistent across files (if you are doing
this with different files this will take some experimentation), makes
all column names lower case, applies the function to get rid of 27 where
it is present and lastly add all the rows to the blank table created in
the last step.

``` r
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
```

Create a rowid and save the dataset as assembled_data.RDS

``` r
df %>%
  rowid_to_column(var = "rowid") %>%
  saveRDS(here("data", "assembled_data.RDS"))
```

The saved dataset should look like this (the following data is all made
up):

``` r
glimpse(fake_ad)
```

    ## Rows: 4
    ## Columns: 28
    ## $ rowid                   <chr> "1", "2", "3", "4"
    ## $ oob                     <chr> "001/07", "002/07", NA, "006/07"
    ## $ ageatpm                 <chr> "45", "67", "81", "90"
    ## $ gender                  <chr> "F", "M", "M", "M"
    ## $ diagnosis               <chr> "Choroidal melanoma", "BDUMP", "Choroidal melanoma", "Choroidal melanom…
    ## $ primarytreatment        <chr> "prxt", "endoresection + plaque", "biopsy", "enucleation"
    ## $ dateofpm                <chr> "01/01/2007", "05/10/2007", "16/07/2007", "15/04/2007"
    ## $ secondarytreatment      <lgl> NA, NA, NA, NA
    ## $ epithelioid             <chr> "Y", NA, "N", NA
    ## $ loops                   <chr> "Y", NA, "Y", NA
    ## $ mitoticcount40hpf       <chr> "6", NA, "1", NA
    ## $ lbd                     <chr> "15.1", "3.5", "18.1", "4.3"
    ## $ uh                      <chr> "2.1", "3.5", "1.1", "2.1"
    ## $ cbi                     <chr> "N", "N", "N", "Y"
    ## $ eospread                <chr> "N", "N", "N", "N"
    ## $ status                  <chr> "A", "A", "D", "D"
    ## $ causeofdeath            <chr> NA, NA, "MM", "Other"
    ## $ dateoflastfollowupdeath <chr> "01/03/2016", "08/09/2011", "09/09/2013", "06/05/2014"
    ## $ genetictest             <chr> "MLPA", NA, "MLPA", NA
    ## $ chr1                    <chr> "N", "N", "N", "N"
    ## $ chr3                    <chr> "L", "N", "N", "N"
    ## $ ch6p                    <chr> "N", "N", "N", "N"
    ## $ ch6q                    <chr> "N", "N", "N", "N"
    ## $ chr8p                   <chr> "N", "N", "N", "N"
    ## $ chr8q                   <chr> "G", "N", "N", "N"
    ## $ nbap1                   <lgl> NA, NA, NA, NA
    ## $ comments                <chr> "No comments", NA, "Primary skin tumour", NA
    ## $ onsetofmm               <chr> NA, NA, "08/09/2013", NA
