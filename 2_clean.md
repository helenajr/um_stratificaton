2_Clean
================
Helena
2023-03-09

## Clean

There is quite a bit of cleaning to do before we get to the analysis.
Cleaning refers to spotting errors and removing them and applying
exclusion criteria, so as far as is possible the dataset is correct and
contains only the appropriate people.

``` r
#Load the libraries and load the data we assembled in 1_assemble
library(tidyverse)
library(lubridate)
library(here)

data <- readRDS(here("data", "assembled_data.RDS"))
```

Correct all identified typos. These are things like someone being 660
years old or having their primary treatment 200 years in the future.
There will always be typos when data are manually entered into a
spreadsheet. Problems like these stand out when doing the exploratory
data analysis in 3_CheckQuality. Where I identified issues these were
checked against the original patient record and then corrected in the
following chunk.

``` r
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
```

The next sections do various bits of tidying up. First some columns are
renamed to make them easier to type.

``` r
data <- data %>%
  rename(followup = dateoflastfollowupdeath,
         mitcount = mitoticcount40hpf)
```

Next the date columns are converted from character to date class and
number columns are made numeric. This is important for calculatiing
intervals and plotting distributions.

``` r
data <- data %>%
  mutate(dateofpm = dmy(dateofpm),
         followup = dmy(followup),
         ageatpm = as.numeric(ageatpm),
         lbd = as.numeric(lbd),
         uh = as.numeric(uh),
         mitcount = as.numeric(mitcount),
         onsetofmm = dmy(onsetofmm))
```

Also there are some quirks specific to this dataset that needed
correction - they are slightly difficult to explain, but feel free to
ask if you need more info. The dataset used death notifications from the
national cancer registry as one source of follow up data and it was
thought that the last update from the registry was received in
“2020-05-01” (so if people were not listed they were assumed alive on
that date), but when this was checked the real date was “2018-09-01”.
Therefore all people listed as “Alive” with last follow up on that date,
needed the date changed. Also some people did not have any follow up
listed (normally overseas patients) which may make it appear as though
their follow-up occurred before their primary treatment (normally the
day before when they arrived at hospital). In these cases I set follow
up to equal the day of their primary treatment.

``` r
data <- data %>%
  mutate(followup = if_else(followup == as_date("2020-05-01") & status != "D", as_date("2018-09-01"), followup),
         followup = if_else(followup < dateofpm, dateofpm, followup))
```

I also made some new columns I needed for data exploration or analysis.
I calculated length of follow up and time to detection of mets in either
years or months (using my favourite operator %–%). I also extracted the
year of primary treatment into a separate column.

``` r
data <- data %>%
  mutate(lfumonths = dateofpm %--% followup / dmonths(1),
         lfuyears = dateofpm %--% followup /dyears(1),
         lfumets = dateofpm %--% onsetofmm / dyears(1),
         dateofpmyear = year(dateofpm))
```

Next, some fields are factors (categories). The primary treatment and
cause of death fields have just been manually enetered as free text and
therefore contain some inconsistencies. Here I combine bunches of
entries that all refer to the same thing into the same category. So we
are just left with the categories which are useful for analysis.

``` r
data <- data %>%
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
         causeofdeath = factor(causeofdeath, levels = c("MM", "Other", "Unknown")))
```

Lastly, tidying up the genetic columns. For our purposes chr3 should
either be “N” (normal), “L” (loss) or NA and chr8q should either be “G”
(gain), “N” (normal) or NA. This sorts that out.

``` r
data <- data %>%
  mutate(chr3 = na_if(chr3, "U"),
         chr3 = if_else(chr3 == "G", "N", chr3),
         chr8q = na_if(chr8q, "U"),
         chr8q = if_else(chr8q == "L", "N", chr8q),
         with3 = if_else(is.na(chr3), "No result", "Chr3 result"))
```

Finally, apply some exclusion criteria. The analysis only includes
people having their primary treatment 2007-2016 and those people who
have the minimum data required for the LUMPOIII website. It excludes
people who already had detectable mets at primary diagnosis and those
who had non-standard primary treatments.

``` r
data <- data %>%
  filter(!is.na(lbd) & !is.na(uh) & !is.na(cbi) &
           !is.na(eospread) & diagnosis == "Choroidal melanoma" & !is.na(ageatpm) & !is.na(gender)) %>%
  filter(dateofpmyear > 2006,
         (is.na(lfumets) | lfumets > 0)) %>%
  filter(!primarytreatment %in% c("ttt", "avastin", "vitrectomy", "observation", "excisionbiopsy",
                                  "treatmentelsewhere")) 
```

``` r
# Save the dataset
saveRDS(data, here("data", "processed_data.RDS"))
```

The dataset should now look something like this (the following data is
made up):

    ## Rows: 3
    ## Columns: 34
    ## $ rowid              <dbl> 1, 3, 4
    ## $ oob                <chr> "001/07", NA, "006/07"
    ## $ ageatpm            <dbl> 45, 81, 90
    ## $ gender             <chr> "F", "M", "M"
    ## $ diagnosis          <chr> "Choroidal melanoma", "Choroidal melanoma", "Choroidal melanoma"
    ## $ primarytreatment   <fct> prxt, biopsy, enucleation
    ## $ dateofpm           <date> 2007-01-01, 2007-07-16, 2007-04-15
    ## $ secondarytreatment <lgl> NA, NA, NA
    ## $ epithelioid        <chr> "Y", "N", NA
    ## $ loops              <chr> "Y", "Y", NA
    ## $ mitcount           <dbl> 6, 1, NA
    ## $ lbd                <dbl> 15.1, 18.1, 4.3
    ## $ uh                 <dbl> 2.1, 1.1, 2.1
    ## $ cbi                <chr> "N", "N", "Y"
    ## $ eospread           <chr> "N", "N", "N"
    ## $ status             <chr> "A", "D", "D"
    ## $ causeofdeath       <fct> NA, MM, Other
    ## $ followup           <date> 2016-03-01, 2013-09-09, 2014-05-06
    ## $ genetictest        <chr> "MLPA", "MLPA", NA
    ## $ chr1               <chr> "N", "N", "N"
    ## $ chr3               <chr> "L", "N", "N"
    ## $ ch6p               <chr> "N", "N", "N"
    ## $ ch6q               <chr> "N", "N", "N"
    ## $ chr8p              <chr> "N", "N", "N"
    ## $ chr8q              <chr> "G", "N", "N"
    ## $ nbap1              <lgl> NA, NA, NA
    ## $ comments           <chr> "No comments", "Primary skin tumour", NA
    ## $ onsetofmm          <date> NA, 2013-09-08, NA
    ## $ lfumonths          <dbl> 109.96304, 73.82341, 84.69815
    ## $ lfuyears           <dbl> 9.163587, 6.151951, 7.058179
    ## $ lfumets            <dbl> NA, 6.149213, NA
    ## $ dateofpmyear       <dbl> 2007, 2007, 2007
    ## $ secondary          <chr> "None", "None", "None"
    ## $ with3              <chr> "Chr3 result", "Chr3 result", "Chr3 result"
