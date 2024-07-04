library(ggplot2)
library(dplyr)
library(tidyr)

rm(list = ls())

phen24 = read.csv('00_raw_data/lomatium_demography/2024_Lomatium_Resurveys.csv')

# Look for NAs in counts
phen24 %>% filter(if_any(starts_with('no.'), ~ is.na(.)))
# none - good

# Look for NAs in anything else
phen24 %>% filter(if_any(everything(), ~ is.na(.)))
# cool - none

head(phen24)

phen24 %>%
  mutate(
    all.umbels = no.buds + no.flowers + no.flp + no.seeding + no.dead,
    liv.umbels = no.buds + no.flowers + no.flp + no.seeding
  ) %>%
  distinct(all.umbels)
# Looks fine with me
                  
str(phen24)
