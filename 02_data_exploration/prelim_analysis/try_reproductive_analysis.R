# Script for testing out analyses of seeds per (reproductive) plant.
# Reading in demo data (with imputed survival) and combined phenology, demo, and
# seed set data (see script 01_data_cleaning_reconcile_demo_phen_seed_umbels.R)

##### Load products

# Read packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(glmmTMB)
# library(lme4)

# Clear namespace
rm(list = ls())

# Read in raw datasets

# Demo (imputed survival, all years)
demo = read.csv('01_data_cleaning/out/demo_imputed_survival.csv')

# Phenology + demography (all matched plants, 2021 and later)
phen.demo = read.csv('01_data_cleaning/out/phen_demo_for_umbel_survival.csv')

# Phenology + demography + seed set
# (not finished yet)


##### Process data as needed

### Demography:

demo.for.umbel.counts = merge(
    x = demo %>%
      # Filter out only surviving plants
      filter(surv) %>%
      # Filter out 2016 plants (size measurements are not reliable)
      filter(Year > 2016) %>%
      # Select columns for year-of data
      select(Year, plantid, Plot, trt, No.umbels),
    y = demo %>%
      # Filter out only surviving plants
      filter(surv) %>%
      # Filter out 2016 plants (size measurements are not reliable)
      filter(Year > 2016) %>%
      # Add column for year matching
      mutate(Year.match = Year + 1) %>% 
      # Select relevant columns
      select(Year.match, plantid, No.umbels, No.leaves, Leaf.length),
    by.x = c('Year', 'plantid'), by.y = c('Year.match', 'plantid'),
    suffixes = c('', '.prev')
  )

head(demo.for.umbel.counts)
nrow(demo.for.umbel.counts) # 3208 observations
length(unique(demo.for.umbel.counts$plantid)) # 842 plants

