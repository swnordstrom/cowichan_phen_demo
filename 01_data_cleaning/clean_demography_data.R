# Script for cleaning and processing raw LoUt demography data for analysis
# initial version focusing only on <YEAR>_Lomatium_Demography_Data.csv files
# (will likely expand beyond that in this script in the future)
# SN - init 10 Oct 2023

library(dplyr)
library(tidyr)

# Considerations:

# Missing values for identifying info (Xcoor, Ycoor, Plot, Tag)
# Duplicate tags (even within plot)

rm(list = ls())

raw.demo.list = dir('00_raw_data/lomatium_demography/') %>%
  grep(pattern = '20\\d\\d\\_Lomatium\\_Demography_Data', x = ., value = TRUE) %>%
  paste0('00_raw_data/lomatium_demography/', .) %>%
  lapply(FUN = read.csv)

# Check for column consistency
sapply(raw.demo.list, ncol)
# some inconsistencies

# Look at column names
lapply(raw.demo.list, names)

# Get unique columns
shared.cols = Reduce(intersect, lapply(raw.demo.list, names))

# Combine common columns into one data frame
raw.demo = raw.demo.list %>%
  lapply(function(df) df[,shared.cols]) %>%
  do.call(what = rbind)

head(raw.demo)

# Create a proc.demo data frame for storing processed demo
proc.demo = raw.demo

##### Check column types

str(proc.demo)
# Stalk height, umbel number, umbel diameter, and year tag are characters...

# Fix stalk height
unique(proc.demo$Stalk_Height) # ah it's a random space in here
proc.demo %>% mutate(Stalk_Height = gsub('\\s', '', Stalk_Height)) %>% pull(Stalk_Height) %>% unique()
proc.demo = propc.demo %>%
  mutate(
    Stalk_Height = gsub('\\s', '', Stalk_Height),
    Stalk_Height = ifelse(grepl('NA', Stalk_Height), NA, Stalk_Height),
    Stalk_Height = as.numeric(Stalk_Height)
  )

# Fix number of umbels
unique(proc.demo$No.umbels) # 7.5 and 0N
raw.demo.list[[3]] %>% filter(No.umbels %in% '7.5') # hmm... no notes
raw.demo.list[[7]] %>% filter(No.umbels %in% '0N') # also no notes
# Looks from this like 0N should be zero (no umbel diameter) and 7.5 can be rounded down to 7
# Test to make sure this chunk of code works 
proc.demo %>%
  mutate(
    numbel = case_when(
      No.umbels %in% '0N' ~ '0',
      No.umbels %in% '7.5' ~ '7',
      .default = No.umbels
    ),
    # Convert to numeric
    numbel = as.numeric(numbel)
  ) %>%
  distinct(numbel, No.umbels)
# looks good!
proc.demo = proc.demo %>%
  mutate(
    No.umbels = case_when(
      No.umbels %in% '0N' ~ '0',
      No.umbels %in% '7.5' ~ '7',
      .default = No.umbels
    ),
    # Convert to numeric
    No.umbels = as.numeric(No.umbels)
  )

# Fix umbel diameters
unique(proc.demo$umbel.diam)
# what is SOS?
proc.demo %>% filter(umbel.diam %in% 'SOS') # ???

# From metadata file:
# All tags with multiple umbles in 2020, 2021, and 2022 have "SOS" (see other
# sheet) for umble diameter in the main demography CSV's. In the multiple umble
# CSV, each umble is a row and the diameter is listed seperately for each one.
# In the future, the "total" umble area can be added together for each tag and
# inserted back into the main demography sheet (need to do this in a script).
# -
# so it seems like this requires merging in the "multiple umble" datasheets

# Look at YrTag
unique(proc.demo$YrTag)
# question mark record seems to be the major offender
proc.demo %>% filter(grepl('\\?', YrTag))
# ahhhhhhh

##### Check for missing info

# Check for missing identifying info
raw.demo %>%
  filter(is.na(Tag) | is.na(Plot) | is.na(Xcoor) | is.na(Ycoor))

# One record - not sure it needs fixing?

