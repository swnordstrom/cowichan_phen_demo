### Script for cleaning/processing seed count data, 2021-2023
### (note: as of init date, 23 Oct 2023, 2023 data has not been entered yet)
###
### init SN 23 Oct 2023

######################################################################
##### Load data and initial assessment 
######################################################################

# Load packages
library(ggplot2)
library(dplyr)
library(tidyr)

# Clear namespace
rm(list = ls())

data.dir = '00_raw_data/lomatium_demography/'

grep('Seed', dir(data.dir), value = TRUE)
# 2023 missing (needs to be entered)
# for now, just focus on 2021-2022 (add 2023 later)

seed.counts.list = paste0(data.dir, grep('Seed', dir(data.dir), value = TRUE)) %>%
  lapply(read.csv)

lapply(seed.counts.list, head)
# different notes...

# 2022 plants - no coordinates listed (hope there are no dupe umbel numbers?)

seed.counts.list[[2]] %>%
  group_by(plot, tag, umble.no) %>%
  filter(n() > 1)

# well... shit.

seed.counts.list[[2]] %>%
  mutate(record.no = 1:nrow(.)) %>%
  group_by(plot, tag, umble.no) %>%
  filter(n() > 1)
# what...?

seed.counts.list[[2]] %>%
  mutate(record.no = 1:nrow(.)) %>%
  filter(tag %in% c(3519, 3149, 3076))
# hmm... going to assume these are the same 

seed.counts.list[[2]] %>%
  group_by(plot, tag) %>%
  mutate(umble.no.fix = ifelse(duplicated(umble.no), max(umble.no) + 1, umble.no)) %>%
  filter(umble.no != umble.no.fix)
# looks good to me

# Well... might as well do these separately...

######################################################################
##### Processing 
######################################################################

### 2021

proc.2021 = seed.counts.list[[1]]

head(proc.2021)

# combine notes
proc.2021 = proc.2021 %>%
  mutate(notes = paste(notes, X, X.1, sep = ';')) %>%
  select(-c(X, X.1))

head(proc.2021)

# Add unique plant ID
proc.2021 = proc.2021 %>% 
  mutate(plantid = paste0(tag, '_', plot, '_', xcoor, ycoor)) %>%
  select(-c(xcoor, ycoor))

# Check for duplicated umbel numbers
proc.2021 %>%
  group_by(plantid, umble.no) %>%
  filter(n() > 1)
# None - good!

str(proc.2021)

# What exactly is in "seed status"?
unique(proc.2021$seed.status)
# seed phenology... maybe keep in mind for later, but remove for now

### 2022

proc.2022 = seed.counts.list[[2]]

head(proc.2022)

# Combine notes column and add unique plant ID

proc.2022 = proc.2022 %>%
  # Combining notes column
  mutate(notes = paste(notes, X, sep = ';')) %>%
  select(-X) %>%
  # Add plant ID (no coords this year)
  mutate(plantid = paste(tag, plot, sep = '_'))

# Fix duplicated umbel counts (using code above)
proc.2022 = proc.2022 %>%
  group_by(plantid) %>%
  mutate(umble.no = ifelse(duplicated(umble.no), max(umble.no) + 1, umble.no)) %>%
  ungroup()

proc.2022 %>% group_by(plantid, umble.no) %>% filter(n() > 1)

head(proc.2022)
str(proc.2022)

### 2023 code can go here


######################################################################
##### Combining
######################################################################

rbind(
  proc.2021 %>% select(year, plantid, plot, tag, umble.no, no.seeds, notes),
  proc.2022 %>% select(year, plantid, plot, tag, umble.no, no.seeds, notes)
) %>%
  # Rename umbel column...
  rename(umbel.no = umble.no) %>%
  # Export
  write.csv(
    file = '01_data_cleaning/out/seed_counts_all_cleaned.csv',
    row.names = FALSE
  )

