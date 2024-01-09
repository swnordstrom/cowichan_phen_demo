library(ggplot2)
library(dplyr)
library(tidyr)

rm(list = ls())

# Processed demo
demo = read.csv('01_data_cleaning/out/demo_imputed_survival.csv') %>%
  # get rid of plants from before phen measurements and non-flowering plants
  filter(Year > 2020, No.umbels > 0 & !is.na(No.umbels)) %>%
  # remove the 'b's from plantids
  # (this will make reconciliation easier)
  mutate(plantid = gsub('b', '', plantid)) %>%
  mutate(finalid = plantid)

# Seed data
seed = read.csv('01_data_cleaning/out/demo_seed_v1.csv') %>%
  mutate(finalid = plantid.seed)

# Phen data
phen = read.csv('01_data_cleaning/out/phenology_all_ind_cleaned.csv') %>%
  mutate(finalid = plantid)

# Umbel data
mumb = read.csv('01_data_cleaning/out/multi_umbel_clean.csv') %>%
  # Turns out, there are 2020 umbel counts
  filter(Year > 2020) %>%
  mutate(finalid = plantid)

### All plantids

all.plantids = rbind(
  # Demo data
  demo %>% 
    distinct(plantid, Year) %>% 
    mutate(source = 'demo', is.in = TRUE),
  # Seed set data
  seed %>% 
    distinct(plantid.seed, year) %>% 
    mutate(source = 'seed', is.in = TRUE) %>% 
    rename(Year = year, plantid = plantid.seed),
  # Phenology data
  phen %>%
    distinct(plantid, year) %>%
    mutate(source = 'phen', is.in = TRUE) %>%
    rename(Year = year),
  # Umbel data
  # NOTE: this dataset only contains *multiple* umbels
  # so I will also add in the demo data for which there is one umbel and umbel
  # measurements
  mumb %>%
    distinct(plantid, Year) %>%
    mutate(source = 'mumb', is.in = TRUE),
  demo %>%
    filter(No.umbels %in% 1) %>%
    distinct(plantid, Year) %>%
    mutate(source = 'mumb', is.in = TRUE)
) %>%
  pivot_wider(
    id_cols = c(plantid, Year), names_from = source, 
    values_from = is.in, values_fill = FALSE
  )

head(all.plantids)
nrow(all.plantids)

all.plantids %>% 
  mutate(n.in = demo + seed + phen + mumb) %>%
  group_by(n.in) %>%
  summarise(n = n())
# lmao, jesus

all.plantids %>% 
  mutate(n.in = demo + seed + phen + mumb) %>%
  group_by(Year, n.in) %>%
  summarise(n = n())

### 2021: only a handful of plants here were included in phen
# so, there shouldn't be any plants from 2021 that are not listed in the
# phenology or seed set data
# I'll remove these plants here

all.plantids = all.plantids %>% filter(Year > 2021 | phen | seed)

### 2022: we don't have coordinates for these plants...
# we'll have to tackle these separately.

####################################################
####################################################
# Start reconciling ################################
####################################################
####################################################

##### Some rules of thumb:
# - it *would* be nice to keep the 'b' markers in here for reconciling all
# datasets... but that'll be tougher to do in code, and more work because it
# will involve editing more datasets
# so instead, I'll just remove the 'b's

##########
### 2022: coordinates not recorded
##########

all.plantids %>%
  filter(Year %in% 2022) %>%
  mutate(plantid = gsub('\\_[0-9]{1,2}[A-Z]$', '', plantid)) %>%
  group_by(plantid) %>%
  summarise(across(c(demo, seed, phen, mumb), any)) %>%
  mutate(n.in = demo + seed + phen + mumb) %>%
  group_by(n.in) %>%
  summarise(n = n())
# okay - much better!

# Plants to edit
all.plantids %>%
  filter(Year %in% 2022) %>%
  mutate(plantid = gsub('\\_[0-9]{1,2}[A-Z]$', '', plantid)) %>%
  group_by(plantid) %>%
  summarise(across(c(demo, seed, phen, mumb), any)) %>%
  filter(!(demo & seed & phen & mumb)) %>%
  print(n = nrow(.))

# 1 2230_4        TRUE  TRUE  FALSE TRUE 
# 2 2238_4        FALSE FALSE TRUE  FALSE

demo %>% filter(grepl('223[08]', plantid))
# no evidence of plant 2238 in demo
phen %>% filter(grepl('223[08]', plantid))
# must have been a misentry
phen = phen %>% mutate(finalid = gsub('2238', '2230', plantid))

# 3 3039_13       TRUE  FALSE FALSE TRUE
demo %>% filter(grepl('3039', plantid))
# note says 'was 3726'
# for purposes of reconciliation, I'll change it back to this
demo = demo %>% mutate(finalid = gsub('3039', '3726', finalid))

# 4 3076_15       FALSE FALSE TRUE  FALSE
phen %>% filter(grepl('3076', plantid))
# weird phen in 2022
# anyway...
demo %>% filter(grepl('3076', plantid)) # tag is not in any demo ids
demo %>% filter(grepl('3076', proc.note))
demo %>% filter(grepl('3076', demo.note))
# ah... tag misentry in phen
phen = phen %>% mutate(finalid = grepl('3076', '3706', finalid))

# 5 3143_9_140.95 TRUE  FALSE FALSE TRUE 
# 6 3151_3        FALSE FALSE TRUE  FALSE
# 7 3164_13       TRUE  TRUE  FALSE TRUE 
# 8 3164_6        TRUE  TRUE  FALSE TRUE 
# 9 3183_9_130.19 TRUE  FALSE FALSE TRUE 
# 10 3296_7        TRUE  TRUE  FALSE TRUE 
# 11 3298_7        FALSE FALSE TRUE  FALSE
# 12 3321_4        TRUE  FALSE FALSE TRUE 
# 13 3363_14       FALSE FALSE TRUE  FALSE
# 14 3367_14       TRUE  TRUE  FALSE TRUE 
# 15 3379_14       FALSE TRUE  TRUE  FALSE
# 16 3389_13       FALSE FALSE FALSE TRUE 
# 17 3456_6        TRUE  FALSE FALSE FALSE
# 18 3469_14       TRUE  TRUE  TRUE  FALSE
# 19 3481_13       TRUE  TRUE  TRUE  FALSE
# 20 3495_14       TRUE  TRUE  FALSE TRUE 
# 21 3542_1        TRUE  TRUE  TRUE  FALSE
# 22 3555_4        TRUE  FALSE FALSE TRUE 
# 23 3558_13       TRUE  FALSE FALSE TRUE 
# 24 3563_4        TRUE  FALSE FALSE TRUE 
# 25 3595_13       FALSE TRUE  FALSE FALSE
# 26 3614_13       FALSE FALSE TRUE  FALSE
# 27 3634_4        TRUE  FALSE FALSE TRUE 
# 28 3651_3        TRUE  TRUE  FALSE TRUE 
# 29 3706_15       TRUE  TRUE  FALSE TRUE 
# 30 3726_13       FALSE TRUE  TRUE  FALSE
# 31 3735_4        TRUE  FALSE FALSE TRUE 
# 32 3750_4        TRUE  FALSE FALSE TRUE 
# 33 3769_13       FALSE TRUE  TRUE  FALSE
# 34 3770_9        TRUE  FALSE FALSE TRUE 
# 35 3782_1        FALSE TRUE  TRUE  FALSE
# 36 3789_6        FALSE FALSE TRUE  FALSE
# 37 3813_15       FALSE FALSE TRUE  FALSE
# 38 3814_1        TRUE  TRUE  TRUE  FALSE
# 39 3831_15       TRUE  TRUE  FALSE TRUE 
# 40 3866_7        TRUE  TRUE  FALSE TRUE 
# 41 3880_13       FALSE FALSE FALSE TRUE 
# 42 3886_7        FALSE FALSE TRUE  FALSE
# 43 5037_7        TRUE  TRUE  FALSE TRUE 
# 44 5040_7        TRUE  TRUE  FALSE TRUE 
# 45 7526_1        FALSE TRUE  FALSE FALSE
# 46 7536_1        TRUE  FALSE TRUE  TRUE 
# 47 7537_13       TRUE  TRUE  TRUE  FALSE
# 48 7587_6        FALSE TRUE  TRUE  TRUE 