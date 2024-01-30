#################################################################################
# Script for combining and reconciling plant ids across datasets for recruitment
# estimation combining together these datasets:
# - demography (2021-2023) using the imputed survival dataset
#   (using this dataset because it's the one with the most processing)
# - seed set data per-umbel
#   (note this was reconciled with an *older* version of demography data)
# - phenology data per-umbel (a cleaned version)
# - the multiple umbel dataset (a cleaned version, with only 2021-2023)
# This script currently is only reconciling plantids among datasets - the
# corrected version is stored in a new `finalid` column.
# Once that's done I'll combine all four datasets into a single CSV for running
# some models of recruitment
# SN - init 9 Jan 2024
# #### Notes:
# - outstanding unreconciled issues, or other issues to be fixed, were flagged
# with the string UNSOLVED - I need to go back to these at some point
# in particular are there missing entries from plot 4, year 2022 in phen+/seed
# - how to handle prematurely dead umbels? ugh.
#################################################################################

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
phen = read.csv('01_data_cleaning/out/phenology_buds_deaths_all.csv') %>%
  # old csv (flowering phen only, fewer plants)
  # read.csv('01_data_cleaning/out/phenology_all_ind_cleaned.csv') %>%
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
phen = phen %>% mutate(finalid = gsub('2238', '2230', finalid))

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
phen = phen %>% mutate(finalid = gsub('3076', '3706', finalid))

# 5 3143_9_140.95 TRUE  FALSE FALSE TRUE 
demo %>% filter(grepl('3143', plantid))
phen %>% filter(grepl('3143', plantid))
seed %>% filter(grepl('3143', plantid.seed))
demo %>% filter(grepl('140\\.95', plantid))
mumb %>% filter(grepl('140\\.95', plantid))
demo = demo %>% mutate(finalid = gsub('140\\.95', '14T', finalid))
mumb = mumb %>% mutate(finalid = gsub('140\\.95', '14T', finalid))

# 6 3151_3        FALSE FALSE TRUE  FALSE
phen %>% filter(grepl('3151', plantid))
demo %>% filter(grepl('3151', plantid)) # nope, diff plant
demo %>% filter(grepl('3151', demo.note))
demo %>% filter(grepl('3151', proc.note)) # nothing here...
demo %>% filter(Year %in% 2022, grepl('\\_3\\_', plantid))
demo %>% filter(grepl('7537', plantid)) # nope, not this one
demo %>% filter(grepl('3651', plantid)) # not this one...
seed %>% filter(year %in% 2022, grepl('\\_3', plantid.seed))
# (solved below!)

# 7 3164_13       TRUE  TRUE  FALSE TRUE 
demo %>% filter(grepl('3164', plantid)) # looks tag misread?
phen %>% filter(grepl('3614', plantid))
demo %>% filter(grepl('3614\\_13', plantid)) # yep - it's this one
phen = phen %>% mutate(finalid = gsub('3614\\_13', '3164_13', finalid))
phen %>% filter(grepl('3614', plantid))
phen %>% filter(grepl('3614', finalid)) # good - fixed

# 8 3164_6        TRUE  TRUE  FALSE TRUE 
phen %>% filter(grepl('3164', plantid))
demo %>% filter(grepl('3164\\_6', plantid)) # some tag chicanery here
phen %>% filter(grepl('3789', plantid)) # and there it is
phen = phen %>% mutate(finalid = gsub('3789\\_6', '3164\\_6', finalid))
phen %>% filter(grepl('3789', finalid)) # fixed
phen %>% filter(grepl('3164\\_6', finalid)) # fixed

# 9 3183_9_130.19 TRUE  FALSE FALSE TRUE 
demo %>% filter(grepl('130\\.19', plantid))
phen %>% filter(grepl('3183', plantid)) # hmm...
phen %>% filter(year %in% 2022, grepl('\\_9$', plantid)) # uh...?
phen %>% filter(grepl('\\_9\\_', plantid)) # super strange...
seed %>% filter(grepl('\\_9\\_', plantid.seed)) # also strange...
# so there are no plants from plot 9 in phen or seed aside from this one?
# weird

# 10 3296_7        TRUE  TRUE  FALSE TRUE 
# 11 3298_7        FALSE FALSE TRUE  FALSE
demo %>% filter(grepl('3296', plantid))
phen %>% filter(grepl('3298', plantid)) # there it is
# change it in all years?
phen = phen %>% mutate(finalid = gsub('3298\\_7', '3296_7', finalid))
phen %>% filter(grepl('3296\\_7', finalid))

# 12 3321_4        TRUE  FALSE FALSE TRUE 
demo %>% filter(grepl('3321', plantid))
phen %>% filter(grepl('3321', plantid))
phen %>% filter(year %in% 2022, grepl('\\_4$', plantid)) # hmm...
phen %>% 
  filter(year %in% 2022, grepl('\\_4$', plantid)) %>%
  group_by(plantid) %>%
  filter(n() == 2) # hmm...
all.plantids %>%
  filter(Year %in% 2022, !demo) %>%
  filter(grepl('3182\\_4', plantid) | grepl('3730\\_4', plantid) | grepl('7545\\_4', plantid))
# demo sheet says no phen for this plant

# 13 3363_14       FALSE FALSE TRUE  FALSE
# 14 3367_14       TRUE  TRUE  FALSE TRUE 
demo %>% filter(grepl('3363', plantid))
demo %>% filter(grepl('3367', plantid))
phen %>% filter(grepl('3363', plantid))
seed %>% filter(grepl('3363', plantid.seed))
phen = phen %>% mutate(finalid = gsub('3363', '3367', finalid))
phen %>% filter(grepl('3363', finalid))
phen %>% filter(grepl('3367', finalid)) # fixed

# 15 3379_14       FALSE TRUE  TRUE  FALSE
phen %>% filter(grepl('3379', plantid))
demo %>% filter(grepl('3379', plantid)) # interesting - here in 2023
seed %>% filter(grepl('3379', plantid.seed))
# ah... listed as having no umbels in 2022...
####### UNSOLVED #######
# phen says stalk was eaten; should do zero seeds...

# 16 3389_13       FALSE FALSE FALSE TRUE 
mumb %>% filter(grepl('3389', plantid))
phen %>% filter(grepl('3481', plantid))
mumb = mumb %>% mutate(finalid = gsub('3389', '3481', finalid))

# 17 3456_6        TRUE  FALSE FALSE FALSE
demo %>% filter(grepl('3456', plantid))
phen %>% filter(grepl('7587', plantid))
all.plantids %>% filter(grepl('7587', plantid))
# change in demo
demo = demo %>% mutate(finalid = gsub('3456', '7587', finalid))
demo %>% filter(grepl('3456', finalid))
demo %>% filter(grepl('7587', finalid))

# 18 3469_14       TRUE  TRUE  TRUE  FALSE
mumb %>% filter(grepl('3469', plantid)) # it's here in 2023
demo %>% filter(grepl('3469', plantid))
# hmm... okay might not need changing?

# 19 3481_13       TRUE  TRUE  TRUE  FALSE # fixed above

# 20 3495_14       TRUE  TRUE  FALSE TRUE 
demo %>% filter(grepl('3495', plantid))
# okay... no notes suggesting other tags...
phen %>% filter(grepl('3495', plantid)) # it's here in 2023
demo %>% filter(Year %in% 2022, grepl('\\_14\\_1A', plantid))
# wow that kinda sucks. nearly identical tags in the same location...
demo %>% filter(grepl('\\_14\\_1A', plantid))
# yeah that sucks
# anyway
# also there's no umbel diameter listed for the 2022 record (damn)
phen %>% filter(grepl('3459', plantid)) # only two umbels are listed here (one clipped)
# okay. I don't think this can be fixed.
phen %>% filter(year %in% 2022, grepl('\\_14', plantid))

# 2X (new in phen - from buds?) 3516_13       FALSE FALSE TRUE  FALSE
phen %>% filter(grepl('3516', plantid))
demo %>% filter(grepl('3516', plantid)) # not in demo
demo %>% filter(grepl('\\_13\\_', plantid))
# demo %>% filter(grepl('3[0-9]16\\_13\\_', plantid), Year %in% 2022) # ooh... bet it's this
# phen %>% filter(grepl('3616', plantid), year %in% 2022) # nope, already has a record...
####### UNSOLVED #######
# there's only one row for this plant in the raw csv
# tag mis-entry?

# 21 3542_1        TRUE  TRUE  TRUE  FALSE
demo %>% filter(grepl('3542', plantid), Year %in% 2022)
demo %>% filter(grepl('3542', plantid)) # no notes about tags...
mumb %>% filter(grepl('\\_1\\_', plantid), Year %in% 2022) %>% distinct(plantid)
mumb %>% filter(grepl('\\_1\\_0I', plantid), Year %in% 2022) # nothing here
mumb %>% filter(grepl('\\_1\\_[01238]I', plantid), Year %in% 2022)
all.plantids %>% filter(Year %in% 2022, grepl('\\_1\\_', plantid), !demo, mumb)
####### UNSOLVED #######
# think it's just... not recorded. or data entry mistake?

# 22 3555_4        TRUE  FALSE FALSE TRUE 
demo %>% filter(grepl('3555', plantid)) # no other tags
phen %>% filter(grepl('3555', plantid)) # here in 2021
seed %>% filter(grepl('3555', plantid.seed)) # here in 2021...
all.plantids %>% filter(grepl('\\_4$', plantid), Year %in% 2022, phen)
# maybe 3656...? no... that would have appeared below...
####### UNSOLVED #######

# 23 3558_13       TRUE  FALSE FALSE TRUE 
demo %>% filter(grepl('3558', plantid)) # no other tags
phen %>% filter(grepl('3558', plantid))
phen %>% filter(grepl('3769', plantid), year %in% 2022) # yep that's it
mumb %>% filter(grepl('3769', plantid)) # sweet
# only one umbel in 2021/2022 so no need to change mumb
demo = demo %>% mutate(finalid = gsub('3558', '3769', finalid))

# 24 3563_4        TRUE  FALSE FALSE TRUE 
demo %>% filter(grepl('3563', plantid)) # no other tags
phen %>% filter(grepl('3563', plantid)) # not here at all...
phen %>% filter(grepl('\\_4$', plantid)) %>% distinct(plantid, .keep_all = TRUE)
seed %>% filter(grepl('3563', plantid.demo)) # nothing here...
seed %>% filter(grepl('\\_4\\_', plantid.demo)) %>% distinct(plantid.demo, .keep_all = TRUE)
# could it be 3656?
demo %>% filter(grepl('3656', plantid), Year %in% 2022) # nope, that's here
####### UNSOLVED #######

# 25 3595_13       FALSE TRUE  FALSE FALSE
seed %>% filter(grepl('3595', plantid.seed)) # ? should be here...
demo %>% filter(grepl('3595', plantid))
# ah... it's not listed as flowering in demo
# although it's not in phen...
phen %>% filter(grepl('3595', plantid)) # here in 2021
phen %>% filter(grepl('\\_13$', plantid)) %>% distinct(plantid, .keep_all = TRUE)
demo %>% filter(grepl('3593', plantid), Year %in% 2022) # nope, not mis-entered as 3593
####### UNSOLVED #######
####### go into demo and add flowering info to 3595_13 in 2022
####### (although it has no phen... maybe I should use budding phen?)

# 26 3614_13       FALSE FALSE TRUE  FALSE # should have been fixed above...
phen %>% filter(grepl('3614', plantid))

# 27 3634_4        TRUE  FALSE FALSE TRUE 
demo %>% filter(grepl('3634', plantid)) # probably not a tag fix... (14 leaves??)
mumb %>% filter(grepl('3634', plantid)) # okay, good
phen %>% filter(grepl('3634', plantid)) # here in 2021
seed %>% filter(grepl('3634', plantid.demo)) # nothing
seed %>% filter(grepl('\\_4\\_', plantid.demo), year %in% 2022)
# ugh... not here?
phen %>% filter(grepl('\\_4$', plantid), year %in% 2022)
####### UNSOLVED #######

# 28 3651_3        TRUE  TRUE  FALSE TRUE 
demo %>% filter(grepl('3651', plantid)) # no other tag
phen %>% filter(grepl('3651', plantid)) # 2023 only
seed %>% filter(grepl('3651', plantid.demo)) # come on man...
all.plantids %>% filter(Year %in% 2022, phen, !seed)
# THAT'S IT! right?
phen %>% filter(grepl('3151', plantid))
demo %>% filter(grepl('3151\\_3', plantid)) # looks promising...
seed %>% filter(grepl('3151\\_3', plantid.seed))
phen = phen %>% mutate(finalid = gsub('3151\\_3', '3651_3', finalid))

# 29 3706_15       TRUE  TRUE  FALSE TRUE  # possibly fixed above
phen %>% filter(grepl('3706', finalid))

# 30 3726_13       FALSE TRUE  TRUE  FALSE
demo %>% filter(grepl('3726', plantid)) # nada
seed %>% filter(grepl('3726', plantid.seed))
demo %>% filter(grepl('\\_13\\_', plantid), Year %in% 2022, No.umbels %in% 1)
demo %>% filter(grepl('3726', proc.note)) # tag was changed
demo %>% filter(grepl('3039', plantid))
demo = demo %>% mutate(finalid = gsub('3039', '3726', finalid))

# 31 3735_4        TRUE  FALSE FALSE TRUE 
demo %>% filter(grepl('3735', plantid))
demo %>% filter(grepl('3735', proc.note))
demo %>% filter(grepl('3735', demo.note))
phen %>% filter(grepl('\\_4$', plantid))
####### UNSOLVED #######

# 32 3750_4        TRUE  FALSE FALSE TRUE 
demo %>% filter(grepl('3750', plantid))
demo %>% filter(grepl('3750', proc.note))
demo %>% filter(grepl('3750', demo.note))
phen %>% filter(grepl('\\_4$', plantid))
# 3730? 7550?
phen %>% filter(grepl('3730', plantid), year %in% 2022) # nope
phen %>% filter(grepl('7550', plantid), year %in% 2022) # nope
####### UNSOLVED #######

# 33 3769_13       FALSE TRUE  TRUE  FALSE
demo %>% filter(grepl('3769', plantid))
demo %>% filter(grepl('3769', demo.note))
demo %>% filter(grepl('3769', proc.note)) # fixed above

# 34 3770_9        TRUE  FALSE FALSE TRUE 
demo %>% filter(grepl('3770', plantid))
phen %>% filter(grepl('3059', plantid))
phen %>% filter(grepl('3770', plantid))
seed %>% filter(grepl('3770', plantid.demo))
seed %>% filter(grepl('3770', plantid.seed))
####### UNSOLVED #######

# 35 3782_1        FALSE TRUE  TRUE  FALSE
demo %>% filter(grepl('3782', plantid)) # it's here in 2023
# ahhh this was that weird pair
demo %>% filter(grepl('3814', plantid))
seed %>% filter(grepl('3782', plantid.seed))
# but which of these is which
# probably requires looking at the original datasheets
####### UNSOLVED #######
# (but solvable!)

# 36 3789_6        FALSE FALSE TRUE  FALSE
phen %>% filter(grepl('3789', plantid)) # already fixed

# 37 3813_15       FALSE FALSE TRUE  FALSE
phen %>% filter(grepl('3813', plantid))
demo %>% filter(grepl('3813', proc.note))
demo %>% filter(grepl('3813', demo.note))
# plant is 3831
demo %>% filter(grepl('3831', plantid))
seed %>% filter(grepl('3831', plantid.seed))
phen %>% filter(grepl('3813\\_15', plantid))
# specifically 3831_15
phen = phen %>% mutate(finalid = gsub('3813\\_15', '3831_15', finalid))
phen %>% filter(grepl('3813\\_15', finalid)) # good
phen %>% filter(grepl('3831\\_15', finalid)) # fixed

# 38 3814_1        TRUE  TRUE  TRUE  FALSE
demo %>% filter(grepl('3814', plantid))
# this is the one from above
# look at the original datasheets
####### UNSOLVED #######
# (but solvable!)

# 39 3831_15       TRUE  TRUE  FALSE TRUE 
phen %>% filter(grepl('3831', plantid)) # fixed above

# 40 3866_7        TRUE  TRUE  FALSE TRUE 
phen %>% filter(grepl('3866', plantid)) # 2021 record only
demo %>% filter(grepl('3866', proc.note))
demo %>% filter(grepl('3866', plantid)) # 3866 (ah.. is below)
phen %>% filter(grepl('3886', plantid))
phen = phen %>% mutate(finalid = gsub('3886', '3866', finalid))
phen %>% filter(grepl('3886', finalid))

# 41 3880_13       FALSE FALSE FALSE TRUE 
mumb %>% filter(grepl('3880', plantid))
# some confusion with 7537
mumb %>% filter(grepl('7537', plantid)) # none here
demo %>% filter(grepl('7537', plantid))
phen %>% filter(grepl('7537', plantid))
mumb = mumb %>% mutate(finalid = gsub('3880', '7537', finalid))

# 42 3886_7        FALSE FALSE TRUE  FALSE # solved above

# 43 5037_7        TRUE  TRUE  FALSE TRUE 
demo %>% filter(grepl('5037', plantid))
# note says no tracking phenology (death)

# 44 5040_7        TRUE  TRUE  FALSE TRUE
demo %>% filter(grepl('5040', plantid))
# note says umbel is dead
seed %>% filter(grepl('5040', plantid.seed))

# 45 7526_1        FALSE TRUE  FALSE FALSE
# 46 7536_1        TRUE  FALSE TRUE  TRUE 
seed %>% filter(grepl('7526', plantid.seed)) # (lol? why did I do this)
phen %>% filter(grepl('7536', plantid))
demo %>% filter(grepl('7536', plantid))
seed = seed %>% mutate(finalid = gsub('7526', '7536', finalid))
seed %>% filter(grepl('7526', finalid))
seed %>% filter(grepl('7536', finalid))

# 47 7537_13       TRUE  TRUE  TRUE  FALSE # solved above

# 48 7587_6        FALSE TRUE  TRUE  TRUE 
demo %>% filter(grepl('7587', proc.note))
demo %>% filter(grepl('7587', plantid))
demo %>% filter(grepl('7587', finalid)) # fixed above


##########
### Re-combine these based on finalids
##########

all.plantids = rbind(
  # Demo data
  demo %>% 
    distinct(finalid, Year) %>% 
    mutate(source = 'demo', is.in = TRUE),
  # Seed set data
  seed %>% 
    distinct(finalid, year) %>% 
    mutate(source = 'seed', is.in = TRUE) %>% 
    rename(Year = year),
  # Phenology data
  phen %>%
    distinct(finalid, year) %>%
    mutate(source = 'phen', is.in = TRUE) %>%
    rename(Year = year),
  # Umbel data
  # NOTE: this dataset only contains *multiple* umbels
  # so I will also add in the demo data for which there is one umbel and umbel
  # measurements
  mumb %>%
    distinct(finalid, Year) %>%
    mutate(source = 'mumb', is.in = TRUE),
  demo %>%
    filter(No.umbels %in% 1) %>%
    distinct(finalid, Year) %>%
    mutate(source = 'mumb', is.in = TRUE)
) %>%
  pivot_wider(
    id_cols = c(finalid, Year), names_from = source, 
    values_from = is.in, values_fill = FALSE
  )

head(all.plantids)
nrow(all.plantids)

all.plantids %>% 
  filter(!(Year %in% 2022)) %>%
  mutate(n.in = demo + seed + phen + mumb) %>%
  group_by(n.in) %>%
  summarise(n = n())
# still quite bad, goodness

all.plantids %>% 
  mutate(n.in = demo + seed + phen + mumb) %>%
  group_by(Year, n.in) %>%
  summarise(n = n())

### 2021: only a handful of plants here were included in seed
# so, there shouldn't be any plants from 2021 that are not listed in the
# seed set data
# I'll remove these plants here

all.plantids = all.plantids %>% filter(Year > 2021 | seed)

####################
# Next - look at 2021
####################

# all.plantids %>%
#   filter(Year %in% 2021) %>%
#   filter(!(demo & mumb & phen & seed)) %>%
#   print(n = nrow(.))

# 1 3135_12_7G   2021 TRUE  TRUE  TRUE  FALSE
demo %>% filter(grepl('3135', plantid))
mumb %>% filter(grepl('3135', plantid))
# notes suggest that one of these two umbels was closed
# so, remove from data?

# 2 3143_9_14T   2021 TRUE  TRUE  TRUE  FALSE
demo %>% filter(grepl('3143', plantid))
mumb %>% filter(grepl('3143', plantid))
# one umbel area is listed - other is closed
# remove?

# 3 3149_14_19A  2021 TRUE  TRUE  FALSE TRUE 
demo %>% filter(grepl('3149', plantid))
phen %>% filter(grepl('\\_14\\_', plantid), year %in% 2021)
# ah... somehow coords got recoreded as '29' instead of '19'
phen = phen %>% mutate(finalid = gsub('3149\\_14\\_29A', '3149_14_19A', finalid))
phen %>% filter(grepl('\\_14\\_', plantid), year %in% 2021)

# 4 3179_6_5E    2021 TRUE  TRUE  FALSE TRUE 
phen %>% filter(grepl('\\_6\\_', finalid), year %in% 2021) %>% distinct(finalid) # that's a lot
demo %>% filter(grepl('3179', plantid))
# in 2022, tag was 3748
phen %>% filter(grepl('3748', plantid)) # and there it is
phen = phen %>% mutate(finalid = gsub('3748', '3179', finalid))
phen %>% filter(grepl('3179', finalid)) # good

# 5 3436_15_17H  2021 TRUE  TRUE  FALSE TRUE 
demo %>% filter(grepl('3436', plantid))
phen %>% filter(grepl('3436', plantid))
phen %>% filter(grepl('3436\\_15\\_17H', plantid))
# mis-entered coordinates
phen = phen %>% mutate(finalid = gsub('3436\\_15\\_12H', '3436_15_17H', finalid))
phen %>% filter(grepl('3436\\_15\\_17H', finalid))

# 6 3453_7_0C    2021 TRUE  TRUE  TRUE  FALSE
demo %>% filter(grepl('3453', plantid)) # no umbels... maybe died? also doesn't say SOS
phen %>% filter(grepl('3453', plantid))
seed %>% filter(grepl('3453', plantid.seed)) # only one is dead
mumb %>% filter(grepl('3453', plantid)) # here in 2022
mumb %>% filter(grepl('\\_7\\_', plantid), Year %in% 2021) # nobody else in plot 7, also not SOS (ugh)
# maybe it was just never entered?
####### UNSOLVED #######

# 7 3567_1_3I    2021 TRUE  TRUE  TRUE  FALSE
demo %>% filter(grepl('3567', plantid)) # only one umbel listed
mumb %>% filter(grepl('3567', plantid)) # yep, here in 2022-23
seed %>% filter(grepl('3567', plantid.seed)) # dead
# guess it'll just get excluded

# 8 3916_6_7C    2021 TRUE  TRUE  FALSE TRUE 
demo %>% filter(grepl('3916', plantid)) # no other tags
phen %>% filter(grepl('\\_6\\_', plantid), year %in% 2021) %>% distinct(plantid)
# bet it's 3918
phen %>% filter(grepl('3918', plantid))
demo %>% filter(grepl('3918', plantid))
phen = phen %>% mutate(finalid = gsub('3918', '3916', finalid))
phen %>% filter(grepl('3916', finalid))

# 9 3978_13_2F   2021 TRUE  TRUE  TRUE  FALSE
demo %>% filter(grepl('3978', plantid)) # only one umbel diameter listed
seed %>% filter(grepl('3978', plantid.seed)) # yes - it's dead

# 10 3365_1_0J    2021 FALSE TRUE  TRUE  FALSE
seed %>% filter(grepl('3365', plantid.seed))
demo %>% filter(grepl('3365\\_1', plantid)) # listed as 
# oh... it's listed as having zero umbels in 2021
phen %>% filter(grepl('3365', plantid))
####### UNSOLVED #######
# check for a data entry error?

# 11 3418_1_18E   2021 FALSE TRUE  TRUE  FALSE
seed %>% filter(grepl('3418', plantid.seed), year %in% 2021)
demo %>% filter(grepl('3418', plantid)) # coords mis-entered
demo = demo %>% mutate(finalid = gsub('3418\\_1\\_18B', '3418_1_18E', finalid))
demo %>% filter(grepl('3418\\_1\\_18E', finalid))

# 12 3122_4_8J    2021 FALSE TRUE  TRUE  FALSE
seed %>% filter(grepl('3122', plantid.seed))
# it's listed as not flowering in 2021 demo (zero umbels)
####### UNSOLVED #######
# check for a data entry error?

# 13 3555_4_11T   2021 FALSE TRUE  TRUE  FALSE
demo %>% filter(grepl('3555', plantid)) # I changed the coord in processing
demo = demo %>% mutate(finalid = gsub('3555\\_4\\_11J', '3555_4_11T', finalid))
demo %>% filter(grepl('3555\\_4\\_11T', finalid))

# 14 3488_5_13J   2021 FALSE TRUE  FALSE FALSE
demo %>% filter(grepl('3488', plantid)) # just a coord issue
mumb %>% filter(grepl('3488', plantid))
phen %>% filter(grepl('3488', plantid)) # hmm... not sure what 2023 will look like
seed %>% filter(grepl('3488', plantid.seed)) # missed in 2023?
demo = demo %>% mutate(finalid = gsub('3488\\_5\\_15I', '3488_5_13J', finalid))
mumb = mumb %>% mutate(finalid = gsub('3488\\_5\\_15I', '3488_5_13J', finalid))

# 15 5045_5_13C   2021 FALSE TRUE  TRUE  FALSE
demo %>% filter(grepl('5045', plantid))
mumb %>% filter(grepl('5045', plantid)) # hmm...
# actually maybe here change seed/phen instead of demo/mumb
seed %>% filter(grepl('5045', plantid.seed))
seed = seed %>% mutate(finalid = gsub('5045\\_5_13C', '5045_5_12B', finalid))
phen %>% filter(grepl('5045', plantid))
phen = phen %>% mutate(finalid = gsub('5045\\_5_13C', '5045_5_12B', finalid))

##########
### Re-combine these based on finalids, again
##########

all.plantids = rbind(
  # Demo data
  demo %>% 
    distinct(finalid, Year) %>% 
    mutate(source = 'demo', is.in = TRUE),
  # Seed set data
  seed %>% 
    distinct(finalid, year) %>% 
    mutate(source = 'seed', is.in = TRUE) %>% 
    rename(Year = year),
  # Phenology data
  phen %>%
    distinct(finalid, year) %>%
    mutate(source = 'phen', is.in = TRUE) %>%
    rename(Year = year),
  # Umbel data
  # NOTE: this dataset only contains *multiple* umbels
  # so I will also add in the demo data for which there is one umbel and umbel
  # measurements
  mumb %>%
    distinct(finalid, Year) %>%
    mutate(source = 'mumb', is.in = TRUE),
  demo %>%
    filter(No.umbels %in% 1) %>%
    distinct(finalid, Year) %>%
    mutate(source = 'mumb', is.in = TRUE)
) %>%
  pivot_wider(
    id_cols = c(finalid, Year), names_from = source, 
    values_from = is.in, values_fill = FALSE
  ) %>%
  filter(Year > 2021 | seed)

head(all.plantids)
nrow(all.plantids)

all.plantids %>% 
  filter(!(Year %in% 2022)) %>%
  mutate(n.in = demo + seed + phen + mumb) %>%
  group_by(n.in) %>%
  summarise(n = n())
# better but not by much

all.plantids %>% 
  mutate(n.in = demo + seed + phen + mumb) %>%
  group_by(Year, n.in) %>%
  summarise(n = n())
# A lot of 2023 to do. Probably coordinate issues (bleh)


####################
# Fix 2023
####################

all.plantids %>%
  filter(Year %in% 2023) %>%
  filter(!(demo & mumb & phen & seed)) %>%
  arrange(finalid) %>%
  print(n = nrow(.))
# that's a lot
# prune out the coords, there's gotta be a clever way to do this... maybe see if
# any tag-plots show up multiple times?

all.plantids %>%
  filter(Year %in% 2023) %>%
  separate(finalid, into = c('tag', 'plot', 'coord'), sep = '_', remove = FALSE) %>%
  select(-coord) %>%
  unite(c(tag, plot), col = planttag, sep = "_") %>%
  group_by(planttag) %>%
  filter(n() > 2)
# only two cases to worry about

all.plantids %>%
  filter(Year %in% 2023) %>%
  separate(finalid, into = c('tag', 'plot', 'coord'), sep = '_', remove = FALSE) %>%
  select(-coord) %>%
  unite(c(tag, plot), col = planttag, sep = "_") %>%
  arrange(planttag) %>% 
  filter(!(demo & seed & phen & mumb))
# jesus that's a lot

all.plantids %>%
  filter(Year %in% 2023) %>%
  separate(finalid, into = c('tag', 'plot', 'coord'), sep = '_', remove = FALSE) %>%
  select(-coord) %>%
  unite(c(tag, plot), col = planttag, sep = "_") %>%
  arrange(planttag) %>% 
  filter(!(demo & seed & phen & mumb)) %>%
  # Get rid of plants where there are two records and the plant-tag combo
  # appears in every data set
  group_by(planttag) %>%
  filter(!(n() %in% 2 & any(demo) & any(seed) & any(phen) & any(mumb))) %>%
  print(n = nrow(.))
# jesus what is up with the seed and phen data in 2023

seed %>% filter(year %in% 2023) %>% distinct(plantid.seed) %>% nrow()
seed %>% 
  filter(year %in% 2023) %>% 
  distinct(plantid.seed, .keep_all = TRUE) %>% 
  group_by(plot) %>%
  summarise(n = n())
phen %>% filter(year %in% 2023) %>% distinct(plantid) %>% nrow()
demo %>% filter(Year %in% 2023, !is.na(No.umbels), No.umbels > 0) %>% nrow()
# ah... ~150 plants in demo but not in phen
# interesting
# So I'll just focus on the plants in phen and seed?

all.plantids %>%
  filter(Year %in% 2023) %>%
  separate(finalid, into = c('tag', 'plot', 'coord'), sep = '_', remove = FALSE) %>%
  select(-coord) %>%
  unite(c(tag, plot), col = planttag, sep = "_") %>%
  arrange(planttag) %>% 
  filter(!(demo & seed & phen & mumb)) %>%
  # Get rid of plants where there are two records and the plant-tag combo
  # appears in every data set
  group_by(planttag) %>%
  filter(!(n() %in% 2 & any(demo) & any(seed) & any(phen) & any(mumb))) %>%
  ungroup() %>%
  arrange(desc(seed), desc(phen), planttag) %>%
  print(n = nrow(.))

# 1 3125_5_2I     3125_5    2023 FALSE TRUE  TRUE  FALSE
demo %>% filter(grepl('3125\\_5', plantid))
phen %>% filter(grepl('3125\\_5', plantid))
seed %>% filter(grepl('3125\\_5', plantid.seed))
phen = phen %>% mutate(finalid = gsub('3125\\_5\\_2I', '3125_5_3I', finalid))
seed = seed %>% mutate(finalid = gsub('3125\\_5\\_2I', '3125_5_3I', finalid))

# 2 3368_15_14B   3368_15   2023 FALSE TRUE  TRUE  FALSE
seed %>% filter(grepl('3368', plantid.seed))
# I checked demo - plant has an all-NA record because it disappeared
# no fix for this

# 3 3391_13_8G    3391_13   2023 FALSE TRUE  TRUE  FALSE
seed %>% filter(grepl('3391\\_13', plantid.seed))
demo %>% filter(grepl('3391\\_13', plantid))
phen = phen %>% mutate(finalid = gsub('3391\\_13\\_8G', '3391_13_7G', finalid))
seed = seed %>% mutate(finalid = gsub('3391\\_13\\_8G', '3391_13_7G', finalid))

# 4 3785_6_15H    3785_6    2023 FALSE TRUE  TRUE  FALSE
seed %>% filter(grepl('3785', plantid.seed))
demo %>% filter(grepl('3785', proc.note)) # lol? it's listed as no plant in 2023...
phen %>% filter(grepl('3785', plantid)) # ?

# 5 5849_7_1E     5849_7    2023 FALSE TRUE  TRUE  FALSE
seed %>% filter(grepl('5849', plantid.seed)) # doesn't match with anything in demo...
demo %>% filter(grepl('5049', plantid))
mumb %>% filter(grepl('5049', plantid))
seed = seed %>% mutate(finalid = gsub('5849', '5049', finalid))
phen = phen %>% mutate(finalid = gsub('5849', '5049', finalid))

# 6 7698_10_9D    7698_10   2023 TRUE  TRUE  TRUE  FALSE
demo %>% filter(grepl('7698', plantid))
mumb %>% filter(grepl('7[0-9]98', plantid))
mumb = mumb %>% mutate(finalid = gsub('7598', '7698', finalid))

# 7 3069_7_1B     3069_7    2023 FALSE TRUE  FALSE FALSE
seed %>% filter(grepl('3069', plantid.seed))
demo %>% filter(grepl('3069', plantid))
phen %>% filter(grepl('3069', plantid)) # it's here now
phen %>% filter(grepl('\\_7\\_[12][ABC]', plantid))
# oh right... umbel is dead.. well even still
seed = seed %>% mutate(finalid = gsub('3069\\_7\\_1B', '3069_7_1C', finalid))

# 8 3297_15_11H   3297_15   2023 FALSE TRUE  FALSE FALSE
seed %>% filter(grepl('3297', plantid.seed))
# bud date (but only bud date) appears in phen
# demo %>% filter(grepl('3297', plantid))

# 9 3454_6_15I    3454_6    2023 FALSE TRUE  FALSE FALSE
seed %>% filter(grepl('3454', plantid.seed)) # ah... missing from phen due to missing records?
demo %>% filter(grepl('3454', plantid)) # hmm...
# oh well, ignore

# 10 3518_3_16J    3518_3    2023 FALSE TRUE  FALSE FALSE
seed %>% filter(grepl('3518', plantid.seed), year %in% 2023)
# wait... data entry error?
# why are there two different plots listed here...
demo %>% filter(grepl('3518', plantid))
# three umbels, just like in seed
# these are all listed as being in plot 6 and not 3
# same with the multi-umbel dataset
# must be a data entry error
seed = seed %>% mutate(finalid = gsub('3518\\_[36]\\_16J', '3518_6_10I', finalid))
phen = phen %>% mutate(finalid = gsub('3518\\_6\\_16J', '3518_6_10I', finalid))

# 11 3590_13_9B    3590_13   2023 FALSE TRUE  FALSE FALSE
seed %>% filter(grepl('3590', plantid.seed))
# listed as 'gone 4 jun'; it's flowering in demo but won't be in phen
phen %>% filter(grepl('3590', plantid)) # it's in phen now
# not sure why it doesn't match with demo though...
# fixed!

# 12 3636_7_3I     3636_7    2023 FALSE TRUE  FALSE FALSE
seed %>% filter(grepl('3636', plantid.seed))
# also listed as "gone" on 24 May
phen %>% filter(grepl('3636', plantid)) #in phen now
demo %>% filter(grepl('3636', plantid))
# hmm... change 3I to 3J in demo and phen
seed = seed %>% mutate(finalid = gsub('3636\\_7\\_3I', '3636_7_3J', finalid))
phen = phen %>% mutate(finalid = gsub('3636\\_7\\_3I', '3636_7_3J', finalid))

# 13 5570_6_18H    5570_6    2023 FALSE TRUE  FALSE FALSE
seed %>% filter(grepl('5570', plantid.seed))
# there are no plants with tags starting with 55, so this is def an entry error
demo %>% filter(grepl('\\_6\\_18[GH]', plantid), Year %in% 2023) # hmm... not liking this
# not listed in processing or demo notes...
####### UNSOLVED #######
# check datasheets

# Bunch of plants - new from budphen
# 12 3010_5_16E    3010_5    2023 FALSE FALSE TRUE  FALSE
phen %>% filter(grepl('3010', plantid), year %in% 2023)
demo %>% filter(grepl('3010', plantid), Year %in% 2023)
seed %>% filter(grepl('3010', plantid.seed)) # no seed in 2023... same for plantid.demo
# I don't think this one is dead though...
# only one record in phen... wonder what happened here. misentry?
phen = phen %>% mutate(finalid = gsub('3010\\_5\\_16E', '3010_5_14G', finalid))
####### UNSOLVED #######
# maybe... at least, go back and see what's up in the raw data

# 13 3065_15_2F    3065_15   2023 TRUE  FALSE TRUE  TRUE 
# ooh... in demo and phen but not in seed...
demo %>% filter(grepl('3065', plantid))
# umbel is listed as dead. so not in seed

# 14 3067_5_17D    3067_5    2023 FALSE FALSE TRUE  FALSE
phen %>% filter(grepl('3067', plantid), year %in% 2023)
demo %>% filter(grepl('3067', plantid))
mumb %>% filter(grepl('3067', plantid))
phen = phen %>% mutate(finalid = gsub('3067\\_5\\_17D', '3067_5_16F', finalid))
seed %>% filter(grepl('3067', plantid.seed)) # still don't know why it's not in seed...

# 15 3068_15_5H    3068_15   2023 FALSE FALSE TRUE  FALSE
# it's in demo
phen = phen %>% mutate(finalid = gsub('3068\\_15\\_5H', '3068_15_6I', finalid))
phen %>% filter(grepl('3068', plantid), year %in% 2023)
seed %>% filter(grepl('3068', plantid.seed)) # here in 2022, not 2023

# 16 3127_5_17B    3127_5    2023 FALSE FALSE TRUE  FALSE
phen = phen %>% mutate(finalid = gsub('3127\\_5\\_17B', '3127_5_16B', finalid))
phen %>% filter(grepl('3127', finalid))
seed %>% filter(grepl('3127', plantid.seed)) # here in 2022, not in 2023...

# 17 3135_5_19C    3135_5    2023 FALSE FALSE TRUE  FALSE
phen = phen %>% mutate(finalid = gsub('3135\\_5\\_19C', '3135_5_18D', finalid))
phen %>% filter(grepl('3135', plantid), year %in% 2023)
seed %>% filter(grepl('3135', plantid.seed), year %in% 2023) # this is a different one
# maybe this is one of those where they just didn't do seed for certain plants

# 18 3144_5_15J    3144_5    2023 FALSE FALSE TRUE  FALSE
phen %>% filter(grepl('3144', plantid))
demo %>% filter(grepl('3144', plantid))
# uh... it's listed as dead in raw demo...
# there's only one record here...
####### UNSOLVED #######

# 19 3153_5_15D    3153_5    2023 FALSE FALSE TRUE  FALSE
phen = phen %>% mutate(finalid = gsub('3153\\_5\\_15D', '3153_5_14D', finalid))
phen %>% filter(grepl('3153', plantid))
demo %>% filter(grepl('3153', plantid), Year %in% 2023)

# 14 7519_5_3D     7519_5    2023 TRUE  TRUE  FALSE TRUE 
seed %>% filter(grepl('7519', plantid.seed))
# okay, I changed the tag at one point...
phen %>% filter(grepl('7512', plantid))
demo %>% filter(grepl('7519', plantid))
demo %>% filter(grepl('7512', plantid))
# ugh... okay going to assume 7512 was misentered as 7519
phen = phen %>% mutate(finalid = gsub('7512\\_5\\_3D', '7519_5_3D', finalid))

# 15 3182_4_17H    3182_4    2023 FALSE FALSE TRUE  FALSE
phen %>% filter(grepl('3182', plantid))
demo %>% filter(grepl('3182', plantid)) # not here in 2023...
# I don't see evidence this anywhere
# oh... listed as NA/NP in demo, disappeared...
####### UNSOLVED #######
# change in demo

# (new from budphen) 21 3187_5_18H    3187_5    2023 FALSE FALSE TRUE  FALSE
phen %>% filter(grepl('3187', plantid), year %in% 2023)
demo %>% filter(grepl('3187\\_5', plantid)) # tag misentry seems likely
# demo %>% filter(grepl('3[0-9]87\\_5', plantid), Year %in% 2023)
# phen %>% filter(grepl('3987', plantid), year %in% 2023) # nope, not this
# demo %>% filter(grepl('31[0-9]7\\_5', plantid), Year %in% 2023)
# phen %>% filter(grepl('31[29]7\\_5', plantid), year %in% 2023) # seems like not this either
# demo %>% filter(grepl('318[0-9]\\_5', plantid), Year %in% 2023) # nothing else here...
# only one record in csv... 
####### UNSOLVED #######

# 16 3192_1_11B    3192_1    2023 TRUE  FALSE TRUE  TRUE 
demo %>% filter(grepl('3192', plantid))
phen %>% filter(grepl('3192', plantid))
# it's listed in seed in 2022
# very weirdly I don't see any plants in 2023 in seed with an x-coord greater than 10
# maybe they were left out from seed? or not entered?
# (seeds were not counted here)

# 17 3193_1_15A    3193_1    2023 FALSE FALSE TRUE  FALSE
phen %>% filter(grepl('3193', plantid)) # maybe a misentry
demo %>% filter(grepl('3193\\_1', plantid))
# yep
seed %>% filter(grepl('3193\\_1', plantid.seed)) # annoying
phen = phen %>% mutate(finalid = gsub('3193\\_1\\_15A', '3193_1_16C', finalid))
# still might not be in seed though (very strange)

# 24 (new from budphen) 3367_14_17I   3367_14   2023 FALSE FALSE TRUE  FALSE
phen = phen %>% mutate(finalid = gsub('3367\\_14\\_17I', '3367_14_15H', finalid))
phen %>% filter(grepl('3367', finalid), year %in% 2023)
seed %>% filter(grepl('336[37]', plantid.seed), year %in% 2023) # doesn't appear in seed...

# 25 (new from budphen) 3417_15_4H    3417_15   2023 FALSE FALSE TRUE  FALSE
phen = phen %>% mutate(finalid = gsub('3417\\_15\\_4H', '3417_15_6I', finalid))
# phen %>% filter(grepl('3417', plantid), year %in% 2023)
# demo %>% filter(grepl('3417', plantid)) # good
seed %>% filter(grepl('3417', plantid.seed)) # only in 2022, not 2023

# 18 3424_1_12G    3424_1    2023 TRUE  FALSE TRUE  TRUE 
phen %>% filter(grepl('3424\\_1\\_', plantid))
demo %>% filter(grepl('3424\\_1\\_', plantid))
mumb %>% filter(grepl('3424', plantid))
seed %>% filter(grepl('3424', plantid.seed))
# ah... again, a missing plant from plot 1, xcoord > 10
# check data entry... I wonder if these were just missed
# (note to past self - yes, they were)

# 19 3432_1_19A    3432_1    2023 TRUE  FALSE TRUE  TRUE 
# seems likely the same as above.,,
phen %>% filter(grepl('3432\\_1\\_', plantid), year %in% 2023)
demo %>% filter(grepl('3432\\_1\\_', plantid))

# 20 3463_14_6A    3463_14   2023 FALSE FALSE TRUE  FALSE
phen %>% filter(grepl('3463', plantid))
demo %>% filter(grepl('3463', plantid)) # hmm... 
# ah, it's listed for some reason as not flowering in 2023?
demo %>% filter(grepl('\\_14\\_5B', plantid), Year %in% 2023) # nope
####### UNSOLVED #######
# but I will change the ID in phen to match the demo id at least
phen = phen %>% mutate(finalid = gsub('3463\\_14\\_6A', '3463_14_5B', finalid))

# (new from budphen) 29 3473_5_18F    3473_5    2023 FALSE FALSE TRUE  FALSE
phen = phen %>% mutate(finalid = gsub('3473\\_5\\_18F', '3473_5_18G', finalid))
# phen %>% filter(grepl('3473', plantid), year %in% 2023)
# demo %>% filter(grepl('3473', plantid))
# seed %>% filter(grepl('3473', plantid.seed)) 2022, not 2023

# (new from budphen) 30 3476_5_11F    3476_5    2023 FALSE FALSE TRUE  FALSE
phen = phen %>% mutate(finalid = gsub('3476\\_5\\_11F', '3476_5_12E', finalid))
# phen %>% filter(grepl('3476', finalid), year %in% 2023)
# demo %>% filter(grepl('3476', plantid))s

# 21 3477_2_5C     3477_2    2023 FALSE FALSE TRUE  FALSE
phen %>% filter(grepl('3477', plantid))
demo %>% filter(grepl('3477', plantid)) # ah... listed as dead,  zero seeds
seed %>% filter(grepl('3477', plantid.seed))
####### UNSOLVED #######
# (add to seed?)

# (from budphen?) 32 3488_5_13J    3488_5    2023 FALSE FALSE TRUE  FALSE
# phen %>% filter(grepl('3488', plantid), year %in% 2023)
# phen = phen %>% mutate(finalid = gsub('3488\\_5\\_13J', '3488_5_15I', finalid))
# demo %>% filter(grepl('3488', plantid))
# it'll get changed later - ignore!

# 22 3519_2_8G     3519_2    2023 FALSE FALSE TRUE  FALSE
demo %>% filter(grepl('3519', plantid)) # also listed as dead (broken stalk)
seed %>% filter(grepl('3519', plantid.seed)) 
####### UNSOLVED #######
# add zero to seed

# (from budphen) 34 3526_14_17C   3526_14   2023 FALSE FALSE TRUE  FALSE
phen = phen %>% mutate(finalid = gsub('3526\\_14\\_17C', '3526_14_14E', finalid))
# phen %>% filter(grepl('3526', plantid))
# demo %>% filter(grepl('3526', plantid))
# seed %>% filter(grepl('3526', plantid.seed))

# 23 3556_14_17I   3556_14   2023 FALSE FALSE TRUE  FALSE
demo %>% filter(grepl('3556', plantid))
phen %>% filter(grepl('3556', plantid))
# super weird
seed %>% filter(grepl('3556', plantid.seed)) # not listed in seed for some reason...
phen = phen %>% mutate(finalid = gsub('3556\\_14\\_17I', '3556_14_15H', finalid))
seed %>% filter(grepl('\\_14\\_', plantid.seed), year %in% 2023) %>% distinct(plantid.seed)
# again, no records from xcoord 10 or above (save for one at x = 10)

# (budphen) 36 3566_5_14I    3566_5    2023 FALSE FALSE TRUE  FALSE
phen = phen %>% mutate(finalid = gsub('3566\\_5\\_14I', '3566_5_13H', finalid))
# phen %>% filter(grepl('3566', finalid))
# demo %>% filter(grepl('3566', finalid))

# (budphen) 37 3659_15_4B    3659_15   2023 FALSE FALSE TRUE  FALSE
# demo %>% filter(grepl('3659', finalid))
# mumb %>% filter(grepl('3659', finalid))
phen = phen %>% mutate(finalid = gsub('3659\\_15\\_4B', '3659_15_4C', finalid))
# phen %>% filter(grepl('3659', finalid))

# (budphen) 38 3663_5_19J    3663_5    2023 FALSE FALSE TRUE  FALSE
phen = phen %>% mutate(finalid = gsub('3663\\_5\\_19J', '3663_5_19I', finalid)) 
# phen %>% filter(grepl('3663', plantid))
# demo %>% filter(grepl('3663', plantid))

# (budphen) 39 3685_5_19B    3685_5    2023 FALSE FALSE TRUE  FALSE
phen = phen %>% mutate(finalid = gsub('3685\\_5\\_19B', '3685_5_18B', finalid))
# phen %>% filter(grepl('3685', plantid))
# demo %>% filter(grepl('3685', plantid))

# 24 3688_7_6G     3688_7    2023 FALSE FALSE TRUE  FALSE
demo %>% filter(grepl('3688', plantid)) # had flower but then disappeared
####### UNSOLVED #######
# add zero to seed

# 25 3697_5_17H    3697_5    2023 FALSE FALSE TRUE  FALSE
demo %>% filter(grepl('3697', plantid)) # umbel is dead
####### UNSOLVED #######
# add zero to seed

# 26 3763_1_12A    3763_1    2023 TRUE  FALSE TRUE  TRUE 
demo %>% filter(grepl('3763', plantid))
# ah... might be missing just because of high x-coord
seed %>% filter(grepl('3763', plantid.seed))
seed %>% filter(grepl('\\_1\\_1[0-9]', plantid.seed), year %in% 2023)

# 27 3807_5_10I    3807_5    2023 FALSE FALSE TRUE  FALSE
demo %>% filter(grepl('3807', plantid)) # dead umbel
####### UNSOLVED #######
# add zero to seed

# (budphen) 44 3833_5_17B    3833_5    2023 FALSE FALSE TRUE  FALSE
phen = phen %>% mutate(finalid = gsub('3833\\_5\\_17B', '3833_5_16B', finalid))
# phen %>% filter(grepl('3833', plantid))
# demo %>% filter(grepl('3833', plantid))

# 28 3850_1_17E    3850_1    2023 TRUE  FALSE TRUE  TRUE 
demo %>% filter(grepl('3850', plantid))
# ah... once again, plot 1, xcoord 17

# (budphen) 46 3869_5_18C    3869_5    2023 FALSE FALSE TRUE  FALSE
phen = phen %>% mutate(finalid = gsub('3869\\_5\\_18C', '3869_5_18E', finalid))
# phen %>% filter(grepl('3869', finalid))
# demo %>% filter(grepl('3869', finalid))

# (budphen) 47 3882_15_1J    3882_15   2023 FALSE FALSE TRUE  FALSE
phen = phen %>% mutate(finalid = gsub('3882\\_15\\_1J', '3882_15_2I', finalid))
# phen %>% filter(grepl('3882', plantid))
# demo %>% filter(grepl('3882', plantid))

# 29 3888_5_18J    3888_5    2023 FALSE FALSE TRUE  FALSE
demo %>% filter(grepl('3888', plantid))
# dead umbel...
####### UNSOLVED #######
# add zero to seed?

# 30 3896_5_10G    3896_5    2023 FALSE FALSE TRUE  FALSE
demo %>% filter(grepl('3896', plantid))
# dead umbel
####### UNSOLVED #######
# add zero to seed?

# 31 3902_1_16I    3902_1    2023 FALSE FALSE TRUE  FALSE
demo %>% filter(grepl('3902', plantid)) # mismatch
phen %>% filter(grepl('3902', plantid))
seed %>% filter(grepl('3902\\_1', plantid.seed))
phen = phen %>% mutate(finalid = gsub('3902\\_1\\_16I', '3902_1_15J', finalid))

# (budphen) 51 5034_5_10H    5034_5    2023 FALSE FALSE TRUE  FALSE
# demo %>% filter(grepl('5034', plantid)) # nothing here...
# ah... demo says NP... shoot...
####### UNSOLVED #######
# check for notes in phen datasheet... oh maybe it died or something
# phen %>% filter(grepl('5034', plantid)) # no deaths here...

# (budphen) 52 5043_5_16D    5043_5    2023 FALSE FALSE TRUE  FALSE
phen = phen %>% mutate(finalid = gsub('5043\\_5\\_16D', '5043_5_15F', finalid))
# phen %>% filter(grepl('5043', plantid))
# demo %>% filter(grepl('5043', plantid))

# 32 5770_6_18H    5770_6    2023 FALSE FALSE TRUE  FALSE
# oh... lol
demo %>% filter(grepl('5[57]70', plantid)) # nowhere in here
seed %>% filter(grepl('5[57]70', plantid.seed))
# okie dokie, here's the tag misread
demo %>% filter(grepl('\\_6\\_1', plantid)) # I don't see anything in here...
phen = phen %>% mutate(finalid = gsub('5770', '5570', finalid))
phen %>% filter(grepl('5570', finalid))

# 33 7512_5_3D     7512_5    2023 FALSE FALSE TRUE  FALSE # fixed above

# 34 7543_5_6D     7543_5    2023 FALSE FALSE TRUE  FALSE
demo %>% filter(grepl('7543', plantid))
# accidental human-induced premature death... so maybe ignore?
# at least, it won't be in seed
phen %>% filter(grepl('7543', plantid))
phen = phen %>% mutate(finalid = gsub('7543\\_5\\_6D', '7543_5_6E', finalid))

# (budphen) 56 7546_5_16E    7546_5    2023 FALSE FALSE TRUE  FALSE
phen = phen %>% mutate(finalid = gsub('7546\\_5\\_16E', '7546_5_14G', finalid))
# phen %>% filter(grepl('7546', plantid))
# demo %>% filter(grepl('7546', plantid))

# (budphen) 57 7555_5_19B    7555_5    2023 FALSE FALSE TRUE  FALSE
phen = phen %>% mutate(finalid = gsub('7555\\_5\\_19B', '7555_5_18B', finalid))
# phen %>% filter(grepl('7555', plantid))
# demo %>% filter(grepl('7555', plantid))

# 35 7572_15_2I    7572_15   2023 FALSE FALSE TRUE  FALSE
demo %>% filter(grepl('7572', plantid)) # not here
demo %>% filter(grepl('\\_15\\_[123][A-Z]', plantid), Year %in% 2023)
# tag must have been modified at some point... not appearing in here though
####### UNSOLVED #######


####################################################
####################################################
# Start combining ##################################
####################################################
####################################################

# Ope I guess I should reconcile the IDs across years
# ugh!

merge(
  demo %>% 
    select(Year, finalid) %>% 
    separate(finalid, into = c("tag", "plot", "coord"), sep = '_', remove = FALSE) %>%
    select(-coord) %>%
    unite(c(tag, plot), col = 'tagplot', sep = '_'),
  seed %>%
    distinct(finalid, .keep_all = TRUE) %>%
    rename(Year = year) %>%
    select(Year, finalid) %>%
    separate(finalid, into = c("tag", "plot", "coord"), sep = '_', remove = FALSE, fill = 'right') %>%
    select(-coord) %>%
    unite(c(tag, plot), col = 'tagplot', sep = '_'),
  by = c('Year', 'tagplot'), suffixes = c('.demo', '.seed')
) %>%
  filter(finalid.demo != finalid.seed, !(Year %in% 2022))
# oh a good question... does any tagplot appear twice in a year?

merge(
  demo %>% 
    select(Year, finalid) %>% 
    separate(finalid, into = c("tag", "plot", "coord"), sep = '_', remove = FALSE) %>%
    select(-coord) %>%
    unite(c(tag, plot), col = 'tagplot', sep = '_'),
  seed %>%
    distinct(finalid, .keep_all = TRUE) %>%
    rename(Year = year) %>%
    select(Year, finalid) %>%
    separate(finalid, into = c("tag", "plot", "coord"), sep = '_', remove = FALSE, fill = 'right') %>%
    select(-coord) %>%
    unite(c(tag, plot), col = 'tagplot', sep = '_'),
  by = c('Year', 'tagplot'), suffixes = c('.demo', '.seed')
) %>%
  group_by(tagplot, Year) %>%
  filter(n() > 1)
# just these two from 2023
# and thanks to the work above we know that the IDs match, so we know which one is in seed

# do any of tag-plot combos appear for different plants in different years?
# ugh... this is difficult to assess
# the number of instances of this will be small... so I won't worry about this now
# (this would be an issue if two different plants got assigned to the same tag-plot label, throwing different plants into the same random intercept)

##### Here - combine the demo and umbel datasets
### I won't change the IDs yet because all of the ID fields are present here

de.um = merge(
  x = demo %>% select(Year, finalid, No.umbels, umbel.diam),
  y = mumb %>%
    group_by(Year, finalid) %>%
    summarise(
      n.umbel.diam = n(),
      n.diam.na    = sum(is.na(umble.diameter)),
      sum.diameter = sum(umble.diameter/2, na.rm = TRUE)
    ),
  by = c("Year", "finalid"), all.x = TRUE
)

head(de.um)
nrow(de.um)
table(de.um$Year) # good

# How many plants don't have an umbel diameter?
de.um %>%
  mutate(
    umbel.status = case_when(
      !(umbel.diam %in% 'SOS') & grepl('[0-9]', umbel.diam) ~ 'umbel.diam',
      !is.na(sum.diameter) & sum.diameter > 0 ~ 'sum.diameter',
      .default = 'other'
    )
  ) %>%
  group_by(Year, umbel.status) %>%
  summarise(n = n())
# hmm... plenty of 'other's

de.um %>%
  mutate(
    umbel.status = case_when(
      !(umbel.diam %in% 'SOS') & grepl('[0-9]', umbel.diam) ~ 'umbel.diam',
      !is.na(sum.diameter) & sum.diameter > 0 ~ 'sum.diameter',
      .default = 'other'
    )
  ) %>%
  filter(umbel.status %in% 'other')
# a couple of cases where no umbel diameter was entered in demo,
# very small number of cases where all of the recorded multi-umbels are NA

# How many plants have an NA-umbel
de.um %>% summarise(p.na.readings = mean(is.na(umbel.diam) | n.diam.na > 0, na.rm = TRUE))
# about 25%... ugh

de.um = de.um %>%
  mutate(
    demo.umbels = No.umbels,
    mumb.umbels = n.umbel.diam,
    diam.umbels = case_when(
      !is.na(sum.diameter) ~ sum.diameter,
      !(is.na(umbel.diam) | umbel.diam %in% 'SOS') ~ as.numeric(umbel.diam)/2,
      .default = NA
    ),
    n.in.measure = case_when(
      !is.na(sum.diameter) ~ n.umbel.diam - n.diam.na,
      !(is.na(umbel.diam) | umbel.diam %in% 'SOS') ~ 1,
      .default = NA
    )
  ) %>%
  select(Year, finalid, demo.umbels, mumb.umbels, diam.umbels, n.in.measure)
# warning message in here...
# not sure what is causing it. 

# oh well
# try to merge in seed set now?

# first, how many NAs are there?
seed %>% filter(is.na(no.seeds)) # three records with NAs here
# annoyingly I can't find the data sheets... so I guess... don't worry about these now
# probably will get removed at some point before analysis anyway

##### Combine seed and phen data
### These also should have consistent labels

# Look for NAs in phen
phen %>% filter(is.na(survey.period), is.na(init.doy))
# none

ph.sd = merge(
  x = phen %>% 
    group_by(Year = year, finalid) %>% 
    summarise(
      mean.doy = mean(init.doy),
      n.phen = n()
    ) %>%
    separate(finalid, into = c("tag", "plot", "coord"), sep = '_', remove = FALSE, fill = 'right') %>%
    select(-coord) %>%
    unite(c(tag, plot), col = 'tagplot', sep = '_'),
  y = seed %>%
    group_by(Year = year, finalid) %>%
    summarise(
      n.seed.counts = sum(!is.na(no.seeds)),
      n.empty = sum(!no.seeds & !is.na(no.seeds)),
      total.seeds = sum(no.seeds),
    ) %>%
    separate(finalid, into = c("tag", "plot", "coord"), sep = '_', remove = FALSE, fill = 'right') %>%
    select(-coord) %>%
    unite(c(tag, plot), col = 'tagplot', sep = '_'),
  by = c('Year', 'tagplot'), all.x = TRUE, all.y = TRUE, suffixes = c('.phen', '.seed')
) %>%
  mutate(
    n.seed.counts = ifelse(is.na(n.seed.counts), 0, n.seed.counts),
    n.phen = ifelse(is.na(n.phen), 0, n.phen)
  )

head(ph.sd)
# how well often do our phenology counts match the seed observations

ph.sd %>%
  group_by(Year, n.phen, n.seed.counts) %>%
  summarise(n = n()) %>%
  ggplot(aes(x = n.phen, y = n.seed.counts, fill = n, label = n)) +
  geom_tile() +
  geom_text(aes(colour = n.phen == n.seed.counts)) +
  scale_fill_gradient(low = 'black', high = 'gray') +
  guides(colour = 'none') +
  facet_wrap(~ Year) +
  theme(panel.background = element_blank())
# okay... surprising number of cases where there are more seed counts than
# phenology observations...

##### Add in seed data now
### But this is the step where we need to reconcile ID codings
# I'll make a new column for the final ID, and assign that tot he 'finalid' if
# it has only one record in demo for that year

de.um = de.um %>%
  separate(finalid, into = c('tag', 'plot', 'coord'), sep = '_', remove = FALSE) %>%
  select(-coord) %>%
  unite(c(tag, plot), col = 'tagplot')

de.um %>% group_by(Year, tagplot) %>% filter(n() > 1)
# only one duplicated ID

dusp = merge(
  x = de.um, y = ph.sd,
  by = c('Year', 'tagplot'),
  all.x = TRUE, all.y = TRUE
)

head(dusp)
nrow(dusp)

# Go in and remove duplicated demo in 2023

dusp %>%
  group_by(Year, tagplot) %>%
  filter(n() > 1)
# two plants here

dusp = dusp %>%
  group_by(Year, tagplot) %>%
  filter(n() == 1 | finalid == finalid.phen) %>%
  ungroup()

# Add in tag and plot info

dusp = dusp %>%
  mutate(plot = gsub('[0-9]{4}\\_', '', tagplot)) %>%
  merge(y = read.csv('00_raw_data/plot_treatments.csv'))

# Evaluate
dusp %>%
  filter(!is.na(total.seeds)) %>%
  ggplot(aes(x = total.seeds, group = interaction(trt, Year))) +
  geom_histogram(aes(fill = trt), position = position_identity(), alpha = 0.5) +
  scale_fill_manual(values = c('black', 'red', 'blue')) +
  facet_wrap(trt ~ Year, scales = 'free_y')
# not super interpretable
# does look like there are a lot of failures in the drought treatments
# distributions don't look super different

# # Export CSV
# write.csv(
#   dusp,
#   file = '01_data_cleaning/out/demo_seed_phen_combined.csv',
#   row.names = FALSE
# )

##### Exporting another version where there is one row per umbel *from seed count data*

ph.sd.by.umbel = merge(
  x = phen %>% 
    group_by(Year = year, finalid) %>% 
    summarise(
      mean.doy = mean(init.doy),
      n.phen = n()
    ) %>%
    separate(finalid, into = c("tag", "plot", "coord"), sep = '_', remove = FALSE, fill = 'right') %>%
    select(-coord) %>%
    unite(c(tag, plot), col = 'tagplot', sep = '_'),
  y = seed %>%
    rename(Year = year) %>%
    select(Year, finalid, no.seeds) %>%
    separate(finalid, into = c("tag", "plot", "coord"), sep = '_', remove = FALSE, fill = 'right') %>%
    select(-coord) %>%
    unite(c(tag, plot), col = 'tagplot', sep = '_'),
  by = c('Year', 'tagplot'), all.x = TRUE, all.y = TRUE, suffixes = c('.phen', '.seed')
)

head(ph.sd.by.umbel)

dusp.u = merge(
  x = de.um, y = ph.sd.by.umbel,
  by = c('Year', 'tagplot'),
  all.x = TRUE, all.y = TRUE
)

dusp.u = dusp.u %>%
  group_by(Year, tagplot) %>%
  # filter(length(unique(finalid)) > 1)
  filter(length(unique(finalid)) == 1 | finalid == finalid.phen) %>%
  ungroup()

head(dusp.u)
nrow(dusp.u)

dusp.u = dusp.u %>%
  mutate(plot = gsub('[0-9]{4}\\_', '', tagplot)) %>%
  merge(y = read.csv('00_raw_data/plot_treatments.csv'))

nrow(dusp.u)

# write.csv(
#   dusp.u,
#   file = '01_data_cleaning/out/demo_seed_phen_by_umbel_combined.csv',
#   row.names = FALSE
# )

####################################################
# For possible analysis ############################
####################################################

# 29 Jan 2024
# Thinking about the following analyses:
# - Number of umbels per plant (possibly also umbel diameter per plant?)
#   - want to test for phen and previous status (size, flowering) in here
#   - would then require demo (2020 and later), phen
# - Umbel probability of surviving/death
#   - testing for phen
#   - also maybe could merge with previou year's status?
# - Seeds per (non-dead) umbel
#   - need seeds, phen, probably demo?

### First: demo
# - want 2020 demo in here too... but didn't do any reconciliation here
# - but I think we can use the finalIDs from processing to merge back in...

demo.prev = read.csv('01_data_cleaning/out/demo_imputed_survival.csv') %>%
  # want to match year column in this data frame to the next year in the
  # processed demo frame
  mutate(year.match = Year + 1) %>%
  # get rid of plants that are before 2019 (we're matching 2021 and later)
  filter(year.match > 2020) %>%
  # remove the 'b's from plantids because I did this in the other df
  mutate(plantid = gsub('b', '', plantid))

head(demo)

# Data set for part 1:
demo.demoprev = merge(
  x = demo %>% select(Year, finalid, plantid, No.umbels, Plot, trt),
  y = demo.prev %>% 
    select(
      plantid, year.match, prev.umbels = No.umbels, 
      prev.leaves = No.leaves, prev.length =  Leaf.length
    ),
  by.x = c('Year', 'plantid'), by.y = c('year.match', 'plantid'),
  all.x = TRUE
)

head(demo.demoprev)

# How many of these plants are missing a previous size?
# Also how many have zero leaves...

demo.demoprev %>%
  mutate(
    case = case_when(
       is.na(prev.leaves) | is.na(prev.length) ~ 'missing.counts',
       !prev.leaves ~ 'dead',
       !is.na(prev.leaves) & !is.na(prev.length) ~ 'has.measures',
       .default = 'other'
    )
  ) %>%
  group_by(case) %>%
  summarise(n = n())
# mostly have measures
# and no dead plants!

demo.demoprev %>%
  mutate(
    case = case_when(
      is.na(prev.leaves) | is.na(prev.length) ~ 'missing.counts',
      !prev.leaves ~ 'dead',
      !is.na(prev.leaves) & !is.na(prev.length) ~ 'has.measures',
      .default = 'other'
    )
  ) %>%
  group_by(Year, case) %>%
  summarise(n = n()) %>%
  pivot_wider(names_from = case, values_from = n)
# Missing a lot of counts from 2020-2021
# (ah... the incomplete sampling in 2020... rats)

# What about umbel counts?
demo.demoprev %>%
  filter(!is.na(prev.leaves) & !is.na(prev.length)) %>%
  group_by(prev.umbels) %>%
  summarise(n = n())
# one NA

demo.demoprev %>% filter(!is.na(prev.leaves) & !is.na(prev.length) & is.na(prev.umbels))
demo.prev %>% filter(grepl('3596', plantid))
# no note here... I'm going to assume it didn't flower.

# Okay well I can live with this I guess.

# Filter out plants with missing measurements, add in 
demo.demoprev.measurements = demo.demoprev %>%
  # Filter out plants with missing records
  filter(!is.na(prev.leaves) & !is.na(prev.length)) %>%
  # Fix missing umbel count (assuming NA is zero)
  mutate(prev.umbels = ifelse(is.na(prev.umbels), 0, prev.umbels)) %>% 
  # Add helpful columns
  mutate(
    # Size in previous year
    prev.size = log(prev.leaves * prev.length),
    # Previously flowered
    prev.flwd = prev.umbels > 0
  )

head(demo.demoprev.measurements)
# Good
nrow(demo.demoprev.measurements)
# 962 observations

# ah... for an analysis of umbel counts would we not just want every year...?

### Second: demo (above) and phen
# - This will be used for probability of umbel death ~ flowering time and other demo
# - so will want to merge in demo.demoprev with other stuff
#   BUT because we don't know if this model will require previous size... let's
#   merge in the original (all plant) dataset rather than the one that has
#   measurements filtered out.
# - Also - want to merge all umbels data frame, correct?

head(phen)

phen.demo = merge(
  x = phen %>%
    group_by(year, finalid) %>%
    mutate(
      n.lost.umbels = sum(varb %in% 'lost'),
      n.phen.umbels  = sum(varb %in% 'new')
    ) %>%
    ungroup() %>%
    # Give me only the umbel budding/appearing dates (not the dead ones)
    filter(varb %in% 'new') %>%
    select(Year = year, plantid, init.doy, finalid, n.lost.umbels, n.phen.umbels) %>%
    # Need to change finalid column to match 2022
    # this is copied from code above
    separate(finalid, into = c("tag", "plot", "coord"), sep = '_', remove = FALSE, fill = 'right') %>%
    select(-coord) %>%
    unite(c(tag, plot), col = 'tagplot', sep = '_'),
  y = demo.demoprev %>% 
    rename(n.demo.umbels = No.umbels) %>%
    separate(finalid, into = c("tag", "plot", "coord"), sep = '_', remove = FALSE, fill = 'right') %>%
    select(-coord) %>%
    unite(c(tag, plot), col = 'tagplot', sep = '_'),
  by = c('Year', 'tagplot'), suffixes = c('.phen', '.demo')
) %>%
  # Change the one previous umbel count to zero
  mutate(
    prev.umbels = ifelse(is.na(prev.umbels) & !is.na(prev.leaves) & !is.na(prev.length), 0, prev.umbels)
  )

head(phen.demo)

# Check for duplicate plantids

phen.demo %>%
  group_by(Year, tagplot) %>%
  filter(length(unique(plantid.demo)) > 1)
# right now (30 jan) only one duplicate plant - neat!

phen.demo = phen.demo %>%
  group_by(Year, tagplot) %>%
  filter(length(unique(plantid.demo)) == 1 | plantid.demo == plantid.phen) %>%
  ungroup() %>%
  select(-contains('finalid'))

nrow(phen.demo)
table(phen.demo$Year)

phen.demo

# # Export this CSV
# write.csv(
#   phen.demo,
#   file = '01_data_cleaning/out/phen_demo_for_umbel_survival.csv',
#   row.names = FALSE
# )

### Third: demo (above), phen, and seeds
# - some thought should go in here...
# - seed counts by individual umbel should be good... but will want to go through
# and remove dead umbels and add non-dead empty umbels
