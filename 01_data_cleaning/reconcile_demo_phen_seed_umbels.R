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
# cool... all of these are missing.
# not sure which it is
####### UNSOLVED #######

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
####### need to change umbel count in demo (is NA, should be 1)

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
mumb = mumb %>% mutate(finalid = gusb('7598', '7698', finalid))

# 7 3069_7_1B     3069_7    2023 FALSE TRUE  FALSE FALSE
seed %>% filter(grepl('3069', plantid.seed))
demo %>% filter(grepl('3069', plantid))
phen %>% filter(grepl('3069', plantid)) # not in phen though...
phen %>% filter(grepl('\\_7\\_[12][ABC]', plantid))
# oh right... umbel is dead.. well even still
seed = seed %>% mutate(finalid = gsub('3069\\_7\\_1B', '3069_7_1C', finalid))

# 8 3297_15_11H   3297_15   2023 FALSE TRUE  FALSE FALSE
# 9 3454_6_15I    3454_6    2023 FALSE TRUE  FALSE FALSE
# 10 3518_3_16J    3518_3    2023 FALSE TRUE  FALSE FALSE
# 11 3590_13_9B    3590_13   2023 FALSE TRUE  FALSE FALSE
# 12 3636_7_3I     3636_7    2023 FALSE TRUE  FALSE FALSE
# 13 5570_6_18H    5570_6    2023 FALSE TRUE  FALSE FALSE
# 14 7519_5_3D     7519_5    2023 TRUE  TRUE  FALSE TRUE 
# 15 3182_4_17H    3182_4    2023 FALSE FALSE TRUE  FALSE
# 16 3192_1_11B    3192_1    2023 TRUE  FALSE TRUE  TRUE 
# 17 3193_1_15A    3193_1    2023 FALSE FALSE TRUE  FALSE
# 18 3424_1_12G    3424_1    2023 TRUE  FALSE TRUE  TRUE 
# 19 3432_1_19A    3432_1    2023 TRUE  FALSE TRUE  TRUE 
# 20 3463_14_6A    3463_14   2023 FALSE FALSE TRUE  FALSE
# 21 3477_2_5C     3477_2    2023 FALSE FALSE TRUE  FALSE
# 22 3519_2_8G     3519_2    2023 FALSE FALSE TRUE  FALSE
# 23 3556_14_17I   3556_14   2023 FALSE FALSE TRUE  FALSE
# 24 3688_7_6G     3688_7    2023 FALSE FALSE TRUE  FALSE
# 25 3697_5_17H    3697_5    2023 FALSE FALSE TRUE  FALSE
# 26 3763_1_12A    3763_1    2023 TRUE  FALSE TRUE  TRUE 
# 27 3807_5_10I    3807_5    2023 FALSE FALSE TRUE  FALSE
# 28 3850_1_17E    3850_1    2023 TRUE  FALSE TRUE  TRUE 
# 29 3888_5_18J    3888_5    2023 FALSE FALSE TRUE  FALSE
# 30 3896_5_10G    3896_5    2023 FALSE FALSE TRUE  FALSE
# 31 3902_1_16I    3902_1    2023 FALSE FALSE TRUE  FALSE
# 32 5770_6_18H    5770_6    2023 FALSE FALSE TRUE  FALSE
# 33 7512_5_3D     7512_5    2023 FALSE FALSE TRUE  FALSE
# 34 7543_5_6D     7543_5    2023 FALSE FALSE TRUE  FALSE
# 35 7572_15_2I    7572_15   2023 FALSE FALSE TRUE  FALSE
