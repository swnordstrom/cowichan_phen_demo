# Script for taking a first stab at reconciling cleaned data
# SN - 24 Oct 2023

library(ggplot2)
library(dplyr)
library(tidyr)

rm(list = ls())

######################################################################
##### Load in processed data
######################################################################

# All demo data
demo = read.csv('01_data_cleaning/out/demo_all_cleaned.csv')

# Multi-umbel data (2021-2023)
mumb = read.csv('01_data_cleaning/out/multi_umbel_clean.csv')

# Seed count data
seed = read.csv('01_data_cleaning/out/seed_counts_all_cleaned.csv')

######################################################################
##### Start combining
######################################################################


#######################
##### Demo + umbel
#######################

head(demo)
head(mumb)

unique(demo$plantid)
unique(mumb$plantid)

nrow(demo)
nrow(mumb)

demo.mumb = merge(
  x = demo %>% 
    # Get relevant columns
    select(Year, No.leaves, Leaf.length, No.umbels, umbel.diam, edited, Note, plantid) %>%
    # Rename notes column to specify source
    rename(
      demo.note = Note,
      demo.edit = edited
    ) %>%
    # This column will help with assessing the merge...
    mutate(in.demo = TRUE),
  y = mumb %>% 
    # Get relevant columns
    select(Year, umble.no, umble.diameter, notes, plantid) %>%
    # Rename column to specify source
    rename(mumb.note = notes) %>%
    # This column will help in merge assessment
    mutate(in.mumb = TRUE),
  by = c('Year', 'plantid'),
  all.x = TRUE, all.y = TRUE
)

nrow(demo.mumb)
# Seems promising...
head(demo.mumb)

###
### Assess merge
###

demo.mumb %>% group_by(in.demo, in.mumb) %>% summarise(n = n())

# The in.mumb = NA cases should all be pre-2020
demo.mumb %>% filter(Year > 2019) %>% filter(is.na(in.mumb))
# (ah... or they should have umbel count 0 or 1 in or after 2020)
demo.mumb %>% 
  filter(Year > 2019, (!is.na(No.umbels) & No.umbels > 1)) %>% 
  filter(is.na(in.mumb))

### 
### Handle these cases: plants that should be in the "mumb" dataset but are not in demo
### (clean them here rather than in original data-cleaning scripts)

missing.mumb = demo.mumb %>% 
  filter(Year > 2019, (!is.na(No.umbels) & No.umbels > 1)) %>% 
  filter(is.na(in.mumb))

missing.mumb %>% select(Year, plantid)

### (2020 data first)

### 1  2020    3563_4_15H
mumb %>% filter(Tag %in% 3563)
# ah... this is one where the coords got changed
mumb = mumb %>%
  mutate(plantid = ifelse(Tag %in% 3563 & Plot %in% 4, '3563_4_15H', plantid))

### 2  2020    3634_4_18E
mumb %>% filter(Tag %in% 3634)
mumb = mumb %>%
  mutate(plantid = ifelse(Tag %in% 3634 & Plot %in% 4, '3634_4_18E', plantid))

### 3  2020    3747_6_18A
mumb %>% filter(Tag %in% 3747)
demo %>% filter(Tag %in% 3747) # 11 umbels??
mumb %>% filter(umble.no > 5) # there is one plant with 11 umbels but it's from 2023
# data entry error? check data sheet
demo %>% filter(No.umbels > 10)
# ah... looks like recording error in 2020
# maybe counting number of flowers instead of number of umbels..?

# how many of the below are simply double digit umbel-plants?

### 4  2020     3749_2_6B
### 5  2020     3751_2_8A
### 6  2020    3784_6_18F
### 7  2020    3785_6_15G
### 8  2020     3786_4_6H
### 9  2020     3787_4_3G
### 10 2020     3788_5_0F
### 11 2020     3790_6_5I
### 12 2020     3791_6_9J

missing.mumb %>% filter(Year %in% 2020, No.umbels > 10) # that's all of them...

### 1  2021    3053_2_14B
mumb %>% filter(Tag %in% 3053)
mumb %>% filter(Plot %in% 2, Year %in% 2021)
# there *is* a plant with two umbels in plot 2 at coords 14B
demo %>% filter(Tag %in% 3164, Plot %in% 2)
demo.mumb %>% filter(Year %in% 2021, grepl('\\_2\\_14B', plantid))
# nope - different plants
# I guess we should assume there is no 3053 in the multi-umbel dataset
demo %>% filter(Tag %in% 3053)
# ah - note in original dataset: "umbels dead"
demo = demo %>%
  mutate(Note = ifelse(Year %in% 2021 & plantid %in% '3053_2_14B', 'umbels dead', Note))

### 2  2021    3135_12_7G
# possible data entry error?

### 3  2021 3143_9_140.95
mumb %>% filter(Tag %in% 3143) # it's just not here in 2021 (but is in other years...)
mumb %>% filter(Plot %in% 9)
# demo has note "14T" - what does this mean...?

### 4  2021    3163_2_11D
mumb %>% filter(Tag %in% 3163) # not here at all
demo %>% filter(Tag %in% 3163)
mumb %>% filter(Plot %in% 2, Year %in% 2021)
# data entry error?

### 5  2021    3352_2_19H
mumb %>% filter(Tag %in% 3352) # here in 2022, not 2021
demo %>% filter(Tag %in% 3352)

### 6  2021    3360_6_12C
mumb %>% filter(Tag %in% 3360) # here in 2022-23
demo %>% filter(Tag %in% 3360)

### 7  2021     3453_7_0C
mumb %>% filter(Tag %in% 3453) # here in 2022, not 2021
demo %>% filter(Tag %in% 3453) # no umbel diameter...
# no note, no umbels listed

### 8  2021     3543_1_2F
mumb %>% filter(Tag %in% 3543) # in 2023, not 2021
demo %>% filter(Tag %in% 3543)

### 9  2021     3567_1_3I
mumb %>% filter(Tag %in% 3567) # present in 2022, 2023
demo %>% filter(Tag %in% 3567)

### 10 2021     3584_2_7F
mumb %>% filter(Tag %in% 3584) # not in mumbel at all
demo %>% filter(Tag %in% 3584)

### 11 2021   3675_13_13I
mumb %>% filter(Tag %in% 3675)
demo %>% filter(Tag %in% 3675)

### 12 2021   3736_15_14B
mumb %>% filter(Tag %in% 3736)
demo %>% filter(Tag %in% 3736)

### 13 2021     3814_1_0J
mumb %>% filter(Tag %in% 3814)
demo %>% filter(Tag %in% 3814)

### 14 2021   3868_15_15A
mumb %>% filter(Tag %in% 3868)
demo %>% filter(Tag %in% 3868)

### 15 2021    3978_13_2F
mumb %>% filter(Tag %in% 3978)
demo %>% filter(Tag %in% 3978)

### 16 2022    3469_14_3A
mumb %>% filter(Tag %in% 3469) # doesn't appear until 2023
demo %>% filter(Tag %in% 3469, Plot %in% 14)
# ah - note in 2022 datasheet says second umbel clipped
demo = demo %>%
  mutate(Note = 
           ifelse(Year %in% 2022 & Tag %in% 3469 & Plot %in% 14,
                  'one umbel clipped',
                  Note)
         )

### 17 2022     3542_1_0I
mumb %>% filter(Tag %in% 3542)
demo %>% filter(Tag %in% 3542) # says "see other sheet"...
mumb %>% filter(Year %in% 2022, Plot %in% 1) %>% distinct(Tag, .keep_all = TRUE)
# was this plant mistaken for 3365?
demo %>% filter(Tag %in% 3365, Year %in% 2022)
# no I don't think so
# hmm... oh well guess it's lost

### 18 2022     3814_1_0J
mumb %>% filter(Tag %in% 3814)
demo %>% filter(Tag %in% 3814) # says "see other sheet"...
mumb %>% filter(Plot %in% 1, Year %in% 2022) %>% distinct(Tag, .keep_all = TRUE)
# maybe never got entered

### 19 2023    3391_13_7G
mumb %>% filter(Tag %in% 3391)
# hmm... address is totally off...
demo %>% filter(Tag %in% 3391)
# okay at some point I went through and changed this in the demo data
# assuming that the 7K plant should be 7G
mumb = mumb %>%
  mutate(plantid = ifelse(plantid %in% '3391_13_7K', '3391_13_7G', plantid))

### 20 2023    7698_10_9D
mumb %>% filter(Tag %in% 7698)
demo %>% filter(Tag %in% 7698) # tag might have been changed
mumb %>% filter(Year %in% 2023, Plot %in% 10) # not seeing anything at x9
# no notes in original demo about tag being replaced... looks like new tag
# assume it's lost I guess

### Next, look at cases where in.demo is missing
demo.mumb %>% filter(is.na(in.demo))
# Seems doable!

### 1 2020  3563_4_15R
demo %>% filter(Tag %in% 3563)
# changed coordinate in demo
mumb = mumb %>%
  mutate(
    plantid = ifelse(Tag %in% 3563 & Plot %in% 4, '3563_4_15H', plantid),
    notes   = ifelse(Tag %in% 3563 & Plot %in% 4, 'plantid corrected', notes)
  )

### 2 2020  3634_4_18O
demo %>% filter(Tag %in% 3634)
# yes - changed coord
mumb = mumb %>%
  mutate(
    plantid = ifelse(Tag %in% 3634 & Plot %in% 4, '3634_4_18E', plantid),
    notes   = ifelse(Tag %in% 3634 & Plot %in% 4, 'plantid corrected', notes)
  )

### 3 2021  3563_4_15R
# fixed this above

### 4 2022 3389_13_12F
demo %>% filter(Tag %in% 3389)
# tag was changed somewhere
demo %>% filter(grepl('3389', Note))
mumb = mumb %>%
  mutate(
    plantid = ifelse(Tag %in% 3389 & Plot %in% 13, '3481_13_12F', plantid),
    notes   = ifelse(Tag %in% 3389 & Plot %in% 13, 'plantid corrected', notes)
  )

mumb %>% filter(Tag %in% 3389)

### 5 2022  3634_4_18O
# fixed above
### 6 2023 3389_13_12F
# fixed above
### 7 2023  3391_13_7K
# fixed above!


#######################
##### Demo + umbel (for real this time)
#######################

demo.mumb = merge(
  x = demo %>% 
    # Get relevant columns
    select(Year, No.leaves, Leaf.length, No.umbels, umbel.diam, edited, Note, plantid) %>%
    # Rename notes column to specify source
    rename(
      demo.note = Note,
      demo.edit = edited
    ) %>%
    # This column will help with assessing the merge...
    mutate(in.demo = TRUE),
  y = mumb %>% 
    # Get relevant columns
    select(Year, umble.no, umble.diameter, notes, plantid) %>%
    # Rename column to specify source
    rename(mumb.note = notes) %>%
    # This column will help in merge assessment
    mutate(in.mumb = TRUE),
  by = c('Year', 'plantid'),
  all.x = TRUE, all.y = TRUE
)

#######################
##### Demo + seed?
#######################

# It looks like only two umbels were observed per plant...

seed %>%
  group_by(year, plantid) %>%
  summarise(n = n()) %>%
  group_by(year, n) %>%
  summarise(n.plants = n())
# oh interesting...

# also interestingly - a lot more seed data in 2022 than 2021?
table(seed$year)
# huh... would be good to enter the 2023 data then I guess!
# might be a lot of data here!

# Look year-by-year

# Just 2021
merge(
  x = mumb %>% 
    filter(Year %in% 2021) %>%
    group_by(plantid) %>% 
    summarise(n1 = n()),
  y = seed %>% 
    filter(year %in% 2021) %>%
    group_by(plantid) %>% 
    filter(n() > 1) %>%
    summarise(n2 = n()),
  all.y = TRUE
)
# this actually looks fine
# but then - why are they nearly all two? should there be more...?

# Will also need to go through at some point and check the 2022 data (missing
# coordinates)

demo.seed.2022 = merge(
  x = seed %>%  filter(year %in% 2022),
  y = demo %>% 
    filter(Year %in% 2022) %>% 
    mutate(tag_plot = paste(Tag, Plot, sep = "_")) %>%
    select(-plantid),
  all.x = TRUE, by.x = "plantid", by.y = "tag_plot"
)

nrow(demo.seed.2022)
demo.seed.2022

demo.seed.2022 %>%
  group_by(plantid) %>%
  filter(any(duplicated(umbel.no)))
# only three cases here!
# and all of them have only one flowering plant! cha ching?

######################################################################
######################################################################

# Should be true that all cases where "umbel.diam" is SOS should have non-NA
# umble.diameter field
with(demo.mumb, table(umbel.diam, umble.diameter, useNA = 'always')) # not useful...
# line below should give an empty data frame
demo.mumb %>% filter(umbel.diam %in% 'SOS' & is.na(umble.diameter)) # damn that's a lot...
demo.mumb %>%
  filter(umbel.diam %in% 'SOS') %>%
  group_by(is.na(umble.diameter)) %>%
  summarise(n = n()) # okay so only a small number have diameter NA...

# Cases where umbel.diam is *not* SOS - do we have anything in umble.diameter?
demo.mumb %>%
  filter(!umbel.diam %in% 'SOS') %>%
  group_by(is.na(umble.diameter)) %>%
  summarise(n = n())
# why... (maybe umbel.diam is NA? merge failure?)