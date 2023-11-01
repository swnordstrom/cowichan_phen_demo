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

##### (assessing data only here)

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
  filter(any(duplicated(umbel.no))) %>%
  select(-c(year, Year, Plot, Tag)) 
# only three cases here!
# and all of them have only one flowering plant! cha ching?

# How many of these are listed as no flowering?

demo.seed.2022 %>% filter(is.na(No.umbels) | No.umbels < 1)
# concerning...
# maybe the better matchup would be to the phen dataset...
# actually I wonder how many of these are due to tags

### What about 2023 data?

demo.seed.2023 = merge(
  x = seed %>% filter(year %in% 2023) %>% mutate(in.seed = TRUE),
  y = demo %>% filter(Year %in% 2023) %>% mutate(in.demo = TRUE),
  all.x = TRUE, by = "plantid"
) %>%
  mutate(across(starts_with('in'), function(x) ifelse(is.na(x), FALSE, x)))

demo.seed.2023
table(demo.seed.2023$in.demo)
# most of these are not in demo...
# probably poor tag matching... ugh

demo.seed.2023 %>%
  distinct(plantid, .keep_all = TRUE) %>%
  group_by(in.demo, in.seed) %>%
  summarise(n = n())

##### Try to merge all of these together

demo.seed = merge(
  x = seed %>% mutate(in.seed = TRUE),
  y = demo %>%
    filter(Year > 2020) %>%
    mutate(flowering.in.demo = !is.na(No.umbels) & No.umbels > 0) %>%
    mutate(coor.demo = paste0(Xcoor, Ycoor)) %>%
    select(-c(Xcoor, Ycoor, YrTag)) %>%
    mutate(in.demo = TRUE),
  by.x = c("year", "plot", "tag"), by.y = c("Year", "Plot", "Tag"),
  all.x = TRUE, all.y = TRUE, suffixes = c('.seed', '.demo')
) %>%
  mutate(across(starts_with('in'), function(x) ifelse(is.na(x), FALSE, x)))

head(demo.seed)
with(demo.seed, table(in.demo, in.seed, year))
# Ugh...

# Okay... what could be going wrong in merging:
# (1) merge duplicates phen and/or demo records
# (2) seed counts for non-flowering demo records
# (3) plain ol' misses...

# (1) merging producing duplicate records
# If this happened, then we should be getting cases where a plantid.demo gets
# matched to >1 unique plantid.seed (and vice versa)

demo.seed %>%
  filter(!is.na(plantid.demo)) %>%
  distinct(year, plantid.demo, plantid.seed, .keep_all = TRUE) %>%
  group_by(year, plantid.demo) %>%
  filter(n() > 1) %>%
  arrange(plantid.demo) %>% distinct(plantid.demo)
# no duplicated demo records?

demo.seed %>%
  filter(!is.na(plantid.seed)) %>%
  distinct(year, plantid.demo, plantid.seed, .keep_all = TRUE) %>%
  group_by(year, plantid.seed) %>%
  filter(n() > 1) %>%
  arrange(plantid.seed) %>%
  select(year, plot, tag, plantid.seed, coor.demo, plantid.demo, flowering.in.demo)
# a couple

# The only one that *needs* to be fixed is the 3125_5 in 2023
# others can be fixed by picking out the "flowering" ones.
# Which brings us to our next issue...

# (2) seed counts for non-flowering demo records
# (feel like there will be a ton of these...)
# Bad record would be: in.seed = TRUE _and_ flowering.in.demo = FALSE

demo.seed %>%
  filter(!is.na(plantid.demo) & !is.na(plantid.seed)) %>%
  distinct(year, plantid.demo, plantid.seed, .keep_all = TRUE) %>%
  group_by(year, flowering.in.demo, in.seed)  %>%
  summarise(n = n()) %>%
  pivot_wider(names_from = year, values_from = n)
# actually this looks okay!  

demo.seed %>%
  filter(!is.na(plantid.demo) & !is.na(plantid.seed)) %>%
  distinct(year, plantid.demo, plantid.seed, .keep_all = TRUE) %>%
  group_by(year, plantid.seed) %>%
  filter(any(!flowering.in.demo)) %>%
  select(year, plot, tag, coor.demo, plantid.seed, No.leaves, Stalk_Height, 
         No.umbels, umbel.no, no.seeds, Note)

# Some of these are just dupes, a couple of them are not...
# Some of them are not though...
# Adds up to maybe ten cases? honestly not bad considering how large the dataset is
# Looks also like only one of these plants has a non-zero number of seeds
# (3365_1_0J) - rest might be dead/eaten

# (3) plain ol' misses...
# all of these plants should have demo records of flowering

demo.seed %>%
  filter(!in.demo | !in.seed) %>%
  filter(flowering.in.demo | is.na(flowering.in.demo)) %>%
  distinct(year, plantid.demo, plantid.seed, .keep_all = TRUE) %>%
  group_by(year) %>%
  summarise(n = n())
# that's a lot..

# in 2021, not every plant was surveyed for seed
# (I think they all were in 2022-2023 though)
demo.seed %>%
  filter(flowering.in.demo | is.na(flowering.in.demo)) %>%
  filter(!in.demo | !in.seed) %>%
  distinct(year, plantid.demo, plantid.seed, .keep_all = TRUE) %>%
  group_by(year, in.demo, in.seed) %>%
  summarise(n = n())

# Maybe should be more concerned with in.seed but not in.demo?
# assume the in.demo but not in.seed just never had measurements taken

demo.seed %>%
  filter(flowering.in.demo | is.na(flowering.in.demo)) %>%
  group_by(year, in.demo, in.seed) %>%
  # group_by(
  #   year, 
  #   in.demo = ifelse(in.demo, 'in.demo', 'not.in.demo'), 
  #   in.seed = ifelse(in.seed, 'in.seed', 'not.in.seed')
  # ) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  mutate(case = case_when(
    in.demo & in.seed ~ 'currently usable',
    in.demo & !in.seed ~ 'probably no seed record',
    !in.demo & in.seed ~ 'seed missed in merge'
    )
  ) %>%
  select(-c(in.demo, in.seed)) %>%
  pivot_wider(names_from = case, values_from = n)

# Okay - have 40-50 seed counts that got missed when merging

##### Try fix cases identified above
# To fix:
# - seed counts not matched to demo records
# - duplicated seed counts (poorly-matched to demo)
# - (maybe also go back and check the flowering records)

### First, look for seed counts not matched to demo records

demo.seed %>%
  filter(flowering.in.demo | is.na(flowering.in.demo)) %>%
  filter(in.seed & !in.demo) %>%
  distinct(year, plantid.seed)

#   year plantid.seed
# 1  2021    3748_6_5E

demo %>% filter(Year %in% 2021 & Tag %in% 3748)
demo %>% filter(Year %in% 2021, Plot %in% 6) %>% arrange(Xcoor)
seed %>% filter(year %in% 2021, tag %in% 3748)
# noted as dead... definitely flowered though...
demo %>% filter(grepl('3748', Note)) # ice cold... not even mentioned in a note...
demo %>% filter(Tag %in% 3748) # I have literally no records of this tag... anywhere
demo.seed %>% filter(year %in% 2021, !in.seed, plot %in% 6, flowering.in.demo)
# could be 3179
demo %>% filter(Tag %in% 3179) # tag replacement... possible here
# yes... note in 2022 says replacement tag 3748 was used in 2021

seed = seed %>%
  mutate(
    plantid = ifelse(tag %in% 3748 & year %in% 2021, gsub('3748', '3179', plantid), plantid),
    notes   = ifelse(tag %in% 3748 & year %in% 2021, '[tag manually edited; 2022 note; was 3748]', notes),
    tag     = ifelse(tag %in% 3748 & year %in% 2021, 3179, tag)
  )

# 2  2021    3928_6_7C
demo %>% filter(Year %in% 2021, Tag %in% 3928) # tag misid or entry
demo.seed %>% filter(year %in% 2021, !in.seed, plot %in% 6, flowering.in.demo)
seed %>% filter(tag %in% 3928)
# one umbel counted in seed, one plant at location 7C in 2021 (also w/ one umbel)
# going to say these are the same
# (datasheet says 3918... assuming this is misread of 3916)
seed = seed %>%
  mutate(
    plantid = ifelse(tag %in% 3928 & year %in% 2021, gsub('3928', '3916', plantid), plantid),
    notes   = ifelse(tag %in% 3928 & year %in% 2021, '[tag manually edited; was 3928]', notes),
    tag     = ifelse(tag %in% 3928 & year %in% 2021, 3916, tag)
  )

# actually going in and looking at demo... yeah there's a tag inconsistency here too
demo = demo %>%
  mutate(
    plantid = ifelse(Tag %in% 3918 & Plot %in% 6, gsub('3918', '3916', plantid), plantid),
    edited = ifelse(Tag %in% 3918 & Plot %in% 6, TRUE, edited),
    Note = ifelse(Tag %in% 3918 & Plot %in% 6, 'tag mis-entered as 3918 in 2016', Note),
    Tag = ifelse(Tag %in% 3918 & Plot %in% 6, 3916, Tag),
  )

# 3  2021    3298_7_9C
demo %>% filter(Tag %in% 3298) # lol
demo.seed %>% filter(year %in% 2021, !in.seed, plot %in% 7, flowering.in.demo)
# seems very likely that it's 3296 (same issue as above...)
seed %>% filter(tag %in% 3298)

seed = seed %>%
  mutate(
    plantid = ifelse(tag %in% 3298 & year %in% 2021, gsub('3298', '3296', plantid), plantid),
    notes   = ifelse(tag %in% 3298 & year %in% 2021, '[tag manually edited; was 3298]', notes),
    tag     = ifelse(tag %in% 3298 & year %in% 2021, 3296, tag)
  )

demo = demo %>%
  mutate(
    plantid = ifelse(Tag %in% 3298 & Plot %in% 7, gsub('3298', '3296', plantid), plantid),
    edited = ifelse(Tag %in% 3298 & Plot %in% 7, TRUE, edited),
    Note = ifelse(Tag %in% 3298 & Plot %in% 7, 'tag mis-entered as 3298 in 2016', Note),
    Tag = ifelse(Tag %in% 3298 & Plot %in% 7, 3296, Tag),
  )

# 4  2021  3363_14_15H
demo %>% filter(Tag %in% 3363) # tag was misentered
seed %>% filter(tag %in% 3363) # mis-entered in both years?
demo %>% filter(grepl('3363', Note)) # tag wasn't manually changed...
demo.seed %>% filter(plot %in% 14, year %in% 2021, !in.seed & flowering.in.demo)
demo.seed %>% filter(plot %in% 14, year %in% 2022, !in.seed & flowering.in.demo)
seed %>% filter(tag %in% 3363) # yes total mis-read/mis-entry
demo %>% filter(Tag %in% 3367) # lol

seed = seed %>%
  mutate(
    plantid = ifelse(tag %in% 3363 & plot %in% 14 & year %in% 2021, gsub('3363', '3367', plantid), plantid),
    notes   = ifelse(tag %in% 3363 & plot %in% 14 & year %in% 2021, '[tag manually edited; was 3363]', notes),
    tag     = ifelse(tag %in% 3363 & plot %in% 14 & year %in% 2021, 3367, tag)
  )

# 5  2022       3782_1
demo %>% filter(Tag %in% 3782)
# note in 2022 says this plant is not new and was mistakenly confused with 3814
demo %>% filter(Tag %in% 3814)
seed %>% filter(tag %in% 3814)
seed %>% filter(tag %in% 3782) # no additional info here
demo.seed %>% filter(year %in% 2022, plot %in% 1, !in.seed, flowering.in.demo)
# hmm... not sure there's a fix for this one 

# 6  2022       5022_5
demo %>% filter(Tag %in% 5022) # haha...
seed %>% filter(tag %in% 5022)
demo %>% filter(Plot %in% 5, Xcoor %in% 8:10, Year %in% 2022)
# lmao... was it 2055? dyscalcula
demo %>% filter(Tag %in% 2055) # uh
demo %>% filter(grepl('^2', Tag))
demo %>% filter(grepl('^5', Tag)) # I am no longer convinced 5022 = 2055
demo.seed %>% filter(plot %in% 5, year %in% 2022, !in.seed, flowering.in.demo)
# okay I'm convinced again

# but are these the same plant? do other records need to be changed?
demo %>% filter(Tag %in% c(5022, 2055)) %>% arrange(Year, Tag)
seed %>% filter(tag %in% c(2055, 5022)) # hmm... okay I think they're the same
# multiple demo records have "tag is 5022" so I'll change the demo

demo = rbind(
  demo %>% filter(!(Tag %in% 2055 & Plot %in% 5)),
  demo %>%
    filter(Tag %in% 2055 & Plot %in% 5) %>%
    mutate(
      # Add new (correct) tag
      Tag = 5022,
      # update plantid (not sure if this will be done elsewhere)
      plantid = gsub('2055', '5022', plantid),
      # Add note
      Note = "tag mis-recorded as 2055",
      edited = TRUE
    )
)


# 7  2022       7536_5
# see above
seed %>% filter(tag %in% 7526)
seed = seed %>%
  mutate(
    plantid = ifelse(tag %in% 7536 & year %in% 2022, gsub('7536', '7526', plantid), plantid),
    notes   = ifelse(tag %in% 7536 & year %in% 2022, '[tag manually edited; was 7536]', notes),
    tag     = ifelse(tag %in% 7536 & year %in% 2022, 7526, tag)
  )

# 8  2022       3354_6
demo %>% filter(Tag %in% 3354, Plot %in% 6) # misentry
seed %>% filter(tag %in% 3354)
demo.seed %>% filter(plot %in% 6, year %in% 2022, !in.seed, flowering.in.demo)
# 3374 seems like a safe bet
seed = seed %>%
  mutate(
    plantid = ifelse(tag %in% 3354 & year %in% 2022 & plot %in% 6, gsub('3354', '3374', plantid), plantid),
    notes   = ifelse(tag %in% 3354 & year %in% 2022 & plot %in% 6, '[tag manually edited; was 3354]', notes),
    tag     = ifelse(tag %in% 3354 & year %in% 2022 & plot %in% 6, 3374, tag)
  )

# 9  2022       3789_6
demo %>% filter(Tag %in% 3789)
demo %>% filter(grepl('3789', Note)) # not a replacement
seed %>% filter(tag %in% 3789)
demo.seed %>% filter(plot %in% 6, year %in% 2022, !in.seed, flowering.in.demo)
# oh... might have been the worst tag-reading job on earth (7->1, 6->8, 4->9)
demo %>% filter(Tag %in% 3164) # plot six...
seed %>% filter(tag %in% 3164) # wait... already fixed?
# oh - no that's for a plant in plot 2 (this one is in plot 6)
seed = seed %>%
  mutate(
    plantid = ifelse(tag %in% 3789 & year %in% 2022 & plot %in% 6, gsub('3789', '3164', plantid), plantid),
    notes   = ifelse(tag %in% 3789 & year %in% 2022 & plot %in% 6, '[tag manually edited; was 3789]', notes),
    tag     = ifelse(tag %in% 3789 & year %in% 2022 & plot %in% 6, 3164, tag)
  )

# 10 2022       3298_7
# see above
seed = seed %>%
  mutate(
    plantid = ifelse(tag %in% 3298 & year %in% 2022, gsub('3298', '3296', plantid), plantid),
    notes   = ifelse(tag %in% 3298 & year %in% 2022, '[tag manually edited; was 3298]', notes),
    tag     = ifelse(tag %in% 3298 & year %in% 2022, 3296, tag)
  )

# 11 2022       3886_7
demo %>% filter(Tag %in% 3886) # hmm... tag change?
seed %>% filter(tag %in% 3886)
# note in 2022 record for 3866 says typo 3886...
demo %>% filter(Tag %in% c(3886, 3866))
# okay - looks like it was simply mis-entered in 3866
# hmm... change the 2016 demo record and the seed data

seed = seed %>%
  mutate(
    plantid = ifelse(tag %in% 3886 & year %in% 2022, gsub('3886', '3866', plantid), plantid),
    notes   = ifelse(tag %in% 3886 & year %in% 2022, '[tag manually edited; was 3886]', notes),
    tag     = ifelse(tag %in% 3886 & year %in% 2022, 3866, tag)
  )

demo = demo %>%
  mutate(
    plantid = ifelse(Tag %in% 3886 & Plot %in% 7, gsub('3886', '3866', plantid), plantid),
    edited = ifelse(Tag %in% 3886 & Plot %in% 7, TRUE, edited),
    Note = ifelse(Tag %in% 3886 & Plot %in% 7, 'tag mis-entered as 3866 in 2016', Note),
    Tag = ifelse(Tag %in% 3886 & Plot %in% 7, 3866, Tag),
  )

# 12 2022      3614_13
demo %>% filter(Tag %in% 3614)
demo %>% filter(grepl('3614', Note))
seed %>% filter(tag %in% 3614)
demo.seed %>% filter(plot %in% 13, year %in% 2022, !in.seed, flowering.in.demo)
# haha... man
# should be 3164? (one of the other 3164s lol)
demo %>% filter(Tag %in% 3164)
# yep... there it is

seed = seed %>%
  mutate(
    plantid = ifelse(tag %in% 3614 & plot %in% 13 & year %in% 2022, gsub('3614', '3164', plantid), plantid),
    notes   = ifelse(tag %in% 3614 & plot %in% 13 & year %in% 2022, '[tag manually edited; was 3614]', notes),
    tag     = ifelse(tag %in% 3614 & plot %in% 13 & year %in% 2022, 3164, tag)
  )

# 13 2022      7537_13
demo %>% filter(Tag %in% 7537)
seed %>% filter(tag %in% 7537)
demo.seed %>% filter(plot %in% 13, year %in% 2022, !in.seed, flowering.in.demo) # ugh...
demo %>% filter(Tag %in% 3880) # there are several...
# will want to edit all demo, no?
# also check to see if old tag appears in seed records...
seed %>% filter(tag %in% 3880, plot %in% 13) # eh change anyway to be safe

demo = rbind(
  demo %>% filter(!(Tag %in% 3880 & Plot %in% 13)),
  demo %>%
    filter(Tag %in% 3880 & Plot %in% 13) %>%
    mutate(
      # Add new (correct) tag
      Tag = 7537,
      # update plantid (not sure if this will be done elsewhere)
      plantid = gsub('3880', '7537', plantid),
      # Add note
      Note = "tag changed to 7537 in 2022",
      edited = TRUE
    )
)

seed = seed %>%
  mutate(
    plantid = ifelse(tag %in% 3880 & plot %in% 13, gsub('3880', '7537', plantid), plantid),
    notes   = ifelse(tag %in% 3880 & plot %in% 13, '[tag manually edited; was 3880]', notes),
    tag     = ifelse(tag %in% 3880 & plot %in% 13, 7537, tag)
  )

# 14 2022      3363_14
demo %>% filter(Tag %in% 3363)
demo %>% filter(Tag %in% 3367, Plot %in% 14) # four umbels in 2022
seed %>% filter(tag %in% c(3363, 3367))

seed = seed %>%
  mutate(
    plantid = ifelse(tag %in% 3363 & plot %in% 14 & year %in% 2022, gsub('3363', '3367', plantid), plantid),
    notes   = ifelse(tag %in% 3363 & plot %in% 14 & year %in% 2022, '[tag manually edited; was 3363]', notes),
    tag     = ifelse(tag %in% 3363 & plot %in% 14 & year %in% 2022, 3367, tag)
  )

# 15 2022      3496_14
demo %>% filter(Tag %in% 3496)
seed %>% filter(tag %in% 3496)
demo.seed %>% filter(plot %in% 14, year %in% 2022, !in.seed, flowering.in.demo)
# it's gotta just be 3495 right?
demo %>% filter(Tag %in% 3495, Plot %in% 14)

seed = seed %>%
  mutate(
    plantid = ifelse(tag %in% 3496 & plot %in% 14 & year %in% 2022, gsub('3496', '3495', plantid), plantid),
    notes   = ifelse(tag %in% 3496 & plot %in% 14 & year %in% 2022, '[tag manually edited; was 3496]', notes),
    tag     = ifelse(tag %in% 3496 & plot %in% 14 & year %in% 2022, 3495, tag)
  )

# 16 2022      3076_15
demo %>% filter(Tag %in% 3076)
seed %>% filter(tag %in% 3076)
demo.seed %>% filter(plot %in% 15, year %in% 2022, !in.seed, flowering.in.demo)
demo %>% filter(Tag %in% 3706, Plot %in% 15)
seed %>% filter(tag %in% 3706)

seed = seed %>%
  mutate(
    plantid = ifelse(tag %in% 3076 & plot %in% 15 & year %in% 2022, gsub('3076', '3706', plantid), plantid),
    notes   = ifelse(tag %in% 3076 & plot %in% 15 & year %in% 2022, '[tag manually edited; was 3076]', notes),
    tag     = ifelse(tag %in% 3076 & plot %in% 15 & year %in% 2022, 3706, tag)
  )

# 17 2022      3813_15
demo %>% filter(Tag %in% 3813)
seed %>% filter(tag %in% 3813)
demo.seed %>% filter(plot %in% 15, year %in% 2022, !in.seed, flowering.in.demo)
demo %>% filter(Tag %in% 3831)
seed %>% filter(tag %in% 3831)

seed = seed %>%
  mutate(
    plantid = ifelse(tag %in% 3813 & plot %in% 15 & year %in% 2022, gsub('3813', '3831', plantid), plantid),
    notes   = ifelse(tag %in% 3813 & plot %in% 15 & year %in% 2022, '[tag manually edited; was 3813]', notes),
    tag     = ifelse(tag %in% 3813 & plot %in% 15 & year %in% 2022, 3831, tag)
  )

# 18 2023   3518_3_16J
demo %>% filter(Tag %in% 3518)
seed %>% filter(tag %in% 3518)
# data entry issue... hopefully fixed

# 19 2023    5022_5_9F
# should have been addressed above...
demo %>% filter(Tag %in% 5022)
seed %>% filter(tag %in% 5022)

# 20 2023   5570_6_18H
demo %>% filter(Tag %in% 5570)
seed %>% filter(tag %in% 5570)
demo.seed %>% filter(plot %in% 5, year %in% 2023, !in.seed, flowering.in.demo)
# dang that's a lot of plants...
# maybe there was just no phen on plants in xcoor 10-19?
demo.seed %>% filter(plot %in% 5, grepl('18', coor.demo), year %in% 2023, !in.seed, flowering.in.demo)
demo.seed %>% filter(plot %in% 5, No.umbels %in% 1, year %in% 2023, !in.seed, flowering.in.demo)

# hmm... might just be lost data (sad)

# 21 2023    5849_7_1E
demo %>% filter(Tag %in% 5849)
seed %>% filter(tag %in% 5849)
demo %>% filter(grepl('5849', Note))
demo.seed %>% filter(plot %in% 7, year %in% 2023, !in.seed, flowering.in.demo)
# also a ton of plants in here
# dude... what's up with that
# replaced tags...?
# anyway I think it's probably 5049 (0 got recorded as an 8?) (no - this has one umbel only)
# (yeah that's def 5049 - I checked the sheet)
# (the other two umbels are dead)

seed = seed %>%
  mutate(
    plantid = ifelse(tag %in% 5849 & plot %in% 15 & year %in% 2022, gsub('5849', '5049', plantid), plantid),
    notes   = ifelse(tag %in% 5849 & plot %in% 15 & year %in% 2022, '[tag manually edited; was 5849]', notes),
    tag     = ifelse(tag %in% 5849 & plot %in% 15 & year %in% 2022, 5049, tag)
  )

# 22 2023   7537_13_4C
demo %>% filter(Tag %in% 7537)
# think I fixed this above.

# 23 2023  3076_15_12H
# see above

seed = seed %>%
  mutate(
    plantid = ifelse(tag %in% 3076 & plot %in% 15 & year %in% 2023, gsub('3076', '3706', plantid), plantid),
    notes   = ifelse(tag %in% 3076 & plot %in% 15 & year %in% 2023, '[tag manually edited; was 3076]', notes),
    tag     = ifelse(tag %in% 3076 & plot %in% 15 & year %in% 2023, 3706, tag)
  )

# 24 2023  3813_15_14C

seed = seed %>%
  mutate(
    plantid = ifelse(tag %in% 3813 & plot %in% 15 & year %in% 2023, gsub('3813', '3831', plantid), plantid),
    notes   = ifelse(tag %in% 3813 & plot %in% 15 & year %in% 2023, '[tag manually edited; was 3813]', notes),
    tag     = ifelse(tag %in% 3813 & plot %in% 15 & year %in% 2023, 3831, tag)
  )

##### 
# Re-try the merge, but with cleaned datasets

demo.seed = merge(
  x = seed %>% mutate(in.seed = TRUE),
  y = demo %>%
    filter(Year > 2020) %>%
    mutate(flowering.in.demo = !is.na(No.umbels) & No.umbels > 0) %>%
    mutate(coor.demo = paste0(Xcoor, Ycoor)) %>%
    select(-c(Xcoor, Ycoor, YrTag)) %>%
    mutate(in.demo = TRUE),
  by.x = c("year", "plot", "tag"), by.y = c("Year", "Plot", "Tag"),
  all.x = TRUE, all.y = TRUE, suffixes = c('.seed', '.demo')
) %>%
  mutate(across(starts_with('in'), function(x) ifelse(is.na(x), FALSE, x)))

head(demo.seed)
with(demo.seed, table(in.demo, in.seed, year))
# Some issues, otherwise okay though.
# (merge with multiple umbels?)

##### 
# (Merge again, but this time without giving us all demo records)

demo.seed = merge(
  x = seed,
  y = demo %>%
    filter(Year > 2020) %>%
    mutate(flowering.in.demo = !is.na(No.umbels) & No.umbels > 0) %>%
    mutate(coor.demo = paste0(Xcoor, Ycoor)) %>%
    select(-c(Xcoor, Ycoor, YrTag)),
  by.x = c("year", "plot", "tag"), by.y = c("Year", "Plot", "Tag"),
  all.x = TRUE, all.y = FALSE, suffixes = c('.seed', '.demo')
)

head(demo.seed)
nrow(demo.seed)

# write.csv(
#   demo.seed,
#   '01_data_cleaning/out/demo_seed_v1.csv',
#   row.names = FALSE
# )

# (ah shoot... also just realized that changes to demo df above would influence
# demo.mumb merge, potentially)

# Write some kind of temp demo data frame too...
# write.csv(
#   demo,
#   '01_data_cleaning/out/demo_all_clean_v2.csv',
#   row.names = FALSE
# )

#######################
##### Demo for survival
#######################

##### Look at NAs in the dataset

apply(demo, 2, function(x) sum(is.na(x)))
# NAs in leaves, stalk height, leaf lengths, number of umbels, umbel diameters...
# some of these are due to non-flowering
# I'd also bet that some of these are due to mortality...

### What's up with the leaf counts?

demo %>% filter(is.na(No.leaves))
# just from these few cases, can see there are some cases where there are umbels listed (argh)

demo %>% filter(is.na(No.leaves)) %>%
  group_by(Year) %>% summarise(n = n())
# occurring throughout, but a ton in 2020

demo %>%
  group_by(leaf.na = is.na(No.leaves), Year) %>% 
  summarise(n = n()) %>%
  ungroup() %>%
  pivot_wider(names_from = leaf.na, values_from = n)
# hmm okay so there definitely were leaf counts being recorded in 2020!

demo %>% filter(Year %in% 2020) # yeah most of these NA records are just NA everywhere...

demo %>%
  filter(Year %in% 2020) %>%
  group_by(Plot, leafno.na = is.na(No.leaves)) %>%
  summarise(n = n()) %>%
  pivot_wider(names_from = leafno.na, values_from = n)
# okay - going to guess these were simply not sampled in 2020
# (oh I do seem to remember there being a note about this somewhere...)

# What about NAs in other years? Are they dead? Or no-checks?

demo %>%
  filter(!Year %in% 2020) %>%
  filter(is.na(No.leaves))
# yeah... not really sure what to do about these

# How many cases are there where there are NAs for leaf size but other data is provided?
demo %>%
  filter(is.na(No.leaves) & (!is.na(Stalk_Height) | !is.na(Leaf.length) | !is.na(No.umbels)))
# well... only six cases... assume they are alive then

### How often do the NAs later have a non-NA record?

demo %>%
  group_by(plantid) %>%
  filter(any(is.na(No.leaves)) & any(!is.na(No.leaves))) %>%
  mutate(min.na = min(Year[is.na(No.leaves)]), max.no.na = max(Year[!is.na(No.leaves) & No.leaves > 0])) %>%
  group_by(is.bad = max.no.na > min.na) %>%
  summarise(n = n())
# oof... lots of baddies it seems  

demo %>%
  group_by(plantid) %>%
  filter(any(is.na(No.leaves)) & any(!is.na(No.leaves))) %>%
  mutate(min.na = min(Year[is.na(No.leaves)]), max.no.na = max(Year[!is.na(No.leaves) & No.leaves > 0])) %>%
  filter(max.no.na > min.na) %>%
  arrange(plantid, Year)
# suggests to me that NAs should just be treated as no-checks

### Okay... how often do we have no-leaf plants appear in a later year?

demo %>%
  group_by(plantid) %>%
  filter(any(!No.leaves & !is.na(No.leaves))) %>%
  filter(Year > min(Year[!No.leaves])) %>%
  mutate(dead.status = any(No.leaves > 0)) %>%
  distinct(plantid, .keep_all = TRUE) %>%
  group_by(dead.status) %>%
  summarise(n = n())
# okay... 120 of 183 of plants ever thought dead were misses!
# (actually might be higher... code above does not account for 2023-only deaths)
# (although it also doesn't account for multiple non-consec zeros)

# (survival should be done with occupancy modeling, methinks)

##### Make a demo table just to see what it would look like

demo.table = demo %>%
  mutate(
    detected = (!is.na(No.leaves) | !is.na(Stalk_Height | !is.na(Leaf.length) | !is.na(No.umbels))),
    obs.alive = detected & (No.leaves > 0),
    no.leaves = ifelse(obs.alive, No.leaves, NA),
    leaf.leng = ifelse(obs.alive, Leaf.length, NA),
    flowering = case_when(
      !detected  ~ NA,
      !obs.alive ~ FALSE,
      No.umbels > 0 ~ TRUE,
      !No.umbels ~ FALSE,
      is.na(No.umbels) ~ FALSE, # assume it would have been recorded if it was flowering
      .default = NA
    ),
    no.umbels = ifelse(No.umbels > 0, No.umbels, NA),
  ) %>%
  select(plantid, Plot, Year, detected, obs.alive, no.leaves, leaf.leng, flowering, no.umbels) %>%
  arrange(Year, Plot, plantid)

# maybe export this for now...

# Write some kind of temp demo data frame too...
write.csv(
  demo.table,
  '01_data_cleaning/out/demo_table.csv',
  row.names = FALSE
)
