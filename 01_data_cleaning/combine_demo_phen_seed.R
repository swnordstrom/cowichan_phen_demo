# Script for reconciling the demography plantids with the seed/phenology plantids.
# Demo: 2016-2024, Phen/seed: 2021-2024
# SN - init 29 June 2024

# TO DO:
# Go back into demo and scrape out notes for which tags were wrong in phen

library(ggplot2)
library(dplyr)
library(tidyr)

rm(list = ls())

######################################################################
##### Load in processed data
######################################################################

# All demo data
demo = read.csv('01_data_cleaning/out/demo_all_cleaned.csv', encoding = 'latin1')

# Seed and phen data
seedphen = read.csv('01_data_cleaning/out/seed_phen_combined_cleaned_all.csv')

######################################################################
##### Assess overlap
######################################################################

overlap.plantids = merge(
  x = demo %>% 
    filter(Year > 2020) %>% 
    select(Year, plantid, Plot, Xcoor, Ycoor, No.leaves, No.umbels, demo.note, proc.note) %>%
    mutate(in.demo = TRUE),
  y = seedphen %>% distinct(year, plantid, plot, in.seed, in.phen),
  by.x = c('Year', 'plantid'), by.y = c('year', 'plantid'), all.x = TRUE, all.y = TRUE
) %>%
  mutate(across(starts_with('in.'), function(x) ifelse(is.na(x), FALSE, x)))

head(overlap.plantids)

overlap.plantids %>%
  group_by(Year, in.demo, in.repr = in.seed | in.phen) %>%
  summarise(n = n()) %>%
  pivot_wider(names_from = Year, values_from = n)

to.fix = overlap.plantids %>%
  # Give me plants that are either
  # - in demo with leaf counts+umbels and missing from phen or seed, or
  # - not in demo
  filter(((in.demo & No.leaves > 0 & No.umbels > 0) & (!in.phen | !in.seed)) | !in.demo) %>%
  arrange(plantid, Year) %>%
  # Filter out plants from 2021 that are missing from seed
  filter(!(Year %in% 2021 & !in.seed)) %>%
  # Filter out plants in demo+phen that were only partially sampled in 2023
  filter(!(Year %in% 2023 & Plot %in% c(1:2, 5:6, 13:14) & Xcoor > 9)) %>%
  filter(!(Year %in% 2023 & Plot %in% 15 & Xcoor < 10)) %>%
  filter(
    !(Year %in% 2023 & Plot %in% 7 & Xcoor > 9) | (plot %in% 7 & Xcoor %in% 15 & Ycoor %in% 'B')
  )

######################################################################
##### Start looking through the list generated above and reconciling
######################################################################

# --- 5022/2205

demo %>% filter(Tag %in% 2055) # lol...
seedphen %>% filter(grepl('5022', plantid))

demo = demo %>%
  mutate(
    edited = ifelse(Tag %in% 2055, TRUE, edited),
    proc.note = ifelse(
      Tag %in% 2055,
      paste0('tag mis-entered as 2055;', proc.note),
      proc.note
    ),
    plantid = gsub('2055', '5022', plantid),
    Tag = ifelse(Tag %in% 2055, 5022, Tag)
  )

# --- 3076

demo %>% filter(Tag %in% c(3076, 3706))

demo = demo %>%
  mutate(
    edited = ifelse(Tag %in% 3706 & Plot %in% 15, TRUE, edited),
    proc.note = ifelse(
      Tag %in% 3706 & Plot %in% 15,
      paste0('tag mis-entered as 3706;', proc.note),
      proc.note
    ),
  plantid = gsub('3706\\_15', '3076_15', plantid),
  Tag = ifelse(Tag %in% 3706 & Plot %in% 15, 3076, Tag)
)

('3076_15' %in% demo$plantid & '3076_15' %in% seedphen$plantid)

# --- 3086

to.fix %>% filter(grepl('\\_14$', plantid), Year %in% 2023)
# Can't find any unexplained plants missing from seed/phen in plot 14 in 2023

# --- 3096

demo %>% filter(Tag %in% 3096)
# demo note in 2024 says maybe outside plot
overlap.plantids %>% filter(plot %in% 13, Year %in% 2023, !in.demo)
# nobody here missing from demo

# --- 3103 (2023)

demo %>% filter(Tag %in% 3103) # ah... in diversity...
overlap.plantids %>% filter(plot %in% 4, Year %in% 2023, !in.demo)
# looks like everyone from plot 14 is accounted for in demo
# assume it was just not included because it was outside the plot

# --- 3108 (2024)

# It's in demo and in seed, just not in phen...
# But it doesn't look like it was excluded in the phen cleaning script...
# on is this one of the weird ones where I had a note about lumping and splitting?
# CHECK RAW DATA

# --- 3125

# We have 3125_5 in phen/seed in 2022, 2023 (maybe elsewhere!)
# Check demo to see which plants were flowering
demo %>% filter(Tag %in% 3125, Year %in% 2022:2023) %>% arrange(Xcoor)

# 2022: 14C was not found, 3I was flowering
# 2023: 14C was flowering but did not in full survey, 3I was in survey
# 2024: plant 14C was recorded in demo as 15D

seedphen = seedphen %>%
  mutate(
    plantid = ifelse(plantid %in% '3125_5' & year %in% 2022:2023, '3125_5.3I', plantid),
    plantid = gsub('3125\\_5\\.15D', '3125_5.14C', plantid)
  )

# Should return false:
('3125_5' %in% seedphen$plantid | '3125_5' %in% demo$plantid)
# Should return true:
('3125_5.14C' %in% seedphen$plantid & '3125_5.3I' %in% seedphen$plantid)
grep('3125', seedphen$plantid, value = TRUE) %>% unique()
# Nice

# --- 3135_15, 2024

# In seed/phen, not in demo
demo %>% filter(Tag %in% 3135)
# in plots 5 and 12, not 15
overlap.plantids %>% filter(Year %in% 2024, Plot %in% 15, in.demo, !in.seed | !in.phen, No.leaves > 0)

# Really not sure what it could be
# Phen+seed were at coords 14J,
# but there's nothing near that coordinate that's missing from phen/seed...

# --- 3143

# One of the weird diversity plot plants from plot 9
seedphen %>% filter(grepl('3143', plantid))
# It's here in 2021 but not 2022-2024
# Probably just missed because it's in the diversity plot.

# In fact, there are *no* plants in plot 9 in 2022-2024

# --- 3164_13

# Demo record says tag is actually 3614
demo %>% filter(Tag %in% 3164, Plot %in% 13)
seedphen %>% filter(grepl('3164\\_13', plantid)) # good - tag 3164 doesn't appear in phen/seed
seedphen %>% filter(grepl('3614\\_13', plantid))

demo = demo %>% 
    mutate(
    edited = ifelse(Tag %in% 3164 & Plot %in% 13, TRUE, edited),
    proc.note = ifelse(
      Tag %in% 3164 & Plot %in% 13,
      paste0('tag mis-entered as 3164;', proc.note),
      proc.note
  ),
  plantid = gsub('3164\\_13', '3614_13', plantid),
  Tag = ifelse(Tag %in% 3164 & Plot %in% 13, 3614, Tag)
)

('3614_13' %in% demo$plantid & '3614_13' %in% seedphen$plantid) # good

# --- 3164_6

# Demo record here says 3164 is old tag, new tag is 3789
demo %>% filter(grepl('3789', demo.note))
seedphen %>% filter(grepl('3789\\_6', plantid))
# Cool
demo = demo %>% 
  mutate(
    edited = ifelse(Tag %in% 3164 & Plot %in% 6, TRUE, edited),
    proc.note = ifelse(
      Tag %in% 3164 & Plot %in% 6,
      paste0('tag switched (old tag 3164);', proc.note),
      proc.note
    ),
    plantid = gsub('3164\\_6', '3789_6', plantid),
    Tag = ifelse(Tag %in% 3164 & Plot %in% 6, 3789, Tag)
  )

('3789_6' %in% demo$plantid & '3789_6' %in% seedphen$plantid) # good

# --- 3185_2

demo %>% filter(Tag %in% 3185, Year %in% 2021:2023) %>% arrange(Xcoor)
# All of these are for the 19B plant

seedphen = seedphen %>% mutate(plantid = gsub('3185\\_2', '3185_2.19B', plantid))

# should be false
('3185_2' %in% seedphen$plantid | '3185_2' %in% demo$plantid)
# should be true
('3185_2.19B' %in% seedphen$plantid)

# --- 3196_13, 2022

# in demo, in seed, not in phen
seedphen %>% filter(grepl('3196', plantid), year %in% 2022)
# excluded from phen data that year due to missing records

# --- 3296/3298

demo %>% filter(Tag %in% 3296) # no other 3296 in here
demo %>% filter(Tag %in% 3298)

# Change in demo
demo = demo %>% 
  mutate(
    edited = ifelse(Tag %in% 3296 & Plot %in% 7, TRUE, edited),
    proc.note = ifelse(
      Tag %in% 3295 & Plot %in% 7,
      paste0('tag misentered as 3296;', proc.note),
      proc.note
    ),
    plantid = gsub('3296\\_7', '3298_7', plantid),
    Tag = ifelse(Tag %in% 3296 & Plot %in% 7, 3298, Tag)
  )

('3298_7' %in% demo$plantid & '3298_7' %in% seedphen$plantid)

# --- 3321_4

demo %>% filter(Tag %in% 3321)
# ah - in diveristy plot...

# --- 3360_14, 2023

demo %>% filter(Tag %in% 3360, Plot %in% 14)
# In seed but not in phen.
overlap.plantids %>% filter(plot %in% 14, in.phen, !in.seed)
# I don't think it's 3463.
# Oh well. Let it be so.

# --- 3360_6, 2024

demo %>% filter(Tag %in% 3360, Plot %in% 6)
# No note in 2024
overlap.plantids %>% filter(plot %in% 6, !in.demo, in.seed | in.phen, Year %in% 2024)
# I don't think it's either of these.
# Maybe just missed in phen?

# --- 3361/5044

demo %>% filter(Tag %in% c(3361, 5044), Plot %in% 5)
seedphen %>% filter(grepl('3361', plantid))
seedphen %>% filter(grepl('5044', plantid))

seedphen = seedphen %>% mutate(plantid = gsub('5044\\_5', '3361_5', plantid))

('3361_5' %in% seedphen$plantid & '3361_5' %in% demo$plantid)
# Below should return false
('5044_5' %in% seedphen$plantid | '5044_5' %in% demo$plantid)

# --- 3362_14, 2023

demo %>% filter(Tag %in% 3362, Plot %in% 14)
seedphen %>% filter(grepl('3362\\_14', plantid))
# nowhere in phen
overlap.plantids %>% filter(plot %in% 14, Year %in% 2023, (in.phen | in.seed) & !in.demo)
# No unaccounted for plants in 14.

# --- 3363

demo %>% filter(Tag %in% 3363)
# These are all in plot 5...
seedphen %>% filter(grepl('3363', plantid))
# 3362 in these years is not flowering or even present...
# 2021 corods in phen: 15H
demo %>% filter(Plot %in% 14, Xcoor %in% 14:16, Ycoor %in% c('G', 'H', 'I'), Year %in% 2021:2022)
# Ah. Demo notes say tag in demo is 3367...
demo %>% filter(Tag %in% 3367, Plot %in% 14) # lmao

demo = demo %>% mutate(
  edited = ifelse(Tag %in% 3367 & Plot %in% 14, TRUE, edited),
  proc.note = ifelse(
    Tag %in% 3367 & Plot %in% 14,
    paste0('tag misentered as 3367;', proc.note),
    proc.note
  ),
  plantid = gsub('3367\\_14', '3363_14', plantid),
  Tag = ifelse(Tag %in% 3367 & Plot %in% 14, 3363, Tag)
)
  

('3363_14' %in% demo$plantid & '3363_14' %in% seedphen$plantid)

# --- 3366_14 in 2023

demo %>% filter(Tag %in% 3366, Plot %in% 14)
seedphen %>% filter(grepl('3366', plantid)) # it's here in 2022
overlap.plantids %>% filter(Plot %in% 14, Year %in% 2023, !in.demo & (in.seed | in.phen))
# no dice... not sure we have this one in data
# (maybe it's in the plot that was insufficiently sampled?)

# --- 3390_13 in 2024

# Note in demo says 'was 3391 in phen'
# So, change seed/phen

seedphen %>% filter(year %in% 2024, grepl('3391\\_13', plantid))
demo %>% filter(Year %in% 2024, plantid %in% '3390_13')

seedphen = seedphen %>% mutate(plantid = gsub('3391\\_13', '3390_13', plantid))

('3390_13' %in% demo$plantid & '3390_13' %in% seedphen$plantid)

# --- 3393_15 in 2024

# Plot 15, 0J, demo comment for 3393 says 'GTP in phen'

seedphen %>% filter(grepl('G[TS]P\\_15', plantid)) # hmm...
overlap.plantids %>% filter(plot %in% 15, Year %in% 2024)
demo %>% filter(Year %in% 2024, Plot %in% 15, Xcoor %in% 0, Ycoor %in% 'J')

# Okay, I am 99.9% sure this is BSP at 0J - GTP in notes must have been a mistake

seedphen %>% filter(grepl('B[ST]P\\_15', plantid))
# Oh... but is it BSP or BTP?
# (I checked - BSP is 0J, BTP is 12G. Sloppy on my part!)

seedphen = seedphen %>% mutate(plantid = gsub('BSP\\_15', '3393_15', plantid))

('3393_15' %in% demo$plantid & '3393_15' %in% seedphen$plantid)
# Should be false below:
('BSP_15' %in% seedphen$plantid | 'BSP_15' %in% seedphen$plantid) # good

# --- 3404_15 in 2024

# Note in demo says this was recorded as 3402 in phen
# Let's just look at the records for 3402/3404 in phen
demo %>% filter(Tag %in% c(3402, 3404), Year %in% 2024, Plot %in% 15)
# Great - 3402 was vegetative, 3404 was flowering. Good stuff!

seedphen %>% filter(grepl('3402\\_15', plantid), year %in% 2024)
# Sweet

seedphen = seedphen %>% 
  mutate(plantid = ifelse(plantid %in% '3402_15' & year %in% 2024, '3404_15', plantid))

# --- 3423_12 in 2024

# Note in demo record for 3423 says it was recorded as 3422 in phen
demo %>% filter(Tag %in% 3422:3423, Plot %in% 12, Year %in% 2024)
# Good - 3422 is actually vegetative
seedphen %>% filter(grepl('342[23]\\_12', plantid), year %in% 2024)

seedphen = seedphen %>%
  mutate(plantid = ifelse(plantid %in% '3422_12' & year %in% 2024, '3423_12', plantid))

# --- 3426 and 3427 in plot 15

# 3426 is in demo and in seed in 2022 but not in phen
# (I checked - it was excluded)
# 3427_15 in 2024 - in demo, not in phen or seed
demo %>% filter(Tag %in% 3427, Plot %in% 15)
seedphen %>% filter(grepl('3427\\_15', plantid)) # here in 2021, 2022
seedphen %>% filter(grepl('342[0-9]\\_15', plantid), year %in% 2024)
# we have seed counts+phen for 3423, 3426
demo %>% filter(Tag %in% 3426, Year %in% 2024, Plot %in% 15)
overlap.plantids %>% filter(plot %in% 15, Year %in% 2024, !in.demo)
# I don't think it's any of these.
# Assume it was just missed in phen?

# --- 3456_6

demo %>% filter(Tag %in% 3456)
# ah - probably in phen/seed as 7587
seedphen %>% filter(grepl('7587', plantid))
# Yup.

seedphen = seedphen %>% mutate(plantid = gsub('7587\\_6', '3456_6', plantid))

('3456_6' %in% seedphen$plantid & '3456_6' %in% demo$plantid)
# Sweet

# --- 3477_2 in 2023

# In demo and phen, but not in seed (weird)
seedphen %>% filter(grepl('3477', plantid))
# no 3477 in raw data...
# hmm... it's listed as just being dead - stalk broke off week after being dead
# Seems fair to call this zero umbels...

# Edit it into phen
seedphen = seedphen %>%
  mutate(
    umbel.no = ifelse(plantid %in% '3477_2' & year %in% 2023, 1, umbel.no),
    no.seeds = ifelse(plantid %in% '3477_2' & year %in% 2023, 0, no.seeds)
  )

# --- 3477_7 (2024)

# In seed/phen, not in demo
demo %>% filter(Tag %in% 3477, Plot %in% 7)
# there is no plant 3477 in plot 7...
overlap.plantids %>% filter(Plot %in% 7, in.demo, !in.seed | !in.phen, Year %in% 2024)

# Oh wow. Plant 3477 is just plant 3447, same exact data. Lol.
seedphen = seedphen %>% filter(!plantid %in% '3477_7')

# --- 3480_14 (2024)

# in phen/seed, not in demo
seedphen %>% filter(plantid %in% '3480_14')
demo %>% filter(Tag %in% 3480)
# no 3480 in plot 14...
# raw phen data says plot 19E
overlap.plantids %>% 
  filter(Plot %in% 14, Year %in% 2024, !(in.phen & in.seed)) %>%
  arrange(desc(Xcoor))
# There's a 3840 in 17F... but, given how messed up the coords are here, that's probably it
seedphen %>% filter(plantid %in% '3840_14') # but, this plant isn't here...

seedphen = seedphen %>% mutate(plantid = gsub('3480\\_14', '3840_14', plantid))

# --- 3482_15

# Probably just missed due to being on plot edge.

# --- 3516_13 (2022)
# in phen only
seedphen %>% filter(plantid %in% '3516_13')
# there is no 3516 in plot 13 (according to demo), only in plot 5
# so this is a tag misread?
# There's only one record for it in phen, and it has coordinate 0A...
overlap.plantids %>% 
  filter(Plot %in% 13, Year %in% 2022, !(in.phen & in.demo)) %>%
  arrange(Xcoor)
# There is no 0A that flowered on this day.

# Going to say this is just a bad record.
# Going to remove this record.
seedphen = seedphen %>% filter(!plantid %in% '3516_13')

# --- 3519_2 (2023)
# in demo and phen but not seed?
seedphen %>% filter(plantid %in% '3519_2')
seedphen %>% filter(year %in% 2023, plot %in% 2, in.seed & !in.phen)
demo %>% filter(Tag %in% 3519, Year %in% 2023)

# Yeah it's just not in phen...
# I guess just add the zero in to phen. It's one dead umbel.
seedphen = seedphen %>%
  mutate(
    umbel.no = ifelse(plantid %in% '3519_2' & year %in% 2023, 1, umbel.no),
    no.seeds = ifelse(plantid %in% '3519_2' & year %in% 2023, 0, no.seeds)
  )

# --- 3671 in 2024
# demo note says that the plant was recorded as tag 3926 in phen

seedphen %>% filter(grepl('3926', plantid))
demo %>% filter(Year %in% 2024, Tag %in% c(3926, 3671))
# Sweet

seedphen = seedphen %>% 
  mutate(
    plantid = ifelse(grepl('3926\\_14', plantid) & year %in% 2024, '3671_14', plantid)
  )

# --- 3687_1 (2024)

demo %>% filter(Tag %in% 3687, Year %in% 2024)
# Hmm... 5I, may be missed because on plot boundary
overlap.plantids %>% filter(Year %in% 2024, plot %in% 1, !in.demo) # nobody missing
# Going to assume this was just missed in phen.

# --- 3599_1 (2024)
# probably missed due to being on edge/outside plot (coord in demo is 0I)
demo %>% filter(Tag %in% 3599)

# --- 3748_6 (2021)
# in phen/seed but not demo
seedphen %>% filter(grepl('3748', plantid))
demo %>% filter(Tag %in% 3748, Plot %in% 6)
# lol there is no plant with this tag...
demo %>% filter(grepl('3748', demo.note))
# new tag for this plant is 3179
seedphen %>% filter(grepl('3179', plantid)) # cool
demo %>% filter(Tag %in% c(3179, 3748), Plot %in% 6)

demo = demo %>%
  mutate(
    edited = ifelse(Tag %in% 3179 & Plot %in% 6, TRUE, edited),
    proc.note = ifelse(
      Tag %in% 3179 & Plot %in% 6,
      'updating tag (old tag 3179 new tag 3748',
      proc.note
    ),
    plantid = gsub('3179\\_6', '3748_6', plantid),
    Tag = ifelse(Tag %in% 3179 & Plot %in% 6, 3748, Tag)
  )

('3748_6' %in% demo$plantid & '3748_6' %in% seedphen$plantid)

# --- 3758_12 (2024)

# Probably missed from phen (19H - might have been outside plot)
demo %>% filter(Tag %in% 3758)
overlap.plantids %>% filter(Year %in% 2024, plot %in% 12, (in.phen | in.seed), !in.demo)
# Nothing unaccounted for in phen
# I think I just missed it in phen.

# --- 3782_1 (2022)
# In phen+seed, not in demo
demo %>% filter(Tag %in% 3782)
# comment in 2023: "not new, separate from 3814"
# so was this 3814 in 2022?
demo %>% filter(Tag %in% c(3814, 3782), Plot %in% 1) %>% arrange(Year, Tag)

# Bleh... split tags but they all have messy records...
# I think in this case it's probably best to just get rid of all records to avoid ambiguity?

# Okay, course of action here:
# - New tag (3782) was placed in 2020
# - Keep tag 3814 for 2017-2018 records, 2023-2024 records
# - Leave proc note saying the plant did not die? (for imputing)
# - Change 2020-2022 records with tag 3814 to 3782 (will catch this particular merge issue)

# demo = demo %>%
#   mutate(
#     edited = ifelse(
#       Tag %in% 3814 & Plot %in% 1 & Year %in% 2020:2022,
#       TRUE,
#       edited
#     ),
#     proc.note = ifelse(
#       Tag %in% 3814 & Plot %in% 1 & Year %in% 2020:2022,
#       paste0(proc.note, '; tag 3782 added in 2020 (old tag: 3814, still in field)'),
#       proc.note
#     ),
#     proc.note = ifelse(
#       Tag %in% 3814 & Plot %in% 1 & Year %in% 2019,
#       paste0(proc.note, 'plant did NOT DIE here - record confusion plant present in 2024'),
#       proc.note
#     ),
#     plantid = ifelse(
#       Tag %in% 3814 & Plot %in% 1 & Year %in% 2020:2022,
#       '3782_1',
#       plantid
#     ),
#     Tag = ifelse(
#       Tag %in% 3814 & Plot %in% 1 & Year %in% 2020:2022,
#       3782,
#       Tag
#     )
#   )

# Check to see what the data for these plants in seed/phen is
seedphen %>% filter(plantid %in% c('3782_1', '3814_1')) %>% arrange(plantid, year)

# 2024: demo note says that the phen records for 3872 are actually for plant 3814
# 2022: 3782 record here can stay as is after making change to demo above
# but... there are also 3814 seed/phen records in 2022?

# okay... we have a seed/phen for both 3814 and 3782 in 2022, but only a demo record for one of them...
# ugh. okay I guess one of these will just not have demo

# --- 3798/3796

### NEED TO REVISE INDIVIDUAL PHEN to split records

# --- 3813/3831

# Notes in demo say tag is actually 3813

demo = demo %>%
  mutate(
    edited = ifelse(Tag %in% 3831 & Plot %in% 15, TRUE, edited),
    proc.note = ifelse(
      Tag %in% 3831 & Plot %in% 15,
      paste0(proc.note, '; tag misentered as 3831'),
      proc.note
    ),
    plantid = gsub('3831\\_15', '3813_15', plantid),
    Tag = ifelse(Tag %in% 3831 & Plot %in% 15, 3813, Tag)
  )

('3813_15' %in% demo$plantid & '3813_15' %in% seedphen$plantid)
# Below should return FALSE
('3831_15' %in% demo$plantid | '3831_15' %in% seedphen$plantid) # good

# --- 3836_5 (2023)
# in demo, not in seed or phen

demo %>% filter(Tag %in% 3836)
overlap.plantids %>% filter(Year %in% 2023, plot %in% 5, !in.demo)
# it's neither of these (3125 or 5022)

# Oh well, assume it was just missed.

# --- 3866/3886

demo %>% filter(Tag %in% 3866)
# multiple notes in here say tag is actually 3886
seedphen %>% filter(grepl('38[68]6', plantid))
# nearly all of these are 3886

# change any instance in either dataset of 3866 to 3886
demo = demo %>%
  mutate(
    edited = ifelse(Tag %in% 3866 & Plot %in% 7, TRUE, edited),
    proc.note = ifelse(
      Tag %in% 3866 & Plot %in% 7,
      paste0(proc.note, '; tag misentered as 3866'),
      proc.note
    ),
    plantid = gsub('3866\\_7', '3886_7', plantid),
    Tag = ifelse(Tag %in% 3866 & Plot %in% 7, 3886, Tag)
  )

seedphen = seedphen %>% mutate(plantid = gsub('3866\\_7', '3886_7', plantid))

('3886_7' %in% demo$plantid & '3886_7' %in% seedphen$plantid)
# below should return false
('3866_7' %in% demo$plantid | '3866_7' %in% seedphen$plantid) # good

# --- 3868_15 (2024)

# wasn't tracked in 2023 (according to demo note)
# coord is 15A in demo so may have been missed in phen?
overlap.plantids %>% filter(plot %in% 15, !in.demo, Year %in% 2024)
# I don't think it's one of these (iut's not one of the picks)
# assume it was missed

# --- 3882 (2024)
# in demo, in phen, not in seed (?)
seedphen %>% filter(grepl('3882', plantid))
# huh. Note in phen dataset says can't find tag.
# Not good data entry on my part.
# Assume it's zero seeds (plant senesced)
seedphen = seedphen %>%
  mutate(
    umbel.no = ifelse(plantid %in% '3882_15' & year %in% 2024, 1, umbel.no),
    no.seeds = ifelse(plantid %in% '3882_15' & year %in% 2024, 0, no.seeds)
  )

# --- 3885_14 (2023)

demo %>% filter(Tag %in% 3885)
seedphen %>% filter(grepl('3885', plantid)) # here in 2022
overlap.plantids %>% filter(plot %in% 14, Year %in% 2023, !in.demo)
# yep... no unaccounted for plants in plot 14
# assume it was missed

# --- 3929
# note in demo says was recorded in phen as 7559
# what does the demo for 7559 look like? would likely have been missed then
demo %>% filter(Tag %in% 7559, Year %in% 2024)
# huh

# Anyway, change tag in phen to 3929
seedphen = seedphen %>%
  mutate(plantid = ifelse(plantid %in% '7559_10' & year %in% 2024, '3929_10', plantid))

# --- 5049_7 (2023)
# in demo, not phen or seed, coord is 2D
seedphen %>% filter(grepl('5049\\_7', plantid), year %in% 2023)
overlap.plantids %>% filter(Year %in% 2023, plot %in% 7, !in.demo)
# aha! 5849
# glad this finally worked for once
seedphen %>% filter(grepl('5849', plantid))
seedphen = seedphen %>% mutate(plantid = gsub('5849', '5049', plantid))

('5049_7' %in% seedphen$plantid | '5049_7' %in% demo$plantid)

# --- 5570
# I looked at this in the phen/seed combining script. It's a Primula tag.
# Coordinate is 18H
# Are t here any demo records at 18H for flowering plants?

# Actually looking at the phen records, this one should probably be excluded

# --- 7521
# 2024: note says missed in phen
# but 2023, it's in seed but not in phen (?)
overlap.plantids %>% filter(plot %in% 15, in.phen, !in.seed, Year %in% 2023)
# no unaccounted for plants (present in phen but not seed)

# (this has been fixed - coordinate was mis-recorded as 1 instead of 15 so it got left out)

# --- 7545

# missed in phen?
# OH... all three umbels were recorded under 7591 for phen/seed
# but, demo has 7591 having two umbels and 7545 having one
# so one of the umbelsassigned to 7591 for phen is actually 7545
# not sure this can be fixed without good notes

# --- 7546 (2024)
# Shoot... I just completely missed this when counting seeds
# Rats!

# --- 7587, but only in 2023?
seedphen %>% filter(grepl('3456', plantid))
# okay - think it's been fixed (was excluded in 2023 because of I misread the
# plot subsampling records)

# --- 7655 (phen/seed), 7665 (demo)
# Same plant? Both in plot 15
demo %>% filter(Tag %in% c(7655, 7665))
# Ah... there was also a 7665 at plot 10
# and there is no 7655 anywhere
# I bet this was a tag mis-entry (yes - data entry mistake - fixing now)
# Ignore the note about BSP/BTP - coordinates don't match

# --- 7681
# demo note says it was recorded as 7692 in phen
demo %>% filter(Tag %in% c(7681, 7692), Year %in% 2024)
seedphen %>% filter(grepl('7681', plantid) | grepl('7692', plantid), year %in% 2024)
# Easy fix in phen
seedphen = seedphen %>% 
  mutate(plantid = ifelse(plantid %in% '7692_10' & year %in% 2024, '7681_10', plantid))

('7681_10' %in% demo$plantid & '7681_10' %in% seedphen$plantid)

# --- BSP in plot 3 (2024)

demo %>% filter(Year %in% 2024, Plot %in% 3, grepl('B[TS]P', demo.note))
# no demo note
# coord is 4F
overlap.plantids %>% filter(Year %in% 2024, plot %in% 3, !in.demo)
# Probably forgot to give it a demo record.

# --- BSP in plot 5 (2024)
# (fixed - was a phen processing issue)

# --- BTP plot 15 (2024)
demo %>% filter(Year %in% 2024, Plot %in% 15, grepl('B[TS]P', demo.note))
# coordinate is 12G
overlap.plantids %>% 
  filter(Year %in% 2024, Plot %in% 15, !(in.phen & in.seed)) %>%
  arrange(Xcoor, Ycoor)
# Not seeing anything.
# Assume I forgot to tag.

# --- RSP plot 14 (2024)
demo %>% filter(Year %in% 2024, Plot %in% 14, grepl('R[ST]P', demo.note))
# Coord is 19E
overlap.plantids %>% filter(Year %in% 2024, Plot %in% 14, Xcoor > 17)
# Assuming I forgot to tag it (bummer)

# --- RSP in plot 15 (2024)
demo %>% filter(Year %in% 2024, Plot %in% 15, grepl('R[ST]P', demo.note))
seedphen %>% filter(grepl('3332', plantid))
# demo note says 3332 was RSP in phen, but seeed+phen 2024 has 3332 instead of rsp
# Plot says 2A, but also, NP by 5/16, so not around for demo

# --- R in plot 7
# records are messy here, going to say let's exclude this
# (I went back and excluded it in phen)

# --- WHTP (white toothpick) in 14
# note says no demot/tag, unsalvagable I think

# --- YTP 13 (2024)
demo %>% filter(grepl('[Yy][ts]p', demo.note))
# no demo note
# eaten, so probably not ready available for demo

# --- The YTP in 5 should probably be ignored (only one record)

######################################################################
##### Merge and export data
######################################################################

# I kind of want to try, for now, saving all data points as one giant df/csv

data.combined = merge(
  x = demo %>% 
    select(
      Year, Plot, No.leaves, Leaf.length, No.umbels, 
      demo.note, proc.note, demo.edited = edited, 
      plantid) %>%
    mutate(in.demo = TRUE),
  y = seedphen %>% select(-plot),
  by.x = c('Year', 'plantid'), by.y = c('year', 'plantid'),
  all.x = TRUE, all.y = TRUE
) %>%
  mutate(across(starts_with('in'), function(x) ifelse(is.na(x), FALSE, x))) %>%
  select(
    Year, plantid, Plot, in.demo, in.phen, in.seed,
    No.leaves, Leaf.length, No.umbels,
    umbel.no, no.seeds, phen.umbels, phen.julis,
    demo.note, proc.note, demo.edited, seed.note = seed.notes
  ) %>%
  arrange(Year, Plot, plantid)

head(data.combined)

write.csv(
  data.combined,
  na = '', row.names = FALSE,
  file = '01_data_cleaning/out/demo_phen_seed_combined.csv'
)
