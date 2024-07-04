# Script for reconciling phenology and seed set data
# This should be straightforward because these two datasets are often based on
# the same plants and stored on the same datasheets.
# SN - init 28 Jun 2024

### --------------------------------------------------------

library(ggplot2)
library(dplyr)
library(tidyr)

rm(list = ls())

### --- Read in data ------------

# Processed phenology data
phen = read.csv('01_data_cleaning/out/phenology_buds_deaths_cleaned.csv')

# Processed seed data
seed = read.csv('01_data_cleaning/out/seed_counts_all_cleaned.csv')

# Processed demo data
# NOT FOR MERGING, but instead for reference
demo = read.csv('01_data_cleaning/out/demo_all_cleaned.csv')

### --- Merge together and compare  ------------

phen.seed = merge(
  x = phen %>% distinct(year, plantid) %>% mutate(in.phen = TRUE),
  y = seed %>% distinct(year, plantid) %>% mutate(in.seed = TRUE),
  by = c('year', 'plantid'),
  all.x = TRUE, all.y = TRUE
) %>%
  mutate(across(c(in.phen, in.seed), function(x) ifelse(is.na(x), FALSE, x)))

head(phen.seed)

phen.seed %>%
  filter(!(in.seed & in.phen)) %>%
  group_by(in.phen, in.seed, year) %>%
  summarise(n = n()) %>%
  pivot_wider(names_from = year, values_from = n)

# in.phen in.seed `2021` `2022` `2023` `2024`
# <lgl>   <lgl>    <int>  <int>  <int>  <int>
# 1 FALSE   TRUE         1     12      7     12
# 2 TRUE    FALSE      166      6     48      6


### --- Start the reconciliation  ------------

### --- 2021 ---------

# Only one plant in 2021 in seed but not in phen
# It may be one of the ~160 plants missing from seed but not phen
# The remaining 2021 plants (in phen, not in seed) are fine - just not included
# in the sample

phen.seed %>% filter(year %in% 2021, in.seed, !in.phen)
# 3928... isn't this actually 3926 or something

phen.seed %>% filter(grepl('392[68]\\_6', plantid))
# hmm...
phen.seed %>% filter(grepl('392[0-9]\\_6', plantid))
phen.seed %>% filter(grepl('392[0-9]', plantid))
# Nope.

seed %>% filter(grepl('3928\\_6', plantid))
phen %>% filter(plot %in% 6, year %in% 2021)
# 3918?
seed %>% filter(grepl('3918', plantid))
# looks plausible

# What is in demo?
demo %>% filter(Plot %in% 6, grepl('39[12]8', plantid))
# only one record...? from 2016...?
demo %>% filter(Plot %in% 6, Xcoor %in% 7:9, Ycoor %in% c('B', 'C', 'D'))
# Yep. Tag was henceforth recorded as 3916.
# Coord in raw phen is the same.

# Change phen AND seed to 3916. This will match demo ID.
seed = seed %>%
  mutate(
    notes = ifelse(
      grepl('39[12][68]\\_6', plantid), 
      paste0(notes, ' fixed tag mistakenly enterd as 3928'), 
      notes
    ),
    tag = ifelse(grepl('39[12][68]\\_6', plantid), 3916, tag),
    plantid = ifelse(grepl('39[12][68]\\_6', plantid), '3916_6', plantid),
  )

phen = phen %>% mutate(plantid = ifelse(grepl('39[12][68]\\_6', plantid), '3916_6', plantid))

('3916_6' %in% seed$plantid & '3916_6' %in% phen$plantid)

### --- 2022 ---------

phen.seed %>% filter(year %in% 2022, in.phen, !in.seed) %>% arrange(plantid) %>% pull(plantid)
# "2238_4"  "3151_3"  "3374_6"  "3495_14" "3516_13" "7526_5" 

phen.seed %>% filter(year %in% 2022, !in.phen, in.seed) %>% arrange(plantid) %>% pull(plantid)
# "2230_4"  "3196_13" "3345_6"  "3354_6"  "3426_15" "3430_14" "3496_14" "3651_3"  "3755_15"
# "3848_15" "3863_15" "7536_5" 

# First one is easy:

# --- 2230/2238
# also not sure why there's a 2-tag at all...)
# what is it in demo?
demo %>% filter(grepl('223[08]', plantid))
# 2230 in demo unequivocally
phen %>% filter(grepl('2238', plantid))

# Change ID in phen
phen = phen %>% mutate(plantid = gsub('2238', '2230', plantid))

('2230_4' %in% phen$plantid & '2230_4' %in% seed$plantid)

# --- 3495/3496
phen %>% filter(grepl('349[56]\\_14', plantid))
seed %>% filter(grepl('349[56]\\_14', plantid))

# Tag is 3495, must have been fat-fingered in entry
seed = seed %>% mutate(plantid = gsub('3496\\_14', '3495_14', plantid))

('3495_14' %in% seed$plantid & '3495_14' %in% phen$plantid)

# --- 7526/7536

# Which is it in demo?
demo %>% filter(Tag %in% c(7536, 7526), Plot %in% 5)
# it's 7526 in demo

seed = seed %>% mutate(plantid = gsub('7536\\_5', '7526_5', plantid))

('7526_5' %in% seed$plantid & '7526_5' %in% phen$plantid)

# --- 3374/3354

seed %>% filter(plantid %in% '3354_6')
phen %>% filter(plantid %in% '3374_6')
demo %>% filter(Tag %in% c(3374, 3354), Plot %in% 6)
# It's 3374.

seed = seed %>% mutate(plantid = gsub('3354\\_6', '3374_6', plantid))

('3374_6' %in% seed$plantid & '3374_6' %in% phen$plantid)

# --- 3151/3651

demo %>% filter(grepl('3[16]51', Tag), Plot %in% 3)
# Tag in demo is 3651

# fix in phen
phen = phen %>% mutate(plantid = gsub('3151\\_3', '3651_3', plantid))

('3651_3' %in% seed$plantid & '3651_3' %in% phen$plantid)

# --- 3516_13... only outstanding plant from seed in plot 13 is 3196

phen %>% filter(grepl('3516\\_13', plantid)) # only one umbel
seed %>% filter(grepl('3196\\_13', plantid)) # two umbels...

demo %>% filter(Tag %in% c(3196, 3516), Plot %in% 13)
# uh... there is no plant 3516 in plot 13, according to demo...

demo %>% filter(grepl('3[0-9]16', Tag), Plot %in% 13)
# 3616?
phen %>% filter(grepl('3616\\_13', plantid))
# no... already here (and in seed..)
demo %>% filter(grepl('35[0-9]6', Tag), Plot %in% 13)
# nothing
demo %>% filter(grepl('351[0-9]', Tag), Plot %in% 13)
# 3515??
seed %>% filter(tag %in% 3515, year %in% 2022) # four umbels...
phen %>% filter(grepl('3515', plantid), year %in% 2022) # four umbels...
# no this is all accounted for...
seed %>% filter(tag %in% 3565)

demo %>% filter(Year %in% 2022, Tag %in% c(3515, 3616, 3565, 3196), Plot %in% 13)
# yugh... and this is the year where phen has no coords...

# Okay well I guess this one remains unresolved...

### --- 2023 ---------

phen.seed %>% filter(year %in% 2023, in.phen, !in.seed) %>% arrange(plantid) %>% pull(plantid)
# [1] "3010_5"    "3065_15"   "3067_5"    "3068_15"   "3127_5"    "3135_5"    "3144_5"    "3153_5"   
# [9] "3182_4"    "3187_5"    "3192_1"    "3193_1"    "3363_14"   "3417_15"   "3424_1"    "3432_1"   
# [17] "3463_14"   "3473_5"    "3476_5"    "3477_2"    "3488_5"    "3519_2"    "3526_14"   "3556_14"  
# [25] "3566_5"    "3659_15"   "3663_5"    "3685_5"    "3688_7"    "3697_5"    "3763_1"    "3807_5"   
# [33] "3833_5"    "3850_1"    "3869_5"    "3882_15"   "3888_5"    "3896_5"    "3902_1"    "5034_5"   
# [41] "5043_5"    "5770_6"    "7512_5.2G" "7512_5.3D" "7543_5"    "7546_5"    "7555_5"    "7572_15" 
# mother of god, what the hell happened here???
# ohhh... I think they only assessed seed set in half of plot 5 so that makes sense...

phen.seed %>% filter(year %in% 2023, !in.phen, in.seed) %>% arrange(plantid) %>% pull(plantid)
# [1] "3338_5"  "3454_6"  "3481_15" "3987_5"  "5570_6"  "7512_5"  "7519_5" 

# Okay, well, first fix this oopsie
# --- 7512/7519
# 7512_5.3D is really 7519, misentered in phen somehow

phen = phen %>%
  mutate(
    plantid = case_match(
      plantid,
      '7512_5.2G' ~ '7512_5',
      '7512_5.3D' ~ '7519_5',
      .default = plantid
    )
  ) # %>%
  # line for checking that this worked...
  # filter(grepl('751[29]\\_5', plantid))

# --- 5570/5770

demo %>% filter(Tag %in% c(5770, 5570)) # don't like that...
seed %>% filter(tag %in% 5570)
# original seed record says coord is 18H
demo %>% filter(Plot %in% 6, Xcoor %in% 17:19, Ycoor %in% c('G', 'H', 'I')) %>% filter(Year %in% 2023)

# Are each of these in seed?
seed %>% filter(year %in% 2023, tag %in% c(3374, 3643, 3767, 3955))
# yes... lol

phen %>%
  filter(plot %in% 6) %>%
  mutate(tag = gsub('\\_6', '', plantid)) %>%
  filter(year %in% 2023, tag %in% c(3374, 3643, 3767, 3955))
# these are all here lol

phen %>%
  filter(plot %in% 6) %>%
  mutate(tag = gsub('\\_6', '', plantid)) %>%
  filter(year %in% 2023, tag %in% c(3374, 3643, 3767, 3955, 5770)) %>% 
  arrange(survey.date)

# Okay. Well, maybe we can't reoncile this with demo.
# But we can still reconcile with each other for phen.

phen = phen %>% mutate(plantid = gsub('5770', '5570', plantid))

('5570_6' %in% phen$plantid & '5570_6' %in% seed$plantid)

# --- 3454_6

phen.seed %>% filter(grepl('\\_6$', plantid), !in.seed, year %in% 2023)
# hmm...
demo %>% filter(Tag %in% 3454, Year %in% 2023)
# okay... tag pulled because it wasn't flowering lol
# wait but it's somehow in seed...? LMAO what?

# Okay I don't know if this is fixable.

# --- 3338

# in raw seed dataset it's coord is 3I
# ahh it was excluded from phenology because of iffy records
# actually these records look okay in raw phenology... going back to fix

# --- 3987

# 3187 in phen?
demo %>% filter(Tag %in% c(3187, 3987), Plot %in% 5)
# hmm... no 3187 in demo (in plot 5)

phen %>% filter(grepl('3187', plantid), year %in% 2023)
# ah... excluded because of a missing week of phen

# 3481_15

demo %>% filter(Tag %in% 3481) # ah this is not hte one that had its tag changed
# ah it was dropped because its stem accidentally was clipped
# so it goes

# --- 2024 ---------

phen.seed %>% filter(year %in% 2024, in.phen, !in.seed) %>% arrange(plantid) %>% pull(plantid)
# "3126_13" "3530_4"  "3796_13" "3882_15" "7546_5"  "BSP_5"

phen.seed %>% filter(year %in% 2024, !in.phen, in.seed) %>% arrange(plantid) %>% pull(plantid)
# [1] "3108_13"    "3176_13"    "3447_7"     "3455_15"    "3515_13"    "3546_5.10D" "3546_5.9G" 
# [8] "3590_13"    "3643_6"     "3730_4"     "3769_13"    "3955_6"    

# 3455_15, 3463, 3955, 3590, 3515 also excluded for from phen during processing

# --- 3108

# --- 3176

phen %>% filter(grepl('31[72]6\\_13', plantid), year %in% 2024)
seed %>% filter(grepl('31[72]6\\_13', plantid), year %in% 2024)
# oh and the tag in demo is 3726?
seed %>% filter(grepl('3726', plantid), year %in% 2024)
# ah... yes... hmm...

demo %>% filter(Tag %in% c(3176, 3126, 3726), Plot %in% 13)
# uh... there is no 3726 in 13
# it's 3126
# hmm...

phen %>% filter(grepl('3726', plantid) | grepl('3126', plantid), year %in% 2024)

# Tag in demo is 3126
# Change it in seed

seed = seed %>%
  mutate(
    tag = ifelse(tag %in% 3176 & plot %in% 13, 3126, tag),
    plantid = gsub('3176\\_13', '3126_13', plantid)
  )

('3126_13' %in% seed$plantid & '3126_13' %in% phen$plantid)

# Are there three or four umbels here?

demo %>% filter(Tag %in% 3126, Year %in% 2024, Plot %in% 13)
# lmao

# OH. 3726 is 3039 in demo (old tag)

phen = phen %>% mutate(plantid = gsub('3726', '3039', plantid))
seed = seed %>% mutate(plantid = gsub('3726', '3039', plantid))

# --- 3447

# Record was crossed out in phen
# (ah - has since been fixed!)

# --- 3730

# There was a typo in one of these.
# Tag in demo is 3730

phen = phen %>% mutate(plantid = gsub('3530\\_4', '3730_4', plantid))

('3730_4' %in% phen$plantid & '3730_4' %in% seed$plantid)

# --- 3769

# 3769/3796

demo %>% filter(Tag %in% c(3769, 3796), Plot %in% 13) # ?
demo %>% filter(grepl('3796', proc.note))
demo %>% filter(grepl('3769', proc.note))
# ah. Okay, in demo, they all have the old tag 3558.

phen = phen %>%
  mutate(
    plantid = gsub('3796\\_13', '3558_13', plantid),
    plantid = gsub('3769\\_13', '3558_13', plantid)
  )

seed = seed %>%
  mutate(
    plantid = gsub('3796\\_13', '3558_13', plantid),
    plantid = gsub('3769\\_13', '3558_13', plantid)
  )

('3558_13' %in% phen$plantid & '3558_13' %in% seed$plantid & '3558_13' %in% demo$plantid)

# --- 3882_15

# this one may have just been left out of seed counts by accident? CHECK RAW DATA.

# --- 7546_5

# it's in seed for 2022 but not for 2024...
# might have just been left out?

# --- BSP_5 

# Left out somehow...


### --- Final merge  ------------

# Our final data frame should look like this:

# - One row for each seed count (i.e., one per umbel)
# - One column with all bud dates for each plant 
# (can use cleverness in separate/unite/grep to separate this out for data)
# - All-to-all merge, allowing us to do tests on each full dataset independently
# - Year (ofc.), plantid, plot, etc

# Step 1: get one row per plant per year of phenology, with all bud days

phen.wide = phen %>%
  # I think we only want new buds, not dead umbels
  filter(varb %in% 'new') %>%
  select(-varb) %>%
  group_by(year, plantid, plot) %>%
  summarise(
    phen.umbels = n(),
    phen.julis = paste0(init.doy, collapse = '; ')
  )

head(phen.wide)

# Now, put combine phen.wide and seed:

seed.phen.combine = merge(
  seed %>% select(-tag) %>% rename(seed.notes = notes) %>% mutate(in.seed = TRUE),
  phen.wide %>% mutate(in.phen = TRUE),
  by = c('year', 'plantid', 'plot'), all.x = TRUE, all.y = TRUE
) %>%
  # turn NAs to FALSEs
  mutate(across(c(in.phen, in.seed), function(x) ifelse(is.na(x), FALSE, x))) %>%
  arrange(year, plot, plantid)

head(seed.phen.combine)

str(seed.phen.combine)

unique(seed.phen.combine$umbel.no)
# ??
# Not sure how one of the umbel numbers became 's', but, let's fix it

seed.phen.combine = seed.phen.combine %>%
  mutate(umbel.no = as.numeric(ifelse(umbel.no %in% 's', 1, umbel.no)))

# Think this is good to export.

write.csv(
  seed.phen.combine,
  row.names = FALSE, na = '',
  '01_data_cleaning/out/seed_phen_combined_cleaned_all.csv'
)

# Last exported 4 Jul 2024
