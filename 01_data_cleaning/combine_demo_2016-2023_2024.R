# Script for combining 2016-2023 demo data with 2024 demo data
# (no seeds yet)
# SN - 4 Jun 2024

library(ggplot2)
library(dplyr)
library(tidyr)

rm(list = ls())

# -------------

d.1623 = read.csv('01_data_cleaning/out/demo_imputed_survival.csv')

head(d.1623)

d.2024 = read.csv('00_raw_data/lomatium_demography/2024_Lomatium_Demography_Data.csv')

head(d.2024)

# -------------

# Prepare 2024 data
# Need to add a plantid column I suppose...

d.2024 = d.2024 %>%
  mutate(plantid = paste0(Tag, '_', Plot, '_', Xcoor, Ycoor))

# -------------

d.compare = merge(
  x = d.1623 %>%
    mutate(plantid = gsub('b', '', plantid)) %>%
    filter(Year %in% 2022:2023) %>% 
    group_by(plantid) %>%
    filter(any(No.leaves > 0)) %>%
    ungroup() %>%
    arrange(desc(Year)) %>%
    distinct(plantid, .keep_all = TRUE) %>% 
    select(plantid, last.old.year = Year, last.old.lfct = No.leaves) %>%
    mutate(in.old = TRUE) %>%
    separate(plantid, into = c('old.tag', 'old.plot', 'old.coords'), sep = '_', remove = FALSE) %>%
    mutate(old.tagplot = paste(old.tag, old.plot, sep = '_')) %>%
    group_by(old.tagplot) %>%
    mutate(n.tags.in.plot = n()) %>%
    ungroup() %>%
    mutate(plantid = ifelse(n.tags.in.plot > 1, plantid, old.tagplot)) %>%
    ungroup(),
  y = d.2024 %>% 
    mutate(new.tagplot = paste(Tag, Plot, sep = '_')) %>%
    group_by(new.tagplot) %>%
    mutate(n.tags.in.plot = n()) %>%
    ungroup() %>%
    mutate(plantid = ifelse(n.tags.in.plot > 1, plantid, new.tagplot)) %>%
    mutate(in.new = TRUE) %>% select(-Year),
  by = 'plantid', all.x = TRUE, all.y = TRUE
) %>%
  mutate(across(starts_with('in.'), function(x) ifelse(is.na(x), FALSE, x))) %>%
  select(-contains('tagplot')) %>%
  select(-contains('n.tags.in.plot'))

head(d.compare)
nrow(d.compare)

d.compare %>%
  group_by(in.old, in.new) %>%
  summarise(n = n())

d.compare %>%
  filter(!(in.old & in.new)) %>%
  arrange(plantid) %>%
  select(
    plantid, in.old, in.new,
    # old data
    old.tag, old.coords, last.old.year, last.old.lfct,
    # new data
    Tag, Xcoor, Ycoor, No.leaves
  ) # %>%
  # Assume (for now?) that if it was seen with no leaves in 2022 or 2023 that it's not 
  # filter(last.old.lfct > 0 | is.na(last.old.lfct))

# --- Fix discrepancies manually ---

d.2024 = d.2024 %>%
  # get rid of empty dupe record
  filter(!(Tag %in% 3180 & Xcoor %in% 17)) %>%
  rename(demo.note = notes) %>%
  mutate(
    proc.note = '',
    edited = FALSE
  ) %>%
  mutate(
    # plant 5022 is weirdly listed as 2055
    proc.note = ifelse(
      grepl('2055\\_5', plantid), 
      paste(proc.note, 'tag changed from 2055 (wrong) to 5022', sep = ';'),
      proc.note
    ),
    edited = ifelse(grepl('2055\\_5', plantid), TRUE, edited),
    plantid = gsub('2055\\_5', '5022_5', plantid),
    # using the two year rule, 3117 in 2024 is a new plant
    proc.note = ifelse(
      grepl('3117\\_13', plantid),
      paste(proc.note, 'tag modified for gap', sep = ';'),
      proc.note
    ),
    edited = ifelse(grepl('3117\\_13', plantid), TRUE, edited),
    plantid = gsub('3117\\_13', '3117b_13', plantid),
    # demo processing uses new tag (3142) instead of old tag (3387)
    proc.note = ifelse(
      grepl('3387\\_13', plantid),
      paste(proc.note, 'tag changed in 2019 (was 3387)', sep = ';'),
      proc.note
    ),
    edited = ifelse(grepl('3387\\_13', plantid), TRUE, edited),
    plantid = gsub('3387\\_13', '3142_13', plantid),
    # demo processing uses new tag (3481) instead of old tag (3389)
    proc.note = ifelse(
      grepl('3389\\_13', plantid),
      paste(proc.note, 'tag changed in 2020 (was 3389)', sep = ';'),
      proc.note
    ),
    edited = ifelse(grepl('3389\\_13', plantid), TRUE, edited),
    plantid = gsub('3389\\_13\\_', '3481_13', plantid),
    # plant with old tag 3396 in 17A had tag replaced with 3989
    proc.note = ifelse(
      grepl('3396\\_13\\_17A', plantid),
      paste(proc.note, 'tag changed based on 2023 note *(was 3396)', sep = ';'),
      proc.note
    ),
    edited = ifelse(grepl('3396\\_13\\_17A', plantid), TRUE, edited),
    plantid = gsub('3396\\_13\\_17A', '3989\\_13\\_17A', plantid),
    # plant with new tag 3769 in field was 3558 in old demo
    proc.note = ifelse(
      grepl('3769\\_13', plantid),
      paste(proc.note, 'tag changed; was 3769', sep = ';'),
      proc.note
    ),
    edited = ifelse(grepl('3769\\_13', plantid), TRUE, edited),
    plantid = gsub('3769\\_13', '3558_13', plantid),
    # plant with new tag 3726 in field was old tag 3039
    proc.note = ifelse(
      grepl('3726\\_13', plantid),
      paste(proc.note, 'tag changed; was 3726', sep = ';'),
      proc.note
    ),
    edited = ifelse(grepl('3726\\_13', plantid), TRUE, edited),
    plantid = gsub('3726\\_13', '3039_13', plantid),
    # plant with new tag 7587 in field was old tag 3456
    proc.note = ifelse(
      grepl('3769\\_13', plantid),
      paste(proc.note, 'edited tag; was 7587', sep = ';'),
      proc.note
    ),
    edited = ifelse(grepl('7587\\_6', plantid), TRUE, edited),
    plantid = gsub('7587\\_6', '3456_6', plantid)
  )

# 3198_1 - tag accidentally pulled (?) in 2023, no new plant around ~16H, assume dead
# 3364_14... I think this is 3838

# --- Add missing columns to 2024 data ---

#  Columns to add:
# - trt (can be merged in)
# - obs.alive (easy to do)
# - surv (will require some thought...)

# Treatment column
d.2024 = merge(
  x = d.2024, y = read.csv('00_raw_data/plot_treatments.csv'),
  by.x = 'Plot', by.y = 'plot'
)

head(d.2024)
# looks good
d.2024 %>% distinct(trt) # no NAs - good

# Observed alive column - just use leaf count > 0
# but there are some cases where it's NA - right?
d.2024 %>% filter(is.na(No.leaves > 0))
# ugh... should combine these at some point I guess...

d.2024 %>% filter(!grepl('NP', demo.note), !No.leaves)
# hmm... two cases here, one of which is for a known-living plant

d.2024 = d.2024 %>%
  mutate(obs.alive = !grepl('NP', demo.note) | grepl('3480\\_2', plantid))

d.2024 %>% group_by(obs.alive) %>% summarise(n = n())
# dang... 222 plants not observed alive this year.

# Now... for survived
# In 2024, we can say true, but otherwise, can't say for certain.
# Probably a good idea to leave this NA for unsen plants
# But we can go back and edit 2023 records


