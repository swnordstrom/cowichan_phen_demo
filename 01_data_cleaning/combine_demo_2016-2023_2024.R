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

d.2024 = d.2024 %>% mutate(obs.alive = !grepl('NP', demo.note) | grepl('3480\\_2', plantid))

d.2024 %>% group_by(obs.alive) %>% summarise(n = n())
# dang... 222 plants not observed alive this year.

# Now... for survived
# In 2024, we can say true, but otherwise, can't say for certain.
# Probably a good idea to leave this NA for unsen plants
# But we can go back and edit 2023 records
# I guess to do this we want to combine the data frames now?

demo.all = rbind(
  d.1623 %>%
    select(
      plantid, Plot, trt, Year, surv, obs.alive, No.leaves, Leaf.length,
      No.umbels, edited, proc.note, demo.note
  ),
  d.2024 %>%
    mutate(
      obs.alive = !grepl('NP', demo.note),
      surv = NA
    ) %>%
    select(
      plantid, Plot, trt, Year, surv, obs.alive, No.leaves, Leaf.length,
      No.umbels, edited, proc.note, demo.note
    )
)

nrow(demo.all)

# Test code to make sure we're doing the survival imputation correctly
# demo.all %>%
#   filter(Year %in% 2023:2024) %>%
#   mutate(tagplotcoord = gsub('b', '', plantid)) %>%
#   group_by(tagplotcoord) %>%
#   # this throws an error...
#   # mutate(
#   #   new.surv = case_when(
#   #     n() < 2 ~ 'one.record',
#   #     diff(obs.alive) > 0 ~ 'impute',
#   #     diff(obs.alive) <= 0 ~ 'stay.same',
#   #     .default = NA
#   #   )
#   # )
#   mutate(rr = n()) %>%
#   ungroup() %>%
#   rbind(
#     . %>% filter(rr < 2) %>% mutate(surv2 = surv),
#     . %>% 
#       filter(rr > 1) %>%
#       group_by(tagplotcoord) %>%
#       mutate(surv2 = diff(obs.alive) > 0) %>%
#       ungroup()
#   )

rbind(
  demo.all %>% 
    filter(Year %in% 2023:2024) %>%
    mutate(tagplotcoord = gsub('b', '', plantid)) %>%
    group_by(tagplotcoord) %>%
    filter(n() < 2) %>%
    ungroup() %>%
    mutate(x = FALSE),
  demo.all %>% 
    filter(Year %in% 2023:2024) %>%
    mutate(tagplotcoord = gsub('b', '', plantid)) %>%
    group_by(tagplotcoord) %>%
    filter(n() > 1) %>%
    filter(diff(obs.alive) <= 0) %>%
    mutate(x = FALSE) %>%
    ungroup(),
  demo.all %>% 
    filter(Year %in% 2023:2024) %>%
    mutate(tagplotcoord = gsub('b', '', plantid)) %>%
    group_by(tagplotcoord) %>%
    filter(n() > 1) %>%
    filter(diff(obs.alive) > 0) %>%
    mutate(x = TRUE) %>%
    ungroup()
)

demo.all = rbind(
  # Separate out all records from pre-2023
  demo.all %>% filter(Year < 2023),
  # 
  # This is SUPER inelegant but the other things I've tried throw annoying
  # recycling errors
  # (each of these requires going through and reconciling 'b' plants)
  # Case 1: only one record for the plant in these two years
  demo.all %>% 
    filter(Year %in% 2023:2024) %>%
    mutate(tagplotcoord = gsub('b', '', plantid)) %>%
    group_by(tagplotcoord) %>%
    filter(n() < 2) %>%
    ungroup() %>%
    select(-tagplotcoord) %>%
    mutate(surv = surv),
  # Case 2: plant has records in both years and does not need correcting
  # (does not go from not observed to observed)
  demo.all %>% 
    filter(Year %in% 2023:2024) %>%
    mutate(tagplotcoord = gsub('b', '', plantid)) %>%
    group_by(tagplotcoord) %>%
    filter(n() > 1) %>%
    # this is the line that tells us if we need to correct anything
    filter(diff(obs.alive) <= 0) %>%
    mutate(
      plantid = ifelse(any(grepl('b', plantid)), grep('b', plantid, value = TRUE), plantid),
      surv = surv
    ) %>%
    ungroup() %>%
    select(-tagplotcoord),
  demo.all %>% 
    filter(Year %in% 2023:2024) %>%
    mutate(tagplotcoord = gsub('b', '', plantid)) %>%
    group_by(tagplotcoord) %>%
    filter(n() > 1) %>%
    # this is the line that tells us if we need to correct anything
    # (obs.alive going fromm FALSE to TRUE)
    filter(diff(obs.alive) > 0) %>%
    mutate(
      plantid = ifelse(any(grepl('b', plantid)), grep('b', plantid, value = TRUE), plantid),
      # just change all of the records to TRUE
      surv = TRUE,
      # leave a note
      proc.note = ifelse(
        Year %in% 2023, 
        paste(proc.note, 'survival imputed based on subsq. record', sep = ';'),
        proc.note
      )
    ) %>%
    ungroup() %>%
    select(-tagplotcoord)
) %>%
  arrange(Plot, plantid, Year)

nrow(demo.all)

# Ah - let's go through now and see if there are any other survival gaps that
# are longer than one year...
# (In which case we'd want to add 'b' to their tags... ugh wish there was a
# smoother way to do this.)

demo.all %>%
  # arrange(plantid, Year) %>%
  filter(surv | obs.alive) %>%
  group_by(plantid) %>%
  filter(n() > 1) %>%
  filter(any(diff(Year) > 1))
# Oh actually I guess we're good?
# (this data frame would be non-empty if there were corrections we needed to make)

# Hmm... okay cool.

# Okay I guess there's that one pair of plants that need to be combined

demo.all %>% filter(grepl('fix', demo.note))
demo.all %>% filter(grepl('3361\\_5', plantid) | grepl('5044\\_5', plantid))
# demo.all %>% filter(grepl('3661', plantid))
# that note should be 3361 not 3661
# oh wow, lmao
# 3361 was not seen in 2020-2022,
# 5044 appeared in 2021-2024
# so what appears to have happened is that the plant was not seen in 2020 (the
# year with sparser sampling) and assumed dead after that, but the plant was in
# fact 3361

demo.all = rbind(
  demo.all %>% filter(!(grepl('3361\\_5', plantid) | grepl('5044\\_5', plantid))),
  demo.all %>%
    filter(grepl('3361\\_5', plantid) | grepl('5044\\_5', plantid)) %>%
    filter(!(grepl('3361\\_5', plantid) & Year > 2020)) %>%
    mutate(
      plantid = '3361_5_13A',
      surv = ifelse(Year %in% 2020, TRUE, surv),
      edited = ifelse(Year >= 2020, TRUE, edited),
      proc.note = ifelse(Year %in% 2020, 'survival imputed based on 2021 record', proc.note),
      proc.note = ifelse(Year > 2020, 'id changed to old tag 3361 (was 5044 which was a dupe)', proc.note)
    )
)

nrow(demo.all)
length(unique(demo.all))

# I guess final step should be to update the surv column.

demo.all = demo.all %>%
  mutate(
    surv = case_when(
      Year %in% 2024 & grepl('NP\\;\\stag\\spull', demo.note) ~ FALSE,
      Year %in% 2024 & !grepl('NP', demo.note) ~ TRUE,
      .default = surv
    )
      # ifelse(Year %in% 2024 & !grepl('NP', demo.note), TRUE, surv)
  )

demo.all %>% group_by(Year, surv) %>% summarise(n = n())
# Looks good to me (for now)

write.csv(
  demo.all,
  na = '',  row.names = FALSE,
  '01_data_cleaning/out/Lomatium_demo_2016-2024.csv'
)
