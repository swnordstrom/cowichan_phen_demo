# Script for cleaning and processing raw Lout "resurvey" (phenology) datasets
# 2021 - 2023
# Objective here: get date for initiation of flowering
# Output is to get a table with unique ID for plant with flowering start date
# SN - init 12 Oct 2023

##### Setup

library(ggplot2)
library(dplyr)
library(tidyr)

rm(list = ls())

# Read in raw demo - one list for each file
raw.demo.list = dir('00_raw_data/lomatium_demography/') %>%
  grep(pattern = 'Resurvey', x = ., value = TRUE) %>%
  paste0('00_raw_data/lomatium_demography/', .) %>%
  lapply(FUN = read.csv)

##### Assess data

length(raw.demo.list)

lapply(raw.demo.list, names)

# 2021 data is different from 2022-2023

# From metadata:
# 2021: capable of pollen transfer in stages: "partial, flat, open, olp"
#   (for start date, can do partial... or flat - ask Jennifer)
# 2022: capable of pollen transfer in stages: "flowers"
#   (can use this for start date - not sure if "pods" is partially open or not)
# 2023: capable of pollen transfer in stages: "flowers"
#   (ibid)

# Also: unique plant identifiers
# 2021: plot, tag, coordinates
raw.demo.list[[1]] %>% group_by(plot, tag, survey.date) %>% filter(n() > 1)
# these have distinct coordinates
# tags are repeated though
raw.demo.list[[1]] %>% distinct(plot, tag) %>% group_by(tag) %>% filter(n() > 1)
# 2022: plot, tag
raw.demo.list[[2]] %>% group_by(plot, tag, survey.date) %>% filter(n() > 1)
# no repeated tags within a plot
# 2023: plot, tag, "coor" (coordinate)
raw.demo.list[[3]] %>% group_by(plot, tag, survey.date) %>% filter(n() > 1)
# have distinct coordinates

# There are notes columns though. I'll take care of these on a case-by-case basis.

##### Start cleaning data

# First initialize separate data frames for processed data

proc21 = raw.demo.list[[1]]
proc22 = raw.demo.list[[2]]
proc23 = raw.demo.list[[3]]

# Add unique plant identifiers ("plantid")
proc21$plantid = with(proc21, paste0(tag, "_", plot, "_", xcoor, ycoor))
proc22$plantid = with(proc22, paste(tag, plot, sep = "_"))
proc23$plantid = with(proc23, paste(tag, plot, coor, sep = "_"))

head(proc21$plantid)
head(proc22$plantid)
head(proc23$plantid)

# Initialize a data frame of plantids to exclude
exclude.plantids = rbind(
  data.frame(plantid = unique(proc21$plantid)),
  data.frame(plantid = unique(proc22$plantid)),
  data.frame(plantid = unique(proc23$plantid))
) %>%
  mutate(exclude = FALSE)

##### Look for missing records

# 2021 data

# Look for cases where there is a gap in records longer than 9 days
proc21$survey.date = as.Date(proc21$survey.date, format = '%m/%d/%Y')

proc21 %>%
  group_by(plot, survey.date) %>%
  summarise(n = n()) %>%
  ggplot(aes(x = survey.date, y = plot)) +
  geom_text(aes(label = n))

proc21 %>%
  group_by(plot, survey.date) %>%
  summarise(n = n()) %>%
  ggplot(aes(x = survey.date, y = plot)) +
  geom_text(aes(label = n)) +
  geom_segment(
    data = data.frame(survey.date = as.Date('2021-04-14') + (0:9)*7) %>%
      mutate(survey.date = ifelse(survey.date %in% as.Date('2021-05-12'),
                                  as.Date('2021-05-11'), as.Date(survey.date)),
             survey.date = ifelse(as.Date(survey.date) %in% as.Date('2021-06-16'),
                                  as.Date('2021-06-14'), as.Date(survey.date))) %>%
      mutate(survey.date = as.Date(survey.date)),
    aes(x = survey.date, xend = survey.date, y = 0, yend = 16),
    linetype = 2
  ) +
  # geom_segment() in here twice to make sure the boundaries set for cut() are
  # right - if they are, the vertical lines should overlap and we the cut()
  # boundaries are correct
  geom_segment(
    data = data.frame(survey.date = min(as.numeric(proc21$survey.date)) + 
                        cumsum(c(0, 7, 7, 7, 6, 8, 7, 7, 7, 5, 9))) %>%
      mutate(survey.date = as.Date(survey.date)),
    aes(x = survey.date, xend = survey.date, y = 0, yend = 16),
    linetype = 2, colour = 'blue'
  )
# Okay I think these boundaries will work

proc21 = proc21 %>%
  mutate(survey.period = cut(
    as.numeric(survey.date),
    breaks = min(as.numeric(survey.date)) -0.5 + cumsum(c(0, 7, 7, 7, 6, 8, 7, 7, 7, 5, 9))
    )
  ) %>%
  mutate(survey.period = as.numeric(factor(survey.period)))

unique(proc21$survey.period) # no NAs
sum(is.na(proc21$survey.period))

# Looking for gaps in records:
# (gaps: missing survey period, but records before/after the missing period)
proc21 %>% 
  group_by(plantid) %>% 
  filter(any(diff(survey.period) > 1)) %>%
  arrange(plantid) %>%
  as.data.frame()

# REMOVE plants if they have a gap that makes the first opening date ambiguous

# exclude 3108_13_12G - 15 day gap between "closed" and "flat"
exclude.plantids = exclude.plantids %>%
  mutate(exclude = case_when(plantid %in% '3108_13_12G' ~ TRUE, .default = exclude))

# exclue 3485_15_3A - 13 day gap between "unfurling" and "partial"
exclude.plantids = exclude.plantids %>%
  mutate(exclude = case_when(plantid %in% '3485_15_3A' ~ TRUE, .default = exclude))

### 2022 data
# look for cases with gap in records longer than 9 days
proc22$survey.date = as.Date(proc22$survey.date, format = '%m/%d/%Y')

# there is a record with date = NA (no other info in record - remove)
proc22 = proc22 %>% filter(!is.na(survey.date))

# Make a plot
proc22 %>%
  group_by(plot, survey.date) %>%
  summarise(n = n()) %>%
  ggplot(aes(x = survey.date, y = plot)) +
  geom_text(aes(label = n)) +
  geom_segment(
    data = data.frame(survey.date = as.Date('2022-04-04') + c(0:11, 13)*7) %>%
      mutate(survey.date = ifelse(as.Date(survey.date) %in% as.Date('2022-04-26'),
                                  as.Date('2022-04-25'), as.Date(survey.date)),
             survey.date = ifelse(as.Date(survey.date) %in% as.Date('2022-05-09'),
                                  as.Date('2022-05-08'), as.Date(survey.date)))  %>%
      mutate(survey.date = as.Date(survey.date)),
    aes(x = survey.date, xend = survey.date, y = 0, yend = 16),
    linetype = 2
  )
# duplication at plot 5 some time in mid-May...

proc22 = proc22 %>%
  # use cut() to assign each survey record to one sampling period
  mutate(survey.period = cut(
    as.numeric(survey.date),
    breaks = data.frame(survey.date = as.Date('2022-04-04') + c(0:11, 13)*7) %>%
      mutate(survey.date = ifelse(as.Date(survey.date) %in% as.Date('2022-04-26'),
                                  as.Date('2022-04-25'), as.Date(survey.date)),
             survey.date = ifelse(as.Date(survey.date) %in% as.Date('2022-05-09'),
                                  as.Date('2022-05-08'), as.Date(survey.date))) %>%
      mutate(survey.date = as.numeric(survey.date) - 0.5) %>%
      pull(survey.date)
    )
  ) %>%
  mutate(survey.period = as.numeric(factor(survey.period)))

# Look for gaps in records
proc22 %>% 
  group_by(plantid) %>%
  filter(any(diff(survey.period) > 1))
# oh... shoot, plots 6 and 15 were missed in the earliest June sampling period
# darn... what to do about this... rats

proc22 %>% 
  filter(!plot %in% c(6, 15)) %>%
  group_by(plantid) %>%
  filter(any(diff(survey.period) > 1))
# One one individual not in the above plots, and it flowered before the gap, so
# it can stay in the dataset

### 2023 data
proc23$survey.date = as.Date(proc23$survey.date, format = '%m/%d/%Y')

# there is a record with date = NA (no other info in record - remove)
proc23 = proc23 %>% filter(!is.na(survey.date))

proc23 %>% 
  group_by(survey.date, plot) %>%
  summarise(n = n()) %>%
  ggplot(aes(x = survey.date, y = plot)) + 
  geom_text(aes(label = n)) +
  geom_segment(
    data = data.frame(survey.date = as.Date('2023-04-18') + c(0:8)*7) %>%
      mutate(survey.date = ifelse(as.Date(survey.date) %in% as.Date('2023-05-16'),
                                  as.Date('2023-05-14'), as.Date(survey.date))) %>%
      mutate(survey.date = as.Date(survey.date)),
    aes(x = survey.date, xend = survey.date, y = 0, yend = 16),
    linetype = 2
  )
# not entirely sure why there are records for plots 6, 13 after the final sampling period...
# also looks like plot 12 was totally missed in sample period 4... ugh

# Add survey period column
proc23 = proc23 %>%
  # use cut() to assign each survey record to one sampling period
  mutate(survey.period = cut(
    as.numeric(survey.date),
    breaks = data.frame(survey.date = as.Date('2023-04-18') + c(0:8)*7) %>%
      mutate(survey.date = ifelse(as.Date(survey.date) %in% as.Date('2023-05-16'),
                                  as.Date('2023-05-14'), as.Date(survey.date))) %>%
      mutate(survey.date = as.numeric(survey.date) - 0.5) %>%
      pull(survey.date)
    )
  ) %>%
  mutate(survey.period = as.numeric(factor(survey.period)))

# Look for gaps in records
proc23 %>%
  group_by(plantid) %>%
  filter(any(diff(survey.period) > 1))
# lots of this from plot 12 of course
proc23 %>%
  filter(!plot %in% 12) %>%
  group_by(plantid) %>%
  filter(any(diff(survey.period) > 1))
# but none other than that


##### Multiple records in the same survey period: why would this happen?

### 2021

proc21 %>%
  # Get the plants (if any) that show up multiple times in the same survey period
  group_by(plantid, survey.period) %>%
  filter(n() > 1) %>%
  arrange(plantid, survey.date)
# as it turns out, these are the only three plants in plot 3
# for some reason, they were surveyed twice...
# (plot 5 was also surveyed over multiple days on 4 may, but it looks like only
# one plant)

# I say, in the case of multiple records for a plant within a sampling period,
# as a rule pick the earliest one and remove the rest.

# Remove these
proc21 = proc21 %>%
  filter(!(plot %in% 3 & as.Date(survey.date) %in% as.Date('2021-05-04')))

### 2022

proc22 %>%
  # Get the plants (if any) that show up multiple times in the same survey period
  group_by(plantid, survey.period) %>%
  filter(n() > 1) %>%
  arrange(plantid, survey.date)
# Not sure why this was done...
# There are 3-4 plants that were not observed the first day, but observed three days later
# oh well! gotta stick with the rule.

# are there any plants surveyed on just the 19th?
proc22 %>% 
  filter(plot %in% 5, survey.period %in% 7) %>% 
  group_by(plantid) %>% 
  summarise(n = n()) %>%
  group_by(n) %>%
  summarise(nn = n())
# everyone was observed twice here...
# okay, good, just take out the 19 May records

proc22 = proc22 %>%
  filter(!(plot %in% 5 & as.Date(survey.date) %in% as.Date('2022-05-19')))

### Check 2023
proc23 %>%
  group_by(plantid, survey.period) %>%
  filter(n() > 1)
# cool!


##### Look into survey periods where only 1-2 plots were surveyed...

### 2021

proc21 %>%
  group_by(survey.period) %>%
  summarise(n.plots = length(unique(plot)))

# hmm...
# one plot surveyed in period 1, six in period 2,
# one plot missing from period 3, one missing from final period
proc21 %>%
  filter(survey.period %in% 1)
# this is only one plant...

proc21 %>% filter(plot %in% 7)
# but there are other plants too...
# looks like inconsistencies in how plants were added (not sure if this is
# inconsistent in this year but seems inconsistent with other years...)

proc21 %>% filter(survey.period %in% 2)
# not seeing anything horrible or alarming here...

# plot 3 is missing from the survey in week 10
proc21 %>% filter(survey.period %in% 9, plot %in% 3)
# uh... assume all done seeding?

### 2022

proc22 %>%
  group_by(survey.period) %>%
  summarise(n.plots = length(unique(plot)))

proc22 %>% filter(survey.period %in% 1)
# maybe all plants in these plots had records back-filled...
# but then why only for a subset of plots?
proc22 %>% filter(survey.period %in% 2:3, plot %in% 10)

# well, okay, this seems fine

### 2023

proc23 %>%
  group_by(survey.period) %>%
  summarise(n.plots = length(unique(plot)))
# one plot missing from period 1, one missing from period 4

# (the one missing from period 4 was covered before - that's pot 12)
proc23 %>% filter(plot %in% 12, survey.period %in% 3:5) %>%
  arrange(tag, survey.date)
# ugh... yep, so many flowering dates might have been missed...

# Plot 3 is the one missing from survey period 1
proc23 %>% filter(plot %in% 3, survey.period < 3)
# hmm... okay I guess this seems fine

# What is up with the two plots that got resurveyed in period 8?
proc23 %>% filter(plot %in% c(6, 13), survey.period %in% 7:8) %>%
  arrange(tag, survey.period)
# maybe these were just the very last plants to be flowering...

##### Plants with no records of flowering at all

proc21 %>%
  group_by(plantid) %>%
  summarise(
    any.fl = any(no.partial > 0 | no.flat > 0 | no.open > 0 | no.olp > 0 | no.seeding > 0)
  ) %>%
  group_by(any.fl) %>%
  summarise(n = n())
# nine NAs... are these just plants that are NAs in the above columns?
# (looking at these... not super closely but it looks like yes,
# and they are NAs because they "stopped checking, dead > 2 weeks")

# proc21 %>%
#   group_by(plantid) %>%
#   summarise(na.fl = all(is.na(contains('no.')))) %>%
#   group_by(any.fl) %>%
#   summarise(n = n())

proc21 %>%
  group_by(plantid) %>%
  filter(is.na(all(no.partial > 0 | no.flat > 0 | no.open > 0 | no.olp > 0 | no.seeding > 0)))
# oh... good catch
# here, the all-NA record is in the middle of the observation period...
# how do I get these out automatically...

proc21 %>%
  group_by(plantid) %>%
  filter(
    any(is.na(no.partial) & is.na(no.flat) & is.na(no.open) & is.na(no.olp) & is.na(no.seeding))
  ) %>%
  arrange(tag, plot, survey.period) %>%
  # Filter out any plant that has any non-NA records after an NA record
  filter(any(survey.period[is.na(no.partial)] < max(survey.period[!is.na(no.partial)])))

# Just this one plant - I say remove it because that missing record is *right*
# where the partial date should be
exclude.plantids = exclude.plantids %>%
  mutate(exclude = ifelse(plantid %in% '3563_4_12S', TRUE, exclude))

###

proc22 %>%
  group_by(plantid) %>%
  summarise(any.fl = any(no.flowers > 0 | no.flp > 0 | no.seeding > 0)) %>%
  group_by(any.fl) %>%
  summarise(n = n())

proc22 %>% group_by(plantid) %>% filter(!any(no.flowers > 0 | no.flp > 0 | no.seeding > 0))
proc22 %>% filter(tag %in% 3516)
proc22 %>% filter(plot %in% 13, survey.period %in% 5)
proc22 %>% filter(plot %in% 13, survey.period %in% 6)  
proc22 %>% filter(tag %in% 3515)
# not sure this needs to be fixed...

proc22 %>%
  group_by(plantid) %>%
  filter(all(is.na(no.flowers) & is.na(no.flp) & is.na(no.seeding)))
# good - nothing for which every record is an NA

proc22 %>%
  group_by(plantid) %>%
  filter(is.na(any(no.flowers > 0 | no.flp > 0 | no.seeding > 0)))
# ah - these are cases where there are NAs, and non-NA records with no
# observations of flowering

# I think this will be okay. I'll just filter out plants that have at least one
# non-zero, non-NA in a relevant field.

### 2023

proc23 %>%
  group_by(plantid) %>%
  summarise(any.fl = any(no.flowers > 0 | no.flp > 0 | no.seeding > 0)) %>%
  group_by(any.fl) %>%
  summarise(n = n())

proc23 %>%
  group_by(plantid) %>%
  filter(is.na(any(no.flowers > 0 | no.flp > 0 | no.seeding > 0))) %>%
  arrange(tag, survey.period)

# ugh... tons of these
# many cases where there are buds listed but nothing but NAs after that
# (why??)
# very few of them have notes, the notes provided aren't super informative
# many go from bud to dead
# a couple of cases with the note "yellow and wilted" - is this flowering?
# (probably good to be consistent and say no)
# also some cases of being eaten

##### Actually do stuff

# Go through each one of these and filter out *only* the plants with a non-zero
# record in a field associated with a reproductive stage.
# Given that there are records where an any() flag fails because of NAs, let's
# first change the NAs to zeros.

### 2021 

flower21 = proc21 %>%
  # Change NAs to zeros
  mutate(
    across(
      contains('no.'), 
      .fn = function(x) ifelse(is.na(x), 0, x)
      )
  ) %>%
  # Now, give me only the plants that have 
  group_by(plantid) %>%
  filter(any(no.partial > 0 | no.flat > 0 | no.open > 0 | no.olp > 0 | no.seeding > 0)) %>%
  ungroup()

# Number of plants (thus far)
length(unique(proc21$plantid))
with(proc21, table(plot, survey.period))

phen21 = flower21 %>%
  group_by(plantid, plot, tag) %>%
  summarise(
    # First survey period where any activity was observed
    fl.init = min(survey.period[no.partial + no.flat + no.open + no.olp + no.seeding > 0]),
    # Final survey period where receptive tissues were observed
    fl.fina = ifelse(
      any(no.partial + no.flat + no.open + no.olp > 0),
      max(survey.period[no.partial + no.flat + no.open + no.olp > 0]),
      NA
    )
  )

head(phen21)  

with(phen21, table(plot, fl.init))

flower21 %>%
  select(-c(no.closed, no.unfurling)) %>%
  mutate(incl.dead = rowSums(across(contains('no.')))) %>%
  group_by(plantid) %>%
  mutate(
    all.incr = diff(c(0, incl.dead)),
    flo.incr = diff(c(0, no.dead + no.eaten))
  ) %>%
  arrange(plantid)
# ah... should assume that if a plant skips from closed/unfurling to dead that it didn't flower

ind.flw.phen21 = flower21 %>%
  arrange(survey.date) %>%
  select(-c(no.closed, no.unfurling)) %>%
  mutate(incl.dead = rowSums(across(contains('no.')))) %>%
  group_by(plantid, plot) %>%
  mutate(
    all.incr = diff(c(0, incl.dead)),
    dea.incr = diff(c(0, no.dead + no.eaten)),
    flo.incr = all.incr - dea.incr
  ) %>%
  filter(flo.incr > 0) %>%
  uncount(flo.incr) %>%
  select(plantid, plot, survey.period)

### 2022

flower22 = proc22 %>%
  # Change NAs to zeros
  mutate(
    across(
      contains('no.'), 
      .fn = function(x) ifelse(is.na(x), 0, x)
    )
  ) %>%
  # Now, give me only the plants that have 
  group_by(plantid) %>%
  filter(any(no.flowers + no.flp + no.seeding + no.broken + no.dead > 0)) %>%
  ungroup()

length(unique(flower22$plantid))
# over 400 plants??

phen22 = flower22 %>%
  group_by(plantid, plot, tag) %>%
  summarise(
    # First survey period where any activity was observed
    fl.init = min(survey.period[no.flowers + no.flp + no.seeding + no.dead > 0]),
    # Final survey period where receptive tissues were observed
    fl.fina = ifelse(
      any(no.flowers + no.flp > 0),
      max(survey.period[no.flowers + no.flp > 0]),
      NA
    )
  )

head(phen22)  

with(phen22, table(plot, fl.init))

ind.flw.phen22 = flower22 %>%
  arrange(survey.date) %>%
  select(-c(no.buds, no.pods)) %>%
  mutate(incl.dead = rowSums(across(contains('no.')))) %>%
  group_by(plantid, plot) %>%
  mutate(
    all.incr = diff(c(0, incl.dead)),
    dea.incr = diff(c(0, no.dead + no.eaten)),
    flo.incr = all.incr - dea.incr
  ) %>%
  filter(flo.incr > 0) %>%
  uncount(flo.incr) %>%
  select(plantid, plot, survey.period)

### 2023

flower23 = proc23 %>%
  # Change NAs to zeros
  mutate(
    across(
      contains('no.'), 
      .fn = function(x) ifelse(is.na(x), 0, x)
    )
  ) %>%
  # Now, give me only the plants that have 
  group_by(plantid) %>%
  filter(any(no.flowers + no.flp + no.seeding + no.broken + no.dead > 0)) %>%
  ungroup()

length(unique(flower23$plantid))
# over 300 plants

phen23 = flower23 %>%
  group_by(plantid, plot, tag) %>%
  summarise(
    # First survey period where any activity was observed
    fl.init = min(survey.period[no.flowers + no.flp + no.seeding + no.dead > 0]),
    # Final survey period where receptive tissues were observed
    fl.fina = ifelse(
      any(no.flowers + no.flp > 0),
      max(survey.period[no.flowers + no.flp > 0]),
      NA
    )
  )

head(phen23)

ind.flw.phen23 = flower23 %>%
  arrange(survey.date) %>%
  select(-c(no.buds, no.pods)) %>%
  mutate(incl.dead = rowSums(across(contains('no.')))) %>%
  group_by(plantid, plot) %>%
  mutate(
    all.incr = diff(c(0, incl.dead)),
    dea.incr = diff(c(0, no.dead + no.eaten)),
    flo.incr = all.incr - dea.incr
  ) %>%
  filter(flo.incr > 0) %>%
  uncount(flo.incr) %>%
  select(plantid, plot, survey.period)

##### Combine all data frames together (remember to add year)

phen.all = rbind(
  phen21 %>% ungroup() %>% mutate(year = 2021),
  phen22 %>% ungroup() %>% mutate(year = 2022),
  phen23 %>% ungroup() %>% mutate(year = 2023)
)

head(phen.all)
nrow(phen.all)
str(phen.all)

# Exclude plants designated for removal above
phen.all = merge(x = phen.all, y = exclude.plantids) %>%
  # Use the "exclude" column to remove these
  filter(!exclude) %>%
  # Remove unnecessary column
  select(-exclude) %>%
  # Arrange columns (for export)
  arrange(year, plot, tag)


# Merge the start days for each year and incorporate actual date in
phen.all = merge(
  x = phen.all,
  y = data.frame(
    year = 2021:2023,
    min.date = c(
      as.numeric(min(proc21$survey.date) - as.Date('2021-01-01')),
      as.numeric(min(proc22$survey.date) - as.Date('2022-01-01')),
      as.numeric(min(proc23$survey.date) - as.Date('2023-01-01'))
    )
  )
) %>%
  mutate(
    init.doy = min.date + (fl.init - 1) * 7,
    fina.doy = min.date + (fl.fina - 1) * 7
  ) %>%
  select(-min.date) %>%
  rename(init.wk = fl.init, fina.wk = fl.fina)

# Do the same for the individual-flowering one

ind.phen.all = rbind(
  # Bind everything together
  ind.flw.phen21 %>% mutate(year = 2021),
  ind.flw.phen22 %>% mutate(year = 2022),
  ind.flw.phen23 %>% mutate(year = 2023)
) %>%
  merge(y = exclude.plantids) %>%
  # Use the "exclude" column to remove these
  filter(!exclude) %>%
  # Remove unnecessary column
  select(-exclude) %>%
  # Arrange columns (for export)
  arrange(year, plot) %>%
  # Merge to add approx. survey date
  merge(
    y = data.frame(
      year = 2021:2023,
      min.date = c(
        as.numeric(min(proc21$survey.date) - as.Date('2021-01-01')),
        as.numeric(min(proc22$survey.date) - as.Date('2022-01-01')),
        as.numeric(min(proc23$survey.date) - as.Date('2023-01-01'))
      )
    )
  ) %>%
  mutate(init.doy = min.date + (survey.period - 1) * 7) %>%
  select(-min.date)

# Export to csvs

write.csv(
  phen.all,
  '01_data_cleaning/out/phenology_all_cleaned.csv',
  row.names = FALSE
)

write.csv(
  ind.phen.all,
  '01_data_cleaning/out/phenology_all_ind_cleaned.csv',
  row.names = FALSE
)


#####
# In future: might want to add eaten or premature death and being eaten to these
