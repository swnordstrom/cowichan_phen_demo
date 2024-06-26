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
raw.demo.list[[4]] %>% group_by(plot, tag, survey.date) %>% filter(n() > 1)
# 3546 in plot 5... some mistake in here (oh this is the 5309 plant)

# There are notes columns though. I'll take care of these on a case-by-case basis.

##### Start cleaning data

# First initialize separate data frames for processed data

proc21 = raw.demo.list[[1]]
proc22 = raw.demo.list[[2]]
proc23 = raw.demo.list[[3]]
proc24 = raw.demo.list[[4]]

# Add unique plant identifiers ("plantid")
proc21$plantid = with(proc21, paste0(tag, "_", plot, "_", xcoor, ycoor))
proc22$plantid = with(proc22, paste(tag, plot, sep = "_"))
proc23$plantid = with(proc23, paste(tag, plot, coor, sep = "_"))
proc24$plantid = with(proc24, paste(tag, plot, coor, sep = "_"))

head(proc21$plantid)
head(proc22$plantid)
head(proc23$plantid)
head(proc24$plantid)

# Initialize a data frame of plantids to exclude
exclude.plantids = rbind(
  data.frame(plantid = unique(proc21$plantid)) %>% mutate(year = 2021),
  data.frame(plantid = unique(proc22$plantid)) %>% mutate(year = 2022),
  data.frame(plantid = unique(proc23$plantid)) %>% mutate(year = 2023),
  data.frame(plantid = unique(proc24$plantid)) %>% mutate(year = 2024)
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
  mutate(exclude = case_when(plantid %in% '3108_13_12G' & year %in% 2021 ~ TRUE, .default = exclude))

# exclue 3485_15_3A - 13 day gap between "unfurling" and "partial"
exclude.plantids = exclude.plantids %>%
  mutate(exclude = case_when(plantid %in% '3485_15_3A' & year %in% 2021  ~ TRUE, .default = exclude))

# Fix mistakes in data entry
proc21 = proc21 %>%
  mutate(
    # number dead for 3127 was misentered on this day
    no.dead = ifelse(tag %in% 3127 & plot %in% 5 & survey.date %in% as.Date('2021-05-27'), 1, no.dead),
    # a seeding plant at 3396 mistakenly entered as dead on this day
    no.dead = ifelse(tag %in% 3396 & plot %in% 13 & survey.date %in% as.Date('2021-05-20'), 0, no.dead),
    no.seeding = ifelse(tag %in% 3396 & plot %in% 13 & survey.date %in% as.Date('2021-05-20'), 1, no.seeding),
    # a dead plant was not listed for 3065 (although probably doesn't matter...)
    no.dead = ifelse(tag %in% 3065 & plot %in% 15 & survey.date %in% as.Date('2021-05-27'), 1, no.dead),
    # seeding umbel not entered for 3567
    no.seeding = ifelse(tag %in% 3567 & plot %in% 1 & survey.date %in% as.Date('2021-5-27'), 2, no.seeding),
    # missing a plant on this day for 3298 (3296 but that is fixed later)...
    no.open = ifelse(tag %in% 3298 & plot %in% 7 & survey.date %in% as.Date('2021-05-20'), 1, no.open)
  )


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

# Fix other data entry mistakes
proc22 = proc22 %>%
  mutate(
    # eaten umbel was mistakenly entered in notes
    no.eaten = ifelse(plantid %in% '3004_14' & survey.date %in% as.Date('2022-05-26'), 1, no.eaten),
    notes = ifelse(plantid %in% '3004_14' & survey.date %in% as.Date('2022-05-26'), '', notes),
    # dead umbel not recorded on this date
    no.dead = ifelse(plantid %in% '3018_13' & survey.date %in% as.Date('2022-06-14'), 1, no.dead),
    # a seeding umbel was mistakenly recorded as dead here
    no.dead = ifelse(plantid %in% '3035_14' & survey.date %in% as.Date('2022-06-08'), 0, no.dead),
    no.seeding = ifelse(plantid %in% '3035_14' & survey.date %in% as.Date('2022-06-08'), 1, no.seeding),
    # mis-entry here - no dead plants
    no.dead = ifelse(plantid %in% '3063_15' & survey.date %in% as.Date('2022-06-06'), 0, no.dead),
    # bunch of mistakes for 3070...
    no.flowers = ifelse(plantid %in% '3070_6' & survey.date %in% as.Date('2022-05-25'), 2, no.flowers),
    no.dead = ifelse(plantid %in% '3070_6' & survey.date %in% as.Date('2022-06-06'), 1, no.dead),
    no.flp = ifelse(plantid %in% '3070_6' & survey.date %in% as.Date('2022-06-06'), 2, no.flp),
    # much confusion in here too
    no.broken = ifelse(plantid %in% '3122_4' & survey.date > as.Date('2022-05-31'), 1, no.broken),
    # this must be a recording mistake - two seeding plants observed before and after
    no.seeding = ifelse(plantid %in% '3122_4' & survey.date %in% as.Date('2022-06-14'), 2, no.seeding),
    # seeding umbel mistakenly entered as dead, also going to propagate forward the eaten plant
    no.seeding = ifelse(plantid %in% '3151_3' & survey.date %in% as.Date('2022-06-07'), 1, no.seeding),
    no.dead = ifelse(plantid %in% '3151_3' & survey.date %in% as.Date('2022-06-07'), 0, no.dead),
    # 3151_1 - going to propagate forward the eaten bud
    no.eaten = ifelse(plantid %in% '3151_1' & survey.date > as.Date('2022-04-25'), 1, no.eaten),
    # mistakenly did not record dead plant
    no.dead = ifelse(plantid %in% '3187_13' & survey.date %in% as.Date('2022-06-14'), 1, no.dead),
    # mistake - plant was eaten, never seeded
    no.seeding = ifelse(plantid %in% '3379_14' & survey.date %in% as.Date('2022-06-28'), 0, no.seeding),
    # did not enter a dead umbel here
    no.dead = ifelse(plantid %in% '3604_13' & survey.date %in% as.Date('2022-06-14'), 1, no.dead),
    # not sure what happened here - def an entry mistake (data sheet says 4D 2S)
    no.dead = ifelse(plantid %in% '3740_15' & survey.date %in% as.Date('2022-06-06'), 2, no.dead),
    no.seeding = ifelse(plantid %in% '3740_15' & survey.date %in% as.Date('2022-06-06'), 4, no.seeding),
    # did not enter a dead umbel here
    no.dead = ifelse(plantid %in% '3978_13' & survey.date %in% as.Date('2022-06-14'), 1, no.dead),
    # didn't enter a dead umbel here
    no.dead = ifelse(plantid %in% '7537_13' & survey.date %in% as.Date('2022-06-14'), 1, no.dead),
    # data entry mistake here - no flp on 5/31 for plant 7591_4
    no.flp = ifelse(plantid %in% '7591_4' & survey.date %in% as.Date('2022-05-31'), 0, no.flp),
    # data entry mistake - there's a question mark by the seeding and it isn't seen in subsq weeks
    no.seeding = ifelse(plantid %in% '3642_13' & survey.date %in% as.Date('2022-05-08'), 0, no.seeding),
    # data entry mistake - two seeding plants, neither dead...
    no.seeding = ifelse(plantid %in% '3185_2' & survey.date %in% as.Date('2022-06-21'), 2, no.seeding),
    no.dead = ifelse(plantid %in% '3185_2' & survey.date %in% as.Date('2022-06-21'), 0, no.dead),
    # data sheet has this as a three/two - makes more sense as a three
    no.flowers = ifelse(plantid %in% '3485_15' & survey.date %in% as.Date('2022-05-16'), 3, no.flowers),
    # this is driving me insane... why did they do this
    no.flp = ifelse(plantid %in% '1411_4' & survey.date %in% as.Date('2022-05-08'), 1, no.flp),
    # mis-entry - 1F, 1FLP, 1S (f missed) - tag is really 2230 but this gets fixed later
    no.flowers = ifelse(plantid %in% '2238_4' & survey.date %in% as.Date('2022-05-31'), 1, no.flowers),
    # mis-entry - didn't enter a dead plant...
    no.dead = ifelse(plantid %in% '3182_4' & survey.date %in% as.Date('2022-05-31'), 1, no.dead),
    # mis-entry - FLP erroneously entered as one
    no.flp = ifelse(plantid %in% '3364_4' & survey.date %in% as.Date('2022-05-31'), 2, no.flp),
    # mis-entry - entered the data from the row above, there is no flowering umbel here this day
    no.flowers = ifelse(plantid %in% '3730_4' & survey.date %in% as.Date('2022-05-31'), 0, no.flowers),
    # mis-entry - entered the data from the row above, there is no seeding umbel on this date 
    # but there is an FLP
    no.flowers = ifelse(plantid %in% '7545_4' & survey.date %in% as.Date('2022-05-31'), 0, no.flowers),
    no.seeding = ifelse(plantid %in% '7545_4' & survey.date %in% as.Date('2022-05-31'), 0, no.seeding),
    no.flp = ifelse(plantid %in% '7545_4' & survey.date %in% as.Date('2022-05-31'), 1, no.flp),
    # mis-entries for two days for this plant - "two" got put in wrong field
    no.flp = ifelse(plantid %in% '7550_4' & survey.date %in% as.Date('2022-05-31'), 1, no.flp),
    no.seeding = ifelse(plantid %in% '7550_4' & survey.date %in% as.Date('2022-06-07'), 2, no.seeding),
    # data recorded as a question mark... come on man that's so not helpful
    # going to assume it's a bud because pollen shedding would be conspicuous
    no.buds = ifelse(plantid %in% '3988_5' & survey.date %in% as.Date('2022-05-16'), 1, no.buds),
    # data is not recorded clearly but it's two FLP not one
    no.flp = ifelse(plantid %in% '5021_5' & survey.date %in% as.Date('2022-05-16'), 2, no.flp),
    # okay... two FLP/Dead (????) plants one week, then one dead+seeding the next week...
    # going to say that it was seeding the week it was missed (but they all died the next week)
    no.seeding = ifelse(plantid %in% '3070_6' & survey.date %in% as.Date('2022-06-14'), 2, no.seeding),
    # dead umbels were not recorded one week (come on man)
    no.dead = ifelse(plantid %in% '3319_6' & survey.date %in% as.Date('2022-06-14'), 2, no.dead),
    # mis-entries - flower got mis-entered as a pod on may 8, then number of
    # flowers got mis-recorded on may 25
    no.flowers = ifelse(plantid %in% '3630_6' & survey.date %in% as.Date('2022-05-08'), 1, no.flowers),
    no.pods = ifelse(plantid %in% '3630_6' & survey.date %in% as.Date('2022-05-08'), 0, no.pods),
    no.flowers = ifelse(plantid %in% '3630_6' & survey.date %in% as.Date('2022-05-25'), 2, no.flowers),
    # mis-entry - two dead umbels mis-recorded as three
    no.dead = ifelse(plantid %in% '3126_13' & survey.date %in% as.Date('2022-06-14'), 2, no.dead),
    # there's a question mark next to this bud record... going to assume it was
    # a mistaken field ID because it isn't seen again
    no.buds = ifelse(plantid %in% '3143_13' & survey.date %in% as.Date('2022-05-08'), 0, no.buds),
    # plant was missed two consecurtive weeks then found eaten... the zeros
    # cause the diff() to assume it's a new bud but it isn't;
    # going to assume it was eaten that first survey it wasn't seen
    no.eaten = ifelse(plantid %in% '3367_13' & survey.date %in% as.Date(c('2022-05-16', '2022-05-25')), 1, no.eaten),
    # mis-entry - one dead umbel mis-recorded as one
    no.seeding = ifelse(plantid %in% '3421_13' & survey.date %in% as.Date('2022-06-24'), 1, no.seeding),
    # mis-entry - didn't record flowering plant for some reason
    no.flowers = ifelse(plantid %in% '3428_13' & survey.date %in% as.Date('2022-05-25'), 1, no.flowers),
    # missing week means seeding plant was not recorded and diff() treats plants as different
    no.seeding = ifelse(plantid %in% '3705_13' & survey.date %in% as.Date('2022-06-14'), 1, no.seeding),
    # 3021... has to be a mistake on 5/16... but what is the truth?
    # god this data is such fucking garbage in some places
    # it's a 3 on the datasheet but it only can be consistent with the other data if it's a 1
    no.flowers = ifelse(plantid %in% '3021_15' & survey.date %in% as.Date('2022-05-16'), 1, no.flowers),
    # 3083... none of this shit makes any sense. god fucking damnit, what were these people thiinking
    # going to assume that the pods are buds
    no.buds = ifelse(plantid %in% '3083_15' & survey.date %in% as.Date('2022-05-16'), 2, no.buds),
    no.pods = ifelse(plantid %in% '3083_15' & survey.date %in% as.Date('2022-05-16'), 2, no.pods),
    # 2 on datasheet mis-entered as a 3
    no.seeding = ifelse(plantid %in% '3423_15' & survey.date %in% as.Date('2022-06-23'), 2, no.seeding),
    # 2 on datasheet mis-entered as a 3
    no.seeding = ifelse(plantid %in% '3892_15' & survey.date %in% as.Date('2022-06-14'), 2, no.seeding)
  )

exclude.plantids = exclude.plantids %>%
  # Plant 3345 in plot 6 - records are a messy disaster here. I don't trust them. Exclude
  mutate(exclude = ifelse(plantid %in% '3345_6' & year %in% 2022, TRUE, exclude)) %>%
  # Plant 3196 - note says "likely missed", data is inconsistent
  mutate(exclude = ifelse(plantid %in% '3196_13' & year %in% 2022, TRUE, exclude)) %>%
  # Plant 3430 - not surveyed one critical week, so possible issues
  mutate(exclude = ifelse(plantid %in% '3430_14' & year %in% 2022, TRUE, exclude)) %>%
  # Plant 3426 - again, just total fucking garbage
  mutate(exclude = ifelse(plantid %in% '3426_15' & year %in% 2022, TRUE, exclude)) %>%
  # Plant 3740 - note says this plant's records include a different plant (3403)
  # but there IS no 3403 in plot 15, only in plot 7
  # what the absolute fuck is going on in this data? christ what a fucking disaster
  mutate(exclude = ifelse(plantid %in% '3403_15' & year %in% 2022, TRUE, exclude)) %>%
  # Plant 3755 - pulled out, records unclear, seems likely missed true bud date
  mutate(exclude = ifelse(plantid %in% '3755_15' & year %in% 2022, TRUE, exclude)) %>%
  # Plant 3848 - missed one week, true flowering date unknown...
  mutate(exclude = ifelse(plantid %in% '3848_15' & year %in% 2022, TRUE, exclude)) %>%
  # Plant 3863 - WHY ARE YOU PUTTING QUESTION MARKS IN HERE - TAKE THE TIME TO
  # FIGURE OUT WHAT THE DATA ACTUALLY IS. I can't impute it because the true
  # flowering date isn't known
  mutate(exclude = ifelse(plantid %in% '3863_15' & year %in% 2022, TRUE, exclude))


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

# Fix data entry mistakes

proc23 = proc23 %>%
  mutate(
    # plant 3400 - data entry mistake - no dead plants on this day
    no.dead = ifelse(tag %in% 3400 & plot %in% 6 & survey.date %in% as.Date('2023-05-15'), 0, no.dead),
    # plant 3518 - day skipped but entries before and after are the same
    no.dead = ifelse(tag %in% 3518 & plot %in% 6 & survey.date %in% as.Date('2023-06-04'), 1, no.dead),
    no.seeding = ifelse(tag %in% 3518 & plot %in% 6 & survey.date %in% as.Date('2023-06-04'), 2, no.seeding),
    # plant 5770 - no record of new dead umbel
    no.dead = ifelse(tag %in% 5770 & plot %in% 6 & survey.date %in% as.Date('2023-06-08'), 0, no.dead),
    # plant 3431 - buds entered as dead
    no.dead = ifelse(tag %in% 3431 & plot %in% 13 & survey.date < as.Date('2023-05-25'), 0, no.dead),
    no.buds = ifelse(tag %in% 3431 & plot %in% 13 & survey.date %in% as.Date('2023-05-04'), 2, no.buds),
    # 3163 - note says eaten but didn't enter it in 'no.eaten' field...
    no.eaten = ifelse(tag %in% 3163 & plot %in% 1 & survey.date %in% as.Date('2023-05-24'), 1, no.eaten),
    # 3005 - why the heck are there four umbels listed here today...?
    # going to assume it flowered on 5/15 because it was in late-stage flowering the next visit
    no.buds = ifelse(tag %in% 3005 & plot %in% 6 & survey.date %in% as.Date('2023-05-15'), 0, no.buds),
    no.flowers = ifelse(tag %in% 3005 & plot %in% 6 & survey.date %in% as.Date('2023-05-15'), 1, no.flowers),
    # 3070 - dead plant not entered
    no.dead = ifelse(tag %in% 3070 & plot %in% 6 & survey.date %in% as.Date('2023-05-15'), 1, no.dead),
    # 3142 - assume can't find because dead
    no.dead = ifelse(tag %in% 3142 & plot %in% 6 & survey.date %in% as.Date('2023-05-15'), 2, no.dead),
    # 3374 - STOP WRITING QUESTION MARKS, THIS IS USELESS, DO YOUR JOB
    # okay. assume the question mark is a dead plant.
    no.dead = ifelse(tag %in% 3374 & plot %in% 6 & survey.date %in% as.Date('2023-05-10'), 1, no.dead),
    # plant 3743 - data entry error... entered data from two rows above...
    no.seeding = ifelse(tag %in% 3743 & plot %in% 6 & survey.date %in% as.Date('2023-06-08'), 0, no.seeding),
    no.eaten = ifelse(tag %in% 3743 & plot %in% 6 & survey.date %in% as.Date('2023-06-08'), 0, no.eaten),
    # plant 3955 - AHHHH THE QUESTION MARKS! DO YOUR JOB AND TAKE ACTUAL DATA
    # okay 1? is probably a dead umbel, then we will add dead records for the other umbels
    no.dead = ifelse(tag %in% 3955 & plot %in% 6 & survey.date %in% as.Date('2023-05-15'), 1, no.dead),
    no.dead = ifelse(tag %in% 3955 & plot %in% 6 & survey.date %in% as.Date('2023-05-25'), 2, no.dead),
    no.dead = ifelse(tag %in% 3955 & plot %in% 6 & survey.date %in% as.Date('2023-06-04'), 3, no.dead),
    # plant 3069 - "yellow and wilted" okay so is it dead or not???? holy shit what are we doing here???
    # yellow means at least it flowered, which we don't have a record for. maybe dying but not dead
    no.flowers = ifelse(tag %in% 3069 & plot %in% 7 & survey.date %in% as.Date('2023-05-15'), 1, no.flowers),
    # plant 5032 - makes sense if two - data sheet has two with someone scratching a 1 over it
    no.flowers = ifelse(tag %in% 5032 & plot %in% 7 & survey.date %in% as.Date('2023-05-04'), 2, no.flowers),
    # plant 5049 (tag erroneously entered as 5849) - stalk broken written in notes, not entered...
    no.broken = ifelse(tag %in% 5849 & plot %in% 7 & survey.date %in% as.Date('2023-05-24'), 2, no.broken),
    # plant 3422 - stem broken, but lost umbels not recorded...
    no.broken = ifelse(tag %in% 3422 & plot %in% 12 & survey.date > as.Date('2023-05-10'), 1, no.broken),
    # plant 3536 - assume that the bud died and was not recorded
    no.dead = ifelse(tag %in% 3536 & plot %in% 12 & survey.date %in% as.Date('2023-05-25'), 1, no.dead),
    # plant 3336 - a three was mis-recorded as a two
    no.seeding = ifelse(tag %in% 3336 & plot %in% 13 & survey.date %in% as.Date('2023-05-25'), 3, no.seeding),
    # plant 3390 - there are definitely three umbels listed on 5/15... assume one died and wasn't recorded
    no.dead = ifelse(tag %in% 3390 & plot %in% 13 & survey.date %in% as.Date('2023-05-25'), 1, no.dead),
    # plant 3978 - a one was mistakenly entered as a three
    no.buds = ifelse(tag %in% 3978 & plot %in% 13 & survey.date %in% as.Date('2023-05-15'), 1, no.buds),
    # plant 3980 - hated question mark, but in context it has to be a bud
    no.buds = ifelse(tag %in% 3980 & plot %in% 13 & survey.date %in% as.Date('2023-05-15'), 1, no.buds),
    # plant 3382 - another question mark... but no evidence of this umbel ever again
    # going to just say that it died
    no.dead = ifelse(tag %in% 3382 & plot %in% 14 & survey.date %in% as.Date('2023-05-15'), 1, no.dead),
    # plant 3354 - bud recorded on datasheet on 15 may... but not listed as dead... ughh....
    # I guess I should add the dead umbel.
    no.dead = ifelse(tag %in% 3354 & plot %in% 15 & survey.date > as.Date('2023-05-15'), 5, no.dead),
    # plant 3377 - 1 was mis-entered as a 2
    no.flp = ifelse(tag %in% 3377 & plot %in% 15 & survey.date %in% as.Date('2023-05-26'), 1, no.flp),
    # plant 3401 - again, a case where a bud was recorded and never seen again
    # assume it died and just wasn't counted? ugh
    no.dead = ifelse(tag %in% 3401 & plot %in% 15 & survey.date %in% as.Date('2023-05-26'), 1, no.dead),
    no.dead = ifelse(tag %in% 3401 & plot %in% 15 & survey.date %in% as.Date('2023-06-04'), 2, no.dead),
    # plant 3482 - another one of these stupid fucking buds. god this data is such shit
    no.dead = ifelse(tag %in% 3482 & plot %in% 15 & survey.date %in% as.Date('2023-05-26'), 1, no.dead),
    no.dead = ifelse(tag %in% 3482 & plot %in% 15 & survey.date %in% as.Date('2023-06-04'), 2, no.dead),
    # plant 3706 (tag mis-entered as 3076, corrected later) - another one of these buds
    # PEOPLE, IF YOU SEE SOMETHING ONE WEEK, CHECK TO SEE IF IT'S THERE THE NEXT WEEK
    # YOU ARE CAUSING SO MANY PROBLEMS
    no.dead = ifelse(tag %in% 3076 & plot %in% 15 & survey.date > as.Date('2023-05-15'), 1, no.dead)
  )

exclude.plantids = exclude.plantids %>%
  # Exclude 3481_15 - stem got clipped, records are hard to parse out...
  mutate(exclude = ifelse(grepl('3481\\_15', plantid) & year %in% 2023, TRUE, exclude)) %>%
  # exclude 3687_1 - missed one week making flowering date uncertain
  mutate(exclude = ifelse(plantid %in% '3867_1' & year %in% 2023, TRUE, exclude)) %>%
  # exclude 3338_5 - STOP WRITING QUESTION MARKS, KEEP TRACK OF OBSERVED UMBELS
  mutate(exclude = ifelse(plantid %in% '3338_5' & year %in% 2023, TRUE, exclude)) %>%
  # exclude 3987_5 - went missing one week, came back with new stuff
  mutate(exclude = ifelse(plantid %in% '3987_5' & year %in% 2023, TRUE, exclude))

### 2024

# Survey periods:
# Just set the survey date equal to the first day of each period
# (scott has a spreadsheet with survey dates for each plot locally stored - ask
# him for data, not put in repo)

proc24 = proc24 %>%
  mutate(
    survey.date = as.Date(survey.date, format = '%m/%d/%y'),
    survey.date = survey.date - as.numeric(
      survey.date %in% as.Date(paste0('2024-', c('04-20', '04-26', '05-02', '05-09', '05-17')))
    )
  )

unique(proc24$survey.date)

# I'll do some data cleaning in here too...

# Plants that had a decrease in umbel count over time
# (many of these won't need fixing because they're dead umbels or lost plants)
# proc24 %>% 
#   group_by(plantid) %>%
#   mutate(n.umbel = no.buds + no.flowers + no.flp + no.seeding + no.dead) %>% 
#   filter(any(diff(n.umbel) < 0)) %>%
#   print(n = nrow(.))
  
proc24 %>%
  mutate(
    # Plant 3168, plot 2: extra umbel was recorded here
    # I think this should have been for plant 3751, which has two umbels in seed set
    # transfer umbels over
    # (just add buds to each date - will preserve bud date which we will use in analysis)
    no.buds = ifelse(grepl('3168\\_2', plantid) & survey.date %in% as.Date('2024-05-01'), 0, no.buds),
    no.buds = ifelse(grepl('3751\\_2', plantid) & survey.date > as.Date('2025-04-30'), 1, no.buds),
    # plant 3546, plot 5: ah this is the plant from above
    # two issues here: a can't find one week, and records for (likely) the same
    # plant under a different tag, etc.
    # first record can be changed by just adding a bud record in the missing
    # week (there was one bud before and after - parsimony, assume it's the
    # same)
    no.buds = ifelse(grepl('3546\\_5', plantid) & survey.date %in% as.Date('2024-04-12'), 1, no.buds),
    # second issue - assume what *was* incorrectly surveyed as 5309 was in fact
    # 3546, which itself was marked as "NT" for the three records where 5309 was
    # found
  ) %>%
  filter(!(grepl('3546\\_5', plantid) & coor %in% '10D' & survey.date > as.Date('2024-04-26'))) %>%
  mutate(
    coor    = ifelse(tag %in% 3546 & coor %in% '9G' & plot %in% 5, '10D', coor),
    plantid = ifelse(grepl('3546\\_5\\_9G', plantid), '3546_5_10D', plantid),
    # RSP/7675 and BSP in plot 5: these are the same plant and I didn't realize it until doing demo
    # wish there was a clever way to combine these... but I will just enter them manually
    # to avoid accidentally making a mistake
    # also - one of these plants goes '1F2B' to '3B' the next weekend
    # I'm going to assume these are the same umbels and they were always buds,
    # because the bud-flower protocol is so finnicky/subjective
    # (found these values by running the following):
    # proc24 %>% 
    #   filter(plot %in% 5, tag %in% c('RSP', 'BSP', 7675)) %>% 
    #   arrange(survey.date) %>% 
    #   group_by(survey.date) %>% 
    #   summarise(across(starts_with('no.'), ~ sum(.)))
    no.buds = case_when(
      grepl('7675\\_5', plantid) & survey.date %in% as.Date('2024-04-19') ~ 2,
      grepl('7675\\_5', plantid) & survey.date %in% as.Date('2024-04-25') ~ 4,
      grepl('7675\\_5', plantid) & survey.date %in% as.Date('2024-05-01') ~ 4,
      .default = no.buds
    ),
    no.flowers = case_when(
      grepl('7675\\_5', plantid) & survey.date %in% as.Date('2024-04-25') ~ 0,
      grepl('7675\\_5', plantid) & survey.date %in% as.Date('2024-05-08') ~ 2,
      .default = no.flowers
    ),
    no.dead = case_when(
      grepl('7675\\_5', plantid) & survey.date %in% as.Date('2024-05-08') ~ 2,
      grepl('7675\\_5', plantid) & survey.date %in% as.Date('2024-05-16') ~ 1,
      .default = no.dead
    ),
    no.seeding = ifelse(grepl('7675\\_5', plantid) & survey.date %in% as.Date('2024-05-16'), 2, no.seeding)
  ) %>%
  filter(!(plantid %in% 'BSP_5_10J')) %>%
  mutate(
    # Plant 3488, plot 5: I think one week I was just looking at the wrong plant
    # (or there was a recording mistake)
    # just assume the record was 3FLP (won't influence analysis whether these are
    # recorded as 3F or 3FLP)
    no.buds = ifelse(grepl('3488\\_5', plantid) & survey.date %in% as.Date('2024-05-08'), 0, no.buds),
    no.flp  = ifelse(grepl('3488\\_5', plantid) & survey.date %in% as.Date('2024-05-08'), 3, no.flp),
    # Plant 3994, plot 7: assume the second bud observed on the first day was
    # missed on the second and was observed on the third, i.e., impute another
    # bud on that day
    no.buds = ifelse(grepl('3994\\_7', plantid) & survey.date %in% as.Date('2024-04-25'), 2, no.buds),
    # Plant 3187, plot 13: likewise impute flowering plant
    # (two Flp observed on the next day, likely that we just missed it)
    no.flowers = ifelse(grepl('3187\\_13', plantid) & survey.date %in% as.Date('2024-05-01'), 2, no.flowers),
    # Plant 3796, plot 13: impute budding plant (records before and after have one bud only)
    no.buds = ifelse(grepl('3796\\_13', plantid) & survey.date %in% as.Date('2024-04-25'), 1, no.buds),
    # Plant 3780, plot 13: because we're not interested in flowering dates, I'm
    # comfortable assigning the missing record here to flower (not that it
    # matters which stage it's in, as long as there is a record there)
    no.flowers = ifelse(grepl('3780\\_13', plantid) & survey.date %in% as.Date('2024-05-01'), 1, no.flowers),
    # Plant 3377, plot 15: impute the missing record (call it a bud)
    no.buds = ifelse(grepl('3377\\_15', plantid) & survey.date %in% as.Date('2024-04-19'), 1, no.buds),
    # Plant 3177 in plot 15: I just really don't trust this strange 2B that appears one day
    # remove them
    no.buds = ifelse(grepl('3177\\_15', plantid) & survey.date %in% as.Date('2024-05-01'), 0, no.buds),
    # Plant 3426 in plot 15: assume that one Flp was missed on this date; makes records consistent
    # (doesn't explain the other bud that appeared and went to Flp in a week but c'est la vie)
    no.flp = ifelse(grepl('3426\\_15', plantid) & survey.date %in% as.Date('2024-05-08'), 2, no.flp),
    # Plant 3479 in plot 15: impute bud
    no.buds = ifelse(grepl('3479\\_15', plantid) & survey.date %in% as.Date('2024-04-25'), 1, no.buds)
  )

exclude.plantids = exclude.plantids %>%
  # (Exclude 3643 and 3599 in plot 6 - these got confused with each other, record messy)
  mutate(
    exclude = ifelse(
      (grepl('3643\\_6', plantid) | grepl('3955\\_6', plantid)) & year %in% 2024, 
      TRUE, 
      exclude
      ),
    # Exclude R(tp) in plot 7 - missed twice so two gaps, just kinda messy
    exclude = ifelse(grepl('R\\_7\\_0I', plantid) & year %in% 2024, TRUE, exclude),
    # Exclude 3515 and 3513 in plot 13 - messy, multiple gaps in each
    exclude = ifelse(grepl('351[35]\\_13', plantid) & year %in% 2024, TRUE, exclude),
    # Exclude 3590 in plot 13 because it was missed on a day that makes
    # assigning bud+flower dates ambiguous
    exclude = ifelse(grepl('3590\\_13', plantid) & year %in% 2024, TRUE, exclude),
    # exclude 3455 in plot 15, I just don't trust that record that has three
    # flowers appearing out of nowhere and then disappearing the next week
    exclude = ifelse(grepl('3455\\_15', plantid) & year %in% 2024, TRUE, exclude),
  )

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

### Check 2024
proc24 %>%
  group_by(plantid, survey.date) %>%
  filter(n() > 1)
# none - good

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

### 2024

# this is fine
# pre-april 12 surveys were done by piper battersby
# 24 may surveys were of only plants with yellow visible (to catch new buds)
# so only a handful are included

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

# # Just this one plant - I say remove it because that missing record is *right*
# # where the partial date should be
# exclude.plantids = exclude.plantids %>%
#   mutate(exclude = ifelse(plantid %in% '3563_4_12S' & year %in% 2021 , TRUE, exclude))
# Actually - going to leave this in here because I'll be doing analysis of bud
# date, not flower date

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

### 2024

proc24 %>%
  group_by(plantid) %>%
  summarise(any.fl = any(no.buds > 0 | no.flowers > 0 | no.flp > 0 | no.seeding > 0)) %>%
  group_by(any.fl) %>%
  summarise(n = n())
# three plants that don't have any records of budding, flowering, etc.

proc24 %>%
  group_by(plantid) %>%
  filter(!any(no.buds > 0 | no.flowers > 0 | no.flp > 0 | no.seeding > 0))

# Oh yes... I know what's going on here.
# These are plants that had records that I crossed out on the datasheet.
# I'll remove them.

proc24 = proc24 %>% filter(!grepl('crossed\\sout', notes))


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
      max(survey.period[no.buds + no.flowers + no.flp > 0]),
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

# Remove some plants that have dead umbels appear very late?
# exclude.plantids = exclude.plantids %>%
#   mutate(exclude = ifelse(plantid %in% c('5037_7', '5040_7', '3495_14'), TRUE, exclude))

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
      max(survey.period[no.buds + no.flowers + no.flp > 0]),
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

### 2024

# (no need to change NAs to zeros because of data entry scheme)
# (also no need to filter out plants that only have budding etc. records because
# of the data entry)

# but to make this consistent with the other years, I will want a survey.period column
# (not sure why I did this...)

proc24 = proc24 %>%
  # set 'date.2024' as julian date but with day 0 being sunday 31 dec. 2024
  mutate(
    date.2024 = as.numeric(survey.date - as.Date('2023-12-31')),
    # convert day to week
    survey.period = date.2024 %/% 7,
    # first observation now will be survey period 1
    survey.period = survey.period - min(survey.period) + 1
  ) %>%
  select(-date.2024)

phen24 = proc24 %>%
  group_by(plantid, plot, tag) %>%
  summarise(
    # First date with observed budding activity
    fl.init = min(survey.period[no.buds + no.flowers + no.flp + no.seeding > 0]),
    fl.fina = ifelse(
      # Hmm.. the flowering flag here may not be as useful here due to the bud
      # classification protocol in 2024
      # also not counting seeding plants would cause an issue no?
      any(no.flowers + no.flp > 0),
      max(survey.period[no.buds + no.flowers + no.flp > 0]),
      NA
    )
  )

# 3341_15... want an algo that handles this


##### Combine all data frames together (remember to add year)

phen.all = rbind(
  phen21 %>% ungroup() %>% mutate(year = 2021),
  phen22 %>% ungroup() %>% mutate(year = 2022),
  phen23 %>% ungroup() %>% mutate(year = 2023),
  phen24 %>% ungroup() %>% mutate(year = 2024)
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
    year = 2021:2024,
    min.date = c(
      as.numeric(min(proc21$survey.date) - as.Date('2020-12-31')),
      as.numeric(min(proc22$survey.date) - as.Date('2021-12-31')),
      as.numeric(min(proc23$survey.date) - as.Date('2022-12-31')),
      as.numeric(min(proc24$survey.date) - as.Date('2023-12-31'))
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

# write.csv(
#   phen.all,
#   '01_data_cleaning/out/phenology_all_cleaned.csv',
#   row.names = FALSE
# )
# 
# write.csv(
#   ind.phen.all,
#   '01_data_cleaning/out/phenology_all_ind_cleaned.csv',
#   row.names = FALSE
# )


#####
# Here: get mean date of *bud* initiation
# Also get number dead/eaten

# 2021 data
buds21 = proc21 %>%
  arrange(survey.date) %>%
  mutate(no.umbels = no.closed + no.unfurling + no.partial + no.flat + no.open + no.olp + no.seeding) %>%
  group_by(plantid, plot) %>%
  mutate(
    live.diff = diff(c(0, no.umbels)),
    lost.diff = diff(c(0, no.dead + no.eaten))
  ) %>%
  filter(!is.na(no.umbels)) %>%
  summarise(
    mean.bud.window = mean(ifelse(live.diff < 0, 0, live.diff) * survey.period),
    n.lost = max(no.eaten + no.dead),
    n.total = max(no.umbels + no.eaten + no.dead)
  ) %>%
  ungroup()

nrow(buds21)

# 2022 data
buds22 = proc22 %>%
  mutate(across(starts_with('no.'), function(x) ifelse(is.na(x), 0, x))) %>%
  arrange(plantid, survey.date) %>%
  # note: NOT counting pods
  mutate(no.umbels = no.buds + no.flowers + no.flp + no.seeding) %>%
  group_by(plantid, plot) %>%
  mutate(
    live.diff = diff(c(0, no.umbels)),
    lost.diff = diff(c(0, no.dead + no.eaten + no.broken))
  ) %>%
  summarise(
    mean.bud.window = mean(ifelse(live.diff < 0, 0, live.diff) * survey.period),
    n.lost = max(no.eaten + no.dead + no.broken),
    n.total = max(no.umbels + no.broken +  no.eaten + no.dead)
  ) %>%
  ungroup()

# 2023 data

buds23 = proc23 %>%
  mutate(across(starts_with('no.'), function(x) ifelse(is.na(x), 0, x))) %>%
  arrange(plantid, survey.date) %>%
  # note: NOT counting pods
  mutate(no.umbels = no.buds + no.flowers + no.flp + no.seeding) %>%
  group_by(plantid, plot) %>%
  mutate(
    live.diff = diff(c(0, no.umbels)),
    lost.diff = diff(c(0, no.dead + no.eaten + no.broken))
  ) %>%
  summarise(
    mean.bud.window = mean(ifelse(live.diff < 0, 0, live.diff) * survey.period),
    n.lost = max(no.eaten + no.dead + no.broken),
    n.total = max(no.umbels + no.broken +  no.eaten + no.dead)
  ) %>%
  ungroup()

buds24 = proc24 %>%
  # I don't think this line is needed but I'll keep it anyway
  mutate(across(starts_with('no.'), function(x) ifelse(is.na(x), 0, x))) %>%
  arrange(plantid, survey.date) %>%
  mutate(no.umbels = no.buds + no.flowers + no.flp + no.seeding) %>%
  group_by(plantid, plot) %>%
  mutate(
    live.diff = diff(c(0, no.umbels)),
    lost.diff = diff(c(0, no.dead))
  ) %>%
  summarise(
    mean.bud.window = mean(ifelse(live.diff < 0, 0, live.diff) * survey.period),
    n.lost = max(no.dead),
    n.total = max(no.umbels)
  ) %>%
  ungroup()

# Bind it all together

buds.all = rbind(
  # Bind everything together
  buds21 %>% mutate(year = 2021),
  buds22 %>% mutate(year = 2022),
  buds23 %>% mutate(year = 2023),
  buds24 %>% mutate(year = 2024)
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
      year = 2021:2024,
      min.date = c(
        as.numeric(min(proc21$survey.date) - as.Date('2020-12-31')),
        as.numeric(min(proc22$survey.date) - as.Date('2021-12-31')),
        as.numeric(min(proc23$survey.date) - as.Date('2022-12-31')),
        as.numeric(min(proc24$survey.date) - as.Date('2023-12-31'))
      )
    )
  ) %>%
  mutate(init.doy = min.date + (mean.bud.window - 1) * 7) %>%
  select(-min.date)

# oh wait I want individual level budding

ind.bud21 = proc21 %>%
  arrange(survey.date) %>%
  # Convert NAs to zeros %>%
  mutate(across(starts_with('no.'), function(x) ifelse(is.na(x), 0, x))) %>%
  # get number of non-dead umbels per plant per day
  mutate(no.umbels = no.closed + no.unfurling + no.partial + no.flat + no.open + no.olp + no.seeding) %>%
  group_by(plantid, plot) %>%
  mutate(
    new = diff(c(0, no.umbels + no.dead + no.eaten)),
    lost = diff(c(0, no.dead + no.eaten))
  ) %>%
  select(plantid, plot, survey.date, survey.period, new, lost) %>%
  pivot_longer(c(new, lost), names_to = 'varb', values_to = 'count') %>%
  # change negative numbers to zeros for uncounting
  mutate(count = ifelse(count < 0, 0, count)) %>%
  uncount(weight = count) %>%
  arrange(plantid)

# looks good to me

ind.bud22 = proc22 %>%
  mutate(across(starts_with('no.'), function(x) ifelse(is.na(x), 0, x))) %>%
  arrange(plantid, survey.date) %>%
  # NOTE: pods removed
  mutate(no.umbels = no.buds + no.flowers + no.flp + no.seeding) %>%
  group_by(plantid, plot) %>%
  mutate(
    new = diff(c(0, no.umbels + no.dead + no.eaten + no.broken)),
    lost = diff(c(0, no.dead + no.eaten + no.broken))
  ) %>%
  select(plantid, plot, survey.date, survey.period, new, lost) %>%
  pivot_longer(c(new, lost), names_to = 'varb', values_to = 'count') %>%
  # change negative numbers to zeros for uncounting
  mutate(count = ifelse(count < 0, 0, count)) %>%
  uncount(weight = count) %>%
  arrange(plantid)

ind.bud23 = proc23 %>%
  mutate(across(starts_with('no.'), function(x) ifelse(is.na(x), 0, x))) %>%
  arrange(plantid, survey.date) %>%
  # NOTE: pods removed
  mutate(no.umbels = no.buds + no.flowers + no.flp + no.seeding) %>%
  group_by(plantid, plot) %>%
  mutate(
    new = diff(c(0, no.umbels + no.dead + no.eaten + no.broken)),
    lost = diff(c(0, no.dead + no.eaten + no.broken))
  ) %>%
  select(plantid, plot, survey.date, survey.period, new, lost) %>%
  pivot_longer(c(new, lost), names_to = 'varb', values_to = 'count') %>%
  # change negative numbers to zeros for uncounting
  mutate(count = ifelse(count < 0, 0, count)) %>%
  uncount(weight = count) %>%
  arrange(plantid)

ind.bud24 = proc24 %>%
  mutate(across(starts_with('no.'), function(x) ifelse(is.na(x), 0, x))) %>%
  arrange(plantid, survey.date) %>%
  mutate(no.umbels = no.buds + no.flowers + no.flp + no.seeding) %>%
  group_by(plantid, plot) %>%
  mutate(
    new = diff(c(0, no.umbels)),
    lost = diff(c(0, no.dead))
  ) %>%
  select(plantid, plot, survey.date, survey.period, new, lost) %>%
  pivot_longer(c(new, lost), names_to = 'varb', values_to = 'count') %>%
  # change negative numbers to zeros for uncounting
  mutate(count = ifelse(count < 0, 0, count)) %>%
  uncount(weight = count) %>%
  arrange(plantid)

ind.buds.all = rbind(
  ind.bud21 %>% mutate(year = 2021),
  ind.bud22 %>% mutate(year = 2022),
  ind.bud23 %>% mutate(year = 2023),
  ind.bud24 %>% mutate(year = 2024)
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
      year = 2021:2024,
      min.date = c(
        as.numeric(min(proc21$survey.date) - as.Date('2020-12-31')),
        as.numeric(min(proc22$survey.date) - as.Date('2021-12-31')),
        as.numeric(min(proc23$survey.date) - as.Date('2022-12-31')),
        as.numeric(min(proc24$survey.date) - as.Date('2023-12-31'))
      )
    )
  ) %>%
  mutate(init.doy = min.date + (survey.period - 1) * 7) %>%
  select(-min.date)

table(ind.buds.all$init.doy)

# write.csv(
#   ind.buds.all,
#   file = '01_data_cleaning/out/phenology_buds_deaths_all.csv',
#   row.names = FALSE
# )

# NOTE:
# 24 jan 2024 - went back and did some edits but *only re-exported the bud dataset*
# I haven't exported the flowering phen dataset with the updated data yet

#-------------------------------------------------------
# Phen overlap on each day?

flw.by.day21 = proc21 %>%
  select(survey.period, plot, tag, no.flat, no.open, no.olp, plantid) %>%
  mutate(across(c(no.open, no.flat, no.olp), function(x) ifelse(is.na(x), 0, x))) %>%
  mutate(no.total = no.open + no.flat + no.olp) %>%
  filter(no.total > 0) %>%
  select(survey.period, plot, tag, no.total, plantid)

flw.by.day22 = proc22 %>%
  select(survey.period, plot, tag, no.flowers, no.flp, plantid) %>%
  mutate(across(c(no.flowers, no.flp), function(x) ifelse(is.na(x), 0, x))) %>%
  mutate(no.total = no.flowers + no.flp) %>%
  filter(no.total > 0) %>%
  select(survey.period, plot, tag, no.total, plantid)

flw.by.day23 = proc23 %>%
  select(survey.period, plot, tag, no.flowers, no.flp, plantid) %>%
  mutate(across(c(no.flowers, no.flp), function(x) ifelse(is.na(x), 0, x))) %>%
  mutate(no.total = no.flowers + no.flp) %>%
  filter(no.total > 0) %>%
  select(survey.period, plot, tag, no.total, plantid)

flw.by.day = rbind(
  flw.by.day21 %>% mutate(year = 2021),
  flw.by.day22 %>% mutate(year = 2022),
  flw.by.day23 %>% mutate(year = 2023)
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

head(flw.by.day)
nrow(flw.by.day)

# write.csv(
#   flw.by.day,
#   row.names = FALSE,
#   file = '01_data_cleaning/out/phen_flowers_open_by_survey.csv'
# )
