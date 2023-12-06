### Script for imputing the survival of missed plants
### and splitting plants with large gaps in records
### These will be used in survival models and models for flowering probability.
### SN - 27 Nov. 2023

library(dplyr)
library(tidyr)

rm(list = ls())

demo = merge(
  x = read.csv('01_data_cleaning/out/demo_postcombine.csv'),
  y = read.csv('00_raw_data/plot_treatments.csv') %>% rename(Plot = plot)
) %>%
  # Get rid of unnecessary columns
  select(-c(Xcoor, Ycoor, YrTag, Tag, Plot)) %>%
  # select(-c(Xcoor, Ycoor, YrTag)) %>%
  arrange(plantid, Year)

head(demo)
nrow(demo)

##### Add "observed alive" column 

demo = demo %>%
  mutate(
    obs.alive = case_when(
      # Plant was observed alive if there was a positive number of umbels
      (No.umbels > 0) & !is.na(No.umbels) ~ TRUE,
      # Plant was definitely not observed alive if recorded as zero leaves
      # (there are three plants that were observed post grazing - zero recorded
      # leaves, but had umbels recorded - these should be caught by above
      # statement)
      (!No.leaves) & !is.na(No.leaves) ~ FALSE,
      # Plant was observed alive if a non-zero number of leaves were recorded
      No.leaves > 0 & !is.na(No.leaves) ~ TRUE,
      # Otherwise, we assume no leaf measurements means nothing there
      is.na(No.leaves) ~ FALSE,
      # Anything not caught by these flags will get NAs
      .default = NA
    )
  )

demo %>% group_by(obs.alive) %>% summarise(n = n())
# no NAs - good

head(demo)

##### Look at dead plants

dead.maybe = demo %>%
  group_by(plantid) %>%
  filter(any(!obs.alive)) %>%
  ungroup()

nrow(dead.maybe)

# How many plants have any records where there are dead records before
# the *latest* alive record
# This will be conservative (plants might be missed in multiple non-contiguous
# stretches)
dead.maybe %>%
  group_by(plantid) %>%
  summarise(bad = any(Year[!obs.alive] < max(Year[obs.alive]))) %>%
  group_by(bad) %>%
  summarise(n = n())
# A lot of cases of missed plants

# Take a look at some of these
dead.maybe %>%
  group_by(plantid) %>%
  filter(any(Year[!obs.alive] < max(Year[obs.alive])))

# Crop out any records before the first observed mortality
dead.maybe %>%
  group_by(plantid) %>%
  filter(any(Year[!obs.alive] < max(Year[obs.alive]))) %>%
  filter(Year >= min(Year[!obs.alive])) %>%
  summarise(gap.size = min(Year[obs.alive]) - min(Year)) %>%
  group_by(gap.size) %>%
  summarise(n = n())

dead.maybe %>%
  group_by(plantid) %>%
  filter(any(Year[!obs.alive] < max(Year[obs.alive]))) %>%
  filter(Year >= min(Year[!obs.alive])) %>%
  summarise(y0 = min(Year), y1 = min(Year[obs.alive])) %>%
  mutate(
    gap.len = y1 - y0,
    in.2020 = y1 >= 2020 & y0 <= 2020
  ) %>%
  group_by(gap.len, in.2020) %>%
  summarise(n = n())
# A lot of these overlap with 2020, where sampling was less thorough.

dead.maybe %>%
  group_by(plantid) %>%
  filter(any(Year[!obs.alive] < max(Year[obs.alive]))) %>%
  filter(Year >= min(Year[!obs.alive])) %>%
  filter(min(Year[obs.alive]) - min(Year) > 1, !(min(Year) < 2020 & min(Year[obs.alive] > 2020)))
# looking at a few of these (not every single one), they look like they could be new plants

# How many of these do we have where the plant is missed more than once, and
# this miss does not include 2020?
dead.maybe %>%
  group_by(plantid) %>%
  filter(any(Year[!obs.alive] < max(Year[obs.alive]))) %>%
  filter(Year >= min(Year[!obs.alive])) %>%
  filter(min(Year[obs.alive]) - min(Year) > 1, !(min(Year) < 2020 & min(Year[obs.alive] > 2020))) %>%
  distinct(plantid) %>%
  nrow()
# 31 plants

# Try to do this in one line...

demo2 = rbind(
  # All alive plants
  demo %>% 
    group_by(plantid) %>% 
    filter(all(obs.alive)) %>%
    ungroup() %>%
    mutate(surv = obs.alive),
  # Dead plants, but without any resurrections
  demo %>%
    group_by(plantid) %>%
    filter(any(!obs.alive)) %>%
    filter(!any(Year[obs.alive] > min(Year[!obs.alive]))) %>%
    ungroup() %>%
    mutate(surv = obs.alive),
  # Dead plants, where gap includes 2020 or is only one year
  # IMPUTE survival here
  demo %>%
    group_by(plantid) %>%
    filter(any(!obs.alive)) %>%
    filter(any(Year[obs.alive] > min(Year[!obs.alive]))) %>%
    mutate(
      # First year observed not alive
      y0 = min(Year[!obs.alive]),
      # First year observed alive after death (resurrected)
      y1 = min(Year[obs.alive & Year > y0]),
      # Gap size (number of years between y0 and y1)
      gapsize = y1 - y0,
      # Numeric for if record gap includes 2020
      in.2020 = as.numeric(y1 > 2020 & y0 <= 2020)
    ) %>%
    ungroup() %>%
    # Separate out plantid into tag, plot, coords
    separate(plantid, into = c('Tag', 'Plot', 'Coord'), sep = '_') %>%
    mutate(
      # add 'b' to tag if gap 
      Tag = ifelse((gapsize > (1 + in.2020)) & (Year >= y1), 
                   paste0(Tag, 'b'), 
                   Tag),
      # Create 'alive' column sensu above
      surv = ifelse((gapsize <= (1 + in.2020) & (Year < y1 & Year >= y0)),
                    TRUE,
                    obs.alive)
    ) %>%
    # Delete old columns
    select(-c(y0, y1, gapsize, in.2020)) %>%
    # Add notes
    mutate(
      proc.note = ifelse(grepl('b', Tag), 
                         paste0(proc.note, '; tag modified for gap'),
                         proc.note),
      proc.note = ifelse(surv != obs.alive,
                         paste0(proc.note, '; survival imputed based on subsq. record'),
                         proc.note),
      edited = ifelse(grepl('b', Tag) | surv != obs.alive, TRUE, edited)
    ) %>%
    unite(Tag, Plot, Coord, col = 'plantid', sep = '_') %>%
    select(c(names(demo), surv))
) %>%
  arrange(plantid, Year)

nrow(demo2)
head(demo2)
length(unique(demo2$plantid))

# Remaining plants with gaps?
demo2 %>% group_by(plantid) %>% filter(any(diff(surv) > 0)) %>% distinct(plantid) %>% nrow()
# Still a considerable number - 

# A while loop would work great for this... but is it worth doing... probably not.
# (I checked - there are no plants with gaps large enough to need tag additions)

### Pass 1
demo2 = rbind(
  # All alive plants
  demo2 %>% 
    group_by(plantid) %>% 
    filter(all(surv)) %>%
    ungroup(),
  # Dead plants, but without any resurrections
  demo2 %>%
    group_by(plantid) %>%
    filter(any(!surv)) %>%
    filter(!any(Year[surv] > min(Year[!surv]))) %>%
    ungroup(),
  # Dead plants, where gap includes 2020 or is only one year
  # IMPUTE survival here
  demo2 %>%
    group_by(plantid) %>%
    filter(any(!surv)) %>%
    filter(any(Year[surv] > min(Year[!surv]))) %>%
    mutate(
      # First year observed not alive
      y0 = min(Year[!surv]),
      # First year observed alive after death (resurrected)
      y1 = min(Year[surv & Year > y0]),
      # Gap size (number of years between y0 and y1)
      gapsize = y1 - y0,
      # Numeric for if record gap includes 2020
      in.2020 = as.numeric(y1 > 2020 & y0 <= 2020)
    ) %>%
    ungroup() %>%
    # Separate out plantid into tag, plot, coords
    mutate(
      # Make notes and add edited flag
      proc.note = ifelse((gapsize <= (1 + in.2020) & (Year < y1 & Year >= y0)),
                         paste0(proc.note, '; survival imputed based on subsq. record'),
                         proc.note),
      edited = ifelse((gapsize <= (1 + in.2020) & (Year < y1 & Year >= y0)), TRUE, edited),
      # Modify 'surv' column to impute survival
      surv = ifelse((gapsize <= (1 + in.2020) & (Year < y1 & Year >= y0)),
                    TRUE,
                    surv)
    ) %>%
    # Delete old columns
    select(-c(y0, y1, gapsize, in.2020)) %>%
    # Add notes
    select(names(demo2))
) %>%
  arrange(plantid, Year)

### Pass 2
demo2 = rbind(
  # All alive plants
  demo2 %>% 
    group_by(plantid) %>% 
    filter(all(surv)) %>%
    ungroup(),
  # Dead plants, but without any resurrections
  demo2 %>%
    group_by(plantid) %>%
    filter(any(!surv)) %>%
    filter(!any(Year[surv] > min(Year[!surv]))) %>%
    ungroup(),
  # Dead plants, where gap includes 2020 or is only one year
  # IMPUTE survival here
  demo2 %>%
    group_by(plantid) %>%
    filter(any(!surv)) %>%
    filter(any(Year[surv] > min(Year[!surv]))) %>%
    mutate(
      # First year observed not alive
      y0 = min(Year[!surv]),
      # First year observed alive after death (resurrected)
      y1 = min(Year[surv & Year > y0]),
      # Gap size (number of years between y0 and y1)
      gapsize = y1 - y0,
      # Numeric for if record gap includes 2020
      in.2020 = as.numeric(y1 > 2020 & y0 <= 2020)
    ) %>%
    ungroup() %>%
    # Separate out plantid into tag, plot, coords
    mutate(
      # Make notes and add edited flag
      proc.note = ifelse((gapsize <= (1 + in.2020) & (Year < y1 & Year >= y0)),
                         paste0(proc.note, '; survival imputed based on subsq. record'),
                         proc.note),
      edited = ifelse((gapsize <= (1 + in.2020) & (Year < y1 & Year >= y0)), TRUE, edited),
      # Modify 'surv' column to impute survival
      surv = ifelse((gapsize <= (1 + in.2020) & (Year < y1 & Year >= y0)),
                    TRUE,
                    surv)
    ) %>%
    # Delete old columns
    select(-c(y0, y1, gapsize, in.2020)) %>%
    # Add notes
    select(names(demo2))
) %>%
  arrange(plantid, Year)

demo2 %>% group_by(plantid) %>% filter(any(diff(surv) > 0)) %>% distinct(plantid) %>% nrow()
# Great - everything has been imputed

# I want to add Plot as a covariate - scrape that out of the plantid

demo2 = demo2 %>%
  separate(col = plantid, into = c("Tag", "Plot", "Coord"), sep = "_", remove = FALSE) %>%
  select(-c(Tag, Coord))

# Now, add in records for the three plants (identified above) that are missing entire records
# (I went through these case-by-case in the combining script)

demo.new = demo2 %>% 
  slice(1:4) %>%
  mutate(across(everything(), function(x) NA)) %>%
  mutate(
    Year = c(2017, 2017, 2019, 2020),
    proc.note = c(
      'plant missed in demo due to tag mishaps; surv imputed',
      'plant missed in demo due to tag mishaps; surv imputed',
      'record missed in data entry; was NA; surv imputed',
      'record missed due to prior data entry mishap; surv imputed'
    ),
    edited = TRUE,
    plantid = c('3585_2_3C', '3844_2_8C', '3699_14_0E', '3699_14_0E'),
    Plot = c(2, 2, 14, 14),
    trt = c("control", "control", "drought", "drought"),
    obs.alive = FALSE,
    surv = TRUE
    )

demo.new

demo2 = rbind(demo2, demo.new) %>% arrange(plantid, Year)

head(demo2)

# Final checks:

# - Any gaps in year (missing records)?
demo2 %>%
  group_by(plantid) %>%
  filter(n() > 1) %>%
  filter(max(diff(Year)) > 1)

# - Any misses in survival?
demo2 %>%
  group_by(plantid) %>%
  filter(n() > 1) %>%
  filter(any(diff(surv) > 0))

# - Any missing plot, treatment, or survival?
demo2 %>% 
  select(Plot, trt, surv) %>%
  apply(2, function(x) sum(is.na(x)))

##### Export this sucker

write.csv(
  demo2,
  row.names = FALSE,
  file = '01_data_cleaning/out/demo_imputed_survival.csv'
)
