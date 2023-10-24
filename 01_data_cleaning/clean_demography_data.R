######################################################################
# Script for cleaning and processing raw LoUt demography data for analysis
# initial version focusing only on <YEAR>_Lomatium_Demography_Data.csv files
# (will likely expand beyond that in this script in the future)
######################################################################
# SN - init 10 Oct 2023
# (24 Oct 2023 - re-ran with 2020 data, added 'edited' col)

library(dplyr)
library(tidyr)

# Considerations:

# Missing values for identifying info (Xcoor, Ycoor, Plot, Tag)
# Duplicate tags (even within plot)

rm(list = ls())

# Read in all documents
# (keep the raw data in list form - one list entry for each data frame)
# (keeping in list because columns differ from year to year)
raw.demo.list = dir('00_raw_data/lomatium_demography/') %>%
  grep(pattern = '20\\d\\d\\_Lomatium\\_Demography_Data', x = ., value = TRUE) %>%
  paste0('00_raw_data/lomatium_demography/', .) %>%
  lapply(FUN = read.csv)

# Check for column consistency
sapply(raw.demo.list, ncol)
# some inconsistencies

# Look at column names
lapply(raw.demo.list, names)

# Get unique columns
# (found this code snippet online)
shared.cols = Reduce(intersect, lapply(raw.demo.list, names))

# Combine common columns into one data frame
raw.demo = raw.demo.list %>%
  lapply(function(df) df[,shared.cols]) %>%
  do.call(what = rbind)

head(raw.demo)

# Create a proc.demo data frame for storing processed demo
proc.demo = raw.demo %>% mutate(Note = NA, edited = FALSE)

######################################################################
##### Need unique identifiers for each plant
######################################################################

# Problem... some tags are recycled, even sometimes in the same plot!
# need to use coordinates
# but note that coordinates are not necessarily identical across years!
# assume here that coordinates will only be off by 1 in either direction from the "true" value...

# Check for missing identifying info
proc.demo %>% filter(is.na(Tag) | is.na(Plot) | is.na(Xcoor) | is.na(Ycoor))
# Not entirely sure this needs fixing...
proc.demo %>% filter(Tag %in% 3044)
# It's in other records as "D"
# Just fix this manually
proc.demo = proc.demo %>%
  mutate(
    Ycoor = ifelse(Tag %in% 3044 & Plot %in% 7, 'D', Ycoor),
    edited = ifelse(Tag %in% 3044 & Plot %in% 7, TRUE, edited),
    Note = ifelse(
      Tag %in% 3044 & Plot %in% 7,
      "one ycoor fixed (was NA)",
      Note
    )
  )

# Now - need to figure out ways to handle off-by-one coordinate changes

proc.demo %>%
  distinct(Plot, Tag, Xcoor, Ycoor) %>%
  group_by(Plot, Tag) %>%
  summarise(n = n()) %>%
  group_by(n) %>%
  summarise(nn = n())
# cool - looks like no tag (w/in a plot) has more than two coordinates listed!

proc.demo %>%
  distinct(Plot, Tag, Xcoor, Ycoor) %>%
  # Do a merge in here to assign the letter-based coordinates to numbers
  # (A = 1, B = 2, ...)
  merge(y = data.frame(Ycoor = LETTERS, Ycoor.num = 1:length(LETTERS))) %>%
  group_by(Plot, Tag) %>%
  filter(n() > 1) %>%
  arrange(Plot, Tag)
# could do something fancy with euclidian distance
# but honestly I think manhattan distance would work just fine here

proc.demo %>%
  distinct(Plot, Tag, Xcoor, Ycoor) %>%
  # Do a merge in here to assign the letter-based coordinates to numbers
  # (A = 1, B = 2, ...)
  merge(y = data.frame(Ycoor = LETTERS, Ycoor.num = 1:length(LETTERS))) %>%
  group_by(Plot, Tag) %>%
  filter(n() > 1) %>%
  mutate(xd = diff(Xcoor), yd = diff(Ycoor.num), manhattan = abs(xd) + abs(yd))

# Seems like...
# - if the distance is 1 or 2, they are the same (only off by 1)
# - there are also some cases where there may have been a data entry mistake? e.g., mistaking C and G
#   (but in this case, we'd expect there to be a change of _only_ Y, not x (or
#   at least not by more than 1))
# Some of these I think should be done manually
# I guess it's also worth noting if the tag shows up multiple times in the same
# year

proc.demo %>%
  group_by(Plot, Tag, Year) %>%
  filter(n() > 1) %>%
  arrange(Plot, Tag, Year)
# okay so plants *can* appear in survey records multiple times in the same year,
# even if they are not listed at two coordinates...

proc.demo %>%
  distinct(Plot, Tag, Year, Xcoor, Ycoor) %>%
  group_by(Plot, Tag, Year) %>%
  filter(n() > 1) %>%
  arrange(Plot, Tag, Year)
# two records in a year means just one record for plant at two different coordinates...

# How many of these duplicate tags within a plot
proc.demo %>%
  group_by(Plot, Tag, Year) %>%
  filter(n() > 1) %>%
  merge(y = data.frame(Ycoor = LETTERS, Ycoord.no = 1:length(LETTERS))) %>%
  group_by(Plot, Tag, Year) %>%
  mutate(xd = abs(diff(Xcoor)), yd = abs(diff(Ycoord.no))) %>%
  mutate(manhatta = xd + yd) %>%
  arrange(Plot, Tag, Year)

# Let's get a list of all possible duplicated cases
proc.demo.dupe = proc.demo %>%
  # Add numeric y coordinate column
  merge(y = data.frame(Ycoor = LETTERS, Ycoor.no = 1:length(LETTERS))) %>%
  group_by(Tag) %>%
  summarise(
    dupe.status = case_when(
      # Only ever appears in one plot at one coordinate
      length(unique(Plot)) == 1 & length(unique(Xcoor)) == 1 & length(unique(Ycoor)) == 1 ~ 
        'Not duplicated',
      # In one plot, with multiple coordinates that are more than 2 off
      length(unique(Plot)) == 1 & (diff(range(Xcoor)) > 2 | diff(range(Ycoor.no)) > 2) ~ 
        'Duplicated in plot, very different coordinates',
      length(unique(Plot)) == 1 & diff(range(Xcoor)) <= 2 & diff(range(Ycoor.no)) <= 2 ~
        'Duplicated in plot, off by one coordinate',
      # Appears in multiple plots
      length(unique(Plot)) > 1 ~ 
        'Duplicate in different plots',
      # catch all for none of the above
      .default = 'other'
    )
  )

# # Multiple coordinates listed with repeat records in same year
# # (these are going to be different plants)
# length(unique(Plot)) > 1 &
#   length(unique(Xcoor)) + length(unique(Ycoor)) > 2*length(unique(Plot)) &
#   any(duplicated(Year)) ~ 'Multiple plants in same year',

proc.demo.dupe %>% group_by(dupe.status) %>% summarise(n = n())
# well... at least we caught all cases

# Guess work through these case by case?

### 3180: 
raw.demo.list %>% lapply(function(df) df %>% filter(Tag %in% 3180))
# two different plants recorded in 2016, but 17B was not seen in 2017
# then two very similar plants observed in 2018, then 18C not seen in 2019,
# then 17B not seen again for duration of study
# and a new distinct plant was observed in plot 13 (...)
# What I think happened:
# two different plants (maybe? not sure) in plot 5
# but one of them might have died, the 2019 record has the wrong coordinates
# (what to do about the 2018 record? ugh)
# The two possibilities for plot 5: it was always only one plant, or it was two but died in year one...

# Here: re-assigning coordinates to the plot 5 year 2019 observations
proc.demo = rbind(
  # Old data:
  proc.demo %>% filter(!(Tag %in% 3180 & Year %in% 2019 & Plot %in% 5)),
  # New data:
  proc.demo %>% 
    filter(Tag %in% 3180 & Year %in% 2019 & Plot %in% 5) %>%
    # Swap coordinates
    # (will assign the "good" record to the proper coordinates)
    mutate(
      Xcoor = rev(Xcoor), 
      Ycoor = rev(Ycoor),
      Note  = "coordinates swapped on 2019 records based on prior records",
      edited = TRUE
    )
)

### 3342: 
raw.demo.list %>% lapply(function(df) df %>% filter(Tag %in% 3342))
# Plot 4 plant is in a diversity plot (does this mean it gets included?)
# In 2016, there were plants observed at 18B and 17B in plot 5
# the one at 17B had "no leaves just an umbel", the one at 18B did not flower
# the one at 18B was never *checked* again, the one at 187B was checked but
# never seen (ever)
# (I suppose we can keep 18B in, it just won't get included in any transitions?)
# (but how can we honestly not think that 18B and 17B are somehow the same plant...)
# (also 17B has no leaf measurements... ahhhhhhhhhhhhh!!!!!)

### 3390
raw.demo.list %>% lapply(function(df) df %>% filter(Tag %in% 3390))
# Different plants

### 3391
raw.demo.list %>% lapply(function(df) df %>% filter(Tag %in% 3391))
# Two different plants (different plots)
# But at some coordinates for one were mis-transcribed

proc.demo = rbind(
  # Old data
  proc.demo %>% filter(!(Tag %in% 3391 & Ycoor %in% 'K')),
  # New data
  proc.demo %>%
    filter(Tag %in% 3391 & Ycoor %in% 'K') %>%
    mutate(
      Ycoor = 'G',
      Note = "y-coordinate changed based on 2023 note (was K)",
      edited = TRUE
    )
)

proc.demo %>% filter(Tag %in% 3391)
# good

### 3396
raw.demo.list %>% lapply(function(df) df %>% filter(Tag %in% 3396))
# some tags were changed or mis-recorded
raw.demo.list %>% lapply(function(df) df %>% filter(Tag %in% 3398))
# the tag that was mis-labeled as 3396 in plot 7 is actually 3398
# (there is also a 3398 in plot 13... but this is a different plant)

proc.demo = rbind(
  # Old data
  proc.demo %>% filter(!(Tag %in% 3396 & Plot %in% 7)),
  # New data
  proc.demo %>%
    filter(Tag %in% 3396 & Plot %in% 7) %>%
    mutate(
      # Correct tag
      Tag = 3398,
      Note = "tag changed based on 2019 note (was 3396)",
      edited = TRUE
    )
)

raw.demo.list %>% lapply(function(df) df %>% filter(Tag %in% 3989))
raw.demo.list %>% lapply(function(df) df %>% filter(Plot %in% 13, Xcoor %in% 17))
# There's a note on the 2023 record saying that 17A is actually tag 3989
# there is no other 3989 in the survey...
# but might be a good idea to make the switch anyway

proc.demo = rbind(
  # Old data
  proc.demo %>% filter(!(Tag %in% 3396 & Plot %in% 13 & Xcoor %in% 17)),
  # New data
  proc.demo %>%
    filter(Tag %in% 3396 & Plot %in% 13 * Xcoor %in% 17) %>%
    mutate(
      # Add new (correct) tag
      Tag = 3989,
      # Add note
      Note = "tag changed based on 2023 note (was 3396)",
      edited = TRUE
    )
)

proc.demo %>% filter(Tag %in% 3396) # good

### Look at the ones with "very different coordinates"

proc.demo.dupe %>% filter(grepl('very different coordinates', dupe.status)) %>% distinct(Tag) %>% arrange(Tag)
# 3090, 3125, 3185, 3356, 3387, 3389, 3555, 3563, 3634

### 3090
raw.demo.list %>% lapply(function(df) df %>% filter(Tag %in% 3090))
# ah this is an obvious one - "F" should be "P"
proc.demo = proc.demo %>% 
  mutate(
    Ycoor = ifelse(Tag %in% 3090 & Plot %in% 4, "P", Ycoor),
    Note  = ifelse(Tag %in% 3090 & Plot %in% 4, "fixed ycoor (was P)", Note),
    edited = ifelse(Tag %in% 3090 & Plot %in% 4, TRUE , edited)
  )

### 3125
raw.demo.list %>% lapply(function(df) df %>% filter(Tag %in% 3125))
# going to say these are in fact different plants - records present in both years

### 3185
raw.demo.list %>% lapply(function(df) df %>% filter(Tag %in% 3185))
# different plants

### 3356
raw.demo.list %>% lapply(function(df) df %>% filter(Tag %in% 3356))
# note in 2019: "no plant; duplicate, changed to 3480" for 16C

proc.demo = rbind(
  # Old data
  proc.demo %>% filter(!(Tag %in% 3356 & Plot %in% 13 & Xcoor %in% 16)),
  # New data
  proc.demo %>%
    filter(Tag %in% 3356 & Plot %in% 13 & Xcoor %in% 16) %>%
    mutate(
      # Add new (correct) tag
      Tag = 3480,
      # Add note
      Note = "tag changed in 2019 (was 3356)",
      edited = TRUE
    )
)

### 3387
raw.demo.list %>% lapply(function(df) df %>% filter(Tag %in% 3387))
# *both* tags changed here!
raw.demo.list %>% lapply(function(df) df %>% filter(Tag %in% 3142)) # plot 6
raw.demo.list %>% lapply(function(df) df %>% filter(Tag %in% 3178)) # plot 5

proc.demo = rbind(
  # Old data
  proc.demo %>% filter(!Tag %in% 3387),
  # New data
  proc.demo %>%
    filter(Tag %in% 3387 & Plot %in% 13 & Xcoor %in% 2) %>%
    mutate(
      # Add new (correct) tag
      Tag = 3178,
      # Add note
      Note = "tag changed in 2019 (was 3387)",
      edited = TRUE
    ),
  proc.demo %>%
    filter(Tag %in% 3387 & Plot %in% 13 & Xcoor %in% 17) %>%
    mutate(
      # Add new (correct) tag
      Tag = 3142,
      # Add note
      Note = "tag changed in 2019 (was 3387)",
      edited = TRUE
    )
)  

### 3555
raw.demo.list %>% lapply(function(df) df %>% filter(Tag %in% 3555))
# J in 2016 should be a T (likely)
proc.demo = proc.demo %>% 
  mutate(
    Ycoor = ifelse(Tag %in% 3555 & Plot %in% 4, 'J', Ycoor),
    Note  = ifelse(Tag %in% 3555 & Plot %in% 4, 'fixed ycoor (was T)', Note),
    edited = ifelse(Tag %in% 3555 & Plot %in% 4, TRUE, edited)
  )

### 3634
raw.demo.list %>% lapply(function(df) df %>% filter(Tag %in% 3634))
# somehow "O" got mis-interpreted as "E"?
proc.demo = proc.demo %>% 
  mutate(
    Ycoor = ifelse(Tag %in% 3634 & Plot %in% 4, 'E', Ycoor),
    Note  = ifelse(Tag %in% 3634 & Plot %in% 4, 'fixed ycoor (was O)', Note),
    edited = ifelse(Tag %in% 3634 & Plot %in% 4, TRUE, edited)
  )

### 3389
raw.demo.list %>% lapply(function(df) df %>% filter(Tag %in% 3389))
raw.demo.list %>% lapply(function(df) df %>% filter(Tag %in% 3481))
# in 2020, tag 3389 in plot 13 (12F) was changed to 3481
# (there is also a 3481 in plot 15 but that is fine)

proc.demo = rbind(
  # Old data
  proc.demo %>% filter(!(Tag %in% 3389 & Plot %in% 13 & Xcoor %in% 12)),
  # New data
  proc.demo %>%
    filter(Tag %in% 3389 & Plot %in% 13 & Xcoor %in% 12) %>%
    mutate(
      # Add new (correct) tag
      Tag = 3481,
      # Add note
      Note = "tag changed in 2020 (was 3389)",
      edited = TRUE
    )
)

### 3563
raw.demo.list %>% lapply(function(df) df %>% filter(Tag %in% 3563)) # proper coord is 15H?
proc.demo = proc.demo %>%
  mutate(
    Ycoor = ifelse(Tag %in% 3563 & Xcoor %in% 15, 'H', Ycoor),
    Note  = ifelse(Tag %in% 3563 & Plot %in% 4, 'fixed ycoor (was R?)', Note),
    edited = ifelse(Tag %in% 3563 & Plot %in% 4, TRUE, edited)
  )

##### Cool!
##### Now, get the off-by-ones

proc.demo.dupe %>% filter(grepl('off by', dupe.status)) %>% distinct(Tag) %>% arrange(Tag)
# 3114, 3142, 3345, 3535, 3554, 3585, 3638, 3699, 3844, 3932

### 3114
raw.demo.list %>% lapply(function(df) df %>% filter(Tag %in% 3114))
# one of these is 3429?
proc.demo %>% filter(Tag %in% 3429)
# hmm...
# the 2023 plant could easily be 3114_8_J... oh well, that's not parsimonious
proc.demo = rbind(
  # Old data
  proc.demo %>% filter(!(Tag %in% 3114 & Plot %in% 13)),
  # New data
  proc.demo %>%
    filter(Tag %in% 3114 & Plot %in% 13 & Ycoor %in% 'H') %>%
    mutate(
      # Add new (correct) tag
      Tag = 3,
      # Add note
      Note = "tag changed in 2021",
      edited = TRUE
    )
)

### 3345
raw.demo.list %>% lapply(function(df) df %>% filter(Tag %in% 3345))
# C/D interchangeable - go with D as it's more commonly used
proc.demo = proc.demo %>% 
  mutate(
    Ycoor = ifelse(Tag %in% 3345 & Plot %in% 6, 'D', Ycoor),
    Note  = ifelse(Tag %in% 3345 & Plot %in% 6, 'fixed ycoor (was C)', Note),
    edited = ifelse(Tag %in% 3345 & Plot %in% 6, TRUE, edited)
  )

### 3535
raw.demo.list %>% lapply(function(df) df %>% filter(Tag %in% 3535))
# 16 was misrecorded as 18
proc.demo = proc.demo %>% 
  mutate(
    Xcoor = ifelse(Tag %in% 3535 & Plot %in% 13, 18, Xcoor),
    Note  = ifelse(Tag %in% 3535 & Plot %in% 13, 'fixed xcoor (was 16)', Note),
    edited = ifelse(Tag %in% 3535 & Plot %in% 13, TRUE, edited)
  )

### 3554
raw.demo.list %>% lapply(function(df) df %>% filter(Tag %in% 3554))
# 6 -> 5
proc.demo = proc.demo %>% 
  mutate(
    Xcoor = ifelse(Tag %in% 3554 & Plot %in% 1, 5, Xcoor),
    Note  = ifelse(Tag %in% 3554 & Plot %in% 1, 'fixed xcoor (was 6)', Note),
    edited = ifelse(Tag %in% 3554 & Plot %in% 1, TRUE, edited)
  )

### 3585
raw.demo.list %>% lapply(function(df) df %>% filter(Tag %in% 3585))
# 2 -> 3
proc.demo = proc.demo %>% 
  mutate(
    Xcoor = ifelse(Tag %in% 3585 & Plot %in% 2, 3, Xcoor),
    Note  = ifelse(Tag %in% 3585 & Plot %in% 2, 'fixed xcoor (was 2)', Note),
    edited = ifelse(Tag %in% 3585 & Plot %in% 2, TRUE, edited)
  )

### 3638
raw.demo.list %>% lapply(function(df) df %>% filter(Tag %in% 3638))
# D -> C
proc.demo = proc.demo %>% 
  mutate(
    Ycoor = ifelse(Tag %in% 3638 & Plot %in% 13, "C", Ycoor),
    Note  = ifelse(Tag %in% 3638 & Plot %in% 13, 'fixed ycoor (was D)', Note),
    edited = ifelse(Tag %in% 3638 & Plot %in% 13, TRUE, edited)
  )

### 3699
raw.demo.list %>% lapply(function(df) df %>% filter(Tag %in% 3699))
# E or F... doesn't matter which
proc.demo %>% filter(Plot %in% 14 & Xcoor %in% 0 & Ycoor %in% c("E", "F"))
# this is the only plant at this location
proc.demo = proc.demo %>% mutate(
  Ycoor = ifelse(Tag %in% 3699 & Plot %in% 14, 'E', Ycoor),
  Note  = ifelse(Tag %in% 3699 & Plot %in% 14, 'fixed ycoor (was F)', Note),
  edited = ifelse(Tag %in% 3699 & Plot %in% 14, TRUE, edited)
)

### 3844
raw.demo.list %>% lapply(function(df) df %>% filter(Tag %in% 3844))
# 7 -> 8
proc.demo = proc.demo %>% mutate(
  Xcoor = ifelse(Tag %in% 3844 & Plot %in% 2, 8, Xcoor),
  Note  = ifelse(Tag %in% 3844 & Plot %in% 2, 'fixed xcoor (was 7)', Note),
  edited = ifelse(Tag %in% 3844 & Plot %in% 2, TRUE, edited)
)

### 3932
raw.demo.list %>% lapply(function(df) df %>% filter(Tag %in% 3932))
# J -> H
proc.demo = proc.demo %>% mutate(
  Ycoor = ifelse(Tag %in% 3932 & Plot %in% 15, 'H', Ycoor),
  Note  = ifelse(Tag %in% 3932 & Plot %in% 15, 'fixed ycoor (was J)', Note),
  edited = ifelse(Tag %in% 3932 & Plot %in% 15, TRUE, edited)
)

unique(proc.demo$Xcoor) # all numeric - good
unique(proc.demo$Ycoor) # not sure what is up with these numbers?
# ALSO: there are spaces in here...

# look at the plants that have numeric Y coordinates
proc.demo %>% filter(grepl('\\.', Ycoor))
# hmm...
proc.demo %>% filter(Plot %in% 9)
# okay... well... guess these should stay in place

# look at the space
proc.demo %>% filter(grepl('\\s', Ycoor))
proc.demo = proc.demo %>% mutate(Ycoor = gsub('\\s', '', Ycoor))
# fixed.

################
##### Assign IDs
################

proc.demo = proc.demo %>% mutate(plantid = paste0(Tag, "_", Plot, "_", Xcoor, Ycoor))

######################################################################
##### Check column types
######################################################################

str(proc.demo)
# Stalk height, umbel number, umbel diameter, and year tag are characters...

# Fix stalk height
unique(proc.demo$Stalk_Height) # ah it's a random space in here
proc.demo %>% mutate(Stalk_Height = gsub('\\s', '', Stalk_Height)) %>% pull(Stalk_Height) %>% unique()
proc.demo = proc.demo %>%
  mutate(
    Stalk_Height = gsub('\\s', '', Stalk_Height),
    Stalk_Height = ifelse(grepl('NA', Stalk_Height), NA, Stalk_Height),
    Stalk_Height = as.numeric(Stalk_Height)
  )

# Fix number of umbels
unique(proc.demo$No.umbels) # 7.5 and 0N
raw.demo.list[[3]] %>% filter(No.umbels %in% '7.5') # hmm... no notes
raw.demo.list[[8]] %>% filter(No.umbels %in% '0N') # also no notes
# Looks from this like 0N should be zero (no umbel diameter) and 7.5 can be rounded down to 7
# Test to make sure this chunk of code works 
proc.demo %>%
  mutate(
    numbel = case_when(
      No.umbels %in% '0N' ~ '0',
      No.umbels %in% '7.5' ~ '7',
      .default = No.umbels
    ),
    # Convert to numeric
    numbel = as.numeric(numbel)
  ) %>%
  distinct(numbel, No.umbels)

# looks good!
proc.demo = proc.demo %>%
  mutate(
    No.umbels = case_when(
      No.umbels %in% '0N' ~ '0',
      No.umbels %in% '7.5' ~ '7',
      .default = No.umbels
    ),
    # Convert to numeric
    No.umbels = as.numeric(No.umbels)
  )

# Fix umbel diameters
unique(proc.demo$umbel.diam)
# SOS - see other sheet - this will require merging in the "multiple umble" data...
# for now we can at least fix the nA to NA
proc.demo = proc.demo %>% mutate(umbel.diam = ifelse(umbel.diam %in% 'nA', NA, umbel.diam))

######################################################################
##### Other cleaning
######################################################################



# There are duplicate records for some reason for plant 3430 in plot 1

raw.demo.list %>% lapply(function(x) x %>% filter(Tag %in% 3430))

# the ones in plot 1 are duplicate records of the same plant
# I say just remove the duplicates

proc.demo = rbind(
  # Records without this tag
  proc.demo %>% filter(!(Tag %in% 3430 & Plot %in% 1)),
  proc.demo %>%
    filter(Tag %in% 3430 & Plot %in% 1) %>%
    arrange(Year, desc(No.leaves)) %>%
    filter(!duplicated(Year)) %>%
    mutate(Note = "dupe records removed", edited = TRUE)
)

# other dupes?

proc.demo %>% group_by(plantid, Year) %>% filter(n() > 1)
# none!

######################################################################
##### Export what we currently have
######################################################################

# Arrange neatly before export!
proc.demo %>%
  arrange(Year, Plot, plantid) %>%
  write.csv(
    '01_data_cleaning/out/demo_all_cleaned.csv', 
    row.names = FALSE
  )

######################################################################
##### Old code 
######################################################################

# YrTag - probably not necessary to reaflly do anything about at this stage.
# # Look at YrTag
# unique(proc.demo$YrTag)
# # question mark record seems to be the major offender
# proc.demo %>% filter(grepl('\\?', YrTag))
# # bummer...
# 
# # Are there any Tags where YrTag is ? that has other records where YrTag is
# # specified?
# proc.demo %>%
#   group_by(Plot, Tag) %>%
#   filter(any(grepl('\\?', YrTag))) %>%
#   arrange(Tag, Plot, YrTag)
# 
# proc.demo %>%
#   group_by(Plot, Tag) %>%
#   filter(any(grepl('\\?', YrTag))) %>%
#   summarise(specified = any(grepl('\\d', YrTag)))
# # Two of the five (both in same plot...)
# # (is YrTag just the first year of observation...?)
# 
# # Try this code chunk out
# proc.demo %>%
#   group_by(Plot, Tag) %>%
#   mutate(yt = ifelse(
#     # If any records have YrTag equal to ? but others with numbers
#     test = any(grepl('\\?', YrTag)) & any(grepl('\\d', YrTag)), 
#     # Set the YrTag equal to the numeric YrTag...
#     yes = unique(grep('\\d', YrTag, value = TRUE)),
#     no = YrTag
#     )
#   ) %>%
#   filter(yt != YrTag) %>%
#   arrange(Plot, Tag)

