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
  lapply(
    function(df) {
      df %>% 
        select(c(all_of(shared.cols), contains('ote'))) %>%
        rename_with(function(x) 'demo.note', contains('ote'))
    }
  ) %>%
  do.call(what = rbind)

head(raw.demo)
tail(raw.demo)

# Create a proc.demo data frame for storing processed demo
proc.demo = raw.demo %>% mutate(proc.note = NA, edited = FALSE)

# Also get in 'notes' column
raw.demo.list %>%
  lapply(function(df) df %>% 
           select(Plot, Tag, Xcoor, Ycoor, contains('ote')) %>%
           rename_with(function(x) 'demo.note', contains('ote')))

######################################################################
##### Fix data entry mistakes
######################################################################

##### 2018 column mistakes
# in 2018 column names are wrong (!)
# leaf length was recorded as "height"
# umble height was recorded as "leaf length"
# but in years where plants did not flower, leaf length = "height"
proc.demo$Leaf.length[proc.demo$Year %in% 2018] = raw.demo.list[[3]]$Height
proc.demo$Stalk_Height[proc.demo$Year %in% 2018] = raw.demo.list[[3]]$Leaf.length

# umbel count of 7.5 is a mis-entry - should be 2 umbels, stalk height 7.5
proc.demo[with(proc.demo, Plot %in% 15 & Tag %in% 3492 & Year %in% 2018), c("Stalk_Height", "No.umbels")] = 
  proc.demo[with(proc.demo, Plot %in% 15 & Tag %in% 3492 & Year %in% 2018), c("No.umbels", "Stalk_Height")]
# (fix umbel heights later...)
proc.demo$edited[with(proc.demo, Plot %in% 15 & Tag %in% 3492 & Year %in% 2018)] = TRUE
proc.demo$proc.note[with(proc.demo, Plot %in% 15 & Tag %in% 3492 & Year %in% 2018)] = 
  "data entry mistake; umbel count and stalk height were swapped"

# I think this alos happened with 3493 in 2018
proc.demo[with(proc.demo, Plot %in% 15 & Tag %in% 3493 & Year %in% 2018), c("Stalk_Height", "No.umbels")] = 
  proc.demo[with(proc.demo, Plot %in% 15 & Tag %in% 3493 & Year %in% 2018), c("No.umbels", "Stalk_Height")]
# (fix umbel heights later...)
proc.demo$edited[with(proc.demo, Plot %in% 15 & Tag %in% 3493 & Year %in% 2018)] = TRUE
proc.demo$proc.note[with(proc.demo, Plot %in% 15 & Tag %in% 3493 & Year %in% 2018)] = 
  "data entry mistake; umbel count and stalk height were swapped"

##### Fix some mistakenly-entered leaf data that I found

proc.demo = proc.demo %>%
  mutate(
    Leaf.length = ifelse(Tag %in% 3687 & Plot %in% 1 & Year %in% 2020, 10.1, Leaf.length),
    edited = ifelse(Tag %in% 3687 & Plot %in% 1 & Year %in% 2020, TRUE, edited),
    proc.note = ifelse(ifelse(Tag %in% 3687 & Plot %in% 1 & Year %in% 2020,
                         "data entry mistake; was entered as 1",
                         proc.note))
  )

### 2018 case - ycoord "O" was mistakenly entered as zero in leaf count field
# (also y-coor was left out in 2017 - add that back in here too)
proc.demo = proc.demo %>%
  mutate(
    Ycoor = ifelse(Tag %in% 3103 & Plot %in% 4 & Year %in% 2017, 'O', Ycoor),
    No.leaves = ifelse(Tag %in% 3103 & Plot %in% 4 & Year %in% 2018, 1, No.leaves),
    edited = ifelse(Tag %in% 3103 & Plot %in% 4 & Year %in% 2018, TRUE, edited),
    proc.note = ifelse(ifelse(Tag %in% 3103 & Plot %in% 4 & Year %in% 2018,
                              "data entry mistake; leaf count was entered as 0",
                              proc.note))
  )

##### At some point someone mistakenly entered decimal leaf-counts
# (I checked the original data - a count was recorded, poorly-scratched out, and
# replaced with another count, and the data-enterer placed both separated by a
# period (?))
# so just take out the first digit and the period

proc.demo %>% filter(grepl('\\.', No.leaves))

proc.demo = proc.demo %>%
  mutate(No.leaves = gsub('\\d\\.', '', as.character(No.leaves)))

##### 2020 umbel counts: there are ten plants where umbel counts were entered
# as double digit (?) - there's one for each number 11-20
# these occur in the last 10 new plants, in order, added to the new plant data
# sheet; no umbel counts are recorded on the data sheet
# for these ten plants there are no umbel counts or stalk heights
# fixeth!

proc.demo = proc.demo %>%
  mutate(
    No.umbels = ifelse(No.umbels %in% as.character(9:20) & Year %in% 2020 & is.na(Stalk_Height), '0', No.umbels)
  )

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
    proc.note = ifelse(
      Tag %in% 3044 & Plot %in% 7,
      "one ycoor fixed (was NA)",
      proc.note
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

### 3142
proc.demo = proc.demo %>%
  mutate(
    edited = ifelse(Tag %in% 31432 & Plot %in% 6 & Ycoor %in% 'C', TRUE, edited),
    proc.note = ifelse(
      Tag %in% 31432 & Plot %in% 6 & Ycoor %in% 'C',
      'coordinate changed (was recorded as C, is actually D)',
      proc.note
    ),
    Ycoor = ifelse(Tag %in% 3142 & Plot %in% 6, 'D', Ycoor),
  )

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
# I'm going to assume it was always one plant and the 2019 record was entered
# for the wrong coordinate.
# 18C is the real plant, 17B is a dupe and should be deleted.

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
      proc.note  = "coordinates swapped in 2019 records based on prior records",
      edited = TRUE
    )
) %>%
  # Remove all records for 17B plant
  filter(!(Tag %in% 3180 & Plot %in% 5 & Xcoor %in% 17))

### 3342: 
raw.demo.list %>% lapply(function(df) df %>% filter(Tag %in% 3342))
# Plot 4 plant is a diversity plant - keep
# Plot 5: assume the 18B plant was a mistake.
# (Only one record for it anyway - it would not matter anyway)
proc.demo = proc.demo %>%
  filter(!(Tag %in% 3342 & Plot %in% 5 & Xcoor %in% 18))

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
      proc.note = "y-coordinate changed based on 2023 note (was K)",
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
      proc.note = "tag changed based on 2019 note (was 3396)",
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
      proc.note = "tag changed based on 2023 note (was 3396)",
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
    proc.note  = ifelse(Tag %in% 3090 & Plot %in% 4, "fixed ycoor (was P)", proc.note),
    edited = ifelse(Tag %in% 3090 & Plot %in% 4, TRUE , edited)
  )

### 3125
raw.demo.list %>% lapply(function(df) df %>% filter(Tag %in% 3125))
# these are two different plants

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
      proc.note = "tag changed in 2019 (was 3356)",
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
      proc.note = "tag changed in 2019 (was 3387)",
      edited = TRUE
    ),
  proc.demo %>%
    filter(Tag %in% 3387 & Plot %in% 13 & Xcoor %in% 17) %>%
    mutate(
      # Add new (correct) tag
      Tag = 3142,
      # Add note
      proc.note = "tag changed in 2019 (was 3387)",
      edited = TRUE
    )
)  

### 3555
raw.demo.list %>% lapply(function(df) df %>% filter(Tag %in% 3555))
# J in 2016 should be a T (likely)
proc.demo = proc.demo %>% 
  mutate(
    Ycoor = ifelse(Tag %in% 3555 & Plot %in% 4, 'J', Ycoor),
    proc.note  = ifelse(Tag %in% 3555 & Plot %in% 4, 'fixed ycoor (was T)', proc.note),
    edited = ifelse(Tag %in% 3555 & Plot %in% 4, TRUE, edited)
  )

### 3634
raw.demo.list %>% lapply(function(df) df %>% filter(Tag %in% 3634))
# somehow "O" got mis-interpreted as "E"?
proc.demo = proc.demo %>% 
  mutate(
    Ycoor = ifelse(Tag %in% 3634 & Plot %in% 4, 'E', Ycoor),
    proc.note  = ifelse(Tag %in% 3634 & Plot %in% 4, 'fixed ycoor (was O)', proc.note),
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
      proc.note = "tag changed in 2020 (was 3389)",
      edited = TRUE
    )
)

### 3563
raw.demo.list %>% lapply(function(df) df %>% filter(Tag %in% 3563)) # proper coord is 15H?
proc.demo = proc.demo %>%
  mutate(
    Ycoor = ifelse(Tag %in% 3563 & Xcoor %in% 15, 'H', Ycoor),
    proc.note  = ifelse(Tag %in% 3563 & Plot %in% 4, 'fixed ycoor (was R?)', proc.note),
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
      Tag = 3114,
      # Add note
      proc.note = "tag changed in 2021",
      edited = TRUE
    )
)

### 3345
raw.demo.list %>% lapply(function(df) df %>% filter(Tag %in% 3345))
# C/D interchangeable - go with D as it's more commonly used
proc.demo = proc.demo %>% 
  mutate(
    Ycoor = ifelse(Tag %in% 3345 & Plot %in% 6, 'D', Ycoor),
    proc.note  = ifelse(Tag %in% 3345 & Plot %in% 6, 'fixed ycoor (was C)', proc.note),
    edited = ifelse(Tag %in% 3345 & Plot %in% 6, TRUE, edited)
  )

### 3535
raw.demo.list %>% lapply(function(df) df %>% filter(Tag %in% 3535))
# 16 was misrecorded as 18
proc.demo = proc.demo %>% 
  mutate(
    Xcoor = ifelse(Tag %in% 3535 & Plot %in% 13, 18, Xcoor),
    proc.note  = ifelse(Tag %in% 3535 & Plot %in% 13, 'fixed xcoor (was 16)', proc.note),
    edited = ifelse(Tag %in% 3535 & Plot %in% 13, TRUE, edited)
  )

### 3554
raw.demo.list %>% lapply(function(df) df %>% filter(Tag %in% 3554))
# 6 -> 5
proc.demo = proc.demo %>% 
  mutate(
    Xcoor = ifelse(Tag %in% 3554 & Plot %in% 1, 5, Xcoor),
    proc.note  = ifelse(Tag %in% 3554 & Plot %in% 1, 'fixed xcoor (was 6)', proc.note),
    edited = ifelse(Tag %in% 3554 & Plot %in% 1, TRUE, edited)
  )

### 3585
raw.demo.list %>% lapply(function(df) df %>% filter(Tag %in% 3585))
# 2 -> 3
proc.demo = proc.demo %>% 
  mutate(
    Xcoor = ifelse(Tag %in% 3585 & Plot %in% 2, 3, Xcoor),
    proc.note  = ifelse(Tag %in% 3585 & Plot %in% 2, 'fixed xcoor (was 2)', proc.note),
    edited = ifelse(Tag %in% 3585 & Plot %in% 2, TRUE, edited)
  )

### 3638
raw.demo.list %>% lapply(function(df) df %>% filter(Tag %in% 3638))
# D -> C
proc.demo = proc.demo %>% 
  mutate(
    Ycoor = ifelse(Tag %in% 3638 & Plot %in% 13, "C", Ycoor),
    proc.note  = ifelse(Tag %in% 3638 & Plot %in% 13, 'fixed ycoor (was D)', proc.note),
    edited = ifelse(Tag %in% 3638 & Plot %in% 13, TRUE, edited)
  )

### 3699
raw.demo.list %>% lapply(function(df) df %>% filter(Tag %in% 3699))
# E or F... doesn't matter which
proc.demo %>% filter(Plot %in% 14 & Xcoor %in% 0 & Ycoor %in% c("E", "F"))
# this is the only plant at this location
proc.demo = proc.demo %>% mutate(
  Ycoor = ifelse(Tag %in% 3699 & Plot %in% 14, 'E', Ycoor),
  proc.note  = ifelse(Tag %in% 3699 & Plot %in% 14, 'fixed ycoor (was F)', proc.note),
  edited = ifelse(Tag %in% 3699 & Plot %in% 14, TRUE, edited)
)

### 3844
raw.demo.list %>% lapply(function(df) df %>% filter(Tag %in% 3844))
# 7 -> 8
proc.demo = proc.demo %>% mutate(
  Xcoor = ifelse(Tag %in% 3844 & Plot %in% 2, 8, Xcoor),
  proc.note  = ifelse(Tag %in% 3844 & Plot %in% 2, 'fixed xcoor (was 7)', proc.note),
  edited = ifelse(Tag %in% 3844 & Plot %in% 2, TRUE, edited)
)

### 3932
raw.demo.list %>% lapply(function(df) df %>% filter(Tag %in% 3932))
# J -> H
proc.demo = proc.demo %>% mutate(
  Ycoor = ifelse(Tag %in% 3932 & Plot %in% 15, 'H', Ycoor),
  proc.note  = ifelse(Tag %in% 3932 & Plot %in% 15, 'fixed ycoor (was J)', proc.note),
  edited = ifelse(Tag %in% 3932 & Plot %in% 15, TRUE, edited)
)

unique(proc.demo$Xcoor) # all numeric - good
unique(proc.demo$Ycoor) # not sure what is up with these numbers?

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
##### Look for plants with with '[Tt]ag' in the notes field
################

# (old code)
# # from here:
# # demo %>% filter(grepl('[Tt]ag', demo.note), np %in% 'none') # np is added col

### 3323/3500
raw.demo %>% filter(Tag %in% c(3323, 3500)) %>% arrange(Plot, Year)
# no 2016 record for 3500 in plot 2

### 3703/3844
raw.demo %>% filter(Tag %in% c(3703, 3844)) %>% arrange(Year, Tag) # 12 leaves???
# actually this looks fine too - two different plants
raw.demo %>% filter(Tag %in% c(3703, 3844)) %>% arrange(Tag, Year)

### 3093
raw.demo %>% filter(Tag %in% 3093) # not a 2017 recruit...
# ('found plant with unlisted tag' - these are not recruits)

### 3090
raw.demo %>% filter(Tag %in% 3090) # ugh...

### 3049
raw.demo %>% filter(Tag %in% 3049)
# not a new recruit...

### 3814/3782
raw.demo %>% filter(Tag %in% c(3814, 3782))
# ugh... which is which?
# 3782 is not new recruit in 2022/2023
# also where is 3782 record for 2022?
raw.demo %>% filter(Plot %in% 1, Xcoor %in% 0, Ycoor %in% 'J')

### 3180
raw.demo %>% filter(Tag %in% 3180) %>% arrange(Plot, Year)
proc.demo %>% filter(Tag %in% 3180) %>% arrange(Plot, Xcoor, Year)

### 3142
raw.demo %>% filter(Tag %in% c(3142, 3792))
proc.demo %>% filter(Tag %in% c(3142, 3792), Plot %in% 6)
# no records for 3792...

### 3179
raw.demo %>% filter(Tag %in% c(3179, 3748), Plot %in% 6)
# no records for 3748...

### 3429
raw.demo %>% filter(Tag %in% c(3, 3429))
raw.demo %>% filter(Plot %in% 13, Xcoor %in% 8, Ycoor %in% 'H')

### 3094
raw.demo %>% filter(Tag %in% 3094)
# woohoo... definitely dead plant!

### 3171
raw.demo %>% filter(Tag %in% c(3171, 3588)) %>% arrange(Year, Tag)
# good demo records!

### 3036
raw.demo %>% filter(Tag %in% c(3036, 3762)) %>% arrange(Year, Tag)
# finally! this routine is useful somehow...

proc.demo = proc.demo %>%
  # Get rid of empty record for 3036 in 2021
  filter(!(Tag %in% 3036 & Year %in% 2021)) %>%
  # Re-assign tag in 3762 record
  mutate(
   proc.note  = ifelse(Tag %in% 3762 & Year %in% 2021, 'manually fixed tag (see demo note)', proc.note),
   edited = ifelse(Tag %in% 3762 & Year %in% 2021, TRUE, edited),
   Tag = ifelse(Tag %in% 3762 & Year %in% 2021, 3036, Tag)
  ) %>%
  # Now, get rid of 3762 records
  filter(!(Tag %in% 3762 & Year > 2021 & Plot %in% 3))
  
### 3915/5051
raw.demo %>% filter(Tag %in% c(3915, 5051)) %>% arrange(Year, Tag)
# okay... same plant
# ugh

proc.demo = rbind(
  # Separate out records of other plants
  proc.demo %>% filter(!(Tag %in% c(3915, 5051) & Plot %in% 7)),
  # Modify the baddies
  proc.demo %>%
    filter(Tag %in% c(3915, 5051) & Plot %in% 7) %>%
    mutate(
      # Add note about editing
      proc.note = ifelse(Year < 2021, 'changed tag manually; was 3195 before 2021', proc.note),
      edited = ifelse(Year < 2021, TRUE, edited),
      # Tag was replaced in 2021, so set prior years' recs to current tag
      Tag = ifelse(Year < 2021, 5051, Tag)
    ) %>%
  # Oh... it looks like records in/after 2021 for 3915 are empty
  # so I can just remove these now...
  filter(!Tag %in% 3915)
)

### 3423/3815
raw.demo %>% filter(Tag %in% c(3423, 3815)) %>% arrange(Plot, Tag, Year)
raw.demo %>% filter(Tag %in% c(3423, 3815), Plot %in% 12) %>% arrange(Year, Tag)
# lol
# okay well I guess no changes to this

### 3119/3388
raw.demo %>% filter(Tag %in% c(3388, 3119)) %>% arrange(Plot, Tag)
raw.demo %>% filter(Tag %in% c(3388, 3119), Plot %in% 13) %>% arrange(Year, Tag)
proc.demo %>% filter(Tag %in% c(3388, 3119), Plot %in% 13) %>% arrange(Year, Tag)
# easy one!
proc.demo = proc.demo %>%
  mutate(
    proc.note = ifelse(Tag %in% 3388 & Plot %in% 13, 'tag manually edited; was 3388', proc.note),
    edited    = ifelse(Tag %in% 3388 & Plot %in% 13, TRUE, edited),
    Tag       = ifelse(Tag %in% 3388 & Plot %in% 13, 3119, Tag)
  )

### 3615
raw.demo %>% filter(Tag %in% c(3615, 3745)) %>% arrange(Plot, Tag)
proc.demo %>% filter(Tag %in% c(3615, 3745), Plot %in% 15) %>% arrange(Year, Tag)
# say 3615 and 3745 are different plants
# assume that if 3745 was alive in 2023, it would have been found
# so, change tag in 2023 record to 3745
proc.demo = proc.demo %>%
  mutate(
    proc.note = ifelse(Tag %in% 3615 & Plot %in% 15 & Year %in% 2023, 
                       'tag manually edited; was 3615', 
                       proc.note),
    edited    = ifelse(Tag %in% 3615 & Plot %in% 15 & Year %in% 2023, TRUE, edited),
    Tag       = ifelse(Tag %in% 3615 & Plot %in% 15 & Year %in% 2023, 3745, Tag)
  )

### 5023
raw.demo %>% filter(Tag %in% c(5023, 3807))
# seems like these are the same plant
proc.demo %>% filter(Tag %in% c(5023, 3807))

proc.demo = rbind(
  # Set aside other plants
  proc.demo %>% filter(!(Tag %in% c(3807, 5023) & Plot %in% 5)),
  # Isolate and fix records for this plant
  proc.demo %>%
    filter(Tag %in% c(3807, 5023), Plot %in% 5) %>%
    # Remove empty records
    filter(!(Tag %in% 3807 & Year %in% 2021)) %>%
    filter(!(Tag %in% 5023 & Year %in% 2022)) %>%
    # One remaining record - 5023 in 2021, change tag to 3807
    mutate(
      proc.note = ifelse(Tag %in% 5023, 'tag and ycoor edited; was 5023-I', proc.note),
      edited    = ifelse(Tag %in% 5023, TRUE, edited),
      Ycoor     = ifelse(Tag %in% 5023, 'H', Ycoor),
      Tag       = ifelse(Tag %in% 5023, 3807, Tag)
    )
    
)

### 3729
raw.demo %>% filter(Tag %in% 3729) # no problems here

### 3773
raw.demo %>% filter(Tag %in% 3773) # no problems here

### 3039
raw.demo %>% filter(Tag %in% 3039)
# ahhhhhhhhh!!!
# actually it looks like 3039 did in fact die in 2020
raw.demo %>% filter(Tag %in% c(3039, 3558, 3726)) %>% arrange(Tag, Year)
# well... 3558 looks like it died
# note suggests to me that 3726 = 3039
proc.demo = rbind(
  # Separate out records for other plants
  proc.demo %>% filter(!(Tag %in% c(3039, 3726) & Plot %in% 13)),
  # Isolate records for 3039/3726
  proc.demo %>%
    filter(Tag %in% c(3039, 3726) & Plot %in% 13) %>%
    # Get rid of records for 3039 after 2020 (when 3726 was first tagged)
    filter(!(Tag %in% 3039 & Year > 2020)) %>%
    mutate(
      proc.note = ifelse(Tag %in% 3726, 'tag manually edited; was 3726', proc.note),
      edited    = ifelse(Tag %in% 3726, TRUE, edited),
      Tag       = ifelse(Tag %in% 3726, 3039, Tag)
    )
)

### 3558
raw.demo %>% filter(Tag %in% c(3558, 3769))
proc.demo %>% filter(Tag %in% c(3558, 3769)) %>% arrange(Year, Tag)
# yep... think these are the same plant
proc.demo = rbind(
  # These records are fine
  proc.demo %>% filter(!(Tag %in% c(3558, 3769) & Plot %in% 13)),
  # Isolate relevant records and fix
  proc.demo %>%
    filter(Tag %in% c(3558, 3769), Plot %in% 13) %>%
    # Gert rid of 3558 records after 2020
    filter(!(Tag %in% 3558 & Year > 2020)) %>%
    # Fix info in remaining records
    mutate(
      proc.note = ifelse(Tag %in% 3769, 'tag changed; was 3769', proc.note),
      edited    = ifelse(Tag %in% 3769, TRUE, edited),
      Ycoor     = ifelse(Tag %in% 3769, 'F', Ycoor),
      Tag       = ifelse(Tag %in% 3769, 3558, Tag)
    )
)

################
##### Look at plants with tag numbers in notes
################

# (list of plants to check from here)
# demo %>%
#   filter(grepl('\\d{4}', demo.note)) %>%
#   filter(!grepl('[Tt]ag', demo.note)) %>%
#   select(-c(plantid, trt)) %>%
#   filter(!edited) %>%
#   arrange(grepl('20[12][01236789]', demo.note), Plot, Tag, Year)

# (going to assume there's no useful info in 'near' or 'next to' notes)

### 3964
raw.demo %>% filter(Tag %in% 3964)
raw.demo %>% filter(Tag %in% c(3964, 3934)) %>% arrange(Year, Tag)
# nah these are separate plants

# 3456
raw.demo %>% filter(Tag %in% c(3456, 7587)) %>% arrange(Year, Tag)
# two leaves in every year... going to assume this is the same plant!
proc.demo %>% filter(Tag %in% c(3456, 7587))
proc.demo = rbind(
  # Records for other plants
  proc.demo %>% filter(!(Tag %in% c(3456, 7587) & Plot %in% 6)),
  # Records for this plant
  proc.demo %>%
    filter(Tag %in% c(3456, 7587)) %>%
    filter(!(Tag %in% 3456 & Year > 2021)) %>%
    mutate(
      proc.note = ifelse(Tag %in% 7587, 'edited tag; was 7587', proc.note),
      edited    = ifelse(Tag %in% 7587, TRUE, edited),
      Ycoor     = ifelse(Tag %in% 7587, 'B', Ycoor),
      Tag       = ifelse(Tag %in% 7587, 3456, Tag)
    )
)

### 3059
proc.demo %>% filter(Tag %in% 3059)
# argh...
proc.demo %>% filter(Tag %in% c(3059, 3731, 3770))
# I mean I guess I should follow the splitting rules... but which of these is it?
# idk, and it would be arbitrary
# maybe just ignore...

### 7511
proc.demo %>% filter(Tag %in% 7511)
proc.demo %>% filter(Tag %in% c(7511, 3417)) %>% arrange(Plot, Year, Tag)
# lol... man
proc.demo %>% filter(Tag %in% c(7511, 3417), Plot %in% 12)
# no - these are in fact different plants based on 2023 record!

### 3393
proc.demo %>% filter(Tag %in% c(3393, 3157), Plot %in% 15)
# hmm... assume these are two different plants
# both seen in 2017, but not again after that
# hmm... oh well

### 3417
raw.demo %>% filter(Tag %in% c(3417, 3471), Plot %in% 14)
proc.demo = proc.demo %>%
  mutate(
    proc.note = ifelse(Tag %in% 3417 & Plot %in% 14, 'changed tag; was mis-entered as 3417', proc.note),
    edited    = ifelse(Tag %in% 3417 & Plot %in% 14, TRUE, edited),
    Tag       = ifelse(Tag %in% 3417 & Plot %in% 14, 3471, Tag)
  )


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
    mutate(proc.note = "dupe records removed", edited = TRUE)
)

# other dupe records?
proc.demo %>% group_by(Plot, Tag, Year) %>% filter(n() > 1)
# Just the ones I already am aware of (in plots 2/5)

# notes about dupes?
proc.demo %>% filter(grepl('[Dd]up[el]', demo.note), !edited)

# Confusing note about tags 3428, 3380, and 7537
proc.demo %>% filter(Tag %in% c(3880, 3428, 7537)) %>% arrange(Tag, Year)
# Plant 3880 got changed to 7537

proc.demo = proc.demo %>%
  mutate(
    edited = ifelse(Tag %in% 3880 & Plot %in% 13, TRUE, edited),
    proc.note = ifelse(
      Tag %in% 3880 & Plot %in% 13,
      'tag changed from 3880 to 7537 in 2022',
      proc.note
    ), 
    Tag = ifelse(Tag %in% 3880 & Plot %in% 13, 7537, Tag)
  )

# Records to fix/combine

proc.demo %>% filter(grepl('fix', demo.note))

# 5044 and 3361 (note mis-entered)
proc.demo %>% filter(Tag %in% c(3361, 5044), Plot %in% 5)
# Plant 3361 was not seen in 2021 or 2022, but new plant 5044 was added
# 5044 is 3361
# (But was recorded with tag 3361 in 2024)

proc.demo = proc.demo %>%
  filter(!(Tag %in% 3361 & Year %in% 2021:2023)) %>%
  filter(!(Tag %in% 5044 & Year %in% 2024)) %>%
  mutate(
    edited = ifelse(Tag %in% 5044 & Plot %in% 5, TRUE, edited),
    proc.note = ifelse(
      Tag %in% 5044 & Plot %in% 5,
      'tag changed from 5044 (old/dupe) to 3361 (still present in field)',
      proc.note
    ),
    Ycoor = ifelse(Tag %in% 5044 & Plot %in% 5, 'A', Ycoor),
    Tag = ifelse(Tag %in% 5044 & Plot %in% 5, 3361, Tag)
  )


################
##### Assign IDs
################

# Plot ID assigned as such:
# Tags that (after processing) appear multiple times in the same plot *in the
# same year* get Tag_plot.coordinate
# Tags that appear only once in a plot (overwhelming majority) get just tag_coord

# These are the plots that are repeated within a plot
proc.demo %>%
  group_by(Plot, Tag) %>%
  filter(any(duplicated(Year))) %>%
  distinct(Tag, Plot)

proc.demo = proc.demo %>% 
  group_by(Plot, Tag) %>%
  mutate(flag = any(duplicated(Year))) %>%
  ungroup() %>%
  mutate(
    plantid = ifelse(
      flag,
      # If duplicates (i.e., multiple tags in same plot)
      paste(paste(Tag, Plot, sep = '_'), paste0(Xcoor, Ycoor), sep = '.'),
      # If no duplicates
      paste(Tag, Plot, sep = '_')
    )
  ) %>%
  select(-flag)

# Tests:

# - Just to make sure: are any plantids recycled within a plot (appearing more than once in a year)

proc.demo %>%
  group_by(plantid, Year) %>%
  filter(any(duplicated(Plot)))

# - Anyone here with different x *and* y coords?

proc.demo %>%
  group_by(plantid) %>%
  filter(any(length(unique(Xcoor)) > 1 & length(unique(Ycoor)) > 1)) %>%
  arrange(plantid)
# None - great

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

