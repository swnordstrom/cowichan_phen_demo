# Script for combining/reconciling Dactylis denography data
# this is a dirty script for making a dataset for Jenn's class, currently
# sn -  jan 2024

library(dplyr)
library(tidyr)

dac16 = read.csv('00_raw_data/dactylis_demography/2016_Dactylis_Demography_Data.csv')
dac17 = read.csv('00_raw_data/dactylis_demography/2017_Dactylis_Demography_Data.csv')

head(dac16)
head(dac17)
nrow(dac16)
nrow(dac17)

all.dac = rbind(dac16, dac17)

# Add plantid
all.dac = all.dac %>%
  mutate(plantid = paste(Plot, paste0(Xcoor, Ycoor), sep = '_'))

head(all.dac)

# Did all of these plants survive?
all.dac %>%
  group_by(Year = paste0('y', Year), No.stems) %>%
  summarise(n = n()) %>%
  pivot_wider(names_from = Year, values_from = n, values_fill = 0) %>%
  arrange(No.stems)
# they're all alive..?

all.dac %>% filter(!No.stems) # lol

all.dac %>% group_by(plantid) %>% filter(max(Year) < 2017)

all.dac = dir('00_raw_data/dactylis_demography/') %>%
  lapply(FUN = function(x) read.csv(paste0('00_raw_data/dactylis_demography/', x))) %>%
  do.call(what = rbind)
# annoying!

all.dac.list = dir('00_raw_data/dactylis_demography/') %>%
  lapply(FUN = function(x) read.csv(paste0('00_raw_data/dactylis_demography/', x)))

lapply(all.dac.list, names)
lapply(all.dac.list, dim) # okay... cool, good, etc.

names(all.dac.list[[4]]) = names(all.dac.list[[1]])

all.dac = do.call(rbind, all.dac.list) %>%
  # Add an ID column
  mutate(plantid = paste(Plot, paste0(Xcoor, Ycoor), sep = '_')) %>%
  # add in treatment data
  merge(read.csv('00_raw_data/plot_treatments.csv') %>% rename(Plot = plot))

head(all.dac)

# Look for missing plants
# all.dac %>%
  # select(Year, Plot, No.stems, plantid) %>%
  
#######

# Maybe just look at 2021-2022 (just picking two later years)
table(all.dac$Year) # okay this is cool, good, etc.

dac20s = all.dac %>% filter(Year %in% 2021:2023)

dac20s

dac20s = dac20s %>% 
  # lots of NPs from 2022... remove these
  filter(!(!No.stems & (Year %in% 2022))) %>%
  # let's mark plants as "new" for recruitment
  group_by(plantid) %>%
  mutate(is.new = Year == min(Year)) %>%
  ungroup() %>%
  filter(Year > 2021) %>%
  # get rid of these 'check to see' plants - probably not dac...
  filter(!grepl('check\\sto\\smake\\ssure\\sthis', Notes)) %>%
  # also these "np for > 2 years"
  filter(!(grepl('np\\sfor\\s', Notes)))

head(dac20s)
nrow(dac20s)
sum(grepl('[Nn]ew', dac20s$Notes))
sum(dac20s$is.new) # ugh

# Add an ID column
dac20s = all.dac %>%
  filter(Year %in% 2022:2023) %>%
  # lots of NPs from 2022... remove these
  filter(!(!No.stems & (Year %in% 2022))) %>%
  # get rid of these 'check to see' plants - probably not dac...
  filter(!grepl('check\\sto\\smake\\ssure\\sthis', Notes)) %>%
  # also these "np for > 2 years"
  filter(!(grepl('np\\sfor\\s', Notes))) %>%
  # add flag for is.seedling
  mutate(is.sling = (grepl('[Nn]ew', Notes) & Year %in% 2022))

dac20s %>% filter(is.sling) # not that many... oh well

dac20s = dac20s %>%
  # fix one mistaken one
  mutate(is.sling = ifelse(grepl('close\\sto\\sa\\snew\\sones', Notes), FALSE, is.sling))

dac20s %>%
  group_by(plantid) %>%
  filter(min(Year) > 2022)
# hmm... let's just call all of these new. lmao

dac20s = dac20s %>%
  group_by(plantid) %>%
  mutate(
    is.sling = ifelse(
      min(Year) > 2022 & Year %in% 2023 & !grepl('\\?', Notes) & No.stems < 10,
      TRUE,
      is.sling
    )
  ) %>%
  ungroup()

nrow(dac20s)
dac20s %>% filter(is.sling)

with(dac20s, table(Year, No.stems))
# holy cow, 80 stems...

dac20s %>%
  mutate(
    stage = case_when(
      !No.stems ~ 'Dead',
      is.sling ~ 'Seedling',
      No.stems %in% 1:8 ~ 'Small',
      No.stems %in% 9:25 ~ 'Medium',
      No.stems > 25 ~ 'Large'
    )
  ) %>%
  group_by(Year, stage) %>%
  summarise(n = n())

# there are NAs?
dac20s %>% filter(is.na(No.stems))
# fix
dac20s = dac20s %>%
  group_by(plantid) %>%
  filter(!any(is.na(No.stems))) %>%
  ungroup()

# any other dead plants that don't have a record?
dac20s %>%
  group_by(plantid) %>%
  filter(max(Year) < 2023)

# hmm... in diversity, so probably just not resurveyed
# oh well! call it dead because this doesn't matter
dac20s %>%
  group_by(plantid) %>%
  mutate(w = 1 + as.numeric(max(Year) < 2023)) %>%
  ungroup() %>%
  uncount(weights = w) %>%
  group_by(plantid) %>%
  mutate(
    across(c(Height, No.stems, No.panicles), function(x) ifelse(duplicated(Year), 0, x)),
    Year = Year + as.numeric(duplicated(Year)),
  ) %>%
  filter(plantid %in% '4_17L')
# sweet

# add dead record using some creativity
dac20s = dac20s %>%
  group_by(plantid) %>%
  mutate(w = 1 + as.numeric(max(Year) < 2023)) %>%
  ungroup() %>%
  uncount(weights = w) %>%
  group_by(plantid) %>%
  mutate(
    across(c(Height, No.stems, No.panicles), function(x) ifelse(duplicated(Year), 0, x)),
    Year = Year + as.numeric(duplicated(Year)),
  )

# Now, add status

dac20s = dac20s %>%
  mutate(
    stage = case_when(
      !No.stems ~ 'Dead',
      is.sling ~ 'Seedling',
      No.stems %in% 1:8 ~ 'Small',
      No.stems %in% 9:25 ~ 'Medium',
      No.stems > 25 ~ 'Large'
    )
  ) %>%
  # Select relevant columns
  select(Year, plantid, Plot, Xcoor, Ycoor, trt, No.stems, No.panicles, stage)

head(dac20s)

dac20s = dac20s %>% arrange(Year, plantid)

head(dac20s)

write.csv(
  dac20s,
  file = '01_data_cleaning/dactylis/out/dactylis-2022-2023.csv',
  row.names = FALSE
)
