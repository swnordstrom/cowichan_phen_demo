library(glmmTMB)
library(dplyr)
library(tidyr)

# Get rid of this super annoying feature
options(dplyr.summarise.inform = FALSE)

rm(list = ls())

# Read in demo data and merge with treatment info
all.data = merge(
  x = read.csv('01_data_cleaning/out/demo_phen_seed_2016-2024_final.csv'),
  y = read.csv('00_raw_data/plot_treatments.csv'),
  by.x = 'Plot', by.y = 'plot'
)

nrow(all.data)
head(all.data)

# Read in phenology means
phen.treatment.means = read.csv('03_construct_kernels/phen_treatment_means.csv')
# Mean that will be used for centering
phen.ctrl.mean = phen.treatment.means$mean.phen[phen.treatment.means$trt %in% 'control']

# We're interested only in plants that are in demo for this analysis
all.demo = all.data %>% 
  filter(in.demo) %>%
  distinct(plantid, Year, .keep_all = TRUE)

nrow(all.demo)

# Get phenology dataset (for testing for growth or survival trade-offs with prior year's phen)
# this is for merging in with plant-level demo datasets, so relevant measure is by plant (not by umbel)
# i.e., need to get an aggregated mean floral emergence date per plant-year
phen.by.plant.for.demo = all.data %>%
  # Give me plants that are in phenology
  filter(in.phen) %>%
  # Extract the bud dates on file
  separate_wider_delim(phen.julis, names = paste0('uu', 1:12), delim = ';', too_few = 'align_start') %>%
  pivot_longer(starts_with('uu'), names_to = 'umbel.number', values_to = 'phen.julian') %>%
  filter(!is.na(phen.julian)) %>%
  # convert to julian date
  mutate(phen.julian = as.numeric(gsub('\\s', '', phen.julian))) %>%
  # give me the mean bud date for each plant
  group_by(plantid, Year) %>%
  summarise(phen.mean = mean(phen.julian)) %>%
  ungroup() %>%
  # Change year to factor
  mutate(Year = as.factor(Year))

# Merge the datasets together:
all.demo = merge(
  all.demo, phen.by.plant.for.demo,
  all.x = TRUE, all.y = FALSE
) %>%
  mutate(phen.c = phen.mean - phen.ctrl.mean)

# Get survival dataset

# Survival dataset:
# (Two versions: a size-dependent one and a size-independent one
# almost surely we will use the size-dependent one for analysis)

demo.surv = merge(
  # Demo in time step t+1 (THIS IS THE CENSUS YEAR DATASET)
  x = all.demo %>% 
    mutate(prev.year = Year - 1) %>%
    rename(surv.year = Year) %>%
    select(Plot, plantid, surv.year, prev.year, No.leaves, Leaf.length, surv, trt),
  # Demo in time step t (THIS CONTAINS THE PRIOR YEAR'S DEMO INFO)
  y = all.demo %>%
    # we are *only* interested in plants alive in time step t
    filter(surv) %>%
    # Select relevant columns
    select(Plot, plantid, Year, No.leaves, Leaf.length, phen.c, trt),
  by.x = c('Plot', 'plantid', 'prev.year', 'trt'),
  by.y = c('Plot', 'plantid', 'Year', 'trt'),
  suffixes = c('', '.pre'),
  all.x = FALSE, all.y = FALSE
) %>%
  # Change years to factors
  mutate(across(contains('year'), as.factor))

head(demo.surv)
# should be less than 1
table(demo.surv$surv, useNA = 'always')
# good

# Want to get a column where we can estimate sizes
demo.surv.sizes = demo.surv %>% 
  filter(
    !is.na(Leaf.length.pre) & !is.na(No.leaves.pre) &
      Leaf.length.pre > 0 & No.leaves.pre > 0
  ) %>%
  # Get rid of 2016 records because the sizes are not reliable
  filter(!(prev.year %in% 2016)) %>%
  # Add size columns
  mutate(size.prev = log(No.leaves.pre * Leaf.length.pre))

nrow(demo.surv.sizes)
table(demo.surv.sizes$surv, useNA = 'always')

# Subset for growth estimation
demo.grow = demo.surv.sizes %>% 
  filter(surv, !is.na(Leaf.length) & !is.na(No.leaves) & Leaf.length > 0 & No.leaves > 0) %>%
  mutate(size.cur = log(Leaf.length * No.leaves))

# Finally: subset surv dataset to not include 2023-2024 surv
demo.surv.sizes = demo.surv.sizes %>% filter(!(surv.year %in% 2024))

# demo.grow = merge(
#   demo.grow, phen.by.plant.for.growth, 
#   by.x = c('prev.year', 'plantid'), by.y = c('Year', 'plantid'),
#   all.x = TRUE, all.y = FALSE
# ) %>%
#   # Center mean around control
#   mutate(phen.c = phen.mean - phen.ctrl.mean)
