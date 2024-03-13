# ------------------------------------------------
# Script for constructing annually-varying kernels (Here: only fitting optimal
# model; this is downstream of the model-fitting process)
# Also assessing spatial and temporal covariation in vital rates (as much as I
# can with what little data is vailable...)
# - SN init 12 Mar 2024
# ------------------------------------------------

# Load packages

library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyr)
library(lme4)
library(glmmTMB)

rm(list = ls())

# ------------------------------------------------
# Load in and set up datasets

# Load in processed demo data
demo = read.csv('01_data_cleaning/out/demo_imputed_survival.csv')

### Survival dataset 

# Create demo dataset for survival
# we'll say that survival in year $t$ means surviving from $t$ to $t+1$
demo.for.surv = merge(
  # x - status in year $t$ - only interested in plants that survived in t
  x = demo %>% 
    filter(surv) %>% 
    select(Year, No.leaves, Leaf.length, No.umbels, demo.note, proc.note, plantid, Plot, trt),
  # y - whether plant was observed alive in year $t+1$
  y = demo %>%
    mutate(p.year = Year - 1) %>%
    select(p.year, plantid, surv, No.leaves, Leaf.length),
  by.x = c('Year', 'plantid'), by.y = c('p.year', 'plantid'),
  # suffixes are t and t+1
  suffixes = c('.t', '.tp1')
) %>%
  # Remove plants iwth no measurements and plants in 2016
  filter(!is.na(No.leaves.t) & !is.na(Leaf.length.t) & No.leaves.t > 0, Year > 2016) %>%
  # add a size column
  mutate(size.t = log(No.leaves.t * Leaf.length.t)) %>%
  # Change year to factor
  mutate(Year = factor(Year))

### Growth rates 

# Get a dataset of demo for growth - subsetting the data frame above
# to get this, we need plants that survived *and* have measurements
demo.for.growth = demo.for.surv %>%
  filter(!is.na(No.leaves.tp1) & !is.na(Leaf.length.tp1) & No.leaves.tp1 > 0) %>%
  # add a size column
  mutate(size.tp1 = log(No.leaves.tp1 * Leaf.length.tp1))

### Probability of flowering

# Demo for flowering dataset
demo.for.flowering = merge(
  x = demo %>%
    # Filter out only surviving plants
    filter(surv) %>%
    # Filter out 2016 plants (size measurements are not reliable)
    filter(Year > 2016) %>%
    # Select columns for year-of data
    select(Year, plantid, Plot, trt, No.umbels, No.leaves, Leaf.length),
  y = demo %>%
    # Filter out only surviving plants
    filter(surv) %>%
    # Filter out 2016 plants (size measurements are not reliable)
    filter(Year > 2016) %>%
    # Add column for year matching
    mutate(Year.match = Year + 1) %>% 
    # Select relevant columns
    select(Year.match, plantid, No.umbels, No.leaves, Leaf.length),
  by.x = c('Year', 'plantid'), by.y = c('Year.match', 'plantid'),
  suffixes = c('', '.prev')
) %>%
  mutate(Year = factor(Year)) %>%
  # Filter out plants with a current demo measurement
  filter(!is.na(No.leaves) & !is.na(Leaf.length) & No.leaves > 0) %>%
  # Add size columns
  mutate(cur.size = log(No.leaves * Leaf.length)) %>%
  select(-c(No.leaves, No.leaves.prev, Leaf.length, Leaf.length.prev)) %>%
  # Replace NAs with zero
  mutate(
    across(contains('umbel'), function(x) ifelse(is.na(x), 0, x)),
    prev.flower = No.umbels.prev > 0,
    flowering = No.umbels > 0
  )

### Distribution of umbels per flowering plant

demo.for.umbel.counts = demo.for.flowering %>%
  filter(flowering) %>%
  select(-flowering)

### Probability of umbel surviving

# Load in phen-demo dataset (for plants matched in both phenology and demography
# data, 2021-2023 so far)
phen.demo = read.csv('01_data_cleaning/out/phen_demo_for_umbel_survival.csv') %>%
  # Change missing umbel counts to zero
  mutate(n.demo.umbels = ifelse(is.na(n.demo.umbels), 0, n.demo.umbels)) %>%
  # Center phenology
  mutate(centered.phen = init.doy - round(mean(init.doy))) %>%
  mutate(Year = factor(Year))

# Get a dataset for umbel success, with mean bud date and number of
# living/dead buds per plant
umbel.fate = phen.demo %>%
  # Get mean budding date for each plant
  group_by(Year, tagplot) %>%
  mutate(mean.doy = mean(init.doy)) %>%
  ungroup() %>%
  # Distinct here because we need only one row per plant
  distinct(Year, tagplot, .keep_all = TRUE) %>%
  # Center the phenology column to help with model convergence
  mutate(centered.phen = mean.doy - round(mean(mean.doy))) %>%
  # Make successful umbel column
  mutate(n.succ.umbels = n.phen.umbels - n.lost.umbels) %>%
  # Year to factor
  mutate(Year = factor(Year))

### Seed set per umbel

# Need to lead in seed-phen-demo dataset, which includes one row per umbel on
# plants appearing in each of seed, phen, and demo (2021 and later) - although
# I'll filter out plants not matched to all datasets just to be safe
seed.phen.demo = read.csv('01_data_cleaning/out/seed_phen_demo_combined.csv') %>%
  # Remove plants that are missing phen
  filter(!is.na(mean.bud.day)) %>%
  # Remove plants that are missing from seed
  filter(!is.na(no.seeds)) %>%
  # Add centered phen column
  mutate(centered.phen = mean.bud.day - round(mean(mean.bud.day))) %>%
  mutate(Year = factor(Year))

# Removing the lost umbels (accounted for in the umbel fate model)
seed.phen.demo = rbind(
  # All umbels with a non-zero number of seeds...
  seed.phen.demo %>% filter(no.seeds > 0),
  # and all umbels with zero seeds, but removing dead umbels from phen
  seed.phen.demo %>%
    filter(!no.seeds) %>%
    group_by(tagplot, Year) %>%
    filter((1:n()) > n.lost.umbel) %>%
    ungroup()
)

### New recruit size distribution

# Subset the demo dataset to include 
# - plants seen for the first time (post-2017)
# - with one leaf only
# - and not flowering

recruits = demo %>%
  # Get plants seen only for the first time
  filter(!duplicated(plantid)) %>%
  # Remove 2016 and 2017 plants
  filter(Year > 2017) %>%
  # Get plants with only one leaf
  filter(No.leaves < 2, !is.na(No.leaves)) %>%
  # Get plants that are not flowering
  filter(!No.umbels | is.na(No.umbels)) %>%
  mutate(Year = factor(Year))

# ------------------------------------------------
# Fit models

### Survival model

surv.mod = glmer(
  formula = surv ~ size.t + (1 | Year) + (1 | Plot / plantid),
  family = 'binomial',
  data = demo.for.surv
)

### Growth model

grow.mod = lmer(
  size.tp1 ~ size.t + trt + (1 | Year) + (1 | Year:trt) + (1 | Plot / plantid),
  data = demo.for.growth,
  control = lmerControl(optimizer = 'bobyqa')
)

### Probability of flowering model

flow.mod = glmer(
  formula = flowering ~ cur.size * trt + (1 | Year) + (1 | Year:trt) + (1 | Plot / plantid),
  family = 'binomial',
  data = demo.for.flowering,
  control = glmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 1e6)) 
)

### Number of umbels model

numb.mod = glmmTMB(
  formula = No.umbels ~ cur.size + trt +  (1 | Year) + (1 | trt:Year) + (1 | Plot / plantid),
  family = 'truncated_poisson',
  data = demo.for.umbel.counts
)

### Umbel survival model

fate.mod = glmer(
  formula = cbind(n.succ.umbels, n.lost.umbels) ~ (1 | Plot) +
    n.phen.umbels + cur.size * trt + trt * Year + Year * centered.phen + trt * centered.phen,
  family = 'binomial',
  data = umbel.fate,
  control = glmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 1e6))
)

### Seeds per umbel model

seed.mod = glmer.nb(
  no.seeds ~ cur.size * trt + trt * Year + centered.phen + (1 | Plot),
  data = seed.phen.demo %>% 
    mutate(
      cur.size = cur.size - mean(cur.size),
      centered.phen = centered.phen / 7
    ),
  control = glmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 1e7))
)

### Recruit size distribution model

recr.mod = lmer(
  log(Leaf.length) ~ (1 | Year) + (1 | Plot),
  data = recruits
)

# ------------------------------------------------
# Assess correlations of random effects

##### Year-varying random effects

surv.yr.rfx = data.frame(int = ranef(surv.mod)$Year) %>%
  mutate(year = row.names(.)) %>%
  rename(int = X.Intercept.)

grow.yr.rfx = data.frame(int = ranef(grow.mod)$Year) %>%
  mutate(year = row.names(.)) %>%
  rename(int = X.Intercept.)

flow.yr.rfx = data.frame(int = ranef(flow.mod)$Year) %>%
  mutate(year = row.names(.)) %>%
  rename(int = X.Intercept.)

numb.yr.rfx = data.frame(int = ranef(numb.mod)$cond$Year) %>%
  mutate(year = row.names(.)) %>%
  rename(int = X.Intercept.)

fate.yr.ffx = data.frame(int = fixef(fate.mod)) %>%
  mutate(term = row.names(.)) %>%
  filter(grepl('^Year20\\d\\d$', term)) %>%
  rbind(data.frame(int = 0, term = 'Year2021'), .) %>%
  mutate(year = gsub('Year', '', term)) %>%
  select(int, year)

seed.yr.ffx = data.frame(int = fixef(seed.mod)) %>%
  mutate(term = row.names(.)) %>%
  filter(grepl('^Year20\\d\\d$', term)) %>%
  rbind(data.frame(int = 0, term = 'Year2021'), .) %>%
  mutate(year = gsub('Year', '', term)) %>%
  select(int, year)

recr.yr.rfx = data.frame(int = ranef(recr.mod)$Year) %>%
  mutate(year = row.names(.)) %>%
  rename(int = X.Intercept.)

all.yr.efx = rbind(
  surv.yr.rfx %>% mutate(rate = 'surv'),
  grow.yr.rfx %>% mutate(rate = 'grow'),
  flow.yr.rfx %>% mutate(rate = 'flow'),
  fate.yr.ffx %>% mutate(rate = 'fate'),
  seed.yr.ffx %>% mutate(rate = 'seed'),
  recr.yr.rfx %>% mutate(rate = 'recr')
)

row.names(all.yr.efx) = NULL

all.yr.efx %>%
  ggplot(aes(x = year, y = int, group = rate, colour = rate)) +
  geom_line() +
  geom_point(size = 3)

long.term.rates = all.yr.efx %>%
  filter(rate %in% c('surv', 'grow', 'flow', 'recr')) %>%
  # Maybe it makes sense to stagger the surv + growth rates...
  # convert year column to numeic
  mutate(year = as.numeric(year)) %>%
  # add one to year for survival and growth
  mutate(year = year + as.numeric(rate %in% c('surv', 'grow')))

long.term.rates %>%
  ggplot(aes(x = year, y = int, colour = rate, int = rate)) +
  geom_point(size = 3) +
  geom_line()
# Looks kinda like positive correlation?

long.term.corrs = long.term.rates %>%
  pivot_wider(names_from = rate, values_from = int) %>%
  select(-year) %>%
  cor()
# obviously very little data to inform these estimates, but...
# survival negatively correlated with growth and flowering,
# growth positively correlated with flowering and recruit size
# recruit size uncorrelated with flowering, weak corr. with survival

##### Spatial random effects

surv.pl.rfx = data.frame(int = ranef(surv.mod)$Plot) %>%
  mutate(plot = row.names(.)) %>%
  rename(int = X.Intercept.)

grow.pl.rfx = data.frame(int = ranef(grow.mod)$Plot) %>%
  mutate(plot = row.names(.)) %>%
  rename(int = X.Intercept.)

flow.pl.rfx = data.frame(int = ranef(flow.mod)$Plot) %>%
  mutate(plot = row.names(.)) %>%
  rename(int = X.Intercept.)

numb.pl.rfx = data.frame(int = ranef(numb.mod)$cond$Plot) %>%
  mutate(plot = row.names(.)) %>%
  rename(int = X.Intercept.)

fate.pl.ffx = data.frame(int = ranef(fate.mod)$Plot) %>%
  mutate(plot = row.names(.)) %>%
  rename(int = X.Intercept.)

seed.pl.ffx = data.frame(int = ranef(seed.mod)$Plot) %>%
  mutate(plot = row.names(.)) %>%
  rename(int = X.Intercept.)

recr.pl.rfx = data.frame(int = ranef(recr.mod)$Plot) %>%
  mutate(plot = row.names(.)) %>%
  rename(int = X.Intercept.)

all.pl.efx = rbind(
  surv.pl.rfx %>% mutate(rate = 'surv'),
  grow.pl.rfx %>% mutate(rate = 'grow'),
  flow.pl.rfx %>% mutate(rate = 'flow'),
  fate.pl.ffx %>% mutate(rate = 'fate'),
  seed.pl.ffx %>% mutate(rate = 'seed'),
  recr.pl.rfx %>% mutate(rate = 'recr')
) %>%
  group_by(rate) %>%
  mutate(across(int, function(x) x/sd(x)))

all.pl.efx %>%
#   pivot_wider(names_from = rate, values_from = int) %>%
  ggplot(aes(x = plot, y = rate, fill = int)) +
  geom_raster() +
  scale_fill_gradient2(low = 'red', high = 'royalblue', mid = 'white', midpoint = 0)

all.pl.efx %>%
  pivot_wider(names_from = rate, values_from = int) %>%
  select(-plot) %>%
  cor() %>%
  round(2)
