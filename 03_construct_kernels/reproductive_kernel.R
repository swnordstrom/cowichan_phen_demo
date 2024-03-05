#-------------------------------------------------------------------
# Script for fitting reproductive/fecundity kernel
# Involves multiple steps:
# - Probability of flowering
# - Number of umbels produced
# - Probability of umbels surviving
# - Seeds per umbel
# - Size distribution of new recruits
# Current plan is to use these to back-estimate seedling recruitment.
# Also doing these comparing fixed and random effect models (for annual variation)
# - SN init 5 Mar 2024
#-------------------------------------------------------------------


#-------------------------------------------------------------------
# Setup

##### Load packages

library(glmmTMB)
library(ggplot2)
library(lme4)
library(dplyr)
library(tidyr)
library(cowplot)

##### Clear namespace
rm(list = ls())

##### Load in and process data

### Flowering probability data (from processed demo, all years)
demo = read.csv('01_data_cleaning/out/demo_imputed_survival.csv')

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
  mutate(Year = factor(Year))

# Make a data frame with size - here, subset to include only plants for which we
# have a prior *and* a current size
demo.for.flowering.both.sizes = demo.for.flowering %>%
  # Filter out only plants with a previous measurement
  filter(!is.na(No.leaves.prev) & !is.na(Leaf.length.prev) & No.leaves.prev > 0) %>%
  # Filter out plants with a current demo measurement
  filter(!is.na(No.leaves) & !is.na(Leaf.length) & No.leaves > 0) %>%
  # Add size columns
  mutate(
    prev.size = log(No.leaves.prev * Leaf.length.prev),
    cur.size = log(No.leaves * Leaf.length)
  ) %>%
  select(-c(No.leaves, No.leaves.prev, Leaf.length, Leaf.length.prev)) %>%
  # Replace NAs with zero
  mutate(
    across(contains('umbel'), function(x) ifelse(is.na(x), 0, x)),
    prev.flower = No.umbels.prev > 0,
    flowering = No.umbels > 0
  )

# Here, subset the dataset to include only plants with a current size
demo.for.flowering = demo.for.flowering %>%
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

### Umbel phenology dataset

phen = merge(
    x = read.csv('01_data_cleaning/out/phenology_buds_deaths_all.csv'), 
    y = phen.demo %>% distinct(tagplot, plantid.phen),
    by.x = 'plantid', by.y = 'plantid.phen'
  ) %>%
  # Filter out only new umbels (i.e., buds)
  filter(varb %in% 'new') %>%
  select(-varb) %>%
  mutate(Year = factor(year)) %>%
  select(-year)
  
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


#-------------------------------------------------------------------
# Fit models for probability of flowering

# First check: which is better, using current size or prior size as a predictor

s_0 = glmer(
  formula = flowering ~ (1 | Plot / plantid) + (1 | Year),
  family = 'binomial',
  data = demo.for.flowering.both.sizes
)

s_cur = glmer(
  formula = flowering ~ cur.size + (1 | Plot / plantid) + (1 | Year),
  family = 'binomial',
  data = demo.for.flowering.both.sizes
)

s_pre = glmer(
  formula = flowering ~ prev.size + (1 | Plot / plantid) + (1 | Year),
  family = 'binomial',
  data = demo.for.flowering.both.sizes
)

AIC(s_0, s_cur, s_pre) %>% mutate(daic = round(AIC - min(AIC), 2))
# Use current size

# Compare fixed and random effects

s_ran = glmer(
  formula = flowering ~ cur.size + (1 | Plot / plantid) + (1 | Year),
  family = 'binomial',
  data = demo.for.flowering,
  control = glmerControl(optimizer = 'bobyqa')
)

s_fix = glmer(
  formula = flowering ~ cur.size + (1 | Plot / plantid) + Year,
  family = 'binomial',
  data = demo.for.flowering,
  control = glmerControl(optimizer = 'bobyqa')
)

AIC(s_fix, s_ran)
# Fixed effects are better (delta AIC ~19)
# Nervous about interactions, but we will see.

# Formulae to test
flow.mod.forms = c(
  'flowering ~ cur.size + Year',
  'flowering ~ cur.size * Year',
  'flowering ~ cur.size + trt + Year',
  'flowering ~ cur.size * trt + Year',
  'flowering ~ cur.size + trt * Year',
  'flowering ~ cur.size * trt + cur.size * Year',
  'flowering ~ cur.size * trt + trt * Year',
  'flowering ~ cur.size * trt + trt * Year + cur.size * Year',
  'flowering ~ cur.size * trt * Year'
) %>%
  paste('+ (1 | Plot / plantid)')

flow.mod.list = vector('list', length = length(flow.mod.forms))

for (j in 1:length(flow.mod.forms)) {
  print(j)
  flow.mod.list[[j]] = glmer(
    formula = flow.mod.forms[j],
    family = 'binomial',
    data = demo.for.flowering,
    control = glmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 1e6)) 
  )
}
# No warnings. Very cool.
# Final model here takes a looooong time to run though

data.frame(
  AIC = sapply(flow.mod.list, AIC) %>% unlist(),
  form = flow.mod.forms
) %>%
  mutate(daic = round(AIC - min(AIC), 2))
# Cool - so best model has all two-way interactions, but no three-way
# Still seems like a ton of parameters estimated.

f_sy_ty_st = glmer(
  formula = 'flowering ~ cur.size * trt + trt * Year + cur.size * Year + (1 | Plot / plantid)',
  family = 'binomial',
  data = demo.for.flowering,
  control = glmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 1e6)) 
)

summary(f_sy_ty_st)

# Try to plot results... this will be very funny.

range(demo.for.flowering$cur.size)

expand.grid(
  Year = factor(2018:2023),
  trt = c('control', 'drought', 'irrigated'),
  cur.size = (7:60)/10
) %>%
  mutate(
    pred = predict(
      f_sy_ty_st,
      newdata = ., type = 'response', 
      re.form = ~ 0, allow.new.levels = TRUE
    )
  ) %>%
  ggplot(aes(x = cur.size, y = pred, colour = trt)) +
  geom_line(linewidth = 1.5) +
  scale_colour_manual(values = c('black', 'red', 'blue')) +
  facet_wrap(~ Year)

# Jeez. That's not pretty.
# Seems like in most (but not all) years drought plants are less likely to flower.

# Will take some thought to think about how to get annual variation...

#-------------------------------------------------------------------
# Fit models for number of umbels per plant

# Here: need glmmTMB for zero-truncated model

# Test for effect of size

n_0 = glmmTMB(
  formula = No.umbels ~ (1 | Plot / plantid) + (1 | Year),
  family = 'truncated_poisson',
  data = demo.for.umbel.counts
)

n_s = glmmTMB(
  formula = No.umbels ~ cur.size +  (1 | Plot / plantid) + (1 | Year),
  family = 'truncated_poisson',
  data = demo.for.umbel.counts
)

AIC(n_s, n_0)
# Huge effect of size, unsurprisingly

# Compare fixed and random effects models

n_fix = glmmTMB(
  formula = No.umbels ~ cur.size + Year + (1 | Plot / plantid),
  family = 'truncated_poisson',
  data = demo.for.umbel.counts
)

AIC(n_fix, n_s)
# Fixed effects it is.

# Using a for loop won't work here... need to fit models separately.

# Treatment effect

n_s_y = glmmTMB(
  formula = No.umbels ~ cur.size + Year + (1 | Plot / plantid),
  family = 'truncated_poisson',
  data = demo.for.umbel.counts
) 

n_sy = glmmTMB(
  formula = No.umbels ~ cur.size * Year + (1 | Plot / plantid),
  family = 'truncated_poisson',
  data = demo.for.umbel.counts
)

AIC(n_sy, n_s_y)
# No size-year interaction (thus far)

n_s_t_y = glmmTMB(
  formula = No.umbels ~ cur.size + trt + Year + (1 | Plot / plantid),
  family = 'truncated_poisson',
  data = demo.for.umbel.counts
)

n_st_y = glmmTMB(
  formula = No.umbels ~ cur.size * trt + Year + (1 | Plot / plantid),
  family = 'truncated_poisson',
  data = demo.for.umbel.counts
)

n_s_ty = glmmTMB(
  formula = No.umbels ~ cur.size + trt * Year + (1 | Plot / plantid),
  family = 'truncated_poisson',
  data = demo.for.umbel.counts
)

AIC(n_s_y, n_s_t_y, n_st_y, n_s_ty) %>% mutate(daic = round(AIC - min(AIC), 2))
# Looks like a treatment-year effect

# Just to be sure, look at size-year effects in this model to
n_sy_ty = glmmTMB(
  formula = No.umbels ~ cur.size * Year + trt * Year + (1 | Plot / plantid),
  family = 'truncated_poisson',
  data = demo.for.umbel.counts
)

AIC(n_sy_ty, n_s_ty)
# Nope.

# Final model includes: size and treatment-year effects.
summary(n_s_ty)
# lmao only one of the year-treatment terms is significant and it's only
# marginal...

expand.grid(
  Year = factor(2018:2023),
  trt = c('control', 'drought', 'irrigated'),
  cur.size = (7:60)/10
) %>%
  mutate(
    pred = predict(
      n_s_ty,
      newdata = ., type = 'response', 
      re.form = ~ 0, allow.new.levels = TRUE
    )
  ) %>%
  ggplot(aes(x = cur.size, y = pred)) +
  geom_point(
    data = demo.for.umbel.counts,
    aes(x = cur.size, y = No.umbels, colour = trt),
    alpha = 0.1, size = 3
  ) +
  geom_line(aes(colour = trt)) +
  scale_colour_manual(values = c('black', 'red', 'blue')) +
  facet_wrap(~ Year) +
  theme(legend.position = 'none')

# Lmao no consistency in the rank ordering.
# Come on man.
# And 12 umbels? Come on man.
# Eight umbels for that one plant in 2018... what the heck man.


#-------------------------------------------------------------------
# Fit models for phen + umbel fate

# NOTE: I got a singularity warning when trying to include a random effect of
# tagplot

# First, look for size effect

u_0 = glmer(
  formula = cbind(n.succ.umbels, n.lost.umbels) ~ 
    (1 | Plot) + Year,
  family = 'binomial',
  data = umbel.fate
)

u_s = glmer(
  formula = cbind(n.succ.umbels, n.lost.umbels) ~ 
    cur.size + (1 | Plot) + Year,
  family = 'binomial',
  data = umbel.fate
)

AIC(u_s, u_0)
# Yes. There is a size effect. (Bleh)

# Test for random effect vs. fixed effect

u_fix =  glmer(
  formula = cbind(n.succ.umbels, n.lost.umbels) ~ 
    cur.size + (1 | Plot) + Year,
  family = 'binomial',
  data = umbel.fate
)

u_ran = glmer(
  formula = cbind(n.succ.umbels, n.lost.umbels) ~ 
    cur.size + (1 | Plot) + (1 | Year),
  family = 'binomial',
  data = umbel.fate
)

AIC(u_fix, u_ran)
# Slight (smaller than expected, ~ 5) advantage to fixed effects model

# Models to test:
# - size*year
# - treatment
# - size*treatment
# - year*treatment
# - phen

fate.mod.forms = c(
  'cbind(n.succ.umbels, n.lost.umbels) ~ cur.size + Year',
  'cbind(n.succ.umbels, n.lost.umbels) ~ cur.size * Year',
  'cbind(n.succ.umbels, n.lost.umbels) ~ cur.size + Year + trt',
  'cbind(n.succ.umbels, n.lost.umbels) ~ cur.size * Year + trt',
  'cbind(n.succ.umbels, n.lost.umbels) ~ cur.size + Year * trt',
  'cbind(n.succ.umbels, n.lost.umbels) ~ cur.size * trt + Year',
  'cbind(n.succ.umbels, n.lost.umbels) ~ cur.size * trt + Year * trt',
  'cbind(n.succ.umbels, n.lost.umbels) ~ cur.size * Year + Year * trt',
  'cbind(n.succ.umbels, n.lost.umbels) ~ cur.size * Year + Year * trt + cur.size * trt',
  'cbind(n.succ.umbels, n.lost.umbels) ~ cur.size * trt * Year'
) %>%
  paste('+ (1 | Plot)')

fate.mod.list = vector('list', length = length(fate.mod.forms))

for (j in 1:length(fate.mod.forms)) {
  print(j)
  fate.mod.list[[j]] = glmer(
    formula = fate.mod.forms[j],
    family = 'binomial',
    data = umbel.fate,
    control = glmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 1e6))
  )
}

data.frame(
  aic = sapply(fate.mod.list, AIC) %>% unlist(),
  form = fate.mod.forms
) %>%
  mutate(daic = round(aic - min(aic), 2)) %>%
  arrange(daic)

# Phew! Nice simple model.

# Now test for phen effects.

fate.phen.mod.forms = c(
  'cbind(n.succ.umbels, n.lost.umbels) ~ cur.size + Year * trt',
  'cbind(n.succ.umbels, n.lost.umbels) ~ cur.size + Year * trt + centered.phen',
  'cbind(n.succ.umbels, n.lost.umbels) ~ cur.size + Year * trt + centered.phen * trt',
  'cbind(n.succ.umbels, n.lost.umbels) ~ cur.size + Year * trt + centered.phen * Year',
  'cbind(n.succ.umbels, n.lost.umbels) ~ cur.size + Year * trt * centered.phen'
) %>%
  paste('+ (1 | Plot)')

fate.phen.mod.list = vector('list', length = length(fate.phen.mod.forms))

for (j in 1:length(fate.phen.mod.forms)) {
  print(j)
  fate.phen.mod.list[[j]] = glmer(
    formula = fate.phen.mod.forms[j],
    family = 'binomial',
    data = umbel.fate,
    control = glmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 1e6))
  )
}

data.frame(
  aic = sapply(fate.phen.mod.list, AIC) %>% unlist(),
  form = fate.phen.mod.forms
) %>%
  mutate(daic = round(aic - min(aic), 2)) %>%
  arrange(daic)

# Okay. Going to take the phen-year interaction only, the three-way support is weak

u_s_yt_py = glmer(
  formula = cbind(n.succ.umbels, n.lost.umbels) ~ (1 | Plot) +
    cur.size + Year * trt + centered.phen * Year,
  family = 'binomial',
  data = umbel.fate,
  control = glmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 1e6))
)

summary(u_s_yt_py)
# okay, definitely some year-treatment effects here

expand.grid(
  cur.size = 3.7,
  trt = c('control', 'drought', 'irrigated'),
  Year = factor(2021:2023),
  centered.phen = -30:30
) %>%
  mutate(
    pred = predict(
      u_s_yt_py,
      newdata = ., type = 'response',
      re.form = ~ 0, allow.new.levels = TRUE
    )
  ) %>%
  ggplot(aes(x = centered.phen, y = pred)) +
  geom_line(aes(colour = trt)) +
  scale_colour_manual(values = c('black', 'red', 'blue')) +
  facet_wrap(~ Year)
# Seems mostly favorable to be in an irrigated plot, drought more likely to lead
# to failures

# What about the effect of size?
expand.grid(
  cur.size = (18:56)/10,
  trt = c('control', 'drought', 'irrigated'),
  Year = factor(2021:2023),
  centered.phen = 0
) %>%
  mutate(
    pred = predict(
      u_s_yt_py,
      newdata = ., type = 'response',
      re.form = ~ 0, allow.new.levels = TRUE
    )
  ) %>%
  ggplot(aes(x = cur.size, y = pred)) +
  geom_line(aes(colour = trt)) +
  scale_colour_manual(values = c('black', 'red', 'blue')) +
  facet_wrap(~ Year)
# Considerable effect of size here!
# Smaller plants much more likely to have umbel death.

# (Note: did not test for phen-size interaction)
# (That's too much stuff!)


#-------------------------------------------------------------------
# Fit models for seeds per umbel

