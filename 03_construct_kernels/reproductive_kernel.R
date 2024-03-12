#-------------------------------------------------------------------
# Script for fitting reproductive/fecundity kernel
# Involves multiple steps:
# - Probability of flowering
# - Number of umbels produced
# - Probability of umbels surviving
# - Seeds per umbel
# - Size distribution of new recruits
# Current plan is to use these to back-estimate seedling recruitment.
# Using rule: fixed effects for models with only three years of data, random
# effects if we have ~6-7 years of data
# (because of this, not testing for year-treatment interactions)
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

# Formulae to test
flow.mod.forms = c(
  'flowering ~ cur.size + (1 | Year)',
  'flowering ~ (cur.size | Year)',
  'flowering ~ cur.size + trt +  (1 | Year)',
  'flowering ~ cur.size * trt +  (1 | Year)',
  'flowering ~ trt +  (cur.size | Year)', # adding trt-year interaction causes singularity
  'flowering ~ cur.size + trt +  (1 | Year) + (1 | Year:trt)',
  'flowering ~ cur.size * trt +  (1 | Year) + (1 | Year:trt)'
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

data.frame(
  AIC = sapply(flow.mod.list, AIC) %>% unlist(),
  form = flow.mod.forms
) %>%
  mutate(daic = round(AIC - min(AIC), 2))
# Of these models, best has size-treatment interaction
# *and* a year-treatment effect

# Best fixed effects model:
# f_sy_ty_st = glmer(
#   formula = 'flowering ~ cur.size * trt + trt * Year + cur.size * Year + (1 | Plot / plantid)',
#   family = 'binomial',
#   data = demo.for.flowering,
#   control = glmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 1e6)) 
# )

# f_st_y = glmer(
#   formula = flowering ~ cur.size * trt +  (1 | Year) + (1 | Plot / plantid),
#   family = 'binomial',
#   data = demo.for.flowering,
#   control = glmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 1e6)) 
# )

f_st_ty = glmer(
  formula = flowering ~ cur.size * trt + (1 | Year) + (1 | Year:trt) + (1 | Plot / plantid),
  family = 'binomial',
  data = demo.for.flowering,
  control = glmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 1e6)) 
)

summary(f_st_ty)
# year-treatment effect is pretty small

# Try to plot results... this will be very funny.

range(demo.for.flowering$cur.size)

expand.grid(
  Year = factor(2018:2023),
  trt = c('control', 'drought', 'irrigated'),
  cur.size = (7:60)/10
) %>%
  mutate(
    pred = predict(
      f_st_ty,
      newdata = ., type = 'response', 
      re.form = ~ (1 | Year), allow.new.levels = TRUE
    )
  ) %>%
  ggplot(aes(x = cur.size, y = pred, colour = trt)) +
  geom_line(linewidth = 1.5) +
  scale_colour_manual(values = c('black', 'red', 'blue')) +
  facet_wrap(~ Year)

# Interesting
# Looks like there is a slightly higher probability of flowering for small plants in drought?
# offset by less flowering for other plants though...

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

# n_fix = glmmTMB(
#   formula = No.umbels ~ cur.size + Year + (1 | Plot / plantid),
#   family = 'truncated_poisson',
#   data = demo.for.umbel.counts
# )

# AIC(n_fix, n_s)
# Fixed effects it is.

# Using a for loop won't work here... need to fit models separately.

# Treatment effect

n_s_y = glmmTMB(
  formula = No.umbels ~ cur.size + (1 | Year) + (1 | Plot / plantid),
  family = 'truncated_poisson',
  data = demo.for.umbel.counts
) 

n_sy = glmmTMB(
  formula = No.umbels ~  (cur.size | Year) + (1 | Plot / plantid),
  family = 'truncated_poisson',
  data = demo.for.umbel.counts
)

AIC(n_sy, n_s_y)
# No size-year interaction (thus far)

n_s_t_y = glmmTMB(
  formula = No.umbels ~ cur.size + trt +  (1 | Year) + (1 | Plot / plantid),
  family = 'truncated_poisson',
  data = demo.for.umbel.counts
)

n_st_y = glmmTMB(
  formula = No.umbels ~ cur.size * trt +  (1 | Year) + (1 | Plot / plantid),
  family = 'truncated_poisson',
  data = demo.for.umbel.counts
)

n_s_ty = glmmTMB(
  formula = No.umbels ~ cur.size + trt +  (1 | Year) + (1 | trt:Year) + (1 | Plot / plantid),
  family = 'truncated_poisson',
  data = demo.for.umbel.counts
)

n_st_ty = glmmTMB(
  formula = No.umbels ~ cur.size * trt +  (1 | Year) + (1 | trt:Year) + (1 | Plot / plantid),
  family = 'truncated_poisson',
  data = demo.for.umbel.counts
)

AIC(n_s_y, n_s_t_y, n_st_y, n_s_ty, n_st_ty) %>% mutate(daic = round(AIC - min(AIC), 2))
# Looks like a treatment-year effect

# Final model includes: size and treatment-year effects.
# summary(n_s_y)
# # lmao only one of the year-treatment terms is significant and it's only
# # marginal...
summary(n_s_ty)
# positive effects of both treatments? on average at least
ranef(n_s_ty)$cond$`trt:Year` %>%
  mutate(param = row.names(.)) %>%
  separate(param, sep = ':', into = c('trt', 'Year')) %>%
  rename(estimate = `(Intercept)`) %>%
  ggplot(aes(x = Year, y = estimate, group = trt, colour = trt)) +
  geom_line() +
  scale_colour_manual(values = c('black', 'red', 'blue'))

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

# 


#-------------------------------------------------------------------
# Fit models for phen + umbel fate

# NOTE: I got a singularity warning when trying to include a random effect of
# tagplot
# Also: using fixed year effects due to low replication

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


# Account for number of umbels produced? (After accounting for size)
# These should be correlated so models may not play nice...
u_n_s = glmer(
  formula = cbind(n.succ.umbels, n.lost.umbels) ~ 
    n.phen.umbels + cur.size + (1 | Plot) + Year,
  family = 'binomial',
  data = umbel.fate,
  control = glmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 1e6))
)

AIC(u_s, u_n_s)
summary(u_n_s)
# Yes. Producing more umbels means your umbels are less likely to succeed.

# Models to test:
# - size*year
# - treatment
# - size*treatment
# - year*treatment
# - phen

fate.mod.forms = c(
  'cbind(n.succ.umbels, n.lost.umbels) ~ n.phen.umbels + cur.size + Year',
  'cbind(n.succ.umbels, n.lost.umbels) ~ n.phen.umbels + cur.size * Year',
  'cbind(n.succ.umbels, n.lost.umbels) ~ n.phen.umbels + cur.size + Year + trt',
  'cbind(n.succ.umbels, n.lost.umbels) ~ n.phen.umbels + cur.size * Year + trt',
  'cbind(n.succ.umbels, n.lost.umbels) ~ n.phen.umbels + cur.size * trt + Year',
  'cbind(n.succ.umbels, n.lost.umbels) ~ n.phen.umbels + cur.size * trt + Year * trt',
  'cbind(n.succ.umbels, n.lost.umbels) ~ n.phen.umbels + cur.size * Year + Year * trt',
  'cbind(n.succ.umbels, n.lost.umbels) ~ n.phen.umbels + cur.size * trt + cur.size * Year',
  'cbind(n.succ.umbels, n.lost.umbels) ~ n.phen.umbels + cur.size * trt + cur.size * Year + Year * trt'
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
  form = fate.mod.forms,
  i = 1:length(fate.mod.forms)
) %>%
  mutate(daic = round(aic - min(aic), 2)) %>%
  arrange(daic)

# Best model here has size-treatment and year-size

# Now test for phen effects.

fate.phen.mod.forms = c(
  'cbind(n.succ.umbels, n.lost.umbels) ~ n.phen.umbels + cur.size * trt + trt * Year',
  'cbind(n.succ.umbels, n.lost.umbels) ~ n.phen.umbels + cur.size * trt + trt * Year + centered.phen',
  'cbind(n.succ.umbels, n.lost.umbels) ~ n.phen.umbels + cur.size * trt + trt * Year + centered.phen * trt',
  'cbind(n.succ.umbels, n.lost.umbels) ~ n.phen.umbels + cur.size * trt + trt * Year + centered.phen * Year',
  'cbind(n.succ.umbels, n.lost.umbels) ~ n.phen.umbels + cur.size * trt + trt * Year + centered.phen * Year + trt * centered.phen'
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
  form = fate.phen.mod.forms,
  i = 1:length(fate.phen.mod.forms)
) %>%
  mutate(daic = round(aic - min(aic), 2)) %>%
  arrange(daic)
# Best model by far has year- and treatment-varying effects of phenology

u_n_st_yt_py_tp = glmer(
  formula = cbind(n.succ.umbels, n.lost.umbels) ~ (1 | Plot) +
    n.phen.umbels + cur.size * trt + trt * Year + Year * centered.phen + trt * centered.phen,
  family = 'binomial',
  data = umbel.fate,
  control = glmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 1e6))
)

u_n_st_yt_py = glmer(
  formula = cbind(n.succ.umbels, n.lost.umbels) ~ (1 | Plot) +
    n.phen.umbels + cur.size * trt + trt * Year + Year * centered.phen,
  family = 'binomial',
  data = umbel.fate,
  control = glmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 1e6))
)

summary(u_n_st_yt_py_tp)
# lots going on in here...

expand.grid(
  cur.size = 3.7,
  Year = factor(2021:2023),
  trt = c('control', 'drought', 'irrigated'),
  centered.phen = -30:30,
  n.phen.umbels = 1
) %>%
  mutate(
    pred = predict(
      u_n_st_yt_py_tp,
      newdata = ., type = 'response',
      re.form = ~ 0, allow.new.levels = TRUE
    ),
    pred2 = predict(
      u_n_st_yt_py,
      newdata = ., type = 'response',
      re.form = ~ 0, allow.new.levels = TRUE
    )
  ) %>%
  pivot_longer(c(pred, pred2), names_to = 'model', values_to = 'estimate') %>%
  ggplot(aes(x = centered.phen, y = estimate, colour = trt)) +
  geom_point(
    data = umbel.fate,
    aes(x = centered.phen, y = (n.succ.umbels/n.phen.umbels), colour = trt),
    inherit.aes = FALSE,
    size = 3, alpha = 0.1
  ) +
  geom_line(aes(linetype = model)) +
  scale_colour_manual(values = c('black', 'red', 'blue')) +
  facet_wrap(trt ~ Year)
# I think the more complicated model is okay here

# Now this is interesting...
# Drought plants *much* more likely to have umbel death in 2021,
# phen effects in drought are less severe in 2022
# and opposite trend in 2023?

# # What about the effect of size?
# expand.grid(
#   cur.size = (18:56)/10,
#   Year = factor(2021:2023),
#   centered.phen = 0
# ) %>%
#   mutate(
#     pred = predict(
#       u_s_y_py,
#       newdata = ., type = 'response',
#       re.form = ~ 0, allow.new.levels = TRUE
#     )
#   ) %>%
#   ggplot(aes(x = cur.size, y = pred)) +
#   geom_line() +
#   scale_colour_manual(values = c('black', 'red', 'blue')) +
#   facet_wrap(~ Year)
# # Considerable effect of size here!
# # Smaller plants much more likely to have umbel death.

# (Note: did not test for phen-size interaction)
# (That's too much stuff!)


#-------------------------------------------------------------------
# Fit models for seeds per umbel

# Singularity issue when trying random effect for plant within plot
# Also doing only fixed effects models here (too few years)

# Test for effect if size
s_0 = glmer.nb(
  formula = no.seeds ~ (1 | Plot) + Year,
  data = seed.phen.demo
)

s_s = glmer.nb(
  formula = no.seeds ~ cur.size + (1 | Plot) + Year,
  data = seed.phen.demo
)

AIC(s_0, s_s)
# Yep - sizeable effect of size

# Models to test:
seed.mod.forms = c(
  'no.seeds ~ cur.size + Year',
  'no.seeds ~ cur.size * Year',
  'no.seeds ~ cur.size + Year + trt',
  'no.seeds ~ cur.size * Year + trt',
  'no.seeds ~ cur.size * trt + Year',
  'no.seeds ~ cur.size + Year * trt',
  'no.seeds ~ cur.size * Year + Year * trt',
  'no.seeds ~ cur.size * trt + cur.size * Year',
  'no.seeds ~ cur.size * trt + cur.size * Year +  Year * trt'
) %>%
  paste('+ (1 | Plot)')

seed.mod.list = vector('list', length = length(seed.mod.forms))

for (j in 1:length(seed.mod.list)) {
  print(j)
  seed.mod.list[[j]] = glmer.nb(
    formula = seed.mod.forms[j],
    data = seed.phen.demo %>% mutate(cur.size = cur.size - mean(cur.size)),
    control = glmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 1e6))
  )
}
# Note: want to include centering in here or else these are slow and don't converge

data.frame(
  aic = sapply(seed.mod.list, AIC) %>% unlist(),
  form = seed.mod.forms,
  i = 1:length(seed.mod.forms)
) %>%
  mutate(daic = round(aic - min(aic), 2)) %>%
  arrange(daic)

# Of course, the complicated model is the best-supported one

s_st_sy_ty = glmer.nb(
  no.seeds ~ cur.size * trt + cur.size * Year + Year * trt + (1 | Plot),
  data = seed.phen.demo %>% mutate(cur.size = cur.size - mean(cur.size)),
  control = glmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 1e6))
)

s_st_sy = glmer.nb(
  no.seeds ~ cur.size * trt + cur.size * Year + (1 | Plot),
  data = seed.phen.demo %>% mutate(cur.size = cur.size - mean(cur.size)),
  control = glmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 1e6))
)

s_st_ty = glmer.nb(
  no.seeds ~ cur.size * trt + trt * Year + (1 | Plot),
  data = seed.phen.demo %>% mutate(cur.size = cur.size - mean(cur.size)),
  control = glmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 1e6))
)


# Cool...
expand.grid(
  Year = factor(2021:2023),
  cur.size = (-15:15)/10,
  trt = c('control', 'drought', 'irrigated')
) %>%
  # mutate(
  #   seed = predict(
  #     object = s_st_sy, newdata = ., type = 'response',
  #     re.form = ~ 0, allow.new.levels = TRUE
  #   )
  # ) %>%
  mutate(
    seed3 = predict(
      object = s_st_sy_ty, newdata = ., type = 'response',
      re.form = ~ 0, allow.new.levels = TRUE,
    ),
    seed2 = predict(
      object = s_st_sy, newdata = ., type = 'response',
      re.form = ~ 0, allow.new.levels = TRUE,
    ),
    seed1 = predict(
      object = s_st_ty, newdata = ., type = 'response',
      re.form = ~ 0, allow.new.levels = TRUE,
    )
  ) %>%
  pivot_longer(c(seed1, seed2, seed3), names_to = 'model', values_to = 'seed') %>%
  mutate(cur.size = cur.size + mean(seed.phen.demo$cur.size)) %>%
  ggplot(aes(x = cur.size, y = seed, colour = trt)) +
  geom_point(
    data = seed.phen.demo,
    aes(x = cur.size, y = no.seeds, colour = trt),
    inherit.aes = FALSE,
    alpha = 0.1, size = 3
  ) +
  geom_line(aes(linetype = model)) +
  scale_colour_manual(values = c('black', 'red', 'blue')) +
  facet_wrap(trt ~ Year) +
  theme(legend.position = 'bottom')
# That 2021-control curve is a lot of extrapolation...
# seed2 is maybe a better fit for drought but it's so little data
# I really prefer the `seed1` model

seed.phen.mod.forms = c(
  'no.seeds ~ cur.size * trt + trt * Year + (1 | Plot)',
  'no.seeds ~ cur.size * trt + trt * Year + centered.phen + (1 | Plot)',
  'no.seeds ~ cur.size * trt + trt * Year + centered.phen * trt + (1 | Plot)',
  'no.seeds ~ cur.size * trt + trt * Year + centered.phen * Year + (1 | Plot)',
  'no.seeds ~ cur.size * trt + trt * Year + centered.phen * trt + centered.phen * Year + (1 | Plot)',
  'no.seeds ~ cur.size * trt + trt * Year + poly(centered.phen, 2) + (1 | Plot)',
  'no.seeds ~ cur.size * trt + trt * Year + poly(centered.phen, 2) * trt + (1 | Plot)'# ,
  # 'no.seeds ~ cur.size * trt + poly(centered.phen, 2) * Year + (1 | Plot)' # does not converge
)

seed.phen.mod.list = vector('list', length = length(seed.phen.mod.forms))

for (j in 1:length(seed.phen.mod.list)) {
  print(j)
  seed.phen.mod.list[[j]] = glmer.nb(
    formula = seed.phen.mod.forms[j],
    data = seed.phen.demo %>% 
      mutate(
        cur.size = cur.size - mean(cur.size),
        centered.phen = centered.phen / 7
      ),
    control = glmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 1e6))
  )
}

data.frame(
  aic = sapply(seed.phen.mod.list, AIC) %>% unlist(),
  form = seed.phen.mod.forms,
  i = 1:length(seed.phen.mod.forms)
) %>%
  mutate(daic = round(aic - min(aic), 2)) %>%
  arrange(daic)
# best model here has a linear effect of phen (nice)
# I'd consider the phen-treatment effect to be weakly supported only

s_st_ty_p = glmer.nb(
  no.seeds ~ cur.size * trt + trt * Year + centered.phen + (1 | Plot),
  data = seed.phen.demo %>% 
    mutate(
      cur.size = cur.size - mean(cur.size),
      centered.phen = centered.phen / 7
    ),
  control = glmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 1e7))
)

expand.grid(
  Year = factor(2021:2023),
  cur.size = c(-1.5, 0, 1.5),
  trt = c('control', 'drought', 'irrigated'),
  centered.phen = (-25:25)/10
) %>%
  mutate(
    seed = predict(
      object = s_st_ty_p, newdata = ., type = 'response',
      re.form = ~ 0, allow.new.levels = TRUE, 
    )
  ) %>%
  mutate(cur.size = factor(cur.size)) %>%
  ggplot(aes(x = centered.phen, y = seed, colour = trt)) +
  geom_point(
    data = seed.phen.demo %>% mutate(centered.phen = centered.phen/7),
    aes(x = centered.phen, y = no.seeds, colour = trt),
    inherit.aes = FALSE,
    alpha = 0.1, size = 3
  ) +
  geom_line(aes(linetype = cur.size)) +
  scale_colour_manual(values = c('black', 'red', 'blue')) +
  scale_linetype_manual(values = c(2, 1, 2)) +
  facet_wrap(trt ~ Year)


#-------------------------------------------------------------------
# Model for distribution of seedling sizes

r_0 = lmer(
  log(Leaf.length) ~ (1 | Plot),
  data = recruits
)

r_fixed = lmer(
  log(Leaf.length) ~ Year + (1 | Plot),
  data = recruits
)

r_ranef = lmer(
  log(Leaf.length) ~ (1 | Year) + (1 | Plot),
  data = recruits
)

AIC(r_0, r_fixed, r_ranef)

# Random effects model is best

r_t_y = lmer(
  log(Leaf.length) ~ trt + (1 | Year) + (1 | Plot),
  data = recruits
)

r_ty = lmer(
  log(Leaf.length) ~ trt + (1 | Year) + (1 | Year:trt) + (1 | Plot),
  data = recruits
)

AIC(r_ty, r_t_y, r_ranef)
# No treatment effect.

summary(r_ranef)
hist(residuals(r_ranef))
# one extreme negartive outlier, otherwise could be worse
qqnorm(residuals(r_ranef))
qqline(residuals(r_ranef))

summary(r_ranef)$sigma

#-------------------------------------------------------------------
# Construct a reproductive kernel
# Start with size only (no phen effects yet)

# # First, get year-coefficient effects for appropriate models
# # (This may not be necessary... might be easier to get estimates on the linear
# # scale?)
# 
# # For umbel fate, mod is u_n_st_yt_py_tp
# umbel.fate.year.terms = with(
#   data = data.frame(t(fixef(u_n_st_yt_py_tp))),
#   data.frame(
#     intercept = c(0, Year2022, Year2023),
#     phen.slope = c(0, Year2022.centered.phen, Year2023.centered.phen),
#     drought.effect = trtdrought + c(0, trtdrought.Year2022, trtdrought.Year2023),
#     irrigat.effect = trtirrigated + c(0, trtirrigated.Year2022, trtirrigated.Year2023)
#   )
# )
# 
# # Plot to see temporal variation in effects:
# umbel.fate.year.terms %>%
#   mutate(Year = 2021:2023) %>%
#   pivot_longer(-Year, names_to = 'term', values_to = 'estimate') %>%
#   ggplot(aes(x = Year, y = estimate, group = term, colour = term)) +
#   geom_line() +
#   geom_point(size = 3)
# # wow - highly variable effects of drought...
# 
# # For seed set, mod is s_st_ty_p
# # code below is probably handling intercept effects the wrong way - above
# # looks closer to correct (but also maybe not necessary?)
# seed.set.year.terms = with(
#   data = data.frame(t(fixef(s_st_ty_p))),
#   data.frame(
#     intercept = c(0, Year2022, Year2023),
#     drought.slope = c(trtdrought, trtdrought.Year2022, trtdrought.Year2023),
#     irrigat.slope = c(trtirrigated, trtirrigated.Year2022, trtirrigated.Year2023)
#   )
# )

# Next, get residual mean and s.d. in size of new recruits
recruit.mu = fixef(r_ranef)[[1]]
recruit.sigma = summary(r_ranef)$sigma

# reprod.kernel = expand.grid(
#   trt = c('control', 'drought', 'irrigated'),
#   cur.size = (5:60)/10,
#   nex.size = (5:60)/10
# ) %>%
#   # Other covariates needed
#   mutate(
#     centered.phen = 0,
#     Year = factor(2021)
#   ) %>%
#   # Model predictions
#   mutate(
#     # Probability of flowering
#     p.flower = predict(
#       f_st_y,
#       newdata = ., type = 'response',
#       re.form = ~ 0, allow.new.levels = TRUE
#     ),
#     # Number of umbels
#     n.phen.umbels = predict(
#       n_s_y,
#       newdata = ., type = 'response',
#       re.form = ~ 0, allow.new.levels = TRUE
#     )
#   ) %>%
#   mutate(
#     # Proportion of umbels surviving
#     p.survive = predict(
#       u_n_s_y_py,
#       newdata = ., type = 'link',
#       re.form = ~ 0, allow.new.levels = TRUE
#     ),
#     p.survive = p.survive + mean(umbel.fate.year.terms$intercept) +
#       mean(umbel.fate.year.terms$slope) * cur.size,
#     p.survive = 1 / (1 + exp(-p.survive))
#   ) %>%
#   rename(n.umbel = n.phen.umbels) %>%
#   # Number of seeds per umbel
#   # first need to convert size
#   mutate(cur.size = cur.size - mean(seed.phen.demo$cur.size)) %>%
#   mutate(
#     n.seeds = predict(
#       s_st_y_p,
#       newdata = ., type = 'link',
#       re.form = ~ 0, allow.new.levels = TRUE
#     ),
#     n.seeds = n.seeds + mean(seed.set.year.terms$intercept),
#     n.seeds = exp(n.seeds),
#     # fix size
#     cur.size = cur.size + mean(seed.phen.demo$cur.size),
#   )

# Here: a way to get the mean effects across years
# For the two models with year-fixed effects, we get the mean effect by
# generating predictions for each year on the linear scale, then average them
# (linear averaging)
# We'll maybe want something else to get the variances...

reprod.kernel = expand.grid(
  # Get combinations of z', z for each treatment and year
  trt = c('control', 'drought', 'irrigated'),
  cur.size = (5:60)/10,
  nex.size = (5:60)/10,
  Year = factor(2021:2023),
  # p.germ = c(.001, .005, .01, .05, .1, .5, 1)
  p.germ = c(.001, .01, .1, 1)
) %>%
  # Other covariates needed for estimating from models
  mutate(
    centered.phen = 0,
  ) %>%
  # Model predictions
  # We'll estimate probability of flowering and umbel phenology on the response
  # scale
  # (although for getting inter-annual estimates we may want the linear scale)
  mutate(
    # Probability of flowering
    p.flower = predict(
      f_st_ty,
      newdata = ., type = 'response',
      re.form = ~ 0, allow.new.levels = TRUE
    ),
    # Number of umbels
    # (using the name n.phen.umbels because that's the predictor used in the the
    # umbel survival model)
    n.phen.umbels = predict(
      n_s_ty,
      newdata = ., type = 'response',
      re.form = ~ 0, allow.new.levels = TRUE
    )
  ) %>%
  # Get linear-scale probability of survival for each year 
  mutate(
    lin.p.survive = predict(
      u_n_st_yt_py,
      newdata = ., type = 'link',
      re.form = ~ 0, allow.new.levels = TRUE
    )
  ) %>%
  # Number of seeds per umbel
  # first need to convert size
  mutate(cur.size = cur.size - mean(seed.phen.demo$cur.size)) %>%
  mutate(
    lin.n.seeds = predict(
      s_st_ty_p,
      newdata = ., type = 'link',
      re.form = ~ 0, allow.new.levels = TRUE
    ),
    # Un-center the size column
    cur.size = cur.size + mean(seed.phen.demo$cur.size),
  ) %>%
  # Rename the numner of umbels column
  rename(n.umbels = n.phen.umbels) %>%
  # Get means of each parameter
  # NOTE: because year is a random effect in the other models, but we generate
  # predictions without year random effects specialized, the estimates are the
  # same across years, so averaging will just give us the same mean estimate
  select(-Year) %>%
  group_by(trt, cur.size, nex.size, centered.phen, p.germ) %>%
  summarise(across(everything(), mean)) %>%
  ungroup() %>%
  # Convert to the prediction scales
  mutate(
    p.survive = 1 / (1 + exp(-lin.p.survive)),
    n.seeds = exp(lin.n.seeds)
  ) %>%
  # Remove columns for linear scale
  select(-c(lin.p.survive, lin.n.seeds)) %>%
  # Get number of seeds produced
  mutate(total.seeds = p.flower * n.umbels * p.survive * n.seeds) %>%
  # Estimate number of germinating seeds (assume this is independent of phen)
  mutate(establishing.seeds = p.germ * total.seeds) %>%
  # Estimate size distribution of seeds
  mutate(
    p.nex.size = 0.1 * dnorm(x = nex.size, mean = recruit.mu, sd = recruit.sigma),
    n.nex.size = p.nex.size * establishing.seeds
  )

head(reprod.kernel)
nrow(reprod.kernel)

reprod.kernel %>%
  filter(p.germ %in% 1) %>%
  pivot_longer(
    cols = c(p.flower, n.umbels, p.survive, n.seeds, total.seeds),
    names_to = 'estim.type',
    values_to = 'estimate'
  ) %>%
  ggplot(aes(x = cur.size, y = estimate)) +
  geom_line(aes(colour = trt)) +
  scale_colour_manual(values = c('black', 'red', 'blue')) +
  # scale_y_log10() +
  facet_wrap(~ estim.type, scale = 'free_y')
# Cool

# Visualization of kernels here... not helpful
reprod.kernel %>%
  ggplot(aes(x = cur.size, y = nex.size)) +
  geom_tile(aes(fill = n.nex.size)) +
  scale_y_reverse() +
  facet_wrap(trt ~ p.germ)

# cool

nrow(reprod.kernel)

# write.csv(
#   reprod.kernel,
#   '03_construct_kernels/out/test_reprod_kernel.csv',
#   row.names = FALSE
# )
