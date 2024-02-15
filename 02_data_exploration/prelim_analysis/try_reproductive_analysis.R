# Script for testing out analyses of seeds per (reproductive) plant.
# Reading in demo data (with imputed survival) and combined phenology, demo, and
# seed set data (see script 01_data_cleaning_reconcile_demo_phen_seed_umbels.R)

##### Load products

# Read packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(glmmTMB)
library(cowplot)
# library(lme4)

# Clear namespace
rm(list = ls())

# Read in raw datasets

# Demo (imputed survival, all years)
demo = read.csv('01_data_cleaning/out/demo_imputed_survival.csv')

# Phenology + demography (all matched plants, 2021 and later)
phen.demo = read.csv('01_data_cleaning/out/phen_demo_for_umbel_survival.csv')

# Phenology + demography + seed set
seed.phen.demo.all = read.csv('01_data_cleaning/out/seed_phen_demo_combined.csv')

# Budding phenology dataset
phen = read.csv('01_data_cleaning/out/phenology_buds_deaths_all.csv')

##### Process data as needed

### Demography:

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
  )

# Assess data

# How many plants are missing umbel counts?
demo.for.flowering %>% filter(is.na(No.umbels))
# ah - these must be due to imputation

# How many plants are missing a *previous* umbel count
demo.for.flowering %>% filter(is.na(No.umbels.prev))
# over 400... missing plants, etc.

# How many plants have a previous size (or do not)
demo.for.flowering %>%
  group_by(has.prev.meas = !is.na(No.leaves.prev) & !is.na(Leaf.length.prev), Year) %>%
  summarise(n = n()) %>%
  pivot_wider(names_from = has.prev.meas, values_from = n)
# I can live with this - still plenty of plants in here...

# How many plants do/do not have current measurements?
demo.for.flowering %>%
  group_by(has.cur.meas = !is.na(No.leaves) & !is.na(Leaf.length), Year) %>%
  summarise(n = n()) %>%
  pivot_wider(names_from = has.cur.meas, values_from = n)

# Make a data frame with size
demo.for.flowering.size.both = demo.for.flowering %>%
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

table(demo.for.flowering.size.both$prev.flower)  
table(demo.for.flowering.size.both$flowering)

# Based on analysis below, probability of flowering best predicted by current
# size
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

head(demo.for.flowering)
nrow(demo.for.flowering) # ~2800 observations
length(unique(demo.for.flowering$plantid)) # almost 850 plants
table(demo.for.flowering$Year) # six years  

# Get a dataset with only flowering plants
# For modeling umbel counts
demo.for.umbel.counts = demo.for.flowering %>%
  filter(flowering) %>%
  select(-flowering)

### Demo + phen

phen.demo.with.size = phen.demo %>%
  # Get rid of plants that have no size measurement (including plants marked as
  # having zero leaves in prior year)
  filter(!is.na(prev.leaves) & !is.na(prev.length)) %>%
  # Change umbel counts to a non-zero quantity
  mutate(n.demo.umbels = ifelse(is.na(n.demo.umbels), 0, n.demo.umbels)) %>%
  # Add a size column
  mutate(prev.size = log(prev.leaves * prev.length)) %>%
  # Center phenology
  mutate(mean.doy = init.doy - round(mean(init.doy)))

head(phen.demo.with.size)

# While I'm here, add a previously-flowering column just in case
phen.demo = phen.demo %>%
  mutate(
    prev.umbels = ifelse(is.na(prev.umbels), 0, prev.umbels),
    prev.flower = prev.umbels > 0
  ) %>%
  mutate(mean.doy = init.doy - round(mean(init.doy)))

### Mean bud date against number of dead umbels
umbel.failure = phen.demo %>%
  group_by(Year, tagplot) %>%
  mutate(mean.doy = mean(init.doy)) %>%
  ungroup() %>%
  distinct(Year, tagplot, .keep_all = TRUE) %>%
  mutate(mean.phen = mean.doy - round(mean(mean.doy)))

umbel.failure.with.size = umbel.failure %>%
  filter(!is.na(prev.length) & !is.na(prev.leaves)) %>%
  mutate(prev.size = log(prev.length * prev.leaves))

# Worth checking at some point if this is truly binomial

### Seed + phen + demo

# Get a restricted dataset
seed.phen.demo = seed.phen.demo.all %>%
  # Remove plants that are missing phen
  filter(!is.na(mean.bud.day)) %>%
  # Remove plants that are missing from seed
  filter(!is.na(no.seeds)) %>%
  # Add centered phen column
  mutate(centered.phen = mean.bud.day - round(mean(mean.bud.day)))

# Removing the lost umbels (which are accounted for in a different model)
seed.phen.demo = rbind(
  # All umbels with a non-zero number of seeds...
  seed.phen.demo %>% filter(no.seeds > 0),
  # and all umbels with zero seeds, but removing 
  seed.phen.demo %>%
    filter(!no.seeds) %>%
    group_by(tagplot, Year) %>%
    # mutate(n = 1:n()) %>%
    filter((1:n()) > n.lost.umbel) %>%
    ungroup()
)

nrow(seed.phen.demo)

seed.phen.demo.size = seed.phen.demo %>%
  filter(!is.na(prev.length) & !is.na(prev.leaves)) %>%
  mutate(prev.size = log(prev.length * prev.leaves))

nrow(seed.phen.demo.size)
seed.phen.demo.size %>% distinct(Year, tagplot) %>% nrow()


### Phen dataset
# - we're reading in the original phen dataset before doing manipulations to correct for IDs
# but otherwise I think the data is the same
# so I'll just merge in IDs from the phen.demo dataset

phen.tagplot = merge(
  x = phen, y = phen.demo %>% distinct(tagplot, plantid.phen),
  by.x = 'plantid', by.y = 'plantid.phen'
) %>%
  # Filter out only new umbels (i.e., buds)
  filter(varb %in% 'new') %>%
  select(-varb)

head(phen.tagplot)
nrow(phen.tagplot)

##### Fit models

### Demography (probability of flowering)

pf_s0 = glmmTMB(
  flowering ~ prev.size + (1 | Year) + (1 | Plot / plantid),
  data = demo.for.flowering.size.both,
  family = 'binomial',
)

pf_s1 = glmmTMB(
  flowering ~ cur.size + (1 | Year) + (1 | Plot / plantid),
  data = demo.for.flowering.size.both,
  family = 'binomial',
)

AIC(pf_s0, pf_s1) %>% mutate(daic = round(AIC - min(AIC), 2))
# Current size is a *much* better predictor than prior size lmao

# (refit on larger dataset)
pf_s = glmmTMB(
  flowering ~ cur.size + (1 | Year) + (1 | Plot / plantid),
  data = demo.for.flowering,
  family = 'binomial',
)

# pf_s_f = glmmTMB(
#   flowering ~ cur.size + prev.flower + (1 | Year) + (1 | Plot / plantid),
#   data = demo.for.flowering,
#   family = 'binomial'
# )
# 
# summary(pf_s_f)
# # plants that flowered last year much more likely to flower - odds increase by 260% (!)
# 
# anova(pf_s_f, pf_s)
# # ah... but this moakes models so much more complicated...
# 
# pf_sy_f = glmmTMB(
#   flowering ~ prev.flower + (cur.size | Year) + (1 | Plot / plantid),
#   data = demo.for.flowering,
#   family = 'binomial'
# )
# 
# summary(pf_sy_f)
# 
# anova(pf_sy_f, pf_s_f) # definitely no evidence of a strong size-year interaction
# 
# pf_s_f_t = glmmTMB(
#   flowering ~ cur.size + prev.flower + trt + (1 | Year) + (1 | Plot / plantid),
#   data = demo.for.flowering,
#   family = 'binomial'
# )
# 
# summary(pf_s_f_t) # doesn't seem like evidence of a treatment effect here
# 
# anova(pf_s_f_t, pf_s_f) # nope - no treatment effect
# # but maybe there is one that depends on size
# 
# pf_st_f = glmmTMB(
#   flowering ~ cur.size * trt + prev.flower + (1 | Year) + (1 | Plot / plantid),
#   data = demo.for.flowering,
#   family = 'binomial'
# )
# 
# anova(pf_st_f, pf_s_f) # interaction not supported
# 
# pf_st_ft = glmmTMB(
#   flowering ~ cur.size * trt + prev.flower * trt + (1 | Year) + (1 | Plot / plantid),
#   data = demo.for.flowering,
#   family = 'binomial'
# )
# 
# anova(pf_st_ft, pf_s_f) 
# # treatment effects are supported, p < 0.02, but delta AIC only barely larger than two...
# 
# pf_s_ft = glmmTMB(
#   flowering ~ cur.size + prev.flower * trt + (1 | Year) + (1 | Plot / plantid),
#   data = demo.for.flowering,
#   family = 'binomial'
# )

# AIC(pf_s_f, pf_st_f, pf_s_ft, pf_st_ft) %>% mutate(daic = round(AIC - min(AIC), 2))

pf_s_t = glmmTMB(
  flowering ~ cur.size + trt + (1 | Year) + (1 | Plot / plantid),
  data = demo.for.flowering,
  family = 'binomial',
)

anova(pf_s_t, pf_s) # no treatment effect

pf_st = glmmTMB(
  flowering ~ cur.size * trt + (1 | Year) + (1 | Plot / plantid),
  data = demo.for.flowering,
  family = 'binomial',
)

anova(pf_st, pf_s) # hmm... it's weak but it is better-performing

AIC(pf_st, pf_s) %>% mutate(daic = round(AIC - min(AIC), 2)) # hmm...

# Final model of these
summary(pf_st)
# hmm.. without mean-centering it's not super easy to discern what's going on
# but... maybe large plants in drought are less likely to flower

### Umbel count models

nu_s = glmmTMB(
  No.umbels ~ cur.size + (1 | Year) + (1 | Plot / plantid),
  data = demo.for.umbel.counts,
  family = 'truncated_poisson'
)

summary(nu_s)

# nu_s_f = glmmTMB(
#   No.umbels ~ prev.size + prev.flower + (1 | Year) + (1 | Plot / plantid),
#   data = demo.for.umbel.counts,
#   family = 'truncated_poisson'
# )
# 
# AIC(nu_s_f, nu_s)
# 
# summary(nu_s_f)
# # Flowering in previous year means more flowers in subsequent year (if flowering)
# 
# # Size and/or previous flowering treatments varying by year
# 
# nu_sy_f = glmmTMB(
#   No.umbels ~  prev.flower + (prev.size | Year) + (1 | Plot / plantid),
#   data = demo.for.umbel.counts,
#   family = 'truncated_poisson'
# )
# 
# nu_s_fy = glmmTMB(
#   No.umbels ~ prev.size + (prev.flower | Year) + (1 | Plot / plantid),
#   data = demo.for.umbel.counts,
#   family = 'truncated_poisson'
# )
# # Convergence issue (this always happens with the categoricals...)
# 
# AIC(nu_sy_f, nu_s_f) # size effects don't vary by year
# 
# # Treatment effects
# 
# nu_s_f_t = glmmTMB(
#   No.umbels ~ prev.size + prev.flower + trt + (1 | Year) + (1 | Plot / plantid),
#   data = demo.for.umbel.counts,
#   family = 'truncated_poisson'
# )
# 
# nu_st_f = glmmTMB(
#   No.umbels ~ prev.size * trt + prev.flower + (1 | Year) + (1 | Plot / plantid),
#   data = demo.for.umbel.counts,
#   family = 'truncated_poisson'
# )
# 
# nu_s_ft = glmmTMB(
#   No.umbels ~ prev.size + prev.flower * trt + (1 | Year) + (1 | Plot / plantid),
#   data = demo.for.umbel.counts,
#   family = 'truncated_poisson'
# )
# 
# nu_st_ft = glmmTMB(
#   No.umbels ~ prev.size * trt + prev.flower * trt + (1 | Year) + (1 | Plot / plantid),
#   data = demo.for.umbel.counts,
#   family = 'truncated_poisson'
# )
# 
# AIC(nu_st_ft, nu_st_f, nu_s_ft, nu_s_f_t, nu_s_f) %>%
#   mutate(daic = round(AIC - min(AIC), 2))
# # No treatment effects at all.

nu_sy = glmmTMB(
  No.umbels ~ (cur.size | Year) + (1 | Plot / plantid),
  data = demo.for.umbel.counts,
  family = 'truncated_poisson'
)

nu_s_t = glmmTMB(
  No.umbels ~ cur.size + trt + (1 | Year) + (1 | Plot / plantid),
  data = demo.for.umbel.counts,
  family = 'truncated_poisson'
)

nu_st = glmmTMB(
  No.umbels ~ cur.size * trt + (1 | Year) + (1 | Plot / plantid),
  data = demo.for.umbel.counts,
  family = 'truncated_poisson'
)

AIC(nu_s, nu_sy, nu_s_t, nu_st) %>% mutate(daic = round(AIC - min(AIC), 2))
# doesn't look like there's any evidence of a year effect here...

summary(nu_s)
# larger plants produce more umbels... seems unsurprising

# NOTE that I had convergence issues with the negative binomial based on this -
# not sure if overdispersion is an issue here

### Umbel budding phenology

phen_0 = glmmTMB(
  formula = mean.doy ~ (1 | Year) + (1 | Plot / tagplot),
  data = phen.demo.with.size,
)

phen_s = glmmTMB(
  formula = mean.doy ~ prev.size + (1 | Year) + (1 | Plot / tagplot),
  data = phen.demo.with.size,
)

anova(phen_s, phen_0)
# no evidence of a size relationship

# test with treatment effects just to be sure

phen_t = glmmTMB(
  formula = mean.doy ~ trt + (1 | Year) + (1 | Plot / tagplot),
  data = phen.demo.with.size,
)

phen_s_t = glmmTMB(
  formula = mean.doy ~ trt + prev.size + (1 | Year) + (1 | Plot / tagplot),
  data = phen.demo.with.size,
)

anova(phen_s_t, phen_t)
# nothing here either
# (I also tried an interaction)
# sweet!

# Refit these with full datasets

phen_0 = glmmTMB(
  formula = mean.doy ~ (1 | Year) + (1 | Plot / tagplot),
  data = phen.demo
)

phen_t = glmmTMB(
  formula = mean.doy ~ trt + (1 | Year) + (1 | Plot / tagplot),
  data = phen.demo,
)

anova(phen_t, phen_0)
# phen effect is here

summary(phen_t)
# as before - drought means earlier flowering on average, irrigation means later

# Just to be sure... is there an effect of previous flowering?
phen_t_f = glmmTMB(
  formula = mean.doy ~ trt + prev.flower + (1 | Year) + (1 | Plot / tagplot),
  data = phen.demo,
)

anova(phen_t_f, phen_t)
AIC(phen_t_f, phen_t)

summary(phen_t_f)
# Interesting. If you flowered in the last year, you're more likely to bud earlier this year?
# (and the mean size of this effect is 1.5 days)
# Why would this be? We could write a story. Is it worth anything?
# throwing this term in also makes a stronger effect of drought...


### Umbel failure
# (Maybe it's smarter to parameterize this as surviving umbels?)

# First try with size

uf_0 = glmmTMB(
  cbind(n.phen.umbels - n.lost.umbels, n.lost.umbels) ~ 1 +
    (1 | Year) + (1 | Plot / tagplot),
  family = 'binomial',
  data = umbel.failure.with.size
)

uf_s = glmmTMB(
  cbind(n.phen.umbels - n.lost.umbels, n.lost.umbels) ~ prev.size +
    (1 | Year) + (1 | Plot / tagplot),
  family = 'binomial',
  data = umbel.failure.with.size
)

anova(uf_s, uf_0)
# okie dokie - no size effects!
# also fit with treatment juuuust to make sure

uf_t = glmmTMB(
  cbind(n.phen.umbels - n.lost.umbels, n.lost.umbels) ~ trt +
    (1 | Year) + (1 | Plot / tagplot),
  family = 'binomial',
  data = umbel.failure.with.size
)

uf_s_t = glmmTMB(
  cbind(n.phen.umbels - n.lost.umbels, n.lost.umbels) ~ prev.size + trt +
    (1 | Year) + (1 | Plot / tagplot),
  family = 'binomial',
  data = umbel.failure.with.size
)

anova(uf_s_t, uf_t)
# Great.

# Fitting without size:

uf_0 = glmmTMB(
  cbind(n.phen.umbels - n.lost.umbels, n.lost.umbels) ~ 1 +
    (1 | Year) + (1 | Plot / tagplot),
  family = 'binomial',
  data = umbel.failure
)

uf_t = glmmTMB(
  cbind(n.phen.umbels - n.lost.umbels, n.lost.umbels) ~ trt +
    (1 | Year) + (1 | Plot / tagplot),
  family = 'binomial',
  data = umbel.failure
)

uf_p = glmmTMB(
  cbind(n.phen.umbels - n.lost.umbels, n.lost.umbels) ~ mean.phen +
    (1 | Year) + (1 | Plot / tagplot),
  family = 'binomial',
  data = umbel.failure
)

uf_p_t = glmmTMB(
  cbind(n.phen.umbels - n.lost.umbels, n.lost.umbels) ~ mean.phen + trt +
    (1 | Year) + (1 | Plot / tagplot),
  family = 'binomial',
  data = umbel.failure
)

uf_p2 = glmmTMB(
  cbind(n.phen.umbels - n.lost.umbels, n.lost.umbels) ~ poly(mean.phen, 2) +
    (1 | Year) + (1 | Plot / tagplot),
  family = 'binomial',
  data = umbel.failure
)

uf_p2_t = glmmTMB(
  cbind(n.phen.umbels - n.lost.umbels, n.lost.umbels) ~ poly(mean.phen, 2) + trt +
    (1 | Year) + (1 | Plot / tagplot),
  family = 'binomial',
  data = umbel.failure
)


AIC(uf_0, uf_p, uf_t, uf_p_t, uf_p2, uf_p2_t) %>% 
  mutate(daic = round(AIC - min(AIC), 2)) %>%
  arrange(daic)
# big ol' phen effect... quadratic or linear
# no treatment effect lmao

summary(uf_p)
# Budding later means slightly higher chance of umbel success

summary(uf_p2)

# Check for a polynomial effect
uf_p2 = glmmTMB(
  cbind(n.phen.umbels - n.lost.umbels, n.lost.umbels) ~ poly(mean.phen, 2) +
    (1 | Year) + (1 | Plot / tagplot),
  family = 'binomial',
  data = umbel.failure
)

anova(uf_p2, uf_p)
# phew! no polynomial effects.

### Seeds per umbel

# First, look for effects of previous size
s_0_size = glmmTMB(
  no.seeds ~ (1 | Plot / tagplot) + (1 | Year),
  data = seed.phen.demo %>%
    filter(!is.na(prev.leaves) & !is.na(prev.length)) %>%
    mutate(prev.size = log(prev.leaves * prev.length)),
  family = 'nbinom2'
)

s_s_size = glmmTMB(
  no.seeds ~ prev.size + (1 | Plot / tagplot) + (1 | Year),
  data = seed.phen.demo %>%
    filter(!is.na(prev.leaves) & !is.na(prev.length)) %>%
    mutate(prev.size = log(prev.leaves * prev.length)),
  family = 'nbinom2'
)

anova(s_0_size, s_s_size)
# god motherfucking damnit
# that's insane.

# Okay.
# Run everything with previous size now I guess.
# (ah... previous size is maybe picking up umbel size)

s_s = glmmTMB(
  no.seeds ~ prev.size + (1 | Plot / tagplot) + (1 | Year),
  data = seed.phen.demo.size,
  family = 'nbinom2'
)

summary(s_s)

s_s_t = glmmTMB(
  no.seeds ~ prev.size + trt + (1 | Plot / tagplot) + (1 | Year),
  data = seed.phen.demo.size,
  family = 'nbinom2'
)

s_st = glmmTMB(
  no.seeds ~ prev.size * trt + (1 | Plot / tagplot) + (1 | Year),
  data = seed.phen.demo.size,
  family = 'nbinom2'
)

AIC(s_st, s_s_t, s_s) %>% mutate(daic = round(AIC - min(AIC), 2))
# no direct treatment effects

# Test for phenology

s_s_p = glmmTMB(
  no.seeds ~ centered.phen + prev.size + (1 | Plot / tagplot) + (1 | Year),
  data = seed.phen.demo.size,
  family = 'nbinom2'
) 

AIC(s_s, s_s_p)
anova(s_s_p, s_s)
summary(s_s_p)
# negative slope - flowering later (higher day) means fewer seeds

s_s_p2 = glmmTMB(
  no.seeds ~ poly(centered.phen, 2) + prev.size + (1 | Plot / tagplot) + (1 | Year),
  data = seed.phen.demo.size,
  family = 'nbinom2'
) 

AIC(s_s, s_s_p, s_s_p2)
summary(s_s_p2)
# definitely not quadratic

# Size-phen interaction...?

s_sp = glmmTMB(
  no.seeds ~ centered.phen * prev.size + (1 | Plot / tagplot) + (1 | Year),
  data = seed.phen.demo.size,
  family = 'nbinom2'
) 

AIC(s_sp, s_s_p) # sweet - no interaction

# Maybe look for treatment-by-phen interactions...

s_s_pt = glmmTMB(
  no.seeds ~ centered.phen * trt + prev.size + (1 | Plot / tagplot) + (1 | Year),
  data = seed.phen.demo.size,
  family = 'nbinom2'
) 

AIC(s_s_pt, s_s_p) # no treatment effects # (also true if we just have fixed effect of phen)

s_s_py = glmmTMB(
  no.seeds ~ prev.size + (1 | Plot / tagplot) + (centered.phen | Year),
  data = seed.phen.demo.size,
  family = 'nbinom2'
) 
# NOTE: I tried size effects varying by year and there was a convergence issue

AIC(s_s_py, s_s_p) # phen effect should not be varying by year

# Best model is s_s_p

summary(s_s_p)
# - flowering on the ~average day is ~33 seeds/umbel
# - one day later flowering causes a ~2% reduction in seed set...
#   - flowering one week later causes a ~10% reduction in seed set
# - effects of body size but these effects are difficult to interpret due to scale
#   (one-unit change is ~14% change to seed set... but 80% of plants are within 2 units...s)

#########################################################
# Figure drafts
#########################################################

#################
# Probability of flowering

pr.flower.data.base = demo.for.flowering %>%
  mutate(flowering.01 = as.numeric(flowering) + ifelse(prev.flower, .05, -.05)) %>%
  ggplot(aes(x = prev.size, y = flowering.01, colour = trt, shape = prev.flower)) +
  geom_point(position = position_jitter(height = 0.01)) +
  scale_colour_manual(values = c('black', 'red', 'blue'), 'treatment') +
  scale_shape_manual(
    values = c(4, 1),
    labels = c('prev. vegetative', 'prev. flowering'),
    ''
  ) +
  facet_wrap(~ Year) +
  labs(x = 'Previous size', y = 'Flowering') +
  theme(
    panel.background = element_blank(),
    legend.position = 'top'
  )

# Model is pf_s_ft

pr.flower.preds = expand.grid(
  prev.size = (10:50)/10,
  prev.flower = c(TRUE, FALSE),
  trt = c('control', 'drought', 'irrigated')
) %>%
  mutate(
    pred.prob = predict(
      object = pf_s_ft,  newdata = ., 
      re.form = ~ 0, allow.new.levels = TRUE, 
      type = 'response'
    )
  )

pr.flower.data.base +
  geom_line(
    data = pr.flower.preds,
    inherit.aes = FALSE,
    aes(
      x = prev.size, y = pred.prob, 
      group = interaction(trt, prev.flower),
      colour = trt, linetype = prev.flower
    )
  ) +
  scale_linetype_manual(
    values = 2:1, 
    labels = c('prev. vegetative', 'prev. flowering'),
    ''
  )


#################
# Number of umbels / plant

# model is nu_s_f

u.count.data.base = demo.for.umbel.counts %>%
  ggplot(aes(x = prev.size, y = No.umbels, colour = trt, shape = prev.flower)) +
  geom_point(position = position_jitter(height = 0.05)) +
  scale_colour_manual(values = c('black', 'red', 'blue'), 'treatment') +
  scale_shape_manual(
    values = c(4, 1),
    labels = c('prev. vegetative', 'prev. flowering'),
    ''
  ) +
  labs(x = 'Previous size', y = 'Number of umbels') +
  facet_wrap(~ Year) +
  theme(
    panel.background = element_blank(),
    legend.position = 'top'
  )

u.count.data.base

u.count.preds = expand.grid(
  prev.size = (10:50)/10,
  prev.flower = c(TRUE, FALSE)
) %>%
  mutate(
    pred.umb = predict(
      object = nu_s_f,  newdata = ., 
      re.form = ~ 0, allow.new.levels = TRUE, 
      type = 'response'
    )
  )

u.count.data.base +
  geom_line(
    data = u.count.preds,
    inherit.aes = FALSE,
    aes(x = prev.size, y = pred.umb, linetype = prev.flower)
  ) +
  scale_linetype_manual(
    values = 2:1, 
    labels = c('prev. vegetative', 'prev. flowering'),
    ''
  ) +
  scale_y_log10()

#################
# Probability of umbels failing

# optimal model is uf_p - only has phenology effect

u.succ.data.base = umbel.failure %>%
  mutate(p.succ.umbel = (n.phen.umbels - n.lost.umbels) / n.phen.umbels) %>%
  ggplot(aes(x = mean.phen, y = p.succ.umbel, colour = trt)) +
  geom_point(aes(size = n.phen.umbels)) +
  scale_colour_manual(values = c('black', 'red', 'blue'), 'treatment') +
  scale_size_continuous('Number of umbels') +
  scale_x_continuous(breaks = 10*(-3:3), labels = c('7 apr', '14 apr', '27 apr', '7 may', '17 may', '27 may', '6 jun')) +
  labs(x = 'Day of budding', y = 'Probability of umbel success') +
  facet_wrap(~ Year) +
  theme(
    panel.background = element_blank(),
    legend.position = 'top'
  )

u.succ.data.base

u.succ.pred = data.frame(
  mean.phen = -25:30
) %>%
  mutate(
    pred.prob = predict(
      object = uf_p,
      newdata = .,
      type = 'response',
      re.form = ~ 0,
      allow.new.levels = TRUE
    )
  )

u.succ.data.base +
  geom_line(
    data = u.succ.pred,
    inherit.aes = FALSE,
    aes(x = mean.phen, y = pred.prob)
  )
# Not amazingly compelling!

#################
# Seeds / reproductive umbel

seed.set.base = seed.phen.demo.size %>%
  ggplot(aes(x = centered.phen, y = no.seeds, colour = trt)) +
  geom_point(size = 3) +
  scale_colour_manual(values = c('black', 'red', 'blue', 'treatment')) +
  labs(x = 'Phenology (centered, in days)', y = 'Number of seeds') +
  facet_wrap(~ Year, ncol = 1) +
  theme(
    panel.background = element_blank(),
    legend.position = 'top'
  )

seed.set.base

seed.set.preds = expand.grid(
  prev.size = c(2.7, 3.8, 4.6),
  centered.phen = -20:30
) %>%
  mutate(
    pred.seed = predict(
      object = s_s_p,
      newdata = .,
      re.form = ~ 0,
      type = 'response',
      allow.new.levels = TRUE
    )
  )

seed.set.base +
  geom_line(
    data = seed.set.preds %>% mutate(prev.size = factor(prev.size)),
    inherit.aes = FALSE,
    aes(
      x = centered.phen, y = pred.seed, 
      group = prev.size,
      linetype = prev.size)
  ) +
  scale_linetype_manual(values = c(2, 1, 2))

#################
# Raw phenology plot
# rows: years + treatments
# y axis: day of year
# (different data input...)

phen.plot.data.base = phen.demo %>%
  mutate(year.dodge = Year + -0.2 * as.numeric(trt %in% 'drought') + 0.2 * as.numeric(trt %in% 'irrigated')) %>%
  ggplot(aes(x = mean.doy, y = year.dodge, colour = trt)) +
  geom_point(position = position_jitter(height = 0.08, width = 2), size = 2) +
  scale_colour_manual(values = c('black', 'red', 'blue')) +
  # scale_y_reverse(breaks = 2020:2023, labels = c('overall', '2021', '2022', '2023')) +
  scale_x_continuous(breaks = 10*(-3:3), labels = c('7 apr', '14 apr', '27 apr', '7 may', '17 may', '27 may', '6 jun')) +
  # date column is centered in data frame such that 0 = may 7 (day 126)
  # use as.Date(10*(-2:3) + 126) to get these
  scale_y_continuous(breaks = 2020:2023, labels = c('overall', '2021', '2022', '2023')) +
  labs(x = 'Day of year (centered)', y = '') +
  theme(
    legend.position = 'none',
    panel.background = element_blank()
  )

phen.plot.data.base

phen.predict = predict(
  object = phen_t,
  newdata = data.frame(trt = c('control', 'drought', 'irrigated')),
  type = 'response',
  se.fit = TRUE,
  re.form = ~ 0,
  allow.new.levels = TRUE
) %>%
  do.call(what = cbind) %>%
  as.data.frame() %>%
  mutate(trt = c('control', 'drought', 'irrigated')) %>%
  mutate(fit.lb = fit - se.fit, fit.ub = fit + se.fit) %>%
  mutate(year.dodge = 2020 + -0.2 * as.numeric(trt %in% 'drought') + 0.2 * as.numeric(trt %in% 'irrigated'))

phen.plot.data.base +
  geom_segment(
    data = phen.predict,
    inherit.aes = FALSE,
    aes(x = fit.lb, xend = fit.ub, y = year.dodge, yend = year.dodge)
  ) +
  geom_point(
    data = phen.predict,
    inherit.aes = FALSE,
    aes(x = fit, y = year.dodge, fill = trt),
    size = 5, shape = 21
  ) +
  scale_fill_manual(values = c('black', 'red', 'blue')) 

#########################################################
# A single figure
#########################################################

# (a) Budding phenology

pan.a = phen.demo %>%
  # Manually dodge y-axis by treatment
  mutate(year.dodge = Year + -0.2 * as.numeric(trt %in% 'drought') + 0.2 * as.numeric(trt %in% 'irrigated')) %>%
  ggplot(aes(x = mean.doy, y = year.dodge, colour = trt)) +
  # Plot raw data
  geom_point(position = position_jitter(height = 0.04, width = 2), size = 1, alpha = 0.25) +
  # Plot model predicted mean bud date
  geom_point(
    data = phen.predict,
    inherit.aes = FALSE,
    aes(x = fit, y = year.dodge, colour = trt),
    size = 5
  ) +
  # Plot model prediction standard errors
  geom_segment(
    data = phen.predict,
    inherit.aes = FALSE,
    aes(x = fit.lb, xend = fit.ub, y = year.dodge, yend = year.dodge)
  ) +
  # Colour/fill aesthetics
  scale_fill_manual(values = c('black', 'red', 'blue'), 'treatment') +
  scale_colour_manual(values = c('black', 'red', 'blue'), 'treatment') +
  scale_x_continuous(breaks = 10*(-3:3), labels = c('7 apr', '14 apr', '27 apr', '7 may', '17 may', '27 may', '6 jun')) +
  # date column is centered in data frame such that 0 = may 7 (day 126)
  # use as.Date(10*(-2:3) + 126) to get these
  scale_y_continuous(breaks = 2020:2023, labels = c('overall', '2021', '2022', '2023')) +
  labs(x = '', y = '') +
  theme(
    legend.position = 'none',
    axis.text.x = element_text(angle = 90),
    panel.background = element_blank()
  )

pan.a

# (b) Umbel success

pan.b = umbel.failure %>%
  mutate(p.succ.umbel = (n.phen.umbels - n.lost.umbels) / n.phen.umbels) %>%
  ggplot(aes(x = mean.phen, y = p.succ.umbel, colour = trt)) +
  geom_point(size = 1, alpha = 0.5) +
  geom_line(
    data = u.succ.pred,
    inherit.aes = FALSE,
    aes(x = mean.phen, y = pred.prob)
  ) +
  scale_colour_manual(values = c('black', 'red', 'blue'), 'treatment') +
  scale_x_continuous(breaks = 10*(-3:3), labels = c('7 apr', '14 apr', '27 apr', '7 may', '17 may', '27 may', '6 jun')) +
  labs(x = '', y = 'Proportion of umbels surviving') +
  theme(
    panel.background = element_blank(),
    axis.text.x = element_text(angle = 90),
    legend.position = 'none'
  )

pan.b

# (c) Seeds per umbel

pan.c = seed.phen.demo.size %>%
  ggplot(aes(x = centered.phen, y = no.seeds, colour = trt)) +
  geom_point(size = 1, alpha = 0.5) +
  geom_line(
    data = seed.set.preds %>% mutate(prev.size = factor(prev.size)),
    inherit.aes = FALSE,
    aes(
      x = centered.phen, y = pred.seed, 
      group = prev.size,
      linetype = prev.size)
  ) +
  scale_linetype_manual(values = c(4, 1, 2)) +
  scale_colour_manual(values = c('black', 'red', 'blue')) +
  scale_x_continuous(breaks = -1 + 10*(-3:3), labels = c('7 apr', '14 apr', '27 apr', '7 may', '17 may', '27 may', '6 jun')) +
  labs(x = '', y = 'Number of seeds in umbel') +
  theme(
    panel.background = element_blank(),
    axis.text.x = element_text(angle = 90),
    legend.position = 'none'
  )

pan.c

# Get single treatment legend


# Use cowplot to put them together

plot_grid(
  get_legend(pan.a + theme(legend.position = 'top')),
  plot_grid(pan.a, pan.b, pan.c, nrow = 1, labels = c('a', 'b', 'c')),
  ncol = 1,
  rel_heights = c(0.1, 1)
) %>%
  save_plot(
    filename = '02_data_exploration/prelim_analysis/fig_phen_repro.png',
    base_height = 5, base_width = 8
  )
