# Script for testing out analyses of seeds per (reproductive) plant.
# Reading in demo data (with imputed survival) and combined phenology, demo, and
# seed set data (see script 01_data_cleaning_reconcile_demo_phen_seed_umbels.R)

##### Load products

# Read packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(glmmTMB)
# library(lme4)

# Clear namespace
rm(list = ls())

# Read in raw datasets

# Demo (imputed survival, all years)
demo = read.csv('01_data_cleaning/out/demo_imputed_survival.csv')

# Phenology + demography (all matched plants, 2021 and later)
phen.demo = read.csv('01_data_cleaning/out/phen_demo_for_umbel_survival.csv')

# Phenology + demography + seed set
# (not finished yet)


##### Process data as needed

### Demography:

demo.for.flowering = merge(
    x = demo %>%
      # Filter out only surviving plants
      filter(surv) %>%
      # Filter out 2016 plants (size measurements are not reliable)
      filter(Year > 2016) %>%
      # Select columns for year-of data
      select(Year, plantid, Plot, trt, No.umbels),
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
  group_by(has.prev.meas = !is.na(No.leaves) & !is.na(Leaf.length), Year) %>%
  summarise(n = n()) %>%
  pivot_wider(names_from = has.prev.meas, values_from = n)
# I can live with this - still plenty of plants in here...

demo.for.flowering = demo.for.flowering %>%
  # Filter out only plants with a previous measurement
  filter(!is.na(No.leaves) & !is.na(Leaf.length) & No.leaves > 0) %>%
  # Add a size column
  mutate(prev.size = log(No.leaves * Leaf.length)) %>%
  select(-c(No.leaves, Leaf.length)) %>%
  # Replace NAs with zero
  mutate(
    across(contains('umbel'), function(x) ifelse(is.na(x), 0, x)),
    prev.flower = No.umbels.prev > 0,
    flowering = No.umbels > 0
  )

table(demo.for.flowering$prev.flower)  
table(demo.for.flowering$flowering)

head(demo.for.flowering)
nrow(demo.for.flowering) # >2700 observations
length(unique(demo.for.flowering$plantid)) # 840 plants
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
  mutate(mean.doy = init.doy - round(mean(init.doy), 2))

head(phen.demo.with.size)

# While I'm here, add a previously-flowering column just in case
phen.demo = phen.demo %>%
  mutate(
    prev.umbels = ifelse(is.na(prev.umbels), 0, prev.umbels),
    prev.flower = prev.umbels > 0
  ) %>%
  mutate(mean.doy = init.doy - round(mean(init.doy), 2))

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

##### Fit models

### Demography (probability of flowering)

pf_s = glmmTMB(
  flowering ~ prev.size + (1 | Year) + (1 | Plot / plantid),
  data = demo.for.flowering,
  family = 'binomial',
)

summary(pf_s)
# Larger plants much more likely to flower

pf_s_f = glmmTMB(
  flowering ~ prev.size + prev.flower + (1 | Year) + (1 | Plot / plantid),
  data = demo.for.flowering,
  family = 'binomial'
)

summary(pf_s_f)
# plants that flowered last year much more likely to flower - odds increase by 230%

anova(pf_s_f, pf_s)

pf_sy_f = glmmTMB(
  flowering ~ prev.flower + (prev.size | Year) + (1 | Plot / plantid),
  data = demo.for.flowering,
  family = 'binomial'
)

summary(pf_sy_f)

anova(pf_sy_f, pf_s_f) # definitely no evidence of a strong size-year interaction


pf_s_f_t = glmmTMB(
  flowering ~ prev.size + prev.flower + trt + (1 | Year) + (1 | Plot / plantid),
  data = demo.for.flowering,
  family = 'binomial'
)

summary(pf_s_f_t) # doesn't seem like evidence of a treatment effect here

anova(pf_s_f_t, pf_s_f) # nope - no treatment effect
# but maybe there is one that depends on size

pf_st_f = glmmTMB(
  flowering ~ prev.size * trt + prev.flower + (1 | Year) + (1 | Plot / plantid),
  data = demo.for.flowering,
  family = 'binomial'
)

anova(pf_st_f, pf_s_f) # interesting... weak support

pf_st_ft = glmmTMB(
  flowering ~ prev.size * trt + prev.flower * trt + (1 | Year) + (1 | Plot / plantid),
  data = demo.for.flowering,
  family = 'binomial'
)

anova(pf_st_ft, pf_s_f) # very well supported

pf_s_ft = glmmTMB(
  flowering ~ prev.size + prev.flower * trt + (1 | Year) + (1 | Plot / plantid),
  data = demo.for.flowering,
  family = 'binomial'
)

AIC(pf_s_f, pf_st_f, pf_s_ft, pf_st_ft) %>% mutate(daic = round(AIC - min(AIC), 2))

summary(pf_s_ft)
# Lol
# Evidence here suggests that:
# - Increasing size increases probability of flowering for all plants
# - Flowering in previous year increases probability of flowering in next year
#   - Effect is much stronger in irrigated plots than other plots
# - For plants that did not flower previously, irrigation may reduce the probability of flowering

### Umbel count models

nu_s = glmmTMB(
  No.umbels ~ prev.size + (1 | Year) + (1 | Plot / plantid),
  data = demo.for.umbel.counts,
  family = 'truncated_poisson'
)

summary(nu_s)

nu_s_f = glmmTMB(
  No.umbels ~ prev.size + prev.flower + (1 | Year) + (1 | Plot / plantid),
  data = demo.for.umbel.counts,
  family = 'truncated_poisson'
)

AIC(nu_s_f, nu_s)

summary(nu_s_f)
# Flowering in previous year means more flowers in subsequent year (if flowering)

# Size and/or previous flowering treatments varying by year

nu_sy_f = glmmTMB(
  No.umbels ~  prev.flower + (prev.size | Year) + (1 | Plot / plantid),
  data = demo.for.umbel.counts,
  family = 'truncated_poisson'
)

nu_s_fy = glmmTMB(
  No.umbels ~ prev.size + (prev.flower | Year) + (1 | Plot / plantid),
  data = demo.for.umbel.counts,
  family = 'truncated_poisson'
)
# Convergence issue (this always happens with the categoricals...)

AIC(nu_sy_f, nu_s_f) # size effects don't vary by year

# Treatment effects

nu_s_f_t = glmmTMB(
  No.umbels ~ prev.size + prev.flower + trt + (1 | Year) + (1 | Plot / plantid),
  data = demo.for.umbel.counts,
  family = 'truncated_poisson'
)

nu_st_f = glmmTMB(
  No.umbels ~ prev.size * trt + prev.flower + (1 | Year) + (1 | Plot / plantid),
  data = demo.for.umbel.counts,
  family = 'truncated_poisson'
)

nu_s_ft = glmmTMB(
  No.umbels ~ prev.size + prev.flower * trt + (1 | Year) + (1 | Plot / plantid),
  data = demo.for.umbel.counts,
  family = 'truncated_poisson'
)

nu_st_ft = glmmTMB(
  No.umbels ~ prev.size * trt + prev.flower * trt + (1 | Year) + (1 | Plot / plantid),
  data = demo.for.umbel.counts,
  family = 'truncated_poisson'
)

AIC(nu_st_ft, nu_st_f, nu_s_ft, nu_s_f_t, nu_s_f) %>%
  mutate(daic = round(AIC - min(AIC), 2))
# No treatment effects at all.

# So final model based on this papears to be:
# - 

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
# (I also tried an interactio)
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
# treatment effect is here, although it appears to be small...

summary(phen_t)
# as before - drought means earlier flowering on average, irrigation means later

# Just to be sure... is there an effect of previous flowering?
phen_t_f = glmmTMB(
  formula = mean.doy ~ trt + prev.flower + (1 | Year) + (1 | Plot / tagplot),
  data = phen.demo,
)

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

AIC(uf_0, uf_p, uf_t, uf_p_t) %>% mutate(daic = round(AIC - min(AIC), 2))
# big ol' phen effect
# no treatment effect lmao

summary(uf_p)
# Budding later means slightly higher chance of umbel success

# Check for a polynomial effect
uf_p2 = glmmTMB(
  cbind(n.phen.umbels - n.lost.umbels, n.lost.umbels) ~ poly(mean.phen, 2) +
    (1 | Year) + (1 | Plot / tagplot),
  family = 'binomial',
  data = umbel.failure
)

anova(uf_p2, uf_p)
# phew! no polynomial effects.

