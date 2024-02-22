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
# size (so we can keep records where the prior year's size is missing)
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
# for modeling umbel counts
demo.for.umbel.counts = demo.for.flowering %>%
  filter(flowering) %>%
  select(-flowering)

### Demo + phen

phen.demo = phen.demo %>%
  # Change missing umbel counts to zero
  mutate(n.demo.umbels = ifelse(is.na(n.demo.umbels), 0, n.demo.umbels)) %>%
  # Center phenology
  mutate(centered.phen = init.doy - round(mean(init.doy)))

head(phen.demo.with.size)

### Mean bud date against number of dead umbels
umbel.failure = phen.demo %>%
  # Get mean budding date for each plant
  group_by(Year, tagplot) %>%
  mutate(mean.doy = mean(init.doy)) %>%
  ungroup() %>%
  # Distinct here because we need only one row per plant
  distinct(Year, tagplot, .keep_all = TRUE) %>%
  # Center the phenology column to help with model convergence
  mutate(centered.phen = mean.doy - round(mean(mean.doy)))

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
  # and all umbels with zero seeds, but removing dead umbels from phen
  seed.phen.demo %>%
    filter(!no.seeds) %>%
    group_by(tagplot, Year) %>%
    filter((1:n()) > n.lost.umbel) %>%
    ungroup()
)

nrow(seed.phen.demo)

nrow(seed.phen.demo)
seed.phen.demo %>% distinct(Year, tagplot) %>% nrow()


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
# Current size is a *much* better predictor than prior size

# Refitting current size model with whole dataset - will be used as null model
# moving forward
pf_s = glmmTMB(
  flowering ~ cur.size + (1 | Year) + (1 | Plot / plantid),
  data = demo.for.flowering,
  family = 'binomial',
)

# Model with treatment intercept
pf_s_t = glmmTMB(
  flowering ~ cur.size + trt + (1 | Year) + (1 | Plot / plantid),
  data = demo.for.flowering,
  family = 'binomial',
)

# Model with size-treatment interaction
pf_st = glmmTMB(
  flowering ~ cur.size * trt + (1 | Year) + (1 | Plot / plantid),
  data = demo.for.flowering,
  family = 'binomial',
)

anova(pf_st, pf_s_t, pf_s)
anova(pf_st, pf_s)
AIC(pf_st, pf_s) %>% mutate(daic = round(AIC - min(AIC), 2))
# It's not unbelievably strong, but according to our cutoffs, the size-treatment
# interaction model performs best.

# Test against a model with size effects varying by year
pf_sy = glmmTMB(
  flowering ~  (cur.size | Year) + (1 | Plot / plantid),
  data = demo.for.flowering,
  family = 'binomial',
)

# Adding treatment to the above model
pf_sy_t = glmmTMB(
  flowering ~ trt + (cur.size | Year) + (1 | Plot / plantid),
  data = demo.for.flowering,
  family = 'binomial',
)

AIC(pf_sy_t, pf_sy, pf_st) %>% mutate(daic = round(AIC - min(AIC), 2))
# no effects of year-varying effects

# Final model of these
summary(pf_st)
# Without mean-centering size effects it's not super easy to discern what's going on
# but... maybe large plants in drought are less likely to flower

### Umbel count models

nu_s = glmmTMB(
  No.umbels ~ cur.size + (1 | Year) + (1 | Plot / plantid),
  data = demo.for.umbel.counts,
  family = 'truncated_poisson'
)

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
# Doesn't look like there's any evidence of a year or treatment effect here...

summary(nu_s)
# larger plants produce more umbels... seems unsurprising

# NOTE that I had convergence issues with the negative binomial based on this -
# not sure if overdispersion is an issue here.

### Umbel budding phenology

phen_0 = glmmTMB(
  formula = centered.phen ~ (1 | Year) + (1 | Plot / tagplot),
  data = phen.demo,
)

phen_s = glmmTMB(
  formula = centered.phen ~ cur.size + (1 | Year) + (1 | Plot / tagplot),
  data = phen.demo,
)

anova(phen_s, phen_0)
# no evidence of a size relationship

# test with treatment effects just to be sure

phen_t = glmmTMB(
  formula = centered.phen ~ trt + (1 | Year) + (1 | Plot / tagplot),
  data = phen.demo,
)

phen_s_t = glmmTMB(
  formula = centered.phen ~ trt + cur.size + (1 | Year) + (1 | Plot / tagplot),
  data = phen.demo,
)

anova(phen_s_t, phen_t)
# nothing here either
# (I also tried an interaction)
# sweet!

# Now, test for a treatment effect on phenology

anova(phen_t, phen_0)

summary(phen_t)
# Drought means buding 3.25 days earlier
# irrigation might mean budding a day later

# I think it's worth doing a residual check here...
hist(residuals(phen_t)) # that's actually fairly normal
qqnorm(residuals(phen_t))
qqline(residuals(phen_t))
# Tails are not great, esp. lower tail, but otherwise actually looks decent

phen.demo %>%
  mutate(residual = residuals(phen_t)) %>%
  ggplot(aes(x = residual, fill = trt)) +
  geom_histogram(position = position_identity(), alpha = 0.25) +
  scale_fill_manual(values = c('black', 'red', 'blue'))

### Umbel failure
# (Maybe it's smarter to parameterize this as surviving umbels?)

# First try with size

uf_0 = glmmTMB(
  cbind(n.phen.umbels - n.lost.umbels, n.lost.umbels) ~ 1 +
    (1 | Year) + (1 | Plot / tagplot),
  family = 'binomial',
  data = umbel.failure
)

uf_s = glmmTMB(
  cbind(n.phen.umbels - n.lost.umbels, n.lost.umbels) ~ cur.size +
    (1 | Year) + (1 | Plot / tagplot),
  family = 'binomial',
  data = umbel.failure
)

anova(uf_s, uf_0)
# there *is* a pronounced size effect here

summary(uf_s)
# larger plants more likely to have umbel success, although effect is not super strong
# 1 unit increase in size means ~48% increase in odds of umbel succeeding

# Evidence of size-effects varying by year?

uf_sy = glmmTMB(
  cbind(n.phen.umbels - n.lost.umbels, n.lost.umbels) ~
    (cur.size | Year) + (1 | Plot / tagplot),
  family = 'binomial',
  data = umbel.failure
)

anova(uf_sy, uf_s)
# No year-varying effects

# Look for treatment effects
# first - treatment standalone
uf_s_t = glmmTMB(
  cbind(n.phen.umbels - n.lost.umbels, n.lost.umbels) ~  cur.size + trt +
    (1 | Year) + (1 | Plot / tagplot),
  family = 'binomial',
  data = umbel.failure
)

# Now - treatment interaction with size
uf_st = glmmTMB(
  cbind(n.phen.umbels - n.lost.umbels, n.lost.umbels) ~  cur.size * trt +
    (1 | Year) + (1 | Plot / tagplot),
  family = 'binomial',
  data = umbel.failure
)

anova(uf_s_t, uf_s)
anova(uf_st, uf_s)
# no evidence of treatment effects

uf_s_p = glmmTMB(
  cbind(n.phen.umbels - n.lost.umbels, n.lost.umbels) ~ centered.phen + cur.size +
    (1 | Year) + (1 | Plot / tagplot),
  family = 'binomial',
  data = umbel.failure
)

uf_sp = glmmTMB(
  cbind(n.phen.umbels - n.lost.umbels, n.lost.umbels) ~ centered.phen * cur.size +
    (1 | Year) + (1 | Plot / tagplot),
  family = 'binomial',
  data = umbel.failure
)

uf_sp_t = glmmTMB(
  cbind(n.phen.umbels - n.lost.umbels, n.lost.umbels) ~ centered.phen * cur.size + trt +
    (1 | Year) + (1 | Plot / tagplot),
  family = 'binomial',
  data = umbel.failure
)

uf_s_p_t = glmmTMB(
  cbind(n.phen.umbels - n.lost.umbels, n.lost.umbels) ~ centered.phen + cur.size + trt +
    (1 | Year) + (1 | Plot / tagplot),
  family = 'binomial',
  data = umbel.failure
)

uf_s_p2 = glmmTMB(
  cbind(n.phen.umbels - n.lost.umbels, n.lost.umbels) ~ poly(centered.phen, 2) + cur.size +
    (1 | Year) + (1 | Plot / tagplot),
  family = 'binomial',
  data = umbel.failure
)

uf_s_p2_t = glmmTMB(
  cbind(n.phen.umbels - n.lost.umbels, n.lost.umbels) ~ poly(centered.phen, 2) + cur.size + trt +
    (1 | Year) + (1 | Plot / tagplot),
  family = 'binomial',
  data = umbel.failure
)

uf_s_py = glmmTMB(
  cbind(n.phen.umbels - n.lost.umbels, n.lost.umbels) ~ cur.size +
    (centered.phen | Year) + (1 | Plot / tagplot),
  family = 'binomial',
  data = umbel.failure
)

uf_s_py_t = glmmTMB(
  cbind(n.phen.umbels - n.lost.umbels, n.lost.umbels) ~ cur.size + trt +
    (centered.phen | Year) + (1 | Plot / tagplot),
  family = 'binomial',
  data = umbel.failure
)

# uf_sy_p = glmmTMB(
#   cbind(n.phen.umbels - n.lost.umbels, n.lost.umbels) ~ mean.phen +
#     (cur.size | Year) + (1 | Plot / tagplot),
#   family = 'binomial',
#   data = umbel.failure
# )
# didn't converge

# uf_sy_py = glmmTMB(
#   cbind(n.phen.umbels - n.lost.umbels, n.lost.umbels) ~
#     (cur.size + mean.phen | Year) + (1 | Plot / tagplot),
#   family = 'binomial',
#   data = umbel.failure
# )
# failed to converge

# uf_s_p2y = glmmTMB(
#   cbind(n.phen.umbels - n.lost.umbels, n.lost.umbels) ~ cur.size +
#     (poly(mean.phen, 2) | Year) + (1 | Plot / tagplot),
#   family = 'binomial',
#   data = umbel.failure
# )
# Quadratic phen varying by year fails to converge

AIC(uf_s, uf_sp, uf_sp_t, uf_s_p, uf_s_p_t, uf_s_p2, uf_s_p2_t, uf_s_py, uf_s_py_t) %>% 
  mutate(daic = round(AIC - min(AIC), 2)) %>%
  arrange(daic)
# Linear phen effect, varying by year

summary(uf_s_py)
# right... this doesn't report mean phen effects...

ranef(uf_s_py)$cond$Year
ranef(uf_s_py)$cond$Year %>% as.data.frame() %>% apply(2, mean)

summary(uf_sp)
# according to this model - larger plants more likely to have umbel success
# budding later means more success for smaller plants
# budding later maybe less impactful for larger plants

# Look at year as fixed effect...
uf_s_py_fixed = glmmTMB(
  cbind(n.phen.umbels - n.lost.umbels, n.lost.umbels) ~ cur.size +
    centered.phen * Year + (1 | Plot / tagplot),
  family = 'binomial',
  data = umbel.failure %>% mutate(Year = factor(Year))
)

AIC(uf_s_py_fixed, uf_s_py)
# huh... works much better as a fixed effect lmao

summary(uf_s_py_fixed)
# effect is strongest in 2021, still apparent in 2022, apparently gone in 2023

### Seeds per umbel

# First, look for effects of previous size
s_0 = glmmTMB(
  no.seeds ~ (1 | Plot / tagplot) + (1 | Year),
  data = seed.phen.demo,
  family = 'nbinom2'
)

s_s = glmmTMB(
  no.seeds ~ cur.size + (1 | Plot / tagplot) + (1 | Year),
  data = seed.phen.demo,
  family = 'nbinom2'
)

anova(s_0, s_s)
# yep - there's a size effect
# (probably picking up umbel size, which increases with plant size)

s_s_t = glmmTMB(
  no.seeds ~ cur.size + trt + (1 | Plot / tagplot) + (1 | Year),
  data = seed.phen.demo,
  family = 'nbinom2'
)

s_st = glmmTMB(
  no.seeds ~ cur.size * trt + (1 | Plot / tagplot) + (1 | Year),
  data = seed.phen.demo,
  family = 'nbinom2'
)

AIC(s_st, s_s_t, s_s) %>% mutate(daic = round(AIC - min(AIC), 2))
# wow - a size-treatment effect here? lol

summary(s_st)

# Test for phenology
s_st_p = glmmTMB(
  no.seeds ~ centered.phen + cur.size * trt + (1 | Plot / tagplot) + (1 | Year),
  data = seed.phen.demo,
  family = 'nbinom2'
) 

anova(s_st_p, s_st)
# yep, phen effect is here
summary(s_st_p)
# - negative phen slope (in controls): fewer seeds with earlier flowering
# - positive intercept for irrigation: more seeds in irrigation (across the board)
# - negative size-irrigation effect: larger plants don't benefit as much from being large

# Polynomial fit?
s_st_p2 = glmmTMB(
  no.seeds ~ poly(centered.phen, 2) + cur.size * trt + (1 | Plot / tagplot) + (1 | Year),
  data = seed.phen.demo,
  family = 'nbinom2'
) 

anova(s_st_p, s_st_p2)
summary(s_st_p2)
# definitely not quadratic

# Treatment against polynomial phenology
s_st_p2t = glmmTMB(
  no.seeds ~ poly(centered.phen, 2) * trt + cur.size * trt + (1 | Plot / tagplot) + (1 | Year),
  data = seed.phen.demo,
  family = 'nbinom2'
) 

# Size-phen interaction?
s_st_sp = glmmTMB(
  no.seeds ~ centered.phen * cur.size + trt * cur.size + (1 | Plot / tagplot) + (1 | Year),
  data = seed.phen.demo,
  family = 'nbinom2'
) 

anova(s_st_p, s_st_sp)
# no size-phen interaction

# Phen-treatment interaction?
s_st_pt = glmmTMB(
  no.seeds ~ centered.phen * trt + trt * cur.size + (1 | Plot / tagplot) + (1 | Year),
  data = seed.phen.demo,
  family = 'nbinom2'
) 

anova(s_st_p, s_st_pt)
# no phen-treatment interaction

# Phen-effects varying by year?
s_st_py = glmmTMB(
  no.seeds ~ cur.size * trt + (1 | Plot / tagplot) + (centered.phen | Year),
  data = seed.phen.demo,
  family = 'nbinom2'
) 

anova(s_st_p, s_st_py)
# definitely no year-varying phen effect

AIC(s_st_p, s_st_sp, s_st_p2, s_st_pt,  s_st_p2t, s_st_py) %>% 
  mutate(daic = round(AIC - min(AIC), 2))
# Best model has a linear effect of phenology, no interactions

# Best model looks like s_st_p
summary(s_st_p)

#########################################################
# Figure drafts
#########################################################

#################
# Probability of flowering

pr.flower.data.base = demo.for.flowering %>%
  mutate(flowering.01 = as.numeric(flowering)) %>%
  ggplot(aes(x = cur.size, y = flowering.01, colour = trt)) +
  geom_point(position = position_jitter(height = 0.01)) +
  scale_colour_manual(values = c('black', 'red', 'blue'), 'treatment') +
  facet_wrap(~ Year) +
  labs(x = 'Size', y = 'Flowering') +
  theme(
    panel.background = element_blank(),
    legend.position = 'top'
  )

# Model is pf_st
pr.flower.preds = expand.grid(
  cur.size = (10:50)/10,
  trt = c('control', 'drought', 'irrigated')
) %>%
  mutate(
    pred.prob = predict(
      object = pf_st,  newdata = ., 
      re.form = ~ 0, allow.new.levels = TRUE, 
      type = 'response'
    )
  )

pr.flower.data.base +
  geom_line(
    data = pr.flower.preds,
    inherit.aes = FALSE,
    aes(
      x = cur.size, y = pred.prob, 
      colour = trt
    )
  )


#################
# Number of umbels / plant

# model is nu_s_f

u.count.data.base = demo.for.umbel.counts %>%
  ggplot(aes(x = cur.size, y = No.umbels, colour = trt)) +
  geom_point(position = position_jitter(height = 0.05)) +
  scale_colour_manual(values = c('black', 'red', 'blue'), 'treatment') +
  labs(x = 'Previous size', y = 'Number of umbels') +
  facet_wrap(~ Year) +
  theme(
    panel.background = element_blank(),
    legend.position = 'top'
  )

u.count.data.base

u.count.preds = data.frame(cur.size = (10:60)/10) %>%
  mutate(
    pred.umb = predict(
      object = nu_s,  newdata = ., 
      re.form = ~ 0, allow.new.levels = TRUE, 
      type = 'response'
    )
  )

u.count.data.base +
  geom_line(
    data = u.count.preds,
    inherit.aes = FALSE,
    aes(x = cur.size, y = pred.umb),
    linewidth = 2
  )

#################
# Probability of umbels failing

# optimal model is uf_s_py, but uf_sp is also plausible

u.succ.data.base = umbel.failure %>%
  mutate(p.succ.umbel = (n.phen.umbels - n.lost.umbels) / n.phen.umbels) %>%
  ggplot(aes(x = centered.phen, y = p.succ.umbel, colour = trt)) +
  geom_point(
    aes(size = n.phen.umbels), 
    size = 2, shape = 21,
    position = position_jitter(height = 1/50)
  ) +
  scale_colour_manual(values = c('black', 'red', 'blue'), 'treatment') +
  scale_size_continuous('Number of umbels') +
  scale_x_continuous(breaks = 10*(-3:3), labels = c('7 apr', '14 apr', '27 apr', '7 may', '17 may', '27 may', '6 jun')) +
  labs(x = 'Day of budding', y = 'Probability of umbel success') +
  facet_wrap(~ Year) +
  theme(
    panel.background = element_blank(),
    axis.text.x = element_text(angle = 90),
    legend.position = 'top'
  )

u.succ.data.base

u.succ.pred.py = expand.grid(
  centered.phen = -25:30,
  size.sd = -2:2
) %>%
  mutate(cur.size = mean(umbel.failure$cur.size) + size.sd * sd(umbel.failure$cur.size)) %>%
  mutate(
    pred.prob = predict(
      object = uf_s_py,
      newdata = .,
      type = 'response',
      re.form = ~ 0,
      allow.new.levels = TRUE
    )
  )

u.succ.data.base +
  geom_line(
    data = u.succ.pred.py,
    inherit.aes = FALSE,
    aes(x = centered.phen, y = pred.prob, group = size.sd, colour = size.sd)
  )
# uh... maybe that slope effect isn't being picked up because it's in a random effect...?
# it definitely looks like it should matter here...

u.succ.pred.sp = expand.grid(
  centered.phen = -25:30,
  size.sd = -2:2
) %>%
  mutate(cur.size = mean(umbel.failure$cur.size) + size.sd * sd(umbel.failure$cur.size)) %>%
  mutate(
    pred.prob = predict(
      object = uf_sp,
      newdata = .,
      type = 'response',
      re.form = ~ 0,
      allow.new.levels = TRUE
    )
  )

u.succ.pred.py.fixed = expand.grid(
  centered.phen = -25:30,
  Year = factor(2021:2023),
  size.sd = -2:2
) %>%
  mutate(cur.size = mean(umbel.failure$cur.size) + size.sd * sd(umbel.failure$cur.size)) %>%
  mutate(
    pred.prob = predict(
      object = uf_s_py_fixed,
      newdata = .,
      type = 'response',
      re.form = ~ 0,
      allow.new.levels = TRUE
    )
  )

u.succ.data.base +
  geom_line(
    data = u.succ.pred.sp,
    inherit.aes = FALSE,
    aes(x = centered.phen, y = pred.prob, group = size.sd, colour = size.sd)
  )
# yeah that's more like it.
# phen matters less for umbel failure of the large plants.
# the 2023 data though... lmao that looks like a bad fit.

u.succ.data.base +
  geom_line(
    data = u.succ.pred.py.fixed %>% filter(!size.sd),
    inherit.aes = FALSE,
    aes(x = centered.phen, y = pred.prob)
  )
# Interesting...

#################
# Seeds / reproductive umbel

seed.set.base = seed.phen.demo %>%
  ggplot(aes(x = centered.phen, y = no.seeds, colour = trt)) +
  geom_point(shape = 21, size = 2) +
  scale_colour_manual(values = c('black', 'red', 'blue'), 'treatment') +
  scale_x_continuous(breaks = 10*(-3:3), labels = c('7 apr', '14 apr', '27 apr', '7 may', '17 may', '27 may', '6 jun')) +
  labs(x = 'Day of budding', y = 'Number of seeds') +
  facet_wrap(~ Year, ncol = 3) +
  theme(
    panel.background = element_blank(),
    axis.text.x = element_text(angle = 90),
    legend.position = 'top'
  )

seed.set.base

seed.set.preds = expand.grid(
  size.sd = 2 * (-1:1),
  # size.sd = 0,
  centered.phen = -20:30,
  trt = c('control', 'drought', 'irrigated')
) %>%
  mutate(cur.size = mean(seed.phen.demo$cur.size) + size.sd * sd(seed.phen.demo$cur.size)) %>%
  mutate(
    pred.seed = predict(
      object = s_st_p,
      newdata = .,
      re.form = ~ 0,
      type = 'response',
      allow.new.levels = TRUE
    )
  )

seed.set.base +
  geom_line(
    data = seed.set.preds %>% filter(!size.sd),
    inherit.aes = FALSE,
    aes(
      x = centered.phen, y = pred.seed, 
      group = interaction(size.sd, trt),
      colour = trt
    ),
    linewidth = 1.2,
  ) +
  scale_linetype_manual(values = c(2, 1, 2))

# Huh...
# why are the irrigated plants not below the drought plants?

#################
# Raw phenology plot
# rows: years + treatments
# y axis: day of year
# (different data input...)

phen.plot.data.base = phen.demo %>%
  mutate(year.dodge = Year + -0.2 * as.numeric(trt %in% 'drought') + 0.2 * as.numeric(trt %in% 'irrigated')) %>%
  ggplot(aes(x = centered.phen, y = year.dodge, colour = trt)) +
  geom_point(position = position_jitter(height = 0.08, width = 2), size = 2) +
  scale_colour_manual(values = c('black', 'red', 'blue')) +
  # scale_y_reverse(breaks = 2020:2023, labels = c('overall', '2021', '2022', '2023')) +
  scale_x_continuous(breaks = 10*(-3:3), labels = c('7 apr', '14 apr', '27 apr', '7 may', '17 may', '27 may', '6 jun')) +
  # date column is centered in data frame such that 0 = may 7 (day 126)
  # use as.Date(10*(-2:3) + 126) to get these
  scale_y_continuous(breaks = 2020:2023, labels = c('overall', '2021', '2022', '2023')) +
  labs(x = 'Day of year', y = '') +
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
  ggplot(aes(x = centered.phen, y = year.dodge, colour = trt)) +
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
  ggplot(aes(x = centered.phen, y = p.succ.umbel, colour = trt)) +
  geom_point(size = 1, alpha = 0.5) +
  geom_line(
    data = u.succ.pred.sp %>% filter(!size.sd),
    inherit.aes = FALSE,
    aes(x = centered.phen, y = pred.prob)
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

pan.c = seed.phen.demo %>%
  ggplot(aes(x = centered.phen, y = no.seeds, colour = trt)) +
  geom_point(size = 1, alpha = 0.5) +
  geom_line(
    data = seed.set.preds %>% filter(!size.sd), # mutate(size.sd = factor(size.sd)),
    inherit.aes = FALSE,
    aes(
      x = centered.phen, y = pred.seed, 
      colour = trt
    )
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
) # %>%
#   save_plot(
#     filename = '02_data_exploration/prelim_analysis/fig_phen_repro.png',
#     base_height = 5, base_width = 8
#   )

#########################################################
# Trying to get a composite mean seed per plant
#########################################################

# NOT adjusted for current size changes

skeleton = expand.grid(
  prev.size = (5:50)/10,
  prev.flower = c(TRUE, FALSE),
  trt = c('control', 'drought', 'irrigated'),
  mean.phen = c(-20:30)
) %>%
  mutate(centered.phen = mean.phen) %>%
  mutate(
    # Probability of flowering
    pred.prob.flower = predict(
      object = pf_s_ft,
      newdata = .,
      allow.new.levels = TRUE,
      re.form = ~ 0,
      type = 'response'
    ),
    # Umbels per successfully flowering plant
    pred.n.umbels = predict(
      object = nu_s_f,
      newdata = .,
      allow.new.levels = TRUE,
      re.form = ~ 0,
      type = 'response'
    ),
    # Probability of succeeding umbel
    pred.p.umbel.succ = predict(
      object = uf_p,
      newdata = .,
      allow.new.levels = TRUE,
      re.form = ~ 0,
      type = 'response'
    ),
    # Seeds per surviving umbel
    pred.seeds = predict(
      object = s_s_p,
      newdata = .,
      allow.new.levels = TRUE,
      re.form = ~ 0,
      type = 'response'
    )
  )
  
skeleton = skeleton %>%
  mutate(reprod.output = pred.prob.flower * pred.n.umbels * pred.p.umbel.succ * pred.seeds)

skeleton %>%
  mutate(prev.flower = factor(prev.flower), reprod.output = log(reprod.output)) %>%
  ggplot(aes(x = mean.phen, y = reprod.output, colour = trt)) +
  geom_line(
    aes(
      group = interaction(prev.size, prev.flower, trt),
      linetype = prev.flower
    )
  )
# useless plot

skeleton %>%
  ggplot(aes(x = prev.size, y = mean.phen)) +
  geom_raster(aes(fill = log(reprod.output))) +
  scale_fill_viridis_c() +
  facet_wrap(prev.flower ~ trt)

# okay... what is obvious from here is that size probably matters more than phen!
# What to make of growth differences, then?
