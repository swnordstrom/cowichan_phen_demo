# Script testing for trade-offs in phenology
# - Effects on subsq year's survival
# - Effects on subsq year's growth
# - Effects on subsq year's flowering, seed set

# ---------------------------------------------
# Setup, data processing, etc.

library(ggplot2)
library(dplyr)
library(tidyr)
library(glmmTMB)

rm(list = ls())

# Read in data
all.data = read.csv('01_data_cleaning/out/demo_phen_seed_2016-2024_final.csv') %>%
  # Merge in treatment levels
  merge(read.csv('00_raw_data/plot_treatments.csv'), by.x = 'Plot', by.y = 'plot') %>%
  # We will only have records in 2021 or later
  filter(Year > 2020) %>%
  # Want only one record per plant
  distinct(plantid, Year, .keep_all = TRUE)

head(all.data)
nrow(all.data)

# Scrape out mean budding dates
phen = all.data %>%
  filter(in.phen) %>%
  separate_wider_delim(phen.julis, names = paste0('uu', 1:12), delim = ';', too_few = 'align_start') %>%
  pivot_longer(starts_with('uu'), names_to = 'umbel.number', values_to = 'phen.julian') %>%
  filter(!is.na(phen.julian)) %>%
  mutate(phen.julian = as.numeric(gsub('\\s', '', phen.julian))) %>%
  group_by(plantid, Year, No.leaves, Leaf.length, phen.umbels) %>%
  summarise(
    phen.mean = mean(phen.julian),
    # analysis of sd - probably will want bessel correction
    phen.sd = sd(phen.julian),
    # absolute first date of budding
    phen.first = min(phen.julian)
  ) %>%
  ungroup()

# Merge phen with next year's data
phen.subsq.demo = merge(
  # Phenology + demo in year t
  phen %>% 
    rename(No.umbels = phen.umbels) %>% 
    mutate(size = log(No.leaves * Leaf.length)) %>%
    select(-c(No.leaves, Leaf.length)),
  # Demo in year t+1
  all.data %>% 
    # Give me only demo records
    filter(in.demo) %>%
    select(Plot, plantid, year.tp1 = Year, No.leaves, Leaf.length, No.umbels, surv, trt) %>%
    mutate(year.t = year.tp1 - 1) %>%
    mutate(size = log(No.leaves * Leaf.length)) %>%
    select(-c(No.leaves, Leaf.length)),
  by.x = c('plantid', 'Year'), by.y = c('plantid', 'year.t'), suffixes = c('.t', '.tp1')
)

head(phen.subsq.demo)

# Get a data frame with growths only (looking for growth trade-offs)
phen.subsq.grow = phen.subsq.demo %>% 
  filter(surv, !is.na(size.t), !is.na(size.tp1)) %>%
  mutate(
    phen.mean.c = phen.mean - round(mean(phen.mean)),
    phen.first.c = phen.first - round(mean(phen.first))
  ) %>%
  group_by(Year = factor(Year)) %>%
  mutate(
    phen.mean.yc = phen.mean - round(mean(phen.mean)),
    phen.first.yc = phen.first - round(mean(phen.first))
  ) %>%
  ungroup()

# Getting rid of 2023-2024 (survival can't be safely estimated)
phen.subsq.surv = phen.subsq.demo %>% 
  filter(Year < 2023) %>%
  mutate(
    phen.mean.c = phen.mean - round(mean(phen.mean)),
    phen.first.c = phen.first - round(mean(phen.first))
  ) %>%
  group_by(Year) %>%
  mutate(
    phen.mean.yc = phen.mean - round(mean(phen.mean)),
    phen.first.yc = phen.first - round(mean(phen.first))
  ) %>%
  ungroup()
  
# ---------------------------------------------
# Do some visuals

# Flowering phenology against survival
phen.subsq.surv %>%
  mutate(surv = as.numeric(surv)) %>%
  ggplot(aes(x = phen.mean, y = surv)) +
  geom_point(
    aes(colour = trt),
    size = 3, position = position_jitter(height = .125), alpha = 0.5
  ) +
  scale_colour_manual(values = c('black', 'red', 'blue'))

phen.subsq.surv %>%
  mutate(surv = as.numeric(surv)) %>%
  ggplot(aes(x = phen.mean, y = surv)) +
  geom_point(
    aes(colour = trt),
    size = 3, position = position_jitter(height = .125), alpha = 0.5
  ) +
  scale_colour_manual(values = c('black', 'red', 'blue')) +
  facet_wrap(~ Year)

phen.subsq.surv %>%
  mutate(surv = as.numeric(surv)) %>%
  ggplot(aes(x = phen.first, y = surv)) +
  geom_point(
    aes(colour = trt),
    size = 3, position = position_jitter(height = .125), alpha = 0.5
  ) +
  scale_colour_manual(values = c('black', 'red', 'blue'))
# Also seems unlikely that anything will happen here.

# I don't think anything interesting will happen here

# Flowering phenology against growth
phen.subsq.grow %>%
  mutate(growth.t = size.tp1 - size.t) %>%
  ggplot(aes(x = phen.mean, y = growth.t)) +
  geom_point(aes(colour = trt), size = 3, alpha = 0.5) +
  scale_colour_manual(values = c('black', 'red', 'blue'))
# That definitely looks like it could be a negative relationship
# (may just be treatment effect, prior size effect, etc.)

phen.subsq.grow %>%
  ggplot(aes(x = size.t, y = size.tp1)) +
  geom_segment(aes(x = 2, xend = 6, y = 2, yend = 6), colour = 'gray', linetype = 2) +
  geom_point(aes(colour = phen.mean), size = 3, alpha = 0.5) +
  scale_colour_viridis_c()

# huh this actually is much closer to the 1-1 line than the typical growth kernel...

phen.subsq.grow %>%
  group_by(Year) %>%
  mutate(phen.yc = phen.mean - mean(phen.mean)) %>%
  ungroup() %>%
  ggplot(aes(x = size.t, y = size.tp1)) +
  geom_segment(aes(x = 2, xend = 6, y = 2, yend = 6), colour = 'gray', linetype = 2) +
  geom_point(aes(colour = phen.yc), size = 3, alpha = 0.5) +
  scale_colour_gradient2(low = 'red', high = 'blue', mid = 'white', midpoint = 0)
# ehh... there might be a slight association with earlier mean flowering

phen.subsq.grow %>%
  ggplot(aes(x = size.t, y = size.tp1)) +
  geom_segment(aes(x = 2, xend = 6, y = 2, yend = 6), colour = 'gray', linetype = 2) +
  geom_point(aes(colour = phen.first), size = 3, alpha = 0.5) +
  scale_colour_viridis_c()

phen.subsq.grow %>%
  mutate(growth.t = size.tp1 - size.t) %>%
  ggplot(aes(x = phen.first, y = growth.t)) +
  geom_point(aes(colour = trt), size = 3, alpha = 0.5) +
  scale_colour_manual(values = c('black', 'red', 'blue')) +
  facet_wrap(~ Year)
# ah this looks slightly different...
# there does look like there's a slight effect in 2021?
# but this could also just be a treatment effect.
# we will soon see...

# ---------------------------------------------
# Fit some models

### Survival models

s_0 = glmmTMB(
  surv ~ (1 | Plot),
  family = 'binomial',
  data = phen.subsq.surv
)

s_p = glmmTMB(
  surv ~ phen.mean.c + (1 | Plot),
  family = 'binomial',
  data = phen.subsq.surv
)

s_p2 = glmmTMB(
  surv ~ poly(phen.mean.c, 2) + (1 | Plot),
  family = 'binomial',
  data = phen.subsq.surv
)

s_py = glmmTMB(
  surv ~ phen.mean.yc + (1 | Plot),
  family = 'binomial',
  data = phen.subsq.surv
)

s_py2 = glmmTMB(
  surv ~ poly(phen.mean.yc, 2) + (1 | Plot),
  family = 'binomial',
  data = phen.subsq.surv
)

AIC(s_0, s_p, s_p2, s_py, s_py2) %>% mutate(daic = round(AIC - min(AIC), 2))
# Yep - no phen effects of mean bud date on survival

s_p = glmmTMB(
  surv ~ phen.first.c + (1 | Plot),
  family = 'binomial',
  data = phen.subsq.surv
)

s_p2 = glmmTMB(
  surv ~ poly(phen.first.c, 2) + (1 | Plot),
  family = 'binomial',
  data = phen.subsq.surv
)

s_py = glmmTMB(
  surv ~ phen.first.yc + (1 | Plot),
  family = 'binomial',
  data = phen.subsq.surv
)

s_py2 = glmmTMB(
  surv ~ poly(phen.first.yc, 2) + (1 | Plot),
  family = 'binomial',
  data = phen.subsq.surv
)

AIC(s_0, s_p, s_p2, s_py, s_py2) %>% mutate(daic = round(AIC - min(AIC), 2))
# same result when considering first bud date

### Growth

g_0a = glmmTMB(
  size.tp1 ~ size.t + Year + (1 | Plot / plantid),
  data = phen.subsq.grow
)

g_0b = glmmTMB(
  size.tp1 ~ size.t * Year + (1 | Plot / plantid),
  data = phen.subsq.grow
)

AIC(g_0a, g_0b)
# Okay - null model will have size-year interaction (different slopes for each year)

# Treatment effects if they are here
g_t = glmmTMB(
  size.tp1 ~ size.t * Year + trt + (1 | Plot / plantid),
  data = phen.subsq.grow
)

g_ty = glmmTMB(
  size.tp1 ~ size.t * Year + trt * Year + (1 | Plot / plantid),
  data = phen.subsq.grow
)

g_t.sy = glmmTMB(
  size.tp1 ~ size.t * Year + trt * size.t + (1 | Plot / plantid),
  data = phen.subsq.grow
)

AIC(g_t, g_ty, g_t.sy)
# Treatment-year effect pretty obvious
# (that AIC on the last one - yikes!)

summary(g_ty)
# Residual plot vs. phen
phen.subsq.grow %>%
  mutate(r = residuals(g_ty)) %>%
  pivot_longer(
    c(phen.mean.c, phen.first.c, phen.mean.yc, phen.first.yc),
    names_to = 'measure', values_to = 'date.c'
  ) %>%
  ggplot(aes(x = date.c, y = r, colour = trt)) +
  geom_point(size = 3, alpha = 0.5) +
  scale_colour_manual(values = c('black', 'red', 'blue')) +
  facet_wrap(~ measure)
# Seems unlikely... we'll see though
# Definitely not seeing any treatment effects though.

# Okay. Now look at phenology

g_m = glmmTMB(
  size.tp1 ~ size.t * Year + trt * Year + phen.mean.c + (1 | Plot / plantid),
  data = phen.subsq.grow
)

g_m2 = glmmTMB(
  size.tp1 ~ size.t * Year + trt * Year + poly(phen.mean.c, 2) + (1 | Plot / plantid),
  data = phen.subsq.grow
)

g_my = glmmTMB(
  size.tp1 ~ size.t * Year + trt * Year + phen.mean.yc + (1 | Plot / plantid),
  data = phen.subsq.grow
)

g_my2 = glmmTMB(
  size.tp1 ~ size.t * Year + trt * Year + poly(phen.mean.yc, 2) + (1 | Plot / plantid),
  data = phen.subsq.grow
)

AIC(g_ty, g_m, g_m2, g_my, g_my2)
# hmm...
# well at least there's no polynomial effect
# the year-centered means and overall means have the same effect.

# let's also try the first bud date just for comparison

g_f = glmmTMB(
  size.tp1 ~ size.t * Year + trt * Year + phen.first.c + (1 | Plot / plantid),
  data = phen.subsq.grow
)

g_f2 = glmmTMB(
  size.tp1 ~ size.t * Year + trt * Year + poly(phen.first.c, 2) + (1 | Plot / plantid),
  data = phen.subsq.grow
)

g_fy = glmmTMB(
  size.tp1 ~ size.t * Year + trt * Year + phen.first.yc + (1 | Plot / plantid),
  data = phen.subsq.grow
)

g_fy2 = glmmTMB(
  size.tp1 ~ size.t * Year + trt * Year + poly(phen.first.yc, 2) + (1 | Plot / plantid),
  data = phen.subsq.grow
)

AIC(g_ty, g_f, g_f2, g_fy, g_fy2) %>% mutate(daic = round(AIC - min(AIC), 2))
# performs just slightly worse

# Okay so there *is* a phen effect

summary(g_m)
# it is negative... but it is incredibly small lmao
# one day of later flying causes a decrease in growth of <.01

# Look at another residual plot
phen.subsq.grow %>%
  mutate(r = residuals(g_m)) %>%
  ggplot(aes(x = trt, y = r)) +
  geom_point(size = 3, position = position_jitter(width = 1/4))
# seems unlikely

g_my = glmmTMB(
  size.tp1 ~ size.t * Year + trt * Year + phen.mean.c * Year + (1 | Plot / plantid),
  data = phen.subsq.grow
)

g_mt = glmmTMB(
  size.tp1 ~ size.t * Year + trt * Year + phen.mean.c * trt + (1 | Plot / plantid),
  data = phen.subsq.grow
)

g_ms = glmmTMB(
  size.tp1 ~ size.t * Year + trt * Year + phen.mean.c * size.t + (1 | Plot / plantid),
  data = phen.subsq.grow
)

AIC(g_my, g_mt, g_ms, g_m)
# phew... these models suck

# Okay so there is a very small growth trade-off to budding later

### AOV-style partitioning of variance in bud date?
# note - only four years... lots of uncertainty in that year-variance estimate...

phen.plot = phen %>% 
  mutate(tagplot = gsub('\\.[0-9]{1,2}[A-Z]', '', plantid)) %>%
  separate(tagplot, into = c('tag', 'plot'), sep = '_') %>%
  merge(read.csv('00_raw_data/plot_treatments.csv'))

v_0 = glmmTMB(
  phen.mean ~ trt + (1 | plot) + (1 | Year),
  data = phen.plot
)

summary(v_0)

v_i = glmmTMB(
  phen.mean ~ trt + (1 | plot / plantid) + (1 | Year),
  data = phen.plot
)

AIC(v_i, v_0)
# I mean... yeah. to be expected

# Breakdown of variance
summary(v_i)
# variance (within a plant) cut down by ~25%
