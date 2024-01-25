# Script for assessing budding phenology and umbel deaths
# SN - init 19 Jan 2024

# Load packages
library(ggplot2)
library(dplyr)
library(tidyr)
# library(lme4) # not working with my updated version of Matrix...
library(glmmTMB)

# Clear namespace
rm(list = ls())

# Read in dataset
buds.all = read.csv('01_data_cleaning/out/phenology_buds_deaths_all.csv') %>%
  # Add in treatment data (because it's useful)
  merge(y = read.csv('00_raw_data/plot_treatments.csv'))

head(buds.all)

# Visualize 
# make plot where one row corresponds to one plant, time is on x-axis

buds.all %>% 
  group_by(plantid) %>%
  mutate(minday = min(init.doy)) %>%
  ungroup() %>%
  arrange(minday, trt, plantid) %>%
  group_by(year) %>%
  # trust me, this works for assigning an integer ID to each plant
  mutate(plant.i = cumsum(!duplicated(plantid))) %>%
  ungroup() %>%
  ggplot(aes(x = init.doy, y = plant.i)) +
  geom_point(aes(colour = trt, shape = varb)) +
  scale_shape_manual(values = c(4, 19), '') +
  scale_colour_manual(values = c('black', 'red', 'blue'), '') +
  facet_wrap(~ year)
# you can kinda see from here 

# Now facet out by treatment
buds.all %>% 
  group_by(plantid) %>%
  mutate(minday = min(init.doy)) %>%
  ungroup() %>%
  arrange(minday, trt, plantid) %>%
  group_by(year, trt) %>%
  mutate(plant.i = cumsum(!duplicated(plantid))) %>%
  ungroup() %>%
  ggplot(aes(x = init.doy, y = plant.i)) +
  geom_point(aes(colour = trt, shape = varb)) +
  scale_shape_manual(values = c(4, 19), '') +
  scale_colour_manual(values = c('black', 'red', 'blue'), '') +
  facet_wrap(year ~ trt)
# eh

# Now... change the plantIDs to be consistent across years

buds.all = buds.all %>%
  separate(plantid, into = c('tag', 'plot', 'coord'), sep = '_', remove = FALSE, fill = 'right') %>%
  select(-coord) %>%
  unite(tagplot, c(tag, plot), remove = FALSE) %>%
  mutate(phen.scaled = (init.doy - round(mean(init.doy))) / 7)

# Subset out new buds only
buds.new = buds.all %>% filter(varb %in% 'new')

# Get a data frame with mean budding day and total number of good/dead umbels
budday.deaths = buds.all %>%
  group_by(tagplot, year, plot, trt) %>%
  summarise(
    n.total.buds = sum(varb %in% 'new'),
    n.dead = sum(varb %in% 'lost'),
    mean.doy = mean(init.doy[varb %in% 'new'])
  ) %>%
  ungroup() %>%
  mutate(phen.scaled = (mean.doy - mean(mean.doy)) / 7)

# Check - are there any plants that have more dead umbels than total ones?
any(budday.deaths$n.dead > budday.deaths$n.total.buds)
# nope = good!

##### Mean date of budding - analysis

# Null model
p_0 = glmmTMB(
  init.doy ~ (1 | plot / tagplot) + (1 | year),
  data = buds.new
)

summary(p_0)

p_t = glmmTMB(
  init.doy ~ trt + (1 | plot / tagplot) + (1 | year),
  data = buds.new
)

summary(p_t)

anova(p_t, p_0)
# well... okay cool, likelihood ratio test suggests an effect of drought on budding phen

predict(
  object = p_t,
  newdata = data.frame(trt = c('control', 'drought', 'irrigated')),
  type = 'response',
  re.form = ~ 0,
  allow.new.levels = TRUE
)

# so ~3 days sooner under drought, ~2 days later under irrigation

##### Date of budding vs. likelihood of bud mortality

budday.deaths %>%
  mutate(p.dead = n.dead / n.total.buds) %>%
  ggplot(aes(x = mean.doy, y = p.dead)) +
  geom_point(aes(colour = trt), size = 3) +
  scale_colour_manual(values = c('black', 'red', 'blue')) +
  facet_wrap(~ year)

d_0 = glmmTMB(
  cbind(n.dead, n.total.buds - n.dead) ~ (1 | plot / tagplot) + (1 | year),
  data = budday.deaths,
  family = 'binomial'
)

summary(d_0)

d_t = glmmTMB(
  cbind(n.dead, n.total.buds - n.dead) ~ trt + (1 | plot / tagplot) + (1 | year),
  data = budday.deaths,
  family = 'binomial'
)

summary(d_t)

d_d = glmmTMB(
  cbind(n.dead, n.total.buds - n.dead) ~ phen.scaled + (1 | plot / tagplot) + (1 | year),
  data = budday.deaths,
  family = 'binomial'
)

summary(d_d)

anova(d_d, d_0) # looks like a significant effect of day
anova(d_t, d_0) # not of treatment though lol

d_d_t = glmmTMB(
  cbind(n.dead, n.total.buds - n.dead) ~ phen.scaled + trt + (1 | plot / tagplot) + (1 | year),
  data = budday.deaths,
  family = 'binomial'
)

anova(d_d_t, d_d)
# no effect of treatment, unless...

d_dt = glmmTMB(
  cbind(n.dead, n.total.buds - n.dead) ~ phen.scaled * trt + (1 | plot / tagplot) + (1 | year),
  data = budday.deaths,
  family = 'binomial'
)

anova(d_dt, d_d)
# nope!
# no treatment effect period

# one more... does doy effect vary by year
d_dy = glmmTMB(
  cbind(n.dead, n.total.buds - n.dead) ~  (1 | plot / tagplot) + (phen.scaled | year),
  data = budday.deaths,
  family = 'binomial'
)

anova(d_dy, d_d)
# whoa... cool, significant effect
AIC(d_dy, d_d, d_0)

# neat.

summary(d_dy)
# huh... okay how do I get the mean effect of phen.scaled out of here lol
ranef(d_dy) # oh hell yeah
ranef(d_dy)$cond$year
# so... negative effect in 2021, weak neutral-ish effects in other years

# In just the model with one overall effects
summary(d_d)
# one *week* of later budding decreases log-odds by -0.17, or, about 16%

