library(ggplot2)
library(dplyr)
library(tidyr)
library(lme4)

rm(list = ls())

all.data = read.csv('01_data_cleaning/out/demo_seed_phen_combined.csv')

head(all.data)
nrow(all.data)

###### Data subsetting

# A dataset where plants with missing seed counts are excluded
seed.no = all.data %>%
  filter(!is.na(total.seeds))

nrow(seed.no)

table(is.na(seed.no$mean.doy)) # 10 missing initiation days
table(is.na(seed.no$demo.umbels)) # 11 missing umbels in demo

# Subset out records that are complete
seed.subs = seed.no %>%
  filter(!is.na(mean.doy) & !is.na(demo.umbels)) %>%
  # center the phen data for convergence
  mutate(mean.doy = (mean.doy - round(mean(mean.doy))) / 7)

nrow(seed.subs)

##### Some crude plots

seed.subs %>%
  ggplot(aes(x = mean.doy, y = total.seeds)) +
  geom_point(aes(colour = trt), position = position_jitter(width = 1/7)) +
  scale_colour_manual(values = c('black', 'red', 'blue')) +
  facet_wrap(~ Year)
# sure looks quadratic to me

seed.subs %>%
  ggplot(aes(x = mean.doy, y = total.seeds)) +
  geom_point(aes(colour = trt), position = position_jitter(width = 1/7)) +
  scale_colour_manual(values = c('black', 'red', 'blue')) +
  facet_wrap(trt ~ Year)
# definitely looks like flowering starts earlier in 2021, maybe a little later
# in 2022

seed.subs %>%
  ggplot(aes(x = demo.umbels, y = total.seeds)) +
  geom_point(aes(colour = trt), position = position_jitter(width = 0.25)) +
  scale_colour_manual(values = c('black', 'red', 'blue')) +
  facet_wrap(trt ~ Year)
# okay well that's a relief at least lmao
# although it does look like some saturation is happening...

seed.subs %>%
  ggplot(aes(x = n.seed.counts, y = total.seeds)) +
  geom_point(aes(colour = trt), position = position_jitter(width = 0.25)) +
  scale_colour_manual(values = c('black', 'red', 'blue')) +
  facet_wrap(trt ~ Year)
# probably should go with the total seed counts

###### Some models

# Null model: number of seeds 
glmer(total.seeds ~ (1 | Year) + (1 | plot), family = 'poisson', data = seed.subs)
# residual devaince, as expected, is *much* larger than df
# overdispersion
# zero-inflated poisson would be good too, maybe in next batch of models

m.0 = glmer.nb(total.seeds ~ (1 | Year) + (1 | plot), data = seed.subs)
# still looking a little overdispersed...

m.t  = glmer.nb(total.seeds ~ trt + (1 | Year) + (1 | plot), data = seed.subs)
m.yt = glmer.nb(total.seeds ~ (trt | Year) + (1 | plot), data = seed.subs)
# aw rats... singularity issue

AIC(m.yt, m.t, m.0)
# wow. zero treatment effect lmao

m.d = glmer.nb(total.seeds ~ mean.doy + (1 | Year) + (1 | plot), data = seed.subs)

AIC(m.d, m.0)
# no effects of flowering day, jesus christ
summary(m.d)

m.dy = glmer.nb(total.seeds ~ (mean.doy | Year) + (1 | plot), data = seed.subs)
# singularity, again
summary(m.dy)
AIC(m.dy) # goodness
