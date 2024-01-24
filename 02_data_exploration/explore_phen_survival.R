# Looking at the subsetted phen and demo data to see if flowering times in 2021
# and 2022 were associated with survival in the next year.
# Also throwing in umbel counts and demo
# SN - init 24 Jan 2024

library(ggplot2)
library(dplyr)
library(tidyr)
library(glmmTMB)

rm(list = ls())

phendemo = read.csv('01_data_cleaning/out/demo_seed_phen_by_umbel_combined.csv')
demodemo = read.csv('01_data_cleaning/out/demo_imputed_survival.csv')

head(phendemo)
head(demodemo)

mergedemo = merge(
  x = phendemo %>% select(Year, Plot = plot, trt, finalid, tagplot, mean.doy, demo.umbels, no.seeds),
  y = demodemo %>% mutate(Year = Year - 1) %>% select(Year, plantid, surv),
  by.x = c('Year', 'finalid'), by.y = c('Year', 'plantid')
) %>%
  # Center phen
  mutate(mean.doy = mean.doy - round(mean(mean.doy, na.rm = TRUE)))

head(mergedemo)
nrow(mergedemo)
table(mergedemo$Year)

apply(mergedemo, 2, function(x) sum(is.na(x)))
# 175 missing seeds (ugh), 31 missing phen

mergedemo.phen = mergedemo %>% filter(!is.na(mean.doy))

mergedemo.phen
table(mergedemo.phen$surv)
# okay - 58 deaths!

##### Visualizaitons

mergedemo.phen %>%
  ggplot(aes(x = mean.doy, y = surv, colour = trt)) +
  geom_point(position = position_jitter(width = 1, height = 0.125)) +
  scale_colour_manual(values = c('black', 'red', 'blue')) +
  facet_wrap(~ Year, nrow = 2)
# interesting

mergedemo.phen %>%
  filter(!is.na(no.seeds)) %>%
  ggplot(aes(x = no.seeds, y = surv, colour = trt)) +
  geom_point(position = position_jitter(width = 0.125, height = 0.125)) +
  scale_colour_manual(values = c('black', 'red', 'blue')) +
  facet_wrap(~ Year, nrow = 2)

mergedemo.phen %>%
  filter(!is.na(demo.umbels)) %>%
  ggplot(aes(x = demo.umbels, y = surv, colour = trt)) +
  geom_point(position = position_jitter(width = 0.125, height = 0.125)) +
  scale_colour_manual(values = c('black', 'red', 'blue')) +
  facet_wrap(~ Year, nrow = 2)

##### Some models

mergedemo.phen = mergedemo.phen %>% mutate(Year = factor(Year))

### Model 1: does timing of flowering influence survival probability

p_0 = glmmTMB(
  formula = surv ~ Year + (1 | Plot),
  family = 'binomial',
  data = mergedemo.phen
)

summary(p_0)
# interestingly, no effect of year...

p_t = glmmTMB(
  formula = surv ~ trt + Year + (1 | Plot),
  family = 'binomial',
  data = mergedemo.phen
)

summary(p_t)
# no treatment effect...
AIC(p_t, p_0)

p_d = glmmTMB(
  formula = surv ~ mean.doy + Year + (1 | Plot),
  family = 'binomial',
  data = mergedemo.phen
)

AIC(p_d, p_0)
# barely

p_yd = glmmTMB(
  formula = surv ~ mean.doy * Year + (1 | Plot),
  family = 'binomial',
  data = mergedemo.phen
)

AIC(p_yd, p_d, p_0) %>% mutate(daic = round(AIC - min(AIC), 2))
  
summary(p_d)

# what happens if I drop year...?

p_d_0 = glmmTMB(
  formula = surv ~ mean.doy + (1 | Plot),
  family = 'binomial',
  data = mergedemo.phen
)

AIC(p_d_0, p_d)
# huh. no year effect.

# and just to make sure there's no treatment effect...

p_d_t = glmmTMB(
  formula = surv ~ mean.doy + trt + (1 | Plot),
  family = 'binomial',
  data = mergedemo.phen
)

p_dt = glmmTMB(
  formula = surv ~ mean.doy * trt + (1 | Plot),
  family = 'binomial',
  data = mergedemo.phen
)

AIC(p_d_t, p_dt, p_d_0) # whoah! interaction effect...

summary(p_dt)
# so drought-phen interaction
# for drought, this looks like this risk is highest for late-flowering plants?

expand.grid(mean.doy = -20:20, trt = c('control', 'drought', 'irrigated')) %>%
  mutate(pred = predict(p_dt, newdata = ., type = 'response', re.form = ~ 0)) %>%
  ggplot(aes(x = mean.doy, y = pred, colour = trt)) +
  geom_line() +
  scale_colour_manual(values = c('black', 'red', 'blue')) +
  scale_y_continuous(limits = 0:1)
# interesting...

### Model 2: what about number of umbels?

mergedemo.phen %>% group_by(isna = is.na(demo.umbels)) %>% summarise(n = n())
# whoa! cool. same dataset...
# could just update our current models

p_u = glmmTMB(
  formula = surv ~ mean.doy * trt + demo.umbels + (1 | Plot),
  family = 'binomial',
  data = mergedemo.phen
)

summary(p_u)
AIC(p_u, p_dt)
# no linear effect of umbels

p_u = glmmTMB(
  formula = surv ~ mean.doy * trt + log(demo.umbels) + (1 | Plot),
  family = 'binomial',
  data = mergedemo.phen
)

summary(p_u)
AIC(p_u, p_dt)
# no effect

p_ut = glmmTMB(
  formula = surv ~ mean.doy * trt + log(demo.umbels) * trt + (1 | Plot),
  family = 'binomial',
  data = mergedemo.phen
)

AIC(p_ut, p_dt)
# ahh look at that effect. there it is.

summary(p_ut)
# whoa... that's nuts
# producing more umbels... oh LOL this is probably just a size thing.
# okay well let's look at it anyway.

expand.grid(mean.doy = -20:20, trt = c('control', 'drought', 'irrigated'), demo.umbels = c(1, 4)) %>%
  mutate(pred = predict(p_ut, newdata = ., type = 'response', re.form = ~ 0)) %>%
  ggplot(aes(x = mean.doy, y = pred, colour = trt)) +
  geom_line(aes(linetype = factor(demo.umbels))) +
  scale_colour_manual(values = c('black', 'red', 'blue')) #+
  # scale_y_continuous(limits = 0:1)
# hmm... more umbels is bad if you're in a drought plot

# Just to be sure let's look at mean phen vs. umbel count

mergedemo.phen %>%
  ggplot(aes(x = demo.umbels, y = mean.doy, colour = trt, shape = surv)) +
  geom_point(size = 3, position = position_jitter(width = 0.25, height = 1)) +
  scale_colour_manual(values = c('black', 'red', 'blue')) +
  scale_shape_manual(values = c(4, 1)) +
  facet_wrap(~ Year)

# yeah the triangle shape does suggest to me that we can't be *too* sure about
# the mean flowering time for multi-umbel plants

