library(ggplot2)
library(dplyr)
library(tidyr)
library(lme4)
library(glmmTMB)
# NOTE had to reinstall glmmTMB from source

rm(list = ls())

all.data = read.csv('01_data_cleaning/out/demo_seed_phen_combined.csv')

head(all.data)
nrow(all.data)

###### Data subsetting

# A dataset where plants with missing seed counts are excluded
seed.no = all.data %>%
  mutate(plot = factor(plot), Year = factor(Year)) %>%
  filter(!is.na(total.seeds))

# nrow(seed.no)

# table(is.na(seed.no$mean.doy)) # 10 missing initiation days
# table(is.na(seed.no$demo.umbels)) # 11 missing umbels in demo

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

seed.subs %>%
  ggplot(aes(x = mean.doy, y = n.seed.counts, colour = trt)) +
  geom_point(position = position_jitter(width = 1/14, height = 0.25)) +
  scale_colour_manual(values = c('black', 'red', 'blue')) +
  facet_wrap(trt ~ Year)
# maybe some slight association here

seed.subs %>%
  ggplot(aes(x = mean.doy, y = total.seeds / n.seed.counts)) +
  geom_point(aes(colour = trt), position = position_jitter(width = 1/7)) +
  scale_colour_manual(values = c('black', 'red', 'blue')) +
  facet_wrap(trt ~ Year)
# definitely looks like flowering starts earlier in 2021, maybe a little later
# in 2022

###### Some models

# Null model: number of seeds 
glmer(total.seeds ~ (1 | Year) + (1 | plot), family = 'poisson', data = seed.subs)
# residual devaince, as expected, is *much* larger than df
# overdispersion
# zero-inflated poisson would be good too, maybe in next batch of models

m.0 = glmer.nb(
  total.seeds ~ (1 | Year) + (1 | plot), 
  data = seed.subs,
)
# still looking a little overdispersed...
# also did not converge...

m.u = glmer.nb(
  total.seeds ~ offset(n.seed.counts) + (1 | Year) + (1 | plot), 
  data = seed.subs
)

AIC(m.u, m.0)
# it is good to include umbel counts

summary(m.u)

m.t  = glmer.nb(
  total.seeds ~ trt + offset(n.seed.counts) + (1 | Year) + (1 | plot), 
  data = seed.subs
)
m.yt = glmer.nb(
  total.seeds ~ offset(n.seed.counts) + (trt | Year) + (1 | plot), 
  data = seed.subs
)
# aw rats... singularity issue

AIC(m.yt, m.t, m.u)
# wow. zero treatment effect

m.d.u = glmer.nb(
  total.seeds ~ mean.doy + offset(n.seed.counts) + (1 | Year) + (1 | plot), 
  data = seed.subs
)

m.d2.u = glmer.nb(
  total.seeds ~ poly(mean.doy, 2) + offset(n.seed.counts) + (1 | Year) + (1 | plot), 
  data = seed.subs
)

AIC(m.d2.u, m.d.u, m.u)
# okay, a linear trend in seeds per umbel

summary(m.d.u)

# But models without adjusting for umbel find a quadratic relationships does
# work

##### Try zero-inflated poisson models?

z_0 = glmmTMB(
  total.seeds ~ (1 | plot) + (1 | Year),
  ziformula = ~ (1 | plot) + (1 | Year),
  family = 'poisson',
  data = seed.subs
)

summary(z_0)
# still looks like it's overdispersed

z_0 = glmmTMB(
  total.seeds ~ (1 | plot) + (1 | Year),
  ziformula = ~ (1 | plot) + (1 | Year),
  family = 'nbinom2',
  data = seed.subs
)

summary(z_0)
# cool?

z_u_0 = glmmTMB(
  total.seeds ~ offset(n.seed.counts) + (1 | plot) + (1 | Year),
  ziformula = ~ (1 | plot) + (1 | Year),
  family = 'nbinom2',
  data = seed.subs
)

AIC(z_u_0, z_0) %>% mutate(daic = AIC - min(AIC))
# oh... the offset is bad now?

z_d_0 = glmmTMB(
  total.seeds ~ mean.doy + (1 | plot) + (1 | Year),
  ziformula = ~ (1 | plot) + (1 | Year),
  family = 'nbinom2',
  data = seed.subs
)

z_d2_0 = glmmTMB(
  total.seeds ~ poly(mean.doy, 2) + (1 | plot) + (1 | Year),
  ziformula = ~ (1 | plot) + (1 | Year),
  family = 'nbinom2',
  data = seed.subs
)

z_0_d = glmmTMB(
  total.seeds ~ (1 | plot) + (1 | Year),
  ziformula = ~ mean.doy + (1 | plot) + (1 | Year),
  family = 'nbinom2',
  data = seed.subs
)

z_0_d2 = glmmTMB(
  total.seeds ~ (1 | plot) + (1 | Year),
  ziformula = ~ poly(mean.doy, 2) + (1 | plot) + (1 | Year),
  family = 'nbinom2',
  data = seed.subs
)

z_d_d = glmmTMB(
  total.seeds ~ mean.doy + (1 | plot) + (1 | Year),
  ziformula = ~ mean.doy + (1 | plot) + (1 | Year),
  family = 'nbinom2',
  data = seed.subs
)

z_d2_d = glmmTMB(
  total.seeds ~ poly(mean.doy, 2) + (1 | plot) + (1 | Year),
  ziformula = ~ mean.doy + (1 | plot) + (1 | Year),
  family = 'nbinom2',
  data = seed.subs
)

z_d_d2 = glmmTMB(
  total.seeds ~ mean.doy + (1 | plot) + (1 | Year),
  ziformula = ~ poly(mean.doy, 2) + (1 | plot) + (1 | Year),
  family = 'nbinom2',
  data = seed.subs
)

z_d2_d2 = glmmTMB(
  total.seeds ~ poly(mean.doy, 2) + (1 | plot) + (1 | Year),
  ziformula = ~ poly(mean.doy, 2) + (1 | plot) + (1 | Year),
  family = 'nbinom2',
  data = seed.subs
)

AIC(z_0, z_0_d, z_0_d2, z_d_0, z_d2_0, z_d_d, z_d2_d, z_d_d2, z_d2_d2) %>%
  mutate(daic = round(AIC - min(AIC), 2))  %>%
  arrange(daic)
# interesting
# probability of zero doesn't seem to vary according to time
# but there is a quadratic effect of phenology

summary(z_d2_0)

### Try treatment effects now

z_d2.t_0 = glmmTMB(
  total.seeds ~ poly(mean.doy, 2) + trt + (1 | plot) + (1 | Year),
  ziformula = ~ (1 | plot) + (1 | Year),
  family = 'nbinom2',
  data = seed.subs
)

z_d2_t = glmmTMB(
  total.seeds ~ poly(mean.doy, 2) + (1 | plot) + (1 | Year),
  ziformula = ~ trt + (1 | plot) + (1 | Year),
  family = 'nbinom2',
  data = seed.subs
)

z_d2.t_t = glmmTMB(
  total.seeds ~ poly(mean.doy, 2) + trt + (1 | plot) + (1 | Year),
  ziformula = ~ trt + (1 | plot) + (1 | Year),
  family = 'nbinom2',
  data = seed.subs
)

AIC(z_d2_0, z_d2.t_0, z_d2_t, z_d2.t_t) %>% 
  mutate(daic = round(AIC - min(AIC), 2))
# no treatment effects (again! lul)

# what happens if I add umbel offset lol

z_u.d2_0 = glmmTMB(
  total.seeds ~ poly(mean.doy, 2) + offset(n.seed.counts) + (1 | plot) + (1 | Year),
  ziformula = ~ (1 | plot) + (1 | Year),
  family = 'nbinom2',
  data = seed.subs
)

AIC(z_u.d2_0, z_d2_0) # lmao why? how does this happen

z_u.d2_0 = glmmTMB(
  total.seeds ~ poly(mean.doy, 2) + n.seed.counts + (1 | plot) + (1 | Year),
  ziformula = ~ (1 | plot) + (1 | Year),
  family = 'nbinom2',
  data = seed.subs
)

summary(z_u.d2_0) # I mean... fucking duh, right?
AIC(z_u.d2_0)


##### Look at umbel failures over time and treatment

glmer(
  formula = cbind(n.empty, n.seed.counts - n.empty) ~ (1 | plot) + (1 | Year),
  family = 'binomial',
  data = seed.subs
)
# glmer is broken
# this sucks.

f_0 = glmmTMB(
  formula = cbind(n.empty, n.seed.counts - n.empty) ~ (1 | plot) + (1 | Year),
  family = 'binomial',
  data = seed.subs
)

summary(f_0)

f_t = glmmTMB(
  formula = cbind(n.empty, n.seed.counts - n.empty) ~ trt + (1 | plot) + (1 | Year),
  family = 'binomial',
  data = seed.subs
)

AIC(f_t, f_0)
# treatment doesn't do a damn thing

f_d = glmmTMB(
  formula = cbind(n.empty, n.seed.counts - n.empty) ~ mean.doy + (1 | plot) + (1 | Year),
  family = 'binomial',
  data = seed.subs
)

f_d2 = glmmTMB(
  formula = cbind(n.empty, n.seed.counts - n.empty) ~ poly(mean.doy, 2) + (1 | plot) + (1 | Year),
  family = 'binomial',
  data = seed.subs
)

AIC(f_d2, f_d, f_0)
# okay. there *is* a relationship here

summary(f_d)
# flowering later *decreases* the odds of having empty umbels
# seems like this effect is considerable

# just to make sure treatment still doesn't belong...

f_d.t = glmmTMB(
  formula = cbind(n.empty, n.seed.counts - n.empty) ~ mean.doy + trt + (1 | plot) + (1 | Year),
  family = 'binomial',
  data = seed.subs
)

f_dt = glmmTMB(
  formula = cbind(n.empty, n.seed.counts - n.empty) ~ mean.doy * trt + (1 | plot) + (1 | Year),
  family = 'binomial',
  data = seed.subs
)

AIC(f_d, f_d.t, f_dt) # yep, really no evidence of treatment effects anywhere in here

# okie dokie.

# uh... hopefully there isn't any effect of the number of umbels, period

f_d.n = glmmTMB(
  formula = cbind(n.empty, n.seed.counts - n.empty) ~ mean.doy + n.seed.counts + (1 | plot) + (1 | Year),
  family = 'binomial',
  data = seed.subs
)

AIC(f_d.n, f_d) # hmm... okay there kinda maybe is but it's very weak
summary(f_d.n) # okay seems okay to me!

expand.grid(mean.doy = (-35:35) / 10) %>%
  mutate(pred = predict(f_d, newdata = ., re.form = ~ 0) %>% (function(x) 1/(1+exp(-x)))) %>%
  ggplot(aes(x = mean.doy, y = pred)) +
  geom_point(
    data = seed.subs %>% mutate(p.empty = n.empty / n.seed.counts),
    aes(x = mean.doy, y = p.empty, colour = trt),
    alpha = 0.25, size = 2,
    position = position_jitter(height = 0.05, width = 0.1)
  ) +
  geom_line() +
  scale_colour_manual(values = c('black', 'red', 'blue'))
# this doesn't look amazing over the data...


# What about hurdle model for the number of non-empty umbels?

seed.subs %>%
  filter(total.seeds > 0) %>%
  ggplot(aes(x = total.seeds)) +
  geom_histogram(aes()) +
  facet_wrap(trt ~ Year)
# ehh... maybe poisson?

seed.subs %>%
  filter(total.seeds > 0) %>%
  ggplot(aes(x = total.seeds / (n.seed.counts - n.empty))) +
  geom_histogram(aes()) +
  facet_wrap(trt ~ Year)
# looks maybe poisson...

#####
##### Looking at only plants with one umbel
#####

one.umbel = seed.subs %>%
  filter(demo.umbels < 2, n.phen < 2, n.seed.counts < 2)

nrow(one.umbel) # that's really not that many, but, oh well
# (of course, this is a non-random sample)
head(one.umbel)

one.umbel %>%
  ggplot(aes(x = total.seeds, fill = trt)) +
  geom_histogram(position = position_identity(), alpha = 0.5) +
  scale_colour_manual(values = c('black', 'red', 'blue')) +
  facet_wrap(trt ~ Year)
# many, many zeros.

one.umbel %>%
  ggplot(aes(x = mean.doy, y = total.seeds, colour = trt)) +
  geom_point(position = position_jitter(width = 2/7), size = 3) +
  scale_colour_manual(values = c('black', 'red', 'blue')) +
  facet_wrap( ~ Year)
# not all that strong...

one.umbel %>%
  ggplot(aes(x = diam.umbels, y = total.seeds, colour = trt)) +
  geom_point(size = 3) +
  scale_colour_manual(values = c('black', 'red', 'blue')) +
  facet_wrap( ~ Year)

one.umbel %>%
  ggplot(aes(x = mean.doy, y = diam.umbels, colour = trt)) +
  geom_point(size = 3) +
  scale_colour_manual(values = c('black', 'red', 'blue')) +
  facet_wrap( ~ Year)
# oh that's interesting...

m_0 = glmmTMB(
  total.seeds ~ (1 | plot) + (1 | Year),
  ziformula = ~ (1 | plot) + (1 | Year),
  family = 'nbinom2',
  data = one.umbel
)

m_d_0 = glmmTMB(
  total.seeds ~ mean.doy + (1 | plot) + (1 | Year),
  ziformula = ~ (1 | plot) + (1 | Year),
  family = 'nbinom2',
  data = one.umbel
)

m_0_d = glmmTMB(
  total.seeds ~ (1 | plot) + (1 | Year),
  ziformula = ~ mean.doy + (1 | plot) + (1 | Year),
  family = 'nbinom2',
  data = one.umbel
)

m_d_d = glmmTMB(
  total.seeds ~ mean.doy + (1 | plot) + (1 | Year),
  ziformula = ~ mean.doy + (1 | plot) + (1 | Year),
  family = 'nbinom2',
  data = one.umbel
)

AIC(m_0, m_d_0, m_0_d, m_d_d) %>% mutate(daic = round(AIC - min(AIC), 2))
# ooh... phen effects in both, but particularly in seeds

summary(m_d_d)
# interesting

# Treatment effects?

m_d.t_d = glmmTMB(
  total.seeds ~ mean.doy + trt + (1 | plot) + (1 | Year),
  ziformula = ~ mean.doy + (1 | plot) + (1 | Year),
  family = 'nbinom2',
  data = one.umbel
)

m_d_d.t = glmmTMB(
  total.seeds ~ mean.doy + (1 | plot) + (1 | Year),
  ziformula = ~ mean.doy + trt + (1 | plot) + (1 | Year),
  family = 'nbinom2',
  data = one.umbel
)

m_d.t_d.t = glmmTMB(
  total.seeds ~ mean.doy + trt + (1 | plot) + (1 | Year),
  ziformula = ~ mean.doy + trt + (1 | plot) + (1 | Year),
  family = 'nbinom2',
  data = one.umbel
)

AIC(m_d_d, m_d.t_d, m_d_d.t, m_d.t_d.t) %>%
  mutate(daic = round(AIC - min(AIC), 2))
# no treatment effects, once again, lmao

# oh I should look into polynomial fits too...

m_d2_d = glmmTMB(
  total.seeds ~ poly(mean.doy, 2) + (1 | plot) + (1 | Year),
  ziformula = ~ mean.doy + (1 | plot) + (1 | Year),
  family = 'nbinom2',
  data = one.umbel
)

m_d_d2 = glmmTMB(
  total.seeds ~ mean.doy + (1 | plot) + (1 | Year),
  ziformula = ~ poly(mean.doy, 2) + (1 | plot) + (1 | Year),
  family = 'nbinom2',
  data = one.umbel
)

m_d2_d2 = glmmTMB(
  total.seeds ~ poly(mean.doy, 2) + (1 | plot) + (1 | Year),
  ziformula = ~ poly(mean.doy, 2) + (1 | plot) + (1 | Year),
  family = 'nbinom2',
  data = one.umbel
)

AIC(m_d_d, m_d2_d, m_d_d2, m_d2_d2) %>%
  mutate(daic = round(AIC - min(AIC), 2))
# no evidence of quadratic fit to either of these.

summary(m_d_d)
# later flowering:
# - increases the mean seed set for non-failing umbels
# - decreases the probability of failing...
# so... it is good to flower early?

predict(m_d_d, type = 'response') # cool
predict(m_d_d, type = 'conditional')
predict(m_d_d, type = 'zprob')

one.umbel.preds = data.frame(mean.doy = (-35:35)/10) %>%
  mutate(
    pred.resp = predict(m_d_d, newdata = ., re.form = ~ 0, type = 'response'),
    pred.cond = predict(m_d_d, newdata = .,  re.form = ~ 0, type = 'conditional'),
    pred.prob = predict(m_d_d, newdata = .,  re.form = ~ 0, type = 'zprob')
  ) %>%
  pivot_longer(-mean.doy, names_to = 'varb', values_to = 'pred')

one.umbel.preds %>%
  ggplot(aes(x = mean.doy, y = pred)) +
  geom_line() +
  facet_wrap(~ varb, scales = 'free_y')
# uhhh... looks suspicitious...

one.umbel.preds %>%
  filter(!varb %in% 'pred.prob') %>%
  ggplot(aes(x = mean.doy, y = pred, group = varb)) +
  geom_point(
    data = one.umbel,
    aes(x = mean.doy, y = total.seeds, fill = !total.seeds),
    inherit.aes = FALSE,
    size = 2, alpha = 0.5, shape = 21
  ) +
  geom_line(aes(colour = varb))
# not super sold on this not being quadratic
# but at the very least it doesn't look better to flower early

# While I'm here... look at random effects
ranef(m_d_d)$cond$plot %>%
  mutate(plot = row.names(.)) %>%
  merge(y = read.csv('00_raw_data/plot_treatments.csv')) %>%
  ggplot(aes(x = `(Intercept)`, fill = trt)) +
  geom_histogram(alpha = 0.5, position = position_identity()) +
  scale_fill_manual(values = c('black', 'red', 'blue'))
# okay, definitely no treatment signal here...
ranef(m_d_d)$zi$plot %>%
  mutate(plot = row.names(.)) %>%
  merge(y = read.csv('00_raw_data/plot_treatments.csv')) %>%
  ggplot(aes(x = `(Intercept)`, fill = trt)) +
  geom_histogram(alpha = 0.5, position = position_identity()) +
  scale_fill_manual(values = c('black', 'red', 'blue'))
# hardly normally distributed... and no treatment effects

# Lesson from this: 
# from *only the one-umbel plants only*, i.e., maybe not generalizable
# but, looks like there is an advantage to flowering later: lower chance of
# reproductive failure and more seeds conditioned on not failing.
# very cool!

#####
##### Next: go back and re-export seed dataset with one row/umbel
##### and run analysis of mean date of flowering ~ seed set
#####
##### also re-run above (one-umbel) but with diamter (log?) offset
#####
