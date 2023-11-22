##### Setup/preparation

library(ggplot2)
library(dplyr)
library(tidyr)
library(lme4)

rm(list = ls())

# Inverse logit wrapper
ilogit = function(x) (1+exp(-x))^(-1)

# Load in a demography table
# demo.table = read.csv('01_data_cleaning/out/demo_table.csv')
demo = read.csv('01_data_cleaning/out/demo_postcombine.csv')

# Merge in plot information
# demo.table = merge(
  # demo.table,
demo = merge(
  demo,
  read.csv('00_raw_data/plot_treatments.csv') %>% rename(Plot = plot)
)

##### For growth, filter out relevant plants
# relevant plants: observed alive in both years, 
# have required size data
# (also - flowered in previous year?)

demo.prev = merge(
  # demo.table,
  # demo.table %>% mutate(Year = Year + 1) %>% select(-c(Plot, trt)),
  demo %>% select(plantid, Plot, trt, Year, No.leaves, Leaf.length),
  demo %>% mutate(Year = Year + 1) %>% select(plantid, Year, No.leaves, Leaf.length, No.umbels),
  by = c('plantid', 'Year'), suffixes = c('', '.prev')
) %>%
  rename(Umbels.prev = No.umbels) %>%
  arrange(plantid, Year, Plot)

head(demo.prev)

# Subset plants we want to put into model
demo.prev.for.mod = demo.prev %>%
  mutate(Year = factor(Year)) %>%
  # # Remove NAs
  filter(
    # !is.na(no.leaves) & !is.na(no.leaves.prev) &
    # !is.na(leaf.leng) & !is.na(leaf.leng.prev)
    !is.na(No.leaves)   & !is.na(No.leaves.prev) &
    !is.na(Leaf.length) & !is.na(Leaf.length.prev)
  ) %>%
  # # Remove records of plants that were not observed in the previous year
  # filter(!is.na(obs.alive.prev) & obs.alive.prev) %>%
  filter(No.leaves > 0 & No.leaves.prev > 0) %>%
  # # Previous umbel counts - update to zero if known to be alive, not flowering
  # # (remove plants that are flowering but have no umbel count, if they exist)
  # filter(!(is.na(no.umbels.prev) & flowering.prev)) %>%
  mutate(Umbels.prev = ifelse(is.na(Umbels.prev), 0, Umbels.prev))

head(demo.prev.for.mod)
table(demo.prev.for.mod$obs.alive, useNA = 'always') # all plants are alive (good)
table(demo.prev.for.mod$detected, useNA = 'always')
table(demo.prev.for.mod$detected.prev, useNA = 'always') # looks good

# Add size variables...

demo.prev.for.mod = demo.prev.for.mod %>%
  mutate(
    # size = log(no.leaves * leaf.leng),
    # size.prev = log(no.leaves.prev * leaf.leng.prev)
    size = log(No.leaves * Leaf.length),
    size.prev = log(No.leaves.prev * Leaf.length.prev),
    flowering.prev = ifelse(Umbels.prev > 0 & !is.na(Umbels.prev),
                            'flowering',
                            'basal')
  )

##### Try doing some plotting

demo.prev.for.mod %>%
  ggplot(aes(x = size)) +
  geom_density(
    aes(group = trt, fill = trt),
    alpha = 0.5, position = 'identity'
  ) +
  facet_wrap(~ Year)
# maybe some compositional effects in 2017, 2020, 2022...

demo.prev.for.mod %>%
  ggplot(aes(x = size.prev, y = size)) +
  geom_point(aes(colour = trt), alpha = 0.5) +
  coord_cartesian() +
  facet_wrap(~ Year) # interesting!

demo.prev.for.mod %>%
  ggplot(aes(x = size.prev, y = size)) +
  geom_point(aes(colour = trt), alpha = 0.5) +
  coord_cartesian() +
  facet_wrap(~ paste(flowering.prev, Year), ncol = 2, dir = 'v') 

demo.prev.for.mod %>%
  ggplot(aes(x = size.prev, y = size)) +
  geom_point(aes(colour = trt, shape = flowering.prev), alpha = 0.5) +
  coord_cartesian() +
  scale_shape_manual(values = c(19, 21)) +
  facet_wrap(~ Year)

##### Okay yeah let's do some modelz

# Null model
g.0 = lmer(
  formula = size ~ size.prev + (1 | Plot),
  data = demo.prev.for.mod
)

summary(g.0) # intercept is very high...?
hist(summary(g.0)$residuals) # kinda have a tail but otherwise okay (but I forget if residuals are ranef-centered or not...)
hist(unlist(ranef(g.0)$Plot))
# hist(unlist(ranef(g.0)$`plantid:Plot`)) # these are normal tho

### Model with year
g.y = lmer(
  formula = size ~ size.prev + (1 | Year) + (1 | Plot),
  data = demo.prev.for.mod
)

summary(g.y)
AIC(g.y, g.0) # yep. year effect is good

### Model with flowering in previous year

g.y.f = lmer(
  formula = size ~ size.prev + flowering.prev + (1 | Year) + (1 | Plot),
  data = demo.prev.for.mod
)

AIC(g.y.f, g.y)
summary(g.y.f) # wait... this effect is positive?

### Flowering in previous year * size interaction

g.y.sf = lmer(
  formula = size ~ size.prev * flowering.prev + (1 | Year) + (1 | Plot),
  data = demo.prev.for.mod
)

summary(g.y.sf) # cool... assuming it's not bunk, this looks like it would make more sense...
AIC(g.y.sf, g.y.f, g.y) # but it is bunk. interesting

### Size-treatment interaction (just to make sure)
g.y.f.t = lmer(
  formula = size ~ size.prev + trt + flowering.prev + (1 | Year) + (1 | Plot),
  data = demo.prev.for.mod
)

AIC(g.y.f.t, g.y.f, g.y)   # no treatment effects here

g.y.f.st = lmer(
  formula = size ~ size.prev * trt + flowering.prev + (1 | Year) + (1 | Plot),
  data = demo.prev.for.mod
)

AIC(g.y.f.st, g.y.f) # no size-treatment effects

### Flowering-treatment interaction
g.y.ft = lmer(
  formula = size ~ size.prev + trt * flowering.prev + (1 | Year) + (1 | Plot),
  data = demo.prev.for.mod
)

AIC(g.y.ft, g.y.f) # no treatment-flowering effects

### Oh I guess I should see if slope varies by year...

g.sy.f = lmer(
  formula = size ~ flowering.prev + (size.prev | Year) + (1 | Plot),
  data = demo.prev.for.mod
)

AIC(g.sy.f, g.y.f) # nope!

### What about effects of prior-year flowering varying by year?

g.yf = lmer(
  formula = size ~ size.prev + (flowering.prev | Year) + (1 | Plot),
  data = demo.prev.for.mod
)

AIC(g.yf, g.y.f) # ooh... super-interesting!! 

summary(g.yf)
hist(unlist(ranef(g.yf)$Plot))
hist(unlist(ranef(g.yf)$Year))
plot(unlist(ranef(g.yf)$Year["(Intercept)"])) # looks a little auto-correlated...
hist(unlist(ranef(g.yf)$Year["flowering.prevflowering"]))
plot(unlist(ranef(g.yf)$Year["flowering.prevflowering"])) # hmm...
# hist(unlist(ranef(g.yf)$`plantid:Plot`)) # big-ol' right tail...

g.yf.s2 = lmer(
  formula = size ~ poly(size.prev, 2) + (flowering.prev | Year) + (1 | Plot),
  data = demo.prev.for.mod
)

AIC(g.yf.s2, g.yf)
# WHAT
# a polynomial effect...?

summary(g.yf.s2)
# what the fuck...

### So far... would appear that the best model is g.yf (flowering effect varying by year, independent of size)
# shall we plot?

data.preds = expand.grid(
  size.prev = seq(.5, 6, by = 0.25),
  Year = factor(2017:2023),
  flowering.prev = c(TRUE, FALSE)
) %>%
  mutate(pred.size = predict(g.yf.s2, newdata = ., re.form = ~ (flowering.prev | Year)))

head(data.preds)  

ggplot() +
  geom_point(
    data = demo.prev.for.mod,
    aes(x = size.prev, y = size, colour = trt, shape = flowering.prev),
    alpha = 0.5
  ) +
  geom_line(
    data = data.preds,
    aes(x = size.prev, y = pred.size, group = flowering.prev, linetype = flowering.prev)
  ) +
  scale_shape_manual(values = c(19, 21)) +
  scale_colour_manual(values = c('black', 'red', 'blue')) +
  facet_wrap(~ Year)

ggplot() +
  geom_point(
    data = demo.prev.for.mod,
    aes(x = size.prev, y = size, colour = trt, shape = flowering.prev),
    alpha = 0.5
  ) +
  geom_line(
    data = data.preds,
    aes(x = size.prev, y = pred.size, group = flowering.prev, linetype = flowering.prev)
  ) +
  scale_shape_manual(values = c(19, 21)) +
  scale_colour_manual(values = c('goldenrod', 'red', 'blue')) +
  facet_wrap(~ paste(Year, flowering.prev), ncol = 4)

# How often is shrinkage happening...?
demo.prev.for.mod %>%
  group_by(Year, trt, flowering.prev = ifelse(flowering.prev, 'fl', 'nfl')) %>%
  summarise(p.grow = mean(size > size.prev)) %>% 
  pivot_wider(names_from = flowering.prev, values_from = p.grow)
  
# whoa... what the fuck

# Okay... look at size components...

##### Should also maybe fit a model with current flowering year...

g.yf0.yf1 = lmer(
  formula = size ~ size.prev + (flowering + flowering.prev | Year) + (1 | Plot / plantid),
  data = demo.prev.for.mod
)

g.yf0.f1 = lmer(
  formula = size ~ size.prev + flowering + (flowering.prev | Year) + (1 | Plot / plantid),
  data = demo.prev.for.mod
)

AIC(g.yf0.yf1, g.yf0.f1, g.yf) # cool

summary(g.yf0.f1)

data.preds = expand.grid(
  size.prev = seq(.5, 6, by = 0.25),
  Year = factor(2017:2023),
  flowering = c(TRUE, FALSE),
  flowering.prev = c(TRUE, FALSE)
) %>%
  mutate(pred.size = predict(g.yf0.f1, newdata = ., re.form = ~ (flowering.prev | Year)))

ggplot() +
  geom_point(
    data = demo.prev.for.mod,
    aes(x = size.prev, y = size, colour = trt, shape = flowering.prev),
    alpha = 0.5
  ) +
  geom_line(
    data = data.preds,
    aes(x = size.prev, y = pred.size, group = interaction(flowering.prev, flowering), linetype = interaction(flowering.prev, flowering))
  ) +
  scale_shape_manual(values = c(19, 21)) +
  scale_colour_manual(values = c('goldenrod', 'red', 'blue')) +
  facet_wrap(~ Year)
