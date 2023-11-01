### Script for playing around with models of probability of flowering
### SN - 31 Oct 2023

### Setup/preparation

library(ggplot2)
library(dplyr)
library(tidyr)

rm(list = ls())

# Inverse logit wrapper
ilogit = function(x) (1+exp(-x))^(-1)

# Load in a demography table
demo.table = read.csv('01_data_cleaning/out/demo_table.csv')

# Merge in plot information
demo.table = merge(
  demo.table,
  read.csv('00_raw_data/plot_treatments.csv') %>% rename(Plot = plot)
)

### Start exploring

demo.table %>%
  filter(obs.alive) %>%
  ggplot(aes(x = no.leaves, y = flowering)) +
  geom_point(position = position_jitter(width = 0.25, height = 0.125), alpha = 0.125) +
  facet_wrap(~ Year)
# hinteresting...

demo.table %>% filter(obs.alive) %>% group_by(flowering) %>% summarise(n = n()) # should be super easy
demo.table %>% filter(obs.alive) %>% group_by(Plot) %>% summarise(p.fl = mean(flowering))

# Fit some models

p.flower.0 = glmer(
  formula = flowering ~ (1 | Plot / plantid),
  data = demo.table %>% filter(obs.alive),
  family = 'binomial'
)

summary(p.flower.0)
# slight overdispersion? not super concerned yet

p.flower.y = glmer(
  formula = flowering ~ (1 | Year) + (1 | Plot / plantid),
  data = demo.table %>% filter(obs.alive) %>% mutate(Year = factor(Year)),
  family = 'binomial'
)

summary(p.flower.y)
# adding a year-level random effect increased the variance in all of the standard effects...

AIC(p.flower.0, p.flower.y) # hmm... leave year out for now...

p.flower.t = glmer(
  formula = flowering ~ trt + (1 | Plot / plantid),
  data = demo.table %>% filter(obs.alive),
  family = 'binomial'
)

summary(p.flower.t) # effects here look pretty small...

AIC(p.flower.t, p.flower.0) # lmao damn wtf

p.flower.ll = glmer(
  formula = flowering ~ log(no.leaves) + (1 | Plot / plantid),
  data = demo.table %>% filter(obs.alive) %>% mutate(Year = factor(Year)),
  family = 'binomial'
)

p.flower.l = glmer(
  formula = flowering ~ no.leaves + (1 | Plot / plantid),
  data = demo.table %>% filter(obs.alive) %>% mutate(Year = factor(Year)),
  family = 'binomial'
)

AIC(p.flower.ll, p.flower.l, p.flower.0) # damn... log is the way to go

p.flower.l = glmer(
  formula = flowering ~ log(no.leaves) + (1 | Plot / plantid),
  data = demo.table %>% filter(obs.alive) %>% mutate(Year = factor(Year)),
  family = 'binomial'
)

summary(p.flower.l) # nice.

# Now try adding in year and trt?

p.flower.l.y = glmer(
  formula = flowering ~ log(no.leaves) + (1 | Year) + (1 | Plot / plantid),
  data = demo.table %>% filter(obs.alive) %>% mutate(Year = factor(Year)),
  family = 'binomial'
)

AIC(p.flower.l.y, p.flower.l) # sweet

p.flower.l.y.t = glmer(
  formula = flowering ~ log(no.leaves) + trt + (1 | Year) + (1 | Plot / plantid),
  data = demo.table %>% filter(obs.alive) %>% mutate(Year = factor(Year)),
  family = 'binomial'
)

AIC(p.flower.l.y.t, p.flower.l.y) # huh... no influence of treatment (maybe it's rolled up into leaf count?)

summary(p.flower.l.y)
hist(unlist(ranef(p.flower.l.y)$Year)) # plausible
plot(unlist(ranef(p.flower.l.y)$Year)) # upward trend...
hist(unlist(ranef(p.flower.l.y)$Plot)) # also plausible, kinda skewed tho

# leaf size too...

p.flower.n.l = glmer(
  formula = flowering ~ log(no.leaves) + log(leaf.leng) + (1 | Plot / plantid),
  data = demo.table %>% filter(obs.alive) %>% mutate(Year = factor(Year)),
  family = 'binomial'
)

AIC(p.flower.n.l, p.flower.l) # aw shit... got some plants with missing leaf lengths
summary(p.flower.n.l)

p.flower.n = glmer(
  formula = flowering ~ log(no.leaves) + (1 | Plot / plantid),
  data = demo.table %>% filter(obs.alive, !is.na(leaf.leng)) %>% mutate(Year = factor(Year)),
  family = 'binomial'
)

p.flower.n.l = glmer(
  formula = flowering ~ log(no.leaves) + log(leaf.leng) + (1 | Plot / plantid),
  data = demo.table %>% filter(obs.alive) %>% mutate(Year = factor(Year)),
  family = 'binomial'
)

AIC(p.flower.n.l, p.flower.n) # some leaf leng too

p.flower.n.l.y = glmer(
  formula = flowering ~ log(no.leaves) + log(leaf.leng) + (1 | Year) + (1 | Plot / plantid),
  data = demo.table %>% filter(obs.alive) %>% mutate(Year = factor(Year)),
  family = 'binomial'
)

AIC(p.flower.n.l.y, p.flower.n.l) # yep - definitely want year in there...

# test for treatment

p.flower.n.l.y.t = glmer(
  formula = flowering ~ log(no.leaves) + log(leaf.leng) + trt + (1 | Year) + (1 | Plot / plantid),
  data = demo.table %>% filter(obs.alive) %>% mutate(Year = factor(Year)),
  family = 'binomial'
)

AIC(p.flower.n.l.y.t, p.flower.n.l.y) # no effect of treatment

summary(p.flower.n.l.y) # wait lmao what? longer longest leaf decreases probability of flowering...?
hist(unlist(ranef(p.flower.n.l.y)$Year)) # looks great!
plot(unlist(ranef(p.flower.n.l.y)$Year)) # upward trend still
hist(unlist(ranef(p.flower.n.l.y)$Plot)) # still skewed..



##### Okay... now let's look at previous-year's effects!
# (oh duh - actually should be looking at prev-year's counts for leaf metrics...)

head(demo.table)

# Use this to check and make sure I did this correctly
demo.table %>% filter(plantid %in% '3658_13_11I') %>% arrange(Year)

demo.prev = merge(
  demo.table,
  demo.table %>% mutate(Year = Year + 1) %>% select(-c(Plot, trt)),
  by = c('plantid', 'Year'), suffixes = c('', '.prev')
)

demo.prev %>% filter(plantid %in% '3658_13_11I') %>% arrange(Year)

head(demo.prev)
nrow(demo.prev) / nrow(demo.table) # still have ~80% of records

# Try a plot of leaf counts...

demo.prev %>%
  filter(!is.na(no.leaves) & !is.na(no.leaves.prev)) %>%
  group_by(Year, no.leaves, no.leaves.prev) %>%
  summarise(n = n()) %>%
  ggplot(aes(x = no.leaves.prev, y = no.leaves)) +
  geom_raster(aes(fill = n)) +
  facet_wrap(~ Year)
# lol not very insightful

demo.prev %>%
  filter(!is.na(no.leaves) & !is.na(no.leaves.prev)) %>%
  group_by(Year, no.leaves, no.leaves.prev) %>%
  summarise(n = n()) %>%
  group_by(Year, no.leaves.prev) %>%
  mutate(p = n / sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x = no.leaves.prev, y = no.leaves)) +
  geom_raster(aes(fill = p)) +
  facet_wrap(~ Year) +
  scale_fill_viridis_c()
# yeah... not super insightful

demo.prev.for.mod = demo.prev %>%
  mutate(Year = factor(Year)) %>%
  # Remove NAs
  filter(
    !is.na(no.leaves) & !is.na(no.leaves.prev) &
    !is.na(leaf.leng) & !is.na(leaf.leng.prev)
  ) %>%
  # Remove records of plants that were not observed in the previous year
  filter(!is.na(obs.alive.prev) & obs.alive.prev) %>%
  # Previous umbel counts - update to zero if known to be alive, not flowering
  # (remove plants that are flowering but have no umbel count, if they exist)
  filter(!(is.na(no.umbels.prev) & flowering.prev)) %>%
  mutate(no.umbels.prev = ifelse(is.na(no.umbels.prev), 0, no.umbels.prev))

nrow(demo.prev.for.mod)
head(demo.prev.for.mod)

sum(is.na(demo.prev.for.mod$flowering))
sum(is.na(demo.prev.for.mod$flowering.prev))

# Guess for now we just want flowering (because survival will require occupancy
# models)

# Terms to test:
# - Year effects
# - Treatment effect
# - Leaf number, previous year
# - Leaf length, previous year
# - Flowering, previous year

##### Start fitting model

# Nullest of models
p.fl.0 = glmer(
  formula = flowering ~ (1 | Plot / plantid),
  data = demo.prev.for.mod,
  family = 'binomial'
)

summary(p.fl.0)

# Include previous year's leaf count
p.fl.n = glmer(
  formula = flowering ~ log(no.leaves.prev) + (1 | Plot / plantid),
  data = demo.prev.for.mod,
  family = 'binomial'
)

AIC(p.fl.n, p.fl.0) # as expected, log of leaf counts should be included

# Now look for longest leaf count
p.fl.n.l = glmer(
  formula = flowering ~ log(no.leaves.prev) + log(leaf.leng.prev) + (1 | Plot / plantid),
  data = demo.prev.for.mod,
  family = 'binomial'
)

AIC(p.fl.n.l, p.fl.n) # yes! makes a big difference!

# Include effect of year
p.fl.n.l.y = glmer(
  formula = flowering ~ log(no.leaves.prev) + log(leaf.leng.prev) + 
    (1 | Year) + (1 | Plot / plantid),
  data = demo.prev.for.mod,
  family = 'binomial'
)

AIC(p.fl.n.l.y, p.fl.n.l) # yes, include year!

# Flowering in previous year
p.fl.n.l.y.f = glmer(
  formula = flowering ~ log(no.leaves.prev) + log(leaf.leng.prev) + flowering.prev +
    (1 | Year) + (1 | Plot / plantid),
  data = demo.prev.for.mod,
  family = 'binomial'
)

p.fl.n.l.y.u = glmer(
  formula = flowering ~ log(no.leaves.prev) + log(leaf.leng.prev) + log(no.umbels.prev + 1) +
    (1 | Year) + (1 | Plot / plantid),
  data = demo.prev.for.mod,
  family = 'binomial'
)

AIC(p.fl.n.l.y.f, p.fl.n.l.y.u, p.fl.n.l.y) 
# yes, include previous year's flowering
# (this is better than umbel counts... which is good for constructing a kernel!)

# Oops haven't done treatment effects yet!

p.fl.n.l.y.f.t = glmer(
  formula = flowering ~ log(no.leaves.prev) + log(leaf.leng.prev) + flowering.prev + trt +
    (1 | Year) + (1 | Plot / plantid),
  data = demo.prev.for.mod,
  family = 'binomial'
)

AIC(p.fl.n.l.y.f.t, p.fl.n.l.y.f) # no discernible treatment effects!

# Look for an interaction between leaves and flowering?

p.fl.l.y.nf = glmer(
  formula = flowering ~ log(no.leaves.prev) * flowering.prev + log(leaf.leng.prev) +
    (1 | Year) + (1 | Plot / plantid),
  data = demo.prev.for.mod,
  family = 'binomial'
)

AIC(p.fl.l.y.nf, p.fl.n.l.y.f) # no discernible flowering x no.leaves interaction

p.fl.n.y.lf = glmer(
  formula = flowering ~ log(no.leaves.prev) + flowering.prev * log(leaf.leng.prev) +
    (1 | Year) + (1 | Plot / plantid),
  data = demo.prev.for.mod,
  family = 'binomial'
)

AIC(p.fl.n.y.lf, p.fl.n.l.y.f) # no discernible flowering x longest-leaf interaction

# Final model would appear to be

p.fl.n.l.y.f = glmer(
  formula = flowering ~ log(no.leaves.prev) + log(leaf.leng.prev) + flowering.prev +
    (1 | Year) + (1 | Plot / plantid),
  data = demo.prev.for.mod,
  family = 'binomial'
)

summary(p.fl.n.l.y.f)
# Positive effects of both previous number of leaves, previous year's leaf length, and previously flowering

# Examine andom effects
# Year effects
hist(unlist(ranef(p.fl.n.l.y.f)$Year)) # oof... does not look normally distributed
plot(unlist(ranef(p.fl.n.l.y.f)$Year)) # if there is something funny happening in 2016, it's not obvious
hist(unlist(ranef(p.fl.n.l.y.f)$Plot)) # these look more normal
hist(unlist(ranef(p.fl.n.l.y.f)$`plantid:Plot`)) # these look normal-ish

# Let's take a look:

demo.prev.for.mod = demo.prev.for.mod #%>%
  # mutate(pred.fl = predict(p.fl.n.l.y.f, newdata = ., re.form = NA) %>% ilogit())

data.scaffold = expand.grid(
  flowering.prev = c(TRUE, FALSE),
  trt = unique(demo.prev.for.mod$trt),
  leaf.leng.prev = 1:30,
  no.leaves.prev = 1:24,
  Year = factor(2017:2023)
)

data.scaffold$pred = predict(
  p.fl.n.l.y.f,
  newdata = data.scaffold,
  re.form = ~ (1 | Year)
) %>% ilogit()

# 
demo.prev.for.mod %>%
  ggplot(aes(x = leaf.leng.prev)) +
  geom_point(
    aes(y = as.numeric(flowering), colour = trt),
    position = position_jitter(width = 0.25, height = 0.125),
    alpha = 0.5
  ) +
  geom_line(
    data = data.scaffold %>% filter(no.leaves.prev %in% c(1, 5)),
    aes(
      y = pred, 
      linetype = flowering.prev,
      group = interaction(flowering.prev, no.leaves.prev)
      ),
    colour = 'gray66'
  ) +
  geom_point(
    data = data.scaffold %>% filter(no.leaves.prev %in% c(1, 5)),
    aes(
      y = pred, 
      shape = factor(no.leaves.prev),
      group = interaction(flowering.prev, no.leaves.prev)
    ),
    colour = 'gray66', size = 1
  ) +
  labs(x = 'Previous year leaf length', y = 'Probability of flowering') +
  scale_shape_discrete('Number of leaves previous year') +
  scale_linetype_discrete('Flowered in previous year') +
  scale_colour_discrete('Treatment') +
  facet_wrap(~ Year) +
  theme(
    legend.position = c(0.8, 0.15),
    legend.direction = 'horizontal',
    panel.background = element_blank()
  )

ggsave('02_data_exploration/figs/flowering_lfln.png', height = 8, width = 8)

demo.prev.for.mod %>%
  ggplot(aes(x = no.leaves.prev)) +
  geom_point(
    aes(y = as.numeric(flowering), colour = trt),
    position = position_jitter(width = 0.25, height = 0.125),
    alpha = 0.5
  ) +
  geom_line(
    data = data.scaffold %>% filter(leaf.leng.prev %in% c(1, 11, 21)),
    aes(
      y = pred, 
      linetype = flowering.prev,
      group = interaction(flowering.prev, leaf.leng.prev)
    ),
    colour = 'gray66'
  ) +
  geom_point(
    data = data.scaffold %>% filter(leaf.leng.prev %in% c(1, 11, 21)),
    aes(
      y = pred, 
      shape = factor(leaf.leng.prev),
      group = interaction(flowering.prev, leaf.leng.prev)
    ),
    colour = 'gray66', size = 1
  ) +
  scale_shape_discrete('Leaf length in previous year') +
  scale_linetype_discrete('Flowered in previous year') +
  scale_colour_discrete('Treatment') +
  labs(x = 'Previous year numbef of leaves', y = 'Probability of flowering') +
  facet_wrap(~ Year) +
  theme(
    legend.position = c(0.8, 0.15),
    legend.direction = 'horizontal',
    panel.background = element_blank()
  )

ggsave('02_data_exploration/figs/flowering_lfct.png', height = 8, width = 8)
