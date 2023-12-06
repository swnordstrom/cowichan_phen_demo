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
  formula = flowering ~ log(no.leaves.prev) + log(leaf.leng.prev) + flowering.prev +
    (trt | Year) + (1 | Plot / plantid),
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

# ggsave('02_data_exploration/figs/flowering_lfln.png', height = 8, width = 8)

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

# ggsave('02_data_exploration/figs/flowering_lfct.png', height = 8, width = 8)

##### Model with one "size" component...

demo.prev.for.mod = demo.prev.for.mod %>%
  mutate(
    size = log(no.leaves * leaf.leng),
    size.prev = log(no.leaves.prev * leaf.leng.prev)
  ) %>%
  # Because it was suggested by Jenn, remove 2017 data (2016 leaf sizes are unreliable)
  filter(!Year %in% 2017)

### Null model as before

p.fl.0 = glmer(
  formula = flowering ~ (1 | Plot / plantid),
  family = 'binomial',
  data = demo.prev.for.mod
)

### Model with current size

p.fl.s1 = glmer(
  formula = flowering ~ size + (1 | Plot / plantid),
  family = 'binomial',
  data = demo.prev.for.mod
)

### Model with prior size

p.fl.s0 = glmer(
  formula = flowering ~ size.prev + (1 | Plot / plantid),
  family = 'binomial',
  data = demo.prev.for.mod
)

### Model with both prior and current size

p.fl.s0.s1 = glmer(
  formula = flowering ~ size + size.prev + (1 | Plot / plantid),
  family = 'binomial',
  data = demo.prev.for.mod
)

### Compare with AIC
AIC(p.fl.0, p.fl.s0, p.fl.s1, p.fl.s0.s1) # both of them! what the heck

summary(p.fl.s0.s1)
# both are positive (greater effect for current size)

### Look for year effects
p.fl.s0.s1.y = glmer(
  formula = flowering ~ size + size.prev + (1 | Year) + (1 | Plot / plantid),
  family = 'binomial',
  data = demo.prev.for.mod
)

AIC(p.fl.s0.s1.y, p.fl.s0.s1) # year effects are good

### Well... do any of the size variables vary by year?

p.fl.ys0.ys1 = glmer(
  formula = flowering ~  (size + size.prev | Year) + (1 | Plot / plantid),
  family = 'binomial',
  data = demo.prev.for.mod
)
# convergence error...

p.fl.ys0.s1 = glmer(
  formula = flowering ~  size + (size.prev | Year) + (1 | Plot / plantid),
  family = 'binomial',
  data = demo.prev.for.mod
)

p.fl.s0.ys1 = glmer(
  formula = flowering ~  size.prev + (size | Year) + (1 | Plot / plantid),
  family = 'binomial',
  data = demo.prev.for.mod
)

AIC(p.fl.s0.ys1, p.fl.ys0.s1, p.fl.s0.s1.y) %>% arrange(AIC)
# no evidence that size effects vary by year

### Look at prior flowering

p.fl.s0.s1.y.f0 = glmer(
  formula = flowering ~ size + size.prev + flowering.prev + (1 | Year) + (1 | Plot / plantid),
  family = 'binomial',
  data = demo.prev.for.mod
)

AIC(p.fl.s0.s1.y.f0, p.fl.s0.s1.y) # yep.
summary(p.fl.s0.s1.y.f0)
# positive effect...

### Prior-flowering varying by year

p.fl.s0.s1.yf0 = glmer(
  formula = flowering ~ size + size.prev + (flowering.prev | Year) + (1 | Plot / plantid),
  family = 'binomial',
  data = demo.prev.for.mod
)

AIC(p.fl.s0.s1.yf0, p.fl.s0.s1.y.f0) %>% arrange(AIC) # nope - fl-effect does not vary by year

### Checking for a treatment effect

p.fl.s0.s1.y.f0.t = glmer(
  formula = flowering ~ size + size.prev + flowering.prev + trt + (1 | Year) + (1 | Plot / plantid),
  family = 'binomial',
  data = demo.prev.for.mod
)

AIC(p.fl.s0.s1.y.f0.t, p.fl.s0.s1.y.f0) %>% arrange(AIC) # nope - no evidence of treatment effect...

### But could treatment have a year-varying effect?

p.fl.s0.s1.f0.yt = glmer(
  formula = flowering ~ size + size.prev + flowering.prev + (trt | Year) + (1 | Plot / plantid),
  family = 'binomial',
  data = demo.prev.for.mod
)

# did not converge

AIC(p.fl.s0.s1.f0.yt, p.fl.s0.s1.y.f0) # it looks slightly better though... oh well.

### Let's plot

summary(p.fl.s0.s1.y.f0)

fl.preds = expand.grid(
  Year = factor(2018:2023),
  size = seq(.5, 6, by = .25),
  size.prev = seq(.5, 6, by = .25),
  flowering.prev = c(TRUE, FALSE)
) %>%
  mutate(fl.pred = ilogit(predict(p.fl.s0.s1.y.f0, newdata = ., re.form = ~ (1 | Year))))

head(fl.preds)

ggplot() +
  geom_tile(
    data = fl.preds %>%
      mutate(flowering.prev = ifelse(flowering.prev, 'Flowered last year', 'Did not flower last year')),
    aes(x = size.prev, y = size, fill = fl.pred)
  ) +
  geom_point(
    data = demo.prev.for.mod %>% 
      group_by(Year, flowering.prev) %>% 
      sample_n(50) %>%
      mutate(
        flowering = ifelse(flowering, 'Flowered', 'Did not flower'),
        flowering.prev = ifelse(flowering.prev, 'Flowered last year', 'Did not flower last year')
      ),
    aes(x = size.prev, y = size, shape = flowering),
    colour = 'black', size = 2
  ) +
  scale_shape_manual(values = c(21, 19), '') +
  scale_fill_stepsn(colours = RColorBrewer::brewer.pal(4, 'Greys'), 'Probability of flowering') +
  facet_wrap(~ paste(Year, flowering.prev, sep = ', ')) +
  labs(x = 'Size in previous year', y = 'Size in current year') +
  theme(
    legend.position = 'bottom',
    panel.background = element_blank()
  )

###### ####################################
###### ####################################
###### ####################################
###### ####################################
###### ####################################
###### ####################################
###### ####################################
###### ####################################
###### ####################################

library(ggplot2)
library(dplyr)
library(tidyr)
library(lme4)

rm(list = ls())

# Inverse logit wrapper
ilogit = function(x) (1+exp(-x))^(-1)

# Read in demography data
demo = read.csv('01_data_cleaning/out/demo_imputed_survival.csv')

head(demo)
nrow(demo)

# Add a "flowering" column
# Check to make sure this works
demo %>%
  mutate(flowering = case_when(
    No.umbels > 0 ~ TRUE,
    !No.umbels ~ FALSE,
    is.na(No.umbels) ~ FALSE,
    .default = NA
  )
) %>%
  distinct(No.umbels, flowering, surv)
# Looks good.

# Add column to data frame
demo = demo %>%
  mutate(
    flowering = case_when(
      No.umbels > 0 ~ TRUE,
      !No.umbels ~ FALSE,
      is.na(No.umbels) ~ FALSE,
      .default = NA
    )
  )

# Attach this year's flowering to prev year's demographic data
demo.flower = merge(
  # x is year of flowering
  # here: need *surviving* plants only
  # also want flowering column
  x = demo %>%
    filter(surv) %>%
    select(Year, plantid, Plot, trt, flowering, demo.note, proc.note, edited),
  # y is previous year's data
  # do this by adding 1 to year (will allow matching)
  # give *only* the plants that have size measurements
  # NOTE: this will cut out some plants from analysis if we don't have their
  # size in the previous year
  y = demo %>%
    mutate(Year = Year + 1) %>%
    filter(!is.na(No.leaves) & !is.na(Leaf.length) & No.leaves > 0) %>%
    mutate(size.prev = log(No.leaves * Leaf.length)) %>%
    select(Year, plantid, size.prev, flowering, demo.note, proc.note, edited),
  all = FALSE, by = c("Year", "plantid"), suffixes = c("", ".prev")
) %>%
  mutate(Plot = factor(Plot))

# Look at this subcase to make sure we got this right
demo %>% filter(grepl('1411', plantid))
demo.flower %>% filter(grepl('1411', plantid))

head(demo.flower)
nrow(demo.flower)

# Read in climate data
clim = read.csv('01_data_cleaning/out/climate_summary.csv')

head(clim)
str(clim)

# Fix last freeze date (convert to numeric)
clim = clim %>%
  mutate(last.freeze = as.numeric(as.Date(last.freeze)) - as.numeric(as.Date('2019-12-31'))) %>%
  # let's standardize some of these
  mutate(across(-Year, ~ (. - mean(., na.rm = TRUE)) / sd(., na.rm = TRUE)))

clim.lag = merge(
  clim,
  clim %>% mutate(Year = Year + 1),
  all.x = TRUE, all.y = TRUE, by = 'Year', suffixes = c('', '.prev')
)

head(clim.lag)
tail(clim.lag)

clim.lag %>% filter(Year %in% 2023)

# Merge these together
# demo.fl.clim = merge(x = demo.flower, y = clim.lag, by = 'Year')

# ACTUALLY for simplicity let's not include lags yet...
demo.fl.clim = merge(x = demo.flower, y = clim, by = 'Year')

nrow(demo.fl.clim)
# Good

demo.fl.clim %>% filter(Year %in% 2022) %>% head()
# Looks right to me, mostly

table(demo.fl.clim$flowering) # wow... it's 50-50
table(demo.fl.clim$flowering.prev)
with(demo.fl.clim, table(flowering, flowering.prev))
# wow... that's a lot on the diagonal
# lots of plants that are repeat-flowerers
table(demo.fl.clim$Year)
# remember - 2016 size data is suspect
# in fact... maybe just take it out? has potential to influence year effects

demo.fl.clim = demo.fl.clim %>% filter(Year > 2017)

#####
##### Before fitting models, do plots
#####

ggplot(demo.fl.clim, aes(x = size.prev, y = flowering)) +
  geom_point(position = position_jitter(height = 0.1), alpha = 0.5) +
  facet_wrap(~ Year)
# certainly looks like it varies by year...

#####
##### Try fitting models
#####

### Null model
f_0 = glmer(
  formula = flowering ~ (1 | Plot / plantid),
  data = demo.fl.clim,
  family = 'binomial'
)

summary(f_0)

### Size
f_s = glmer(
  formula = flowering ~ size.prev + (1 | Plot / plantid),
  data = demo.fl.clim,
  family = 'binomial'
) 

summary(f_s)
AIC(f_0, f_s)

### Size + flowering
f_s.f = glmer(
  formula = flowering ~ size.prev + flowering.prev + (1 | Plot / plantid),
  data = demo.fl.clim,
  family = 'binomial'
) 
# argh.... singularity
# should I remove the plant-level random effect?

f_s.f = glmer(
  formula = flowering ~ size.prev + flowering.prev + (1 | Plot),
  data = demo.fl.clim,
  family = 'binomial'
) 
# yeah that does it
# argh... guess we're flying without an individual-level random effect!

summary(f_s.f)
AIC(f_0, f_s, f_s.f)

### Size + flowering + year
f_s.f.y = glmer(
  formula = flowering ~ size.prev + flowering.prev + (1 | Year) + (1 | Plot),
  data = demo.fl.clim,
  family = 'binomial'
) 

summary(f_s.f.y)
AIC(f_s.f.y, f_s.f)
# Oh yes - keep year effects

### Treatment effects?
f_s.f.y.t = glmer(
  formula = flowering ~ size.prev + flowering.prev + trt + (1 | Year) + (1 | Plot),
  data = demo.fl.clim,
  family = 'binomial'
) 

summary(f_s.f.y.t) # haha
AIC(f_s.f.y, f_s.f.y.t) # lmao

### Now... I have so many climate variables...
# let's just try a couple

f_s.f.freeze = glmer(
  formula = flowering ~ size.prev + flowering.prev + last.freeze + (1 | Plot),
  data = demo.fl.clim,
  family = 'binomial'
)

summary(f_s.f.freeze) # ooh that looks significant
# (but in the opposite way I'd expect... freezing later means increased flowering?)

AIC(f_s.f.freeze, f_s.f.y) # whoa but the AIC is HORRIBLE

f_s.f.summertemp = glmer(
  formula = flowering ~ size.prev + flowering.prev + summer_meanTemp + (1 | Plot),
  data = demo.fl.clim,
  family = 'binomial'
)

summary(f_s.f.summertemp) # also looks positive... warmer previous summer means more flowering?
AIC(f_s.f.summertemp, f_s.f.y) # AIC here is better, but still...
# (I mean... will a single term *ever* do as well as a year random effect?)

f_s.f.earlyPrec = glmer(
  formula = flowering ~ size.prev + flowering.prev + early_prec.sum + (1 | Plot),
  data = demo.fl.clim,
  family = 'binomial'
)

summary(f_s.f.earlyPrec) # once again... positive association
# (wetter early season means... more flowering)
# (why are all of these positive...)
AIC(f_s.f.earlyPrec, f_s.f.y) # still not great

# Okay... if we try fitting all combos here, I'm not convinced that any will
# do nearly as well with AIC
# how would we know which of the AICs was the worth keeping?
#
# although it looks like glmnet doesn't work with random effects...
# (rats)
# but... glmmLasso does!

library(glmmLasso)

l_s.f = glmmLasso(
  fix = flowering ~ size.prev + as.factor(flowering.prev) +
    early_prec.sum + early_minnTemp + last.freeze + 
    grow_meanTemp + summer_meanTemp + winter_meanTemp,
  rnd = list(Plot = ~ 1),
  data = demo.fl.clim,
  family = binomial(),
  lambda = 50,
  control = list(
    index = c(NA, NA, 1, 2, 3, 4, 5, 6),
    print.iter = TRUE
  )
)
# that's a lot...

l_s.f$coefficients
# wow... did not expect that lmao

# test.lambdas = c(0.1, 0.5, 1, 5, 10, 50, 100, 500)
test.lambdas = 10^((-5:8)/3)

mod.output = vector('list', length = length(test.lambdas))

for (i in 1:length(test.lambdas)) {
  mod.fit = glmmLasso(
    fix = flowering ~ size.prev + as.factor(flowering.prev) +
      early_prec.sum + early_minnTemp + last.freeze + 
      grow_meanTemp + summer_meanTemp + winter_meanTemp,
    rnd = list(Plot = ~ 1),
    data = demo.fl.clim,
    family = binomial(),
    lambda = test.lambdas[i],
    control = list(index = c(NA, NA, 1, 2, 3, 4, 5, 6))
  )
   mod.output[[i]] = mod.fit 
}

lapply(mod.output, function(x) x$coefficients) %>%
  do.call(what = rbind) %>%
  as.data.frame() %>%
  select(-`(Intercept)`) %>%
  mutate(lambda = test.lambdas) %>%
  pivot_longer(-lambda, names_to = 'coef', values_to = 'coef.est') %>%
  ggplot(aes(x = log(lambda, base = 10), y = coef.est, group = coef)) +
  geom_line(aes(colour = coef %in% c('(Intercept)', 'size.prev', 'as.factor(flowering.prev)TRUE')), size = 1.5) +
  geom_point(size = 1.5) +
  scale_x_reverse(breaks = (-1:3), labels = c(0.1, 1, 10, 100, 1000)) +
  labs(x = 'lambda (penalization strength)', y = 'coefficient estimate') +
  theme(legend.position = 'none', panel.background = element_blank())

cbind(lambda = test.lambdas, bic = sapply(mod.output, function(x) x$bic)) %>% plot(type = 'l')
# small lambda - BIC barely influenced
# but it looks like there's an inflection at ~ lambda = 100

mod.output[[which(test.lambdas == 100)]]
# these are all non-zero...
# although two are effectively zero?

sapply(mod.output, function(x) x$bic)
sapply(mod.output, function(x) x$bic) %>% plot(type = 'l')

# Okay... not nearly as helpful as I had been hoping!

# I guess cross validation would probably be ideal here...

##### Try the Harmony approach...

ranef.corrs = data.frame(
  Year = row.names(ranef(f_s.f.y.t)$Year),
  mod.ranef = unlist(ranef(f_s.f.y.t)$Year)
) %>%
  merge(y = clim)

# feel like there's a way to do pairwise correlations with purrr...
# actually I guess just a correlation matrix will work lol
ranef.corrs[,-1] %>%
  cor() %>%
  round(2)
# oh... summer and winter mean temperatures and are correlated...
# so is growing season mean temp and winter mean temp
# bleh
# Best correlates with the year-level random effects are summer temp, early
# season precip, and early season minimum temperature
# (I don't like min temp as a predictor though! ugh.)

ranef.corrs %>%
  pivot_longer(-c(Year, mod.ranef), names_to = 'coef', values_to = 'coef.val') %>%
  ggplot(aes(x = mod.ranef, y = coef.val, colour = coef)) +
  geom_point(size = 3) +
  facet_wrap(~ coef)
# lmao with this few data points anything is possible

f_s.f.earlyPrec = glmer(
  formula = flowering ~ size.prev + flowering.prev + early_prec.sum + (1 | Plot),
  data = demo.fl.clim,
  family = 'binomial'
)

f_s.f.earlytemp = glmer(
  formula = flowering ~ size.prev + flowering.prev + early_minnTemp + (1 | Plot),
  data = demo.fl.clim,
  family = 'binomial'
)

f_s.f.summertemp = glmer(
  formula = flowering ~ size.prev + flowering.prev + summer_meanTemp + (1 | Plot),
  data = demo.fl.clim,
  family = 'binomial'
)

AIC(f_s.f.earlyPrec, f_s.f.earlytemp, f_s.f.summertemp, f_s.f.y)
# yeah... lol.
# seems like way possible overfitting though
