# Script for exploring "recruit" data
# SN 11 Dec 2023

library(ggplot2)
library(dplyr)
library(tidyr)
library(lme4)

rm(list = ls())

# Read in demo data
demo = read.csv('01_data_cleaning/out/demo_imputed_survival.csv') %>%
  arrange(plantid, Year) %>%
  mutate(
    flowering = case_when(
      No.umbels > 0 ~ TRUE,
      !No.umbels ~ FALSE,
      is.na(No.umbels) ~ FALSE,
      .default = NA
    )
  )

head(demo)

# How many of plants were seen previously ("reseen") and how many are only being
# seen for the first time?
# Break this down by year
# Also - how many of these plants (in either group) are flowering?
demo %>%
  filter(surv) %>%
  mutate(is.reseen = ifelse(duplicated(plantid), 'reseen', 'firsty')) %>%
  group_by(Year, is.reseen) %>%
  summarise(n = n(), pfl = mean(flowering)) %>%
  pivot_wider(names_from = is.reseen, values_from = c(n, pfl), values_fill = 0)
  
# Okay... based on this, I don't trust that the 2017 new plants are new recruits
# There are a lot of them, and most of those seen in 2017 flowered
# likewise incomplete sampling in 2020 means those shouldn't be analyzed
# going to assume 2018-2019, 2021-2023 are accurate

# Get only those plants being seen for the first time
croots = demo %>%
  filter(surv) %>%
  filter(!duplicated(plantid)) %>%
  # get rid of plants in 2017
  filter(Year > 2017)

head(croots)
nrow(croots)

# Distribution of sizes over time
croots %>%
  filter(!is.na(Leaf.length) & !is.na(No.leaves)) %>%
  mutate(size = log(Leaf.length * No.leaves)) %>%
  # removing 2020
  filter(!Year %in% 2020) %>%
  mutate(Year = factor(Year)) %>%
  ggplot(aes(x = size, group = Year, fill = Year)) +
  geom_histogram(position = position_identity(), alpha = 0.5, binwidth = 0.25) +
  scale_fill_brewer(palette = 'Set2')
# well... size peak varies from year to year. so that's interesting at least.

# Sizes by treatment over time
# (leave 2020 in for now...)
croots %>%
  filter(!is.na(Leaf.length) & !is.na(No.leaves)) %>%
  mutate(size = log(Leaf.length * No.leaves)) %>%
  ggplot(aes(x = Year, y = size, colour = trt)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.5)) +
  scale_colour_manual(values = c('black', 'red', 'blue'))
# interesting.
# some years seem to have treatment effects - looks often like this is irrigated
# plots having smaller recruits

# Look at leaf sizes, etc. by 
croots %>%
  mutate(Year = factor(Year)) %>%
  filter(!is.na(No.leaves) & !is.na(Leaf.length)) %>%
  ggplot(aes(x = No.leaves, y = Leaf.length, fill = Year)) +
  geom_point(
    size = 3, shape = 21,
    position = position_jitterdodge(dodge.width = .75, jitter.width = .4)
  )
# six leaves?
# some of these leaves are effing massive!!!

# Looking at the USDA docs on congeners (specifically, L. nudicaule), and some
# photos, I suppose I could imagine multi-leaved plants. But also, 30 cm??
# Seriously? That's a foot-long leaf on a supposedly-newborn (or never before
# seen) plant.

# How many of these are plants with a 'b' in the tag, i.e., plants with tags
# re-assigned after a long gap in between sightings?

croots %>%
  mutate(
    Year = factor(Year),
    btag = grepl('\\db', plantid)
  ) %>%
  filter(!is.na(No.leaves) & !is.na(Leaf.length)) %>%
  ggplot(aes(x = No.leaves, y = Leaf.length, fill = Year, alpha = btag)) +
  geom_point(
    size = 3, shape = 21,
    position = position_jitterdodge(dodge.width = .75, jitter.width = .4)
  ) +
  scale_alpha_manual(values = c(0.05, 1))
# So some, but not that many!

# What does this look like if we (more or less) ignore those plants?
croots %>%
  mutate(
    Year = factor(Year),
    btag = grepl('\\db', plantid)
  ) %>%
  filter(!is.na(No.leaves) & !is.na(Leaf.length)) %>%
  ggplot(aes(x = No.leaves, y = Leaf.length, fill = Year, alpha = btag)) +
  geom_point(
    size = 3, shape = 21,
    position = position_jitterdodge(dodge.width = .75, jitter.width = .4)
  ) +
  scale_alpha_manual(values = c(1, 0.05))
# Still plenty of plants with >20cm leaves

# Not sure how much I trust this size data!

#####
#####
#####

# What about "recruits" per plot vs. number of flowering plants in that plot?
# (need to do this in two steps... for some reason)
flw.rct.same.year = demo %>%
  mutate(first = !duplicated(plantid)) %>%
  filter(surv) %>%
  group_by(Year, Plot, trt) %>%
  summarise(
    n.rct = sum(first),
    n.flw = sum(flowering)
  ) %>%
  ungroup()

flw.rct = merge(
    x = flw.rct.same.year %>% select(Year, Plot, trt, n.rct),
    # add 1 to year so that flowering merges with year recruit was seen
    y = flw.rct.same.year %>% select(Year, Plot, trt, n.flw) %>% mutate(Year = Year + 1)
  ) %>%
  filter(Year > 2017, !Year %in% 2020) %>%
  arrange(Year, Plot) %>%
  mutate(across(c(Year, Plot), factor))

head(flw.rct)

# Check to make sure this is right
flw.rct.same.year %>% filter(Plot %in% 2)
flw.rct %>% filter(Plot %in% 2)
# in 2018, plot 2 had 5 new plants
# in 2017, plot 2 had 7 flowering plants
# good

flw.rct.plot = flw.rct %>%
  ggplot(aes(x = n.flw, y = n.rct, colour = trt)) +
  geom_point(
    size = 3,
    position = position_jitter(width = .25, height = .25)
  ) +
  scale_colour_manual(values = c('black', 'red', 'blue'))

flw.rct.plot
# no obvious pattern based on this plot...

# What woudl a glmer say?
flw.rct.mod = glmer(
  n.rct ~ n.flw + (1 | Plot) + (1 | Year),
  family = 'poisson',
  data = flw.rct
)

summary(flw.rct.mod)
# LMAO A NEGATIVE EFFECT

flw.rct.trt.mod = glmer(
  n.rct ~ trt * n.flw + (1 | Plot) + (1 | Year),
  family = 'poisson',
  data = flw.rct
)

summary(flw.rct.trt.mod)
# Unbelievable.

flw.rct.trt.mod.out = flw.rct %>%
  mutate(
    trt.pred = predict(flw.rct.trt.mod, newdata = ., re.form = ~ 0),
    plt.pred = predict(flw.rct.trt.mod, newdata = ., re.form = ~ (1 | Plot))
)

flw.rct.plot +
  geom_line(
    data = flw.rct.trt.mod.out,
    aes(y = trt.pred, colour = trt)
  ) +
  geom_line(
    data = flw.rct.trt.mod.out,
    aes(y = plt.pred, colour = trt, group = Plot),
    linetype = 2, linewidth = 0.5
  )

rct.per.flw.mod = glmer(
  n.rct ~ trt + (1 | Plot) + (1 | Year),
  offset = n.flw,
  family = 'poisson',
  data = flw.rct
)
# (non-positive definite VcV)
# sweet. awesome. amazing.

flw.rct %>% filter(!n.flw)
# three plots with zero flowering plants, two of which have new plants...
flw.rct.same.year %>% filter(Plot %in% c(5, 10)) %>% arrange(Plot, Year)
# looks okay to me...

# Let's look at plot 10, year 2020: no flowering plants. sampling?
demo %>% filter(Year %in% 2020, Plot %in% 10)
# there are just not that many plants in this plot
demo %>% filter(Year %in% 2017, Plot %in% 10)
# then how are we still getting recruits there... lol

###
###
### 

# Looking now at what happens if you say only one-leaved plants are new recruits
head(croots)

croots1 = croots %>% filter(No.leaves < 2)

nrow(croots1) # 309 plants
table(croots1$Year) # year-to-year variation
with(croots1, table(trt, Year))
# Looks like control has the least... lmao.

# Oh... how many are flowering?
croots1 %>% filter(flowering) # ~20 plants in here that flowered
# Are plants that flower larger than plants that don't?
croots1 %>%
  ggplot(aes(x = flowering, y = Leaf.length, colour = trt)) +
  geom_point(size = 4, position = position_jitter(width = 0.4)) +
  scale_colour_manual(values = c('black', 'red', 'blue'))
# Slightly larger but not by much
# (I guess there's also an NA in here? lol)
# Looks like the drought plants are huge, wow

croots1 = croots1 %>% filter(!flowering, !is.na(flowering), !is.na(Leaf.length))

nrow(croots1) # 285 plants
table(croots1$Year)
with(croots1, table(trt, Year))

# Visualize distribution for time, treatment
croots1 %>%
  ggplot(aes(x = Year, y = Leaf.length)) +
  geom_point(
    aes(colour = trt, group = trt), 
    position = position_dodge(width = 0.25), 
    size = 3, alpha = 0.25) +
  scale_colour_manual(values = c('black', 'red', 'blue'))
# 2022 is much larger than the rest
# there's def year-to-year variation here

# Are these lognormally distributed?

croots1 %>%
  mutate(log.len = log(Leaf.length)) %>%
  ggplot(aes(x = log.len, group = trt, fill = trt)) +
  geom_histogram(alpha = 0.5, position = position_identity()) +
  scale_fill_manual(values = c('black', 'red', 'blue'))
# hmm... long left tail

croots1 %>%
  ggplot(aes(x = Leaf.length, group = trt, fill = trt)) +
  geom_histogram(alpha = 0.5, position = position_identity()) +
  scale_fill_manual(values = c('black', 'red', 'blue'))
# This certainly isn't normally distributed
# Oh well.

# Some models

croots1 = croots1 %>% mutate(log.len = log(Leaf.length))

s_0 = lmer(
  formula = log.len ~ (1 | Plot),
  data = croots1
)

s_y = lmer(
  formula = log.len ~ (1 | Year) + (1 | Plot),
  data = croots1
)

AIC(s_y, s_0)
# Year-to-year variation, very cool

s_y_t = lmer(
  formula = log.len ~ trt + (1 | Year) + (1 | Plot),
  data = croots1
)

# Just in case - no year effect, treatment effect
s_t = lmer(
  formula = log.len ~ trt + (1 | Plot),
  data = croots1
)

AIC(s_y, s_t, s_y_t) %>% mutate(daic = round(AIC - min(AIC), 2))
# no treatment effect - very cool!

qqnorm(residuals(s_y))
qqline(residuals(s_y))
# left tail, not great, something weird happening upper right
# oh well

summary(s_y)

croots1 %>%
  mutate(resid = residuals(s_y)) %>%
  ggplot(aes(x = Year, y = resid, colour = trt)) +
  geom_point(position = position_jitter(width = 0.25), size = 3) +
  scale_colour_manual(values = c('black', 'red', 'blue'))
# 2022 is a little small but otherwise looks fine to me
