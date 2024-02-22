# Is thre density dependence in seed production?
# Here: combining the seed set data with the demo data to look at number of
# neighbors in plot vs. number of seeds produced
# this is a quick and dirty script rather than a formal analysis (so far)

library(ggplot2)
library(dplyr)
library(tidyr)
library(glmmTMB)

rm(list = ls())

#--------------------------------------------------------------------------
# Overall plant density, no phenology

seed.phen.demo = read.csv('01_data_cleaning/out/demo_seed_phen_by_umbel_combined.csv') %>%
  # give me only the records where we have seed counts
  filter(!is.na(no.seeds))

demo = read.csv('01_data_cleaning/out/demo_imputed_survival.csv') %>%
  # Give me only the plants in 2021-2023
  filter(Year > 2020, No.umbels > 0, !is.na(No.umbels)) %>%
  # give me plot/coords
  separate(col = plantid, remove = FALSE, into = c('tag', 'plot', 'coord'), sep = '_') %>%
  select(-plot) %>%
  # give me only alphabetical coordinates and then extract them
  # this step removes five records
  filter(grepl('[A-Z]{1}$', coord)) %>%
  # ah... the above all removes plants from plot 9... just remove these
  filter(!Plot %in% 9) %>%
  mutate(x = as.numeric(gsub('[A-Z]', '', coord)), ylet = gsub('[0-9]', '', coord)) %>%
  merge(y = data.frame(ylet = LETTERS, y = 1:length(LETTERS)))

unique(demo$y) # sweet.

# To answer questions about density of flowering plants and its effects on seed
# set... we need/want
# - Seeds (obviously)
# - Number of flowering plants in each coordinate (!)
#   - this requires coords for each plant

plants.per.plot = demo %>%
  group_by(Year, Plot, x, y) %>%
  summarise(n = n())

plants.per.plot %>%
  ggplot(aes(x = x, y = y, fill = n)) +
  geom_tile() +
  facet_grid(Plot ~ Year)

plants.doubled = merge(
  x = plants.per.plot, y = plants.per.plot,
  by = c('Year', 'Plot'),
  suffixes = c('0', '1')
) %>%
  filter(abs(x0 - x1) < 3 & abs(y0 - y1) < 3) %>%
  group_by(Year, Plot, x0, y0) %>%
  summarise(n = sum(n1)) %>%
  # remove boundary plants
  filter(x0 > 0, x0 < 19, y0 > 1, y0 < 16) %>%
  rename(neighbors = n) %>%
  ungroup()

head(plants.doubled)
nrow(plants.doubled)

# Do slightly more post-processing on seed for merging

seed.phen.demo = seed.phen.demo %>%
  filter(!plot %in% 9) %>%
  separate(col = finalid, remove = FALSE, into = c('tag', 'Plot', 'coord'), sep = '_') %>%
  select(-c(tag, Plot)) %>%
  filter(grepl('[A-Z]{1}$', coord)) %>%
  mutate(x = as.numeric(gsub('[A-Z]', '', coord)), ylet = gsub('[0-9]', '', coord)) %>%
  merge(y = data.frame(ylet = LETTERS, y = 1:length(LETTERS))) %>%
  select(-ylet)


seed.phen.combined = merge(
  x = seed.phen.demo,
  y = plants.doubled,
  by.x = c('Year', 'plot', 'x', 'y'),
  by.y = c('Year', 'Plot', 'x0', 'y0'),
)

head(seed.phen.combined)
nrow(seed.phen.combined)

# Oh maybe I should add size
seed.demo.size = merge(
  seed.phen.combined, 
  demo %>% mutate(tagplot = paste(tag, Plot, sep = '_')) %>% select(tagplot, Year, No.leaves, Leaf.length, plantid),
  by = c('Year', 'tagplot')
)
# Some plants get lost (bummer)
# Maybe some duplicates?

# Check duplicates
# seed.demo.size %>% group_by(Year, tagplot) %>% filter(length(unique(plantid)) > 1) #%>% View()
# It's just one - the 3125 plant
seed.demo.size = seed.demo.size %>%
  group_by(Year, tagplot) %>%
  filter(length(unique(plantid)) == 1 | finalid == plantid) %>%
  ungroup() %>%
  select(-plantid)

# Now, add sizes
seed.demo.size = seed.demo.size %>%
  # loses eight records
  filter(!is.na(No.leaves) & !is.na(Leaf.length)) %>%
  mutate(cur.size = log(No.leaves * Leaf.length)) %>%
  select(-c(No.leaves, Leaf.length))

##### Plots

seed.demo.size %>%
  ggplot(aes(x = neighbors, y = no.seeds, colour = trt)) +
  geom_point(position = position_jitter(width = 0.25)) +
  scale_colour_manual(values = c('black', 'red', 'blue')) +
  scale_x_log10() +
  facet_wrap(~ Year)
# well... there's something going on here...
# wait this looks maybe like negative density dependence??

seed.demo.size %>%
  ggplot(aes(x = neighbors, y = no.seeds, colour = trt)) +
  geom_point(position = position_jitter(width = 0.5)) +
  scale_colour_manual(values = c('black', 'red', 'blue')) +
  scale_x_log10() +
  facet_grid(plot ~ Year) +
  theme(legend.position = 'none')
# oh looks like there really isn't any effect when faceted out like this...

##### Models

s_0 = glmmTMB(
  no.seeds ~ cur.size + (1 | Year) + (1 | plot / tagplot),
  family = 'nbinom2',
  ziformula = ~ (1 | plot / tagplot) + (1 | Year),
  data = seed.demo.size
)
# test says no treatment effect

s_n = glmmTMB(
  no.seeds ~ log(neighbors) + cur.size + (1 | Year) + (1 | plot / tagplot),
  family = 'nbinom2',
  ziformula = ~ (1 | plot / tagplot) + (1 | Year),
  data = seed.demo.size
)

s_nn = glmmTMB(
  no.seeds ~ log(neighbors) + cur.size + (1 | Year) + (1 | plot / tagplot),
  family = 'nbinom2',
  ziformula = ~ log(neighbors) + (1 | plot / tagplot) + (1 | Year),
  data = seed.demo.size
)

anova(s_n, s_0)
# p = 0.23

anova(s_nn, s_0)
# p = 0.42

### Conclusion: no evidence of influence of neighbors
# taking out some of the random effects from the zero inflation term doesn't
# hurt the model, but doing this still doesn't produce a signficant effect of
# neighbor density

#--------------------------------------------------------------------------
# Overall plant density, no phenology

phen.dates.open = read.csv('01_data_cleaning/out/phen_flowers_open_by_survey.csv') %>%
  # Remove plot 9
  filter(!plot %in% 9) %>%
  rename(Year = year)

# To get number of umbels open per day per coordinate,
# I need to merge in coords from demo...
phen.dates.open = merge(
  x = phen.dates.open,
  y = demo %>% select(tag, Plot, Year, x, y, plantid),
  by.x = c('Year', 'tag', 'plot'), 
  by.y = c('Year', 'tag', 'Plot'),
  suffixes = c('', '.demo')
) %>%
  group_by(Year, tag, plot, survey.period) %>%
  filter(length(unique(plantid)) == 1 | (plantid == plantid.demo)) %>%
  ungroup() %>%
  select(-plantid.demo)

nrow(phen.dates.open)

# Okay - now summarise by coordinate
# I want to get the number of plants (not umbels, but plants) at each coordinate on each day
phen.dates.summ = phen.dates.open %>%
  group_by(Year, plot, x, y, survey.period) %>%
  summarise(n.plants = n())

# Combine now?
phen.dates.neighbors = merge(
  x = phen.dates.summ, y = phen.dates.summ,
  by = c('Year', 'plot', 'survey.period'),
  suffixes = c('0', '1')
) %>%
  filter(abs(x0 - x1) < 3 & abs(y0 - y1) < 3) %>%
  group_by(Year, plot, survey.period, x0, y0) %>%
  summarise(n = sum(n.plants1)) %>%
  # remove boundary plants
  filter(x0 > 0, x0 < 19, y0 > 1, y0 < 16) %>%
  rename(neighbors = n) %>%
  ungroup()

head(phen.dates.neighbors)
head(phen.dates.neighbors)

# Now... merge the neighbor list back in with the flower list

phen.dates.combined = merge(
  phen.dates.open,
  phen.dates.neighbors,
  by.x = c('Year', 'plot', 'survey.period', 'x', 'y'),
  by.y = c('Year', 'plot', 'survey.period', 'x0', 'y0')
)

nrow(phen.dates.open)
nrow(phen.dates.neighbors)
nrow(phen.dates.combined) 
# ah - some of these are taken out because they're on the boundary

# Last step - average over periods observed to get one neighbor density per
# plant
phen.dates.combined = phen.dates.combined %>%
  group_by(Year, plot, tag, plantid, x, y) %>%
  summarise(
    n.periods = n(),
    mean.neighbors = mean(neighbors)
  )

nrow(phen.dates.combined) 
# damn... not that many

# Oh right... don't have seed set in here yet lol

seed.phen.dates = merge(
  seed.demo.size %>% select(Year, tagplot, plot, x, y, no.seeds, cur.size, trt),
  phen.dates.combined
)

nrow(seed.phen.dates)
# cool. nice round number

##### Plot

seed.phen.dates %>%
  ggplot(aes(x = mean.neighbors, y = no.seeds, colour = trt)) +
  geom_point(position = position_jitter(width = 0.5)) +
  scale_colour_manual(values = c('black', 'red', 'blue')) +
  scale_x_log10() +
  facet_grid(plot ~ Year)
# not super compelling

##### Model

s_0 = glmmTMB(
  formula = no.seeds ~ cur.size + (1 | plot / plantid) + (1 | Year),
  ziformula = ~ (1 | plot / plantid) + (1 | Year),
  family = 'nbinom2',
  data = seed.phen.dates
)

s_n = glmmTMB(
  formula = no.seeds ~ cur.size + log(mean.neighbors, base = 2) + (1 | plot / plantid) + (1 | Year),
  ziformula = ~ (1 | plot / plantid) + (1 | Year),
  family = 'nbinom2',
  data = seed.phen.dates
)

anova(s_n, s_0)
# Wow. Shit. Fuck. Goddamnit. How?

summary(s_n)
# negative density dependence... More neighbors means fewer seeds?

s_n_t = glmmTMB(
  formula = no.seeds ~ cur.size + trt + log(mean.neighbors, base = 2) + (1 | plot / plantid) + (1 | Year),
  ziformula = ~ (1 | plot / plantid) + (1 | Year),
  family = 'nbinom2',
  data = seed.phen.dates
)

anova(s_n_t, s_n) # no treatment effect

s_n_p = glmmTMB(
  formula = no.seeds ~ cur.size + n.periods + log(mean.neighbors, base = 2) + (1 | plot / plantid) + (1 | Year),
  ziformula = ~ (1 | plot / plantid) + (1 | Year),
  family = 'nbinom2',
  data = seed.phen.dates
)

# Okay...
# - coefficient is ~ -0.1, so a doubling of the number of neighbors means about a 10% reduction in seed set
# cool (not cool)
