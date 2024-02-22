# Is thre density dependence in seed production?
# Here: combining the seed set data with the demo data to look at number of
# neighbors in plot vs. number of seeds produced
# this is a quick and dirty script rather than a formal analysis (so far)

library(ggplot2)
library(dplyr)
library(tidyr)
library(glmmTMB)

rm(list = ls())

seed.phen.demo = read.csv('01_data_cleaning/out/demo_seed_phen_by_umbel_combined.csv') %>%
  # give me only the records where we have seed counts
  filter(!is.na(no.seeds))

head(phen.seed.demo)

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
  summarise(n = n()) %>%
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
seed.phen.combined %>%
  filter(!is.na(No.leaves) & !is.na(Leaf.length)) %>% nrow()
# oh somehow I don't have the measurements in here
# ehhhhhh fuck it, let's not worry about it for now

# Plot

seed.phen.combined %>%
  ggplot(aes(x = neighbors, y = no.seeds, colour = trt)) +
  geom_point(position = position_jitter(width = 0.5)) +
  scale_colour_manual(values = c('black', 'red', 'blue')) +
  facet_wrap(~ Year)
# well... there's something going on here...
# wait wtf this looks maybe like negative density dependence??


