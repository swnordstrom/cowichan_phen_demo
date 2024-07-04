# Script for finalizing reconciliation of seed and demography data
# This script reconciles 2024 data ONLY
# It reads in a version of the seed dataset (2021-2024) where 2021-2023 have
# been reconciled with demo (see script currently named
# `01_data_cleaning/combine_demo_seed_mumb_data_2021-2023.R`) with the final-ish
# version of demo (including 2024)
# sn - 26 jun 2024 init

# Packages
library(ggplot2)
library(dplyr)
library(tidyr)

# Clear namespace
rm(list = ls())

# Read in seed data
seed = read.csv('01_data_cleaning/out/seed_reconciled_2021-2023.csv')

# Read in phen data
phen = read.csv('01_data_cleaning/out/phenology_buds_deaths_all.csv')

# Read in demo data
demo = read.csv('01_data_cleaning/out/Lomatium_demo_2016-2024.csv')

# Merge together phen and seed

phen.seed.plants = merge(
  x = seed %>% distinct(plantid, year) %>% mutate(in.seed = TRUE), 
  y = phen %>% distinct(plantid, year) %>% mutate(in.phen = TRUE), 
  by = c('plantid', 'year'), all.x = TRUE, all.y = TRUE
) %>%
  mutate(across(c(in.seed, in.phen), function(x) ifelse(is.na(x), FALSE, x)))

phen.seed.plants

phen.seed.plants %>%
  group_by(year, in.seed, in.phen) %>%
  summarise(n = n())

phen.seed.plants %>%
  filter(!(in.seed & in.phen)) %>%
  # filter(!(!i))
  filter(year > 2021) %>%
  arrange(year, plantid)

# Subset out 2024 data in each dataset
# Make tagplot column (coordinates unlikely to match)
# Then merge?

merged.2024 = merge(
  x = seed %>% 
    filter(year %in% 2024) %>% 
    select(-year) %>%
    mutate(tagplot = paste(tag, plot, sep = '_')) %>%
    mutate(in.seed = TRUE),
  y = demo %>% 
    filter(Year %in% 2024) %>%
    separate(plantid, into = c('tag', 'plot', 'coord'), sep = '_', remove = FALSE) %>%
    mutate(tagplot = paste(gsub('b', '', tag), plot, sep = '_')) %>%
    select(-c(tag, plot, Plot, trt, coord, Year)) %>%
    mutate(in.demo = TRUE),
  by = 'tagplot', all.x = TRUE, all.y = TRUE, suffixes = c('.seed', '.demo'),
)

merged.2024 %>%
  group_by(in.demo, in.seed) %>%
  summarise(n = n())
# Okay - ~250 cases where it is in demo but not seed
# (most of these are probably fine)
# and ~50 cases where it's in seed but not demo

merged.2024 %>% 
  filter(is.na(in.demo)) %>% 
  distinct(tagplot, plantid.seed, plot, tag) %>%
  arrange(plot)
# ah. Some of these are toothpicks that we just won't be able to merge.

# Go through these one by one?

# # tagplot plantid.seed plot  tag
# 1    BSP_3     BSP_3_4F    3  BSP
# 7    YTP_5    YTP_5_11E    5  YTP
# 8   3789_6   3789_6_11A    6 3789
# 9   7587_6   7587_6_12A    6 7587
# 10  3298_7   3298_7_10B    7 3298
# 11  3477_7   3477_7_14D    7 3477
# 12  3886_7   3886_7_10B    7 3886
# 13     R_7       R_7_0I    7    R
# 14 3176_13   3176_13_0E   13 3176
# 15 3180_13  3180_13_16G   13 3180
# 16 3481_13  3481_13_12F   13 3481
# 17 3726_13   3726_13_0D   13 3726
# 18 3769_13   3769_13_1F   13 3769
# 19  YTP_13    YTP_13_6J   13  YTP
# 20 3480_14  3480_14_19E   14 3480
# 21  RSP_14   RSP_14_19E   14  RSP
# 22 WHTP_14   WHTP_14_0E   14 WHTP
# 23 3076_15  3076_15_13G   15 3076
# 24 3135_15  3135_15_14J   15 3135
# 25 3813_15  3813_15_14C   15 3813
# 26 7655_15   7655_15_4J   15 7655
# 27  BSP_15    BSP_15_0J   15  BSP
# 28  BTP_15   BTP_15_12G   15  BTP
# 29  RSP_15    RSP_15_2A   15  RSP

