# =============================================================================== =
# Script for performing bootstrapping on phenology-treatment effects
# Reads in data demo+phen data, fits models on both observed and bootstrapped
# datasets
# These bootstrapped model coefficients are then used to get bootstrapped
# treatment means in an identical manner used to get kernel bootstrappped entries
# Doing this in a separate script allows us to read in and use the *same*
# treatment means within bootstrap replicates across subkernels
# This exports a csv i.e. so it only needs to be run once.
# 22 Jan 2025
# ===============================================================================

# ====== Setup ==================================================================

# Load packages
library(parallel)
library(ggplot2) # used for aux plot but not necessary
library(glmmTMB)
library(dplyr)
library(tidyr)

# Get rid of this super annoying feature
options(dplyr.summarise.inform = FALSE)

# Clear the namespace
rm(list = ls())

# Read in demo data and merge with treatment info
all.data = merge(
  x = read.csv('01_data_cleaning/out/demo_phen_seed_2016-2024_final.csv'),
  y = read.csv('00_raw_data/plot_treatments.csv'),
  by.x = 'Plot', by.y = 'plot'
)

# Umbel-level budding phenology data
# (used only in bootstrapping)
phen = all.data %>% 
  filter(in.phen) %>%
  # Split out the bud dates for bud date models; the most umbels seen in a
  # plant is 12, so use separate() to kick these out and then pivot_long to get
  # one row per umbel
  # (first - need to get one row per plant - do a distinct())
  distinct(Year, plantid, .keep_all = TRUE) %>%
  separate_wider_delim(phen.julis, names = paste0('uu', 1:12), delim = ';', too_few = 'align_start') %>%
  pivot_longer(starts_with('uu'), names_to = 'umbel.number', values_to = 'phen.julian') %>%
  filter(!is.na(phen.julian)) %>%
  mutate(phen.julian = as.numeric(gsub('\\s', '', phen.julian))) %>%
  mutate(Year = factor(Year))

# ====== Bootstrap models ====================================================
# Here: resampling the observed dataset, fitting models, storing *model
# coefficients*

# Set seed for reproducibility
set.seed(9908847)

# Fit original model
# This is needed for centering the bootstrapped coefficients
d_t = glmmTMB(
  phen.julian ~ trt + Year + (1 | Plot / plantid),
  data = phen
)

# Define number of bootstraps
n.straps = 100

# Data resampling and model fitting
# phen.boot output df has one row per bootstrap, column for each model
# coefficient
phen.boot = phen %>%
  # Add an 'obs.no' for distinguishing umbels
  group_by(plantid, Year) %>%
  mutate(obs.no = 1:n()) %>%
  ungroup() %>%
  uncount(weights = n.straps) %>%
  # Label these entries with a `samp` (sample) column to delineate different
  # bootstrap samples
  group_by(plantid, Year, obs.no) %>%
  mutate(samp = 1:n()) %>%
  # Perform the resampling, preserving plot and survival structure
  group_by(Plot, Year, samp) %>%
  sample_n(size = n(), replace = TRUE) %>%
  ungroup() %>%  
  # Split the dataset by each sample and re-fit the umbel success/seed model
  split(.$samp) %>%
  mclapply(
    function(df) {
      glmmTMB(
        phen.julian ~ trt + Year + (1 | Plot / plantid),
        data = df
      ) %>%
        # extract model parameters
        (function(mod) mod$fit$par)
    },
    mc.cores = 6
  ) %>%
  # Combine into single data frame
  do.call(rbind, .) %>%
  data.frame() %>%
  mutate(mod = 'phen', i = 1:n.straps) %>%
  select(mod, i, everything())

# Center bootstrapped coefficients
phen.boot[,-(1:2)] = phen.boot[-(1:2)] + matrix(
  (d_t$fit$par - colMeans(phen.boot[,-(1:2)])),
  nrow = n.straps, ncol = length(d_t$fit$par), byrow = TRUE
)


# ====== Bootstrap means =====================================================
# Taking the bootstrapped coefficients from above and using them to estimate
# treatment means

# List to store outputs
phen.list.out = vector('list', length = n.straps)

# Initiate a scaffold to estimate estimates over
phen.backbone = expand.grid(
  trt = c('control', 'drought', 'irrigated'),
  Year = factor(2021:2024)
)

# Do the predictions
for (i in 1:n.straps) {
  phen.list.out[[i]] = phen.backbone %>%
    mutate(
      pred.phen = predict(
        d_t, newdata = phen.backbone, allow.new.levels = TRUE, re.form = ~ 0,
        newparams = phen.boot[i, -(1:2)]
      ) 
    ) %>%
    group_by(trt) %>%
    summarise(mean.phen = mean(pred.phen)) %>%
    mutate(boot = i)
  
  print(i)
}

# Combine list into data frame
phen.boot.trt = do.call(rbind, phen.list.out) %>% mutate(boot = as.numeric(boot))

# Visualization (just to make sure everything works properly)
phen.boot.trt %>% 
  ggplot(aes(x = trt, y = mean.phen, colour = trt)) + 
  geom_point(position = position_jitter(width = 0.5), size = 3) +
  scale_colour_manual(values = c('black', 'red', 'blue'))

# Export file
phen.boot.trt %>% 
  # Pivot to wider form to save space
  pivot_wider(names_from = trt, values_from = mean.phen) %>%
  write.csv(
    row.names = FALSE, na = '',
    file = '03_construct_kernels/out/phenology_bootstrapped_means.csv'
  )

