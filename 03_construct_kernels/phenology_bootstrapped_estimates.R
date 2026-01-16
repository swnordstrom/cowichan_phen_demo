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
library(purrr)

cat('Bootstrapping phenology means...\n')

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

# Get phenology at umbel level
phen.by.umbel = all.data %>% 
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

# However, for our purposes we want to get the mean umbel emergence date by plant

phen.by.plant = phen.by.umbel %>%
  group_by(Plot, plantid, trt, Year) %>%
  summarise(phen.julian = mean(phen.julian), n.umbel = n()) %>%
  mutate(Year = factor(Year)) %>%
  ungroup()

# ====== Bootstrap models ====================================================
# Here: resampling the observed dataset, fitting models, storing *model
# coefficients*

# Fit original model
# This is needed for centering the bootstrapped coefficients
d_t = glmmTMB(
  phen.julian ~ trt + Year + (1 | Plot / plantid),
  data = phen.by.plant
)

# Define number of bootstraps
n.straps = 500

# Set seed for reproducibility
set.seed(9908847)

# Data resampling and model fitting
# phen.boot output df has one row per bootstrap, column for each model
# coefficient
phen.boot = phen.by.plant %>%
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
  lapply(
      function(df) {
        glmmTMB(
          phen.julian ~ trt + Year + (1 | Plot / plantid),
          data = df
        ) %>%
          # extract model parameters
          (function(mod) mod$fit$par)
      }
  ) %>%
  # mclapply(
  #   function(df) {
  #     glmmTMB(
  #       phen.julian ~ trt + Year + (1 | Plot / plantid),
  #       data = df
  #     ) %>%
  #       # extract model parameters
  #       (function(mod) mod$fit$par)
  #   },
  #   mc.cores = 6
  # ) %>%
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
# for (i in 1:n.straps) {
#   phen.list.out[[i]] = phen.backbone %>%
#     mutate(
#       pred.phen = predict(
#         d_t, newdata = phen.backbone, allow.new.levels = TRUE, re.form = ~ 0,
#         newparams = phen.boot[i, -(1:2)]
#       ) 
#     ) %>%
#     group_by(trt) %>%
#     summarise(mean.phen = mean(pred.phen)) %>%
#     mutate(boot = i)
#   
#   print(i)
# }
phen.list.out = map(
  1:n.straps,
  \(i) phen.backbone %>%
      mutate(
        pred.phen = predict(
          d_t, newdata = phen.backbone, allow.new.levels = TRUE, re.form = ~ 0,
          newparams = phen.boot[i, -(1:2)]
        )
      ) %>%
      group_by(trt) %>%
      summarise(mean.phen = mean(pred.phen)) %>%
      mutate(boot = i),
  .progress = FALSE
)


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

cat('Done.\n')

# ====== Bootstrap means in interaction model ==============================
# Repeating the bootstrapping procedure on a model with an interaction
# *for visualization purposes*.

# # Fit original model
# # This is needed for centering the bootstrapped coefficients
# d_ty = glmmTMB(
#   phen.julian ~ trt * Year + (1 | Plot / plantid),
#   data = phen.by.plant
# )
# 
# set.seed(8823104)
# 
# phen.boot = phen.by.plant %>%
#   # Add an 'obs.no' for distinguishing umbels
#   group_by(plantid, Year) %>%
#   mutate(obs.no = 1:n()) %>%
#   ungroup() %>%
#   uncount(weights = n.straps) %>%
#   # Label these entries with a `samp` (sample) column to delineate different
#   # bootstrap samples
#   group_by(plantid, Year, obs.no) %>%
#   mutate(samp = 1:n()) %>%
#   # Perform the resampling, preserving plot and survival structure
#   group_by(Plot, Year, samp) %>%
#   sample_n(size = n(), replace = TRUE) %>%
#   ungroup() %>%  
#   # Split the dataset by each sample and re-fit the umbel success/seed model
#   split(.$samp) %>%
#   lapply(
#     function(df) {
#       glmmTMB(
#         phen.julian ~ trt * Year + (1 | Plot / plantid),
#         data = df
#       ) %>%
#         # extract model parameters
#         (function(mod) mod$fit$par)
#     }
#   ) %>%
#   # mclapply(
#   #   function(df) {
#   #     glmmTMB(
#   #       phen.julian ~ trt + Year + (1 | Plot / plantid),
#   #       data = df
#   #     ) %>%
#   #       # extract model parameters
#   #       (function(mod) mod$fit$par)
#   #   },
#   #   mc.cores = 6
#   # ) %>%
#   # Combine into single data frame
#   do.call(rbind, .) %>%
#   data.frame() %>%
#   mutate(mod = 'phen', i = 1:n.straps) %>%
#   select(mod, i, everything())
# 
# # Center bootstrapped coefficients
# phen.boot[,-(1:2)] = phen.boot[-(1:2)] + matrix(
#   (d_ty$fit$par - colMeans(phen.boot[,-(1:2)])),
#   nrow = n.straps, ncol = length(d_t$fit$par), byrow = TRUE
# )
# 
# # Get predicted means from the interaction model
# # NOTE: not aggregating by year here (because want annual differences)
# phen.list.out = map(
#   1:n.straps,
#   \(i) phen.backbone %>%
#     mutate(
#       pred.phen = predict(
#         d_ty, newdata = phen.backbone, allow.new.levels = TRUE, re.form = ~ 0,
#         newparams = phen.boot[i, -(1:2)]
#       )
#     ) %>%
#     mutate(boot = i),
#   .progress = FALSE
# )
# 
# Get treatment effects (by year)
# phen.interaction.df = do.call(rbind, phen.list.out) %>%
#   pivot_wider(names_from = trt, values_from = pred.phen) %>%
#   mutate(d.c = drought - control, i.c = irrigated - control) %>%
#   select(-c(control, drought, irrigated)) %>%
#   pivot_longer(c(d.c, i.c), names_to = 'contrast', values_to = 'phen.diff')
# 
# phen.interaction.df %>%
#   mutate(
#     Year = as.numeric(as.character(Year)),
#     contrast = paste(ifelse(contrast %in% 'd.c', 'drought', 'irrigated'), 'vs. controls')
#   ) %>%
#   ggplot(aes(x = Year, y = phen.diff, colour = contrast)) +
#   annotate('segment', x = 2020.5, xend = 2024.5, y = 0, yend = 0, linetype = 2, colour = 'gray77') +
#   geom_point(size = 3, alpha = 0.125, position = position_jitter(width = 0.25)) +
#   scale_colour_manual(values = c('goldenrod1', 'dodgerblue')) +
#   scale_y_continuous(breaks = (-4:3) * 2) +
#   labs(y = 'Phenological shift (days)') +
#   theme(
#     panel.background = element_blank(),
#     legend.position = 'none'
#   ) +
#   facet_wrap(~ contrast)
# 
# ggsave('04_analysis/figures/fig_supp_trt_year.png', width = 5, height = 3)
# 
# empirical p-values for each year
# phen.interaction.df %>%
#   group_by(Year, contrast) %>%
#   mutate(bootstrap.mean = mean(phen.diff)) %>%
#   summarise(p = mean(sign(phen.diff) != sign(bootstrap.mean)))
# 