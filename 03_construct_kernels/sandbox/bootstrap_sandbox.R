# --- Setup ---------------------------------------------------------

library(ggplot2)
library(dplyr)
library(tidyr)
# library(purrr)
library(parallel)
library(glmmTMB)

rm(list = ls())

all.data = merge(
  x = read.csv('01_data_cleaning/out/demo_phen_seed_2016-2024_final.csv'),
  y = read.csv('00_raw_data/plot_treatments.csv'),
  by.x = 'Plot', by.y = 'plot'
)

nrow(all.data)
head(all.data)

all.demo = all.data %>% 
  filter(in.demo) %>%
  distinct(plantid, Year, .keep_all = TRUE)

nrow(all.demo)

# Get survival dataset

# Survival dataset:
# (Two versions: a size-dependent one and a size-independent one
# almost surely we will use the size-dependent one for analysis)

demo.surv = merge(
  # Demo in time step t+1
  x = all.demo %>% 
    mutate(prev.year = Year - 1) %>%
    rename(surv.year = Year) %>%
    select(Plot, plantid, surv.year, prev.year, No.leaves, Leaf.length, surv, trt),
  # Demo in time step t
  y = all.demo %>%
    # we are *only* interested in plants alive in time step t
    filter(surv) %>%
    # Select relevant columns
    select(Plot, plantid, Year, No.leaves, Leaf.length, trt),
  by.x = c('Plot', 'plantid', 'prev.year', 'trt'),
  by.y = c('Plot', 'plantid', 'Year', 'trt'),
  suffixes = c('', '.pre'),
  all.x = FALSE, all.y = FALSE
)

head(demo.surv)
# should be less than 1
table(demo.surv$surv, useNA = 'always')
# good

# Want to get a column where we can estimate sizes
demo.surv.sizes = demo.surv %>% 
  filter(
    !is.na(Leaf.length.pre) & !is.na(No.leaves.pre) &
      Leaf.length.pre > 0 & No.leaves.pre > 0
  ) %>%
  # Get rid of 2016 records because the sizes are not reliable
  filter(prev.year > 2016) %>%
  # Add size columns
  mutate(size.prev = log(No.leaves.pre * Leaf.length.pre))

nrow(demo.surv.sizes)
table(demo.surv.sizes$surv, useNA = 'always')

# Subset for growth estimation
demo.grow = demo.surv.sizes %>% 
  filter(surv, !is.na(Leaf.length) & !is.na(No.leaves) & Leaf.length > 0 & No.leaves > 0) %>%
  mutate(size.cur = log(Leaf.length * No.leaves))

# Finally: subset surv dataset to not include 2023-2024 surv
demo.surv.sizes = demo.surv.sizes %>% filter(surv.year < 2024)

# --- Boots ---------------------------------------------------------

set.seed(11225)

# Number of bootstraps
n.straps = 500

# Set up the bootstraps

surv.boots = demo.surv.sizes %>%
  sample_n(size = n.straps * nrow(.), replace = TRUE) %>%
  mutate(samp = rep_len(1:n.straps, length.out = nrow(.))) %>%
  split(.$samp) %>%
  mclapply(function(df) {
    glmmTMB(
      formula = surv ~ size.prev + (1 | Plot),
      family = 'binomial',
      data = df
    ) %>%
      (function(mod) mod$fit$par)
    },
    mc.cores = 6
  ) %>%
  do.call(rbind, .) %>%
  data.frame() %>%
  mutate(mod = 'surv', i = 1:n.straps) %>%
  select(mod, i, everything())

grow.boots = demo.grow %>%
  sample_n(size = n.straps * nrow(.), replace = TRUE) %>%
  mutate(samp = rep_len(1:n.straps, length.out = nrow(.))) %>%
  split(.$samp) %>%
  mclapply(
    function(df) {
      glmmTMB(
        size.cur ~ size.prev + size.prev * trt + (1 | prev.year) + (1 | prev.year:trt) + (1 | Plot / plantid),
        data = df
      ) %>%
        (function(mod) c(mod$fit$par, sigma = summary(mod)$sigma))
      },
    mc.cores = 6
  ) %>%
  # lapply(function(mod) c(mod$fit$par, sigma = summary(mod)$sigma)) %>%
  do.call(rbind, .) %>%
  data.frame() %>%
  mutate(mod = 'grow', i = 1:n.straps) %>%
  select(mod, i, sigma, everything())

# --- Original models -------------------------------------------------------

# Survival model
s_s = glmmTMB(
  formula = surv ~ size.prev + (1 | Plot),
  family = 'binomial',
  data = demo.surv.sizes
)

# Growth model
g_st.ty = glmmTMB(
  size.cur ~ size.prev + size.prev * trt + (1 | prev.year) + (1 | prev.year:trt) + (1 | Plot / plantid),
  data = demo.grow
)

# Get bootstrapped kernels

bootstrap.backbone = expand.grid(
  size.prev = (5:60)/10,
  size.cur = (5:60)/10,
  trt = c('control', 'drought', 'irrigated')
)

kernel.list = vector('list', length = n.straps)

for (i in 1:n.straps) {
  
  kernel.list[[i]] = bootstrap.backbone %>%
    # Predicted survival
    mutate(
      pred.surv = predict(
        newdata = .,
        object = s_s, type = 'response',
        newparams = surv.boots[i, -(1:2)],
        re.form = ~ 0, allow.new.levels = TRUE
      )
    ) %>%
    # Predicted growth
    mutate(
      pred.grow.mean = predict(
        newdata = .,
        object = g_st.ty, type = 'response',
        newparams = grow.boots[i, -(1:3)],
        re.form = ~ 0, allow.new.levels = TRUE
      )
    ) %>%
    # Predicted distribution of sizes in next time step
    mutate(p.grow.size = 0.1 * dnorm(size.cur, mean = pred.grow.mean, sd = grow.boots$sigma[i])) %>%
    # Combine all together to get overall size distribution in next time step
    mutate(p.size.cur = pred.surv * p.grow.size) %>%
    select(-c(pred.surv, pred.grow.mean, p.grow.size)) %>%
    mutate(boot = paste0('k', i))
  
  if (!(i %% 10)) print(i)
  
}

do.call(rbind, kernel.list) %>%
  pivot_wider(names_from = boot, values_from = p.size.cur) %>%
  write.csv(
    '03_construct_kernels/out/deterministic_growsurv_bootstrap.csv',
    row.names = FALSE, na = ''
  )
