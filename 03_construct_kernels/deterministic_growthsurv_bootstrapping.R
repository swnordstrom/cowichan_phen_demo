# Script for generating bootstrapped survival-growth subkernels, including
# perturbations
# Outputs:
# - 100 bootstrapped growth+surv sub-kernels
# - 100 bootstrapped perturbed growth+surv sub-kernels
# - 100 parameter differences used to generate sub-kernels

# --- Setup ---------------------------------------------------------

library(ggplot2)
library(dplyr)
library(tidyr)
library(parallel)
library(glmmTMB)

rm(list = ls())

# # Data preparation 
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

# --- Original models -------------------------------------------------------
# Need to run original models 

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

# --- Boots ---------------------------------------------------------

set.seed(11225)

# Number of bootstraps
n.straps = 100

# Generate bootstrapped estimates

surv.boots = demo.surv.sizes %>%
  # Do bootstrapped resampling
  # copy the data frame (each row duplicated)
  uncount(weights = n.straps) %>%
  # label these entries with a `samp` (sample) column to delineate different
  # bootstrap samples
  group_by(plantid, surv.year) %>%
  mutate(samp = 1:n.straps) %>%
  # Perform the resampling, preserving plot and year structure
  group_by(Plot, surv.year, samp) %>%
  sample_n(size = n(), replace = TRUE) %>%
  ungroup() %>%
  # Split the dataset by each sample and re-fit the survival model
  split(.$samp) %>%
  mclapply(
    function(df) {
      glmmTMB(
        formula = surv ~ size.prev + (1 | Plot),
        family = 'binomial',
        data = df
      ) %>%
        # Collect model parameters (for making predictions)
        (function(mod) mod$fit$par)
    },
    mc.cores = 6
  ) %>%
  # Combine these together into a single data frame
  do.call(rbind, .) %>%
  data.frame() %>%
  mutate(mod = 'surv', i = 1:n.straps) %>%
  select(mod, i, everything())

grow.boots = demo.grow %>%
  # Do bootstrapped resampling
  # copy the data frame (each row duplicated)
  uncount(weights = n.straps) %>%
  # label these entries with a `samp` (sample) column to delineate different
  # bootstrap samples
  group_by(plantid, surv.year) %>%
  mutate(samp = 1:n.straps) %>%
  # Perform the resampling, preserving plot and year structure
  group_by(Plot, surv.year, samp) %>%
  sample_n(size = n(), replace = TRUE) %>%
  ungroup() %>%
  # Split the dataset by each sample and re-fit the growth model
  split(.$samp) %>%
  mclapply(
    function(df) {
      glmmTMB(
        size.cur ~ size.prev + size.prev * trt + (1 | prev.year) + (1 | prev.year:trt) + (1 | Plot / plantid),
        data = df
      ) %>%
        # Collect model parameters (for making predictions)
        (function(mod) mod$fit$par)
    },
    mc.cores = 6
  ) %>%
  # Combine these together into a single data frame
  do.call(rbind, .) %>%
  data.frame() %>%
  mutate(mod = 'grow', i = 1:n.straps) %>%
  select(mod, i, everything())

# # Un-biasing the bootstrapped estimates
# Manually shifting individual columns of the bootstrap to have the mean of
# boostrapped samples be identical to the parameter values of the true model

# Survival bootstraps
surv.boots[,-(1:2)] = surv.boots[,-(1:2)] + matrix(
  (s_s$fit$par - colMeans(surv.boots[,-(1:2)])), 
  nrow = n.straps, ncol = length(s_s$fit$par), byrow = TRUE,
)

# Growth bootstraps
grow.boots[,-(1:2)] = grow.boots[,-(1:2)] + matrix(
  (g_st.ty$fit$par - colMeans(grow.boots[,-(1:2)])), 
  nrow = n.straps, ncol = length(g_st.ty$fit$par), byrow = TRUE,
)

# very small differences, all numerical rounding
mean((colMeans(surv.boots[,-(1:2)]) - s_s$fit$par)^2)
mean((colMeans(grow.boots[,-(1:2)]) - g_st.ty$fit$par)^2)

# --- Get bootstrapped kernels ------------------------------------

# Data frame for making predictions at each level
bootstrap.backbone = expand.grid(
  size.prev = (5:60)/10,
  size.cur = (5:60)/10,
  trt = c('control', 'drought', 'irrigated')
)

# List for storing each bootstrapped sample in
kernel.list = vector('list', length = n.straps)

# Do kernel estimation on each bootstrapped set of parameters

for (i in 1:n.straps) {
  
  kernel.list[[i]] = bootstrap.backbone %>%
    # Predicted survival
    mutate(
      pred.surv = predict(
        newdata = .,
        object = s_s, type = 'response',
        # survival prediction made with ith bootstrap parameter set
        newparams = surv.boots[i, -(1:2)],
        re.form = ~ 0, allow.new.levels = TRUE
      )
    ) %>%
    # Predicted growth
    mutate(
      pred.grow.mean = predict(
        newdata = .,
        object = g_st.ty, type = 'response',
        # growth prediction made with ith bootstrap parameter set
        newparams = grow.boots[i, -(1:2)],
        re.form = ~ 0, allow.new.levels = TRUE
      )
    ) %>%
    # Predicted distribution of sizes in next time step
    # note: 'betad' parameter is the log of the model's estimated residual variance term
    mutate(p.grow.size = 0.1 * dnorm(size.cur, mean = pred.grow.mean, sd = sqrt(exp(grow.boots$betad[i])))) %>%
    # Combine all together to get overall size distribution in next time step
    mutate(p.size.cur = pred.surv * p.grow.size) %>%
    select(-c(pred.surv, pred.grow.mean, p.grow.size)) %>%
    mutate(boot = paste0('b', i))
  
  if (!(i %% 10)) print(i)
  
}

# Combine kernels, convert to wide form, and export

do.call(rbind, kernel.list) %>%
  pivot_wider(names_from = boot, values_from = p.size.cur) %>%
  write.csv(
    '03_construct_kernels/out/deterministic_growsurv_bootstrap.csv',
    row.names = FALSE, na = ''
  )

# ------
# ------ Perturbation bootstrapping
# ------ (repeat above procedure but with perturbations at each vital rate)
# ------

# Parameters to perturb:
# - Survival: nothing (no treatment effects in any terms)
# - Growth:
#   - Intercept
#   - Slope

bootstrap.backbone = expand.grid(
  size.prev = (5:60)/10,
  size.cur = (5:60)/10,
  trt = c('control', 'drought', 'irrigated')
)

# Perturbation amount
delta = 0.0001

# Get a list for outputs
gs.pert.boot = vector('list', n.straps)

# Loop through each boot sample
for (i in 1:n.straps) {
  
  # Set up list for outputs
  this.boot = vector('list', 2)
  
  # 1: Growth model intercept
  
  this.boot[[1]] = bootstrap.backbone %>%
    # Predicted survival
    mutate(
      pred.surv = predict(
        newdata = .,
        object = s_s, type = 'response',
        newparams = surv.boots[i,-(1:2)],
        re.form = ~ 0, allow.new.levels = TRUE
      )
    ) %>%
    # Predicted growth
    mutate(
      pred.grow.mean = predict(
        newdata = .,
        object = g_st.ty, type = 'response',
        # Apply the perturbation
        newparams = grow.boots[i, -(1:2)] %>%
          (function(x) {
            x[1] <- x[1] + delta
            return(x)
          }),
        re.form = ~ 0, allow.new.levels = TRUE
      )
    ) %>%
    mutate(p.grow.size = 0.1 * dnorm(size.cur, pred.grow.mean, sd = sqrt(exp(grow.boots$betad[i])))) %>%
    # Combine all together to get overall size distribution in next time step
    mutate(p.size.cur = pred.surv * p.grow.size) %>%
    select(-c(pred.surv, pred.grow.mean, p.grow.size)) %>%
    mutate(perturb.param = 'grow.int')
  
  # 2: Growth model slope
  
  this.boot[[2]] = bootstrap.backbone %>%
    # Predicted survival
    mutate(
      pred.surv = predict(
        newdata = .,
        object = s_s, type = 'response',
        newparams = surv.boots[i,-(1:2)],
        re.form = ~ 0, allow.new.levels = TRUE
      )
    ) %>%
    # Predicted growth
    mutate(
      pred.grow.mean = predict(
        newdata = .,
        object = g_st.ty, type = 'response',
        # apply the perturbation
        newparams = grow.boots[i,-(1:2)] %>%
          (function(x) {
            x[2] <- x[2] + delta
            return(x)
          }),
        re.form = ~ 0, allow.new.levels = TRUE
      )
    ) %>%
    mutate(p.grow.size = 0.1 * dnorm(size.cur, pred.grow.mean, sd = sqrt(exp(grow.boots$betad[i])))) %>%
    # Combine all together to get overall size distribution in next time step
    mutate(p.size.cur = pred.surv * p.grow.size) %>%
    select(-c(pred.surv, pred.grow.mean, p.grow.size)) %>%
    mutate(perturb.param = 'grow.slope')
  
  gs.pert.boot[[i]] = do.call(rbind, this.boot) %>% mutate(boot = paste0('b', i))
  
  print(i)
  
}

# Export perturbed boostrap (first converting to wide format to make the file
# smaller)

do.call(rbind, gs.pert.boot) %>%
  pivot_wider(names_from = boot, values_from = p.size.cur) %>%
  write.csv(
    '03_construct_kernels/out/deterministic_growsurv_perturb_bootstraps.csv',
    row.names = FALSE, na = ''
  )


# Export perturbed parameters (for vital rate differences in LTRE)

cbind(
  boot = 1:n.straps,
  grow.int_control = unlist(grow.boots[,-(1:2)][1]),
  grow.int_drought = unlist(grow.boots[,-(1:2)][1] + grow.boots[,-(1:2)][3]),
  grow.int_irrigated = unlist(grow.boots[,-(1:2)][1] + grow.boots[,-(1:2)][4]),
  grow.slope_control = unlist(grow.boots[,-(1:2)][2]),
  grow.slope_drought = unlist(grow.boots[,-(1:2)][2] + grow.boots[,-(1:2)][5]),
  grow.slope_irrigated = unlist(grow.boots[,-(1:2)][2] + grow.boots[,-(1:2)][6])
) %>%
  write.csv(
    file = '03_construct_kernels/out/growsurv_bootstrapped_perturbed_params.csv',
    row.names = FALSE, na = ''
  )
