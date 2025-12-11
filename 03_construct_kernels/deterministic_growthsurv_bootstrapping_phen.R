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
library(purrr)

rm(list = ls())

source('03_construct_kernels/prepare_demo_data_growsurv.R')

# Read in bootstrapped phen estimates
# done in script `phenology_bootstrapped_estimares.R`
phen.boot.trt = read.csv('03_construct_kernels/out/phenology_bootstrapped_means.csv') %>%
  # Convert to long data frame
  pivot_longer(-boot, names_to = 'trt', values_to = 'mean.phen')

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

# Growth-phen model
g_phen = glmmTMB(
  size.cur ~ size.prev * prev.year + trt * prev.year + phen.c + (1 | Plot / plantid),
  data = demo.grow %>% filter(!is.na(phen.c))
)

# Needed aux effects

# --- Boots ---------------------------------------------------------

cat('Bootstrapping parameters for growth+survival subkernels...\n')

set.seed(11225)

# Number of bootstraps
n.straps = 500

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
  map(
    function(df) {
      glmmTMB(
        formula = surv ~ size.prev + (1 | Plot),
        family = 'binomial',
        data = df
      ) %>%
        # Collect model parameters (for making predictions)
        (function(mod) mod$fit$par)
    },
    .progress = FALSE
  ) %>%
  # mclapply(
  #   function(df) {
  #     glmmTMB(
  #       formula = surv ~ size.prev + (1 | Plot),
  #       family = 'binomial',
  #       data = df
  #     ) %>%
  #       # Collect model parameters (for making predictions)
  #       (function(mod) mod$fit$par)
  #   },
  #   mc.cores = 6
  # ) %>%
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
  map(
    function(df) {
      glmmTMB(
        size.cur ~ size.prev + size.prev * trt + (1 | prev.year) + (1 | prev.year:trt) + (1 | Plot / plantid),
        data = df
      ) %>%
        # Collect model parameters (for making predictions)
        (function(mod) mod$fit$par)
    },
    .progress = FALSE
  ) %>%
  # mclapply(
  #   function(df) {
  #     glmmTMB(
  #       size.cur ~ size.prev + size.prev * trt + (1 | prev.year) + (1 | prev.year:trt) + (1 | Plot / plantid),
  #       data = df
  #     ) %>%
  #       # Collect model parameters (for making predictions)
  #       (function(mod) mod$fit$par)
  #   },
  #   mc.cores = 6
  # ) %>%
  # Combine these together into a single data frame
  do.call(rbind, .) %>%
  data.frame() %>%
  mutate(mod = 'grow', i = 1:n.straps) %>%
  select(mod, i, everything())

# Bootstrapped growth with phenology model
phen.effect.boots = demo.grow %>% 
  # Subset to only plants with phenology records
  filter(!is.na(phen.c)) %>%
  # Do bootstrapped resampling as before:
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
  map(
    function(df) {
      glmmTMB(
        size.cur ~ size.prev * prev.year + trt * prev.year + phen.c + (1 | Plot / plantid),
        data = df
      ) %>%
        # Collect *only model parameters 7 and 14* 
        # 7 is the phen effect (phen intercept)
        # 14 is the log of the sqrt of the residual variance
        (function(mod) mod$fit$par[c(7, 14)])
    },
    .progress = FALSE
  ) %>%
  # mclapply(
  #   function(df) {
  #     glmmTMB(
  #       size.cur ~ size.prev * prev.year + trt * prev.year + phen.c + (1 | Plot / plantid),
  #       data = df
  #     ) %>%
  #       # Collect *only model parameters 7 and 14* 
  #       # 7 is the phen effect (phen intercept)
  #       # 14 is the log of the sqrt of the residual variance
  #       (function(mod) mod$fit$par[c(7, 14)])
  #   },
  #   mc.cores = 6
  # ) %>%
  # Combine these together into a single data frame
  do.call(rbind, .) %>%
  data.frame() %>%
  mutate(mod = 'grow', i = 1:n.straps) %>%
  select(mod, i, everything())

### Un-biasing the bootstrapped estimates
# Manually shifting individual columns of the bootstrap to have the mean of
# boostrapped samples be identical to the parameter values of the true model

# Survival bootstraps
surv.boots[,-(1:2)] = surv.boots[,-(1:2)] + matrix(
  (s_s$fit$par - colMeans(surv.boots[,-(1:2)])), 
  nrow = n.straps, ncol = length(s_s$fit$par), byrow = TRUE,
)

# Growth bootstraps, all data/vegetative plants
grow.boots[,-(1:2)] = grow.boots[,-(1:2)] + matrix(
  (g_st.ty$fit$par - colMeans(grow.boots[,-(1:2)])), 
  nrow = n.straps, ncol = length(g_st.ty$fit$par), byrow = TRUE,
)

# Growth-phen per-day effect
phen.effect.boots[,-(1:2)] = phen.effect.boots[,-(1:2)] + matrix(
  (g_phen$fit$par[c(7, 14)] - colMeans(phen.effect.boots[,-(1:2)])), 
  nrow = n.straps, ncol = length(g_phen$fit$par[c(7, 14)]), byrow = TRUE,
)

# # very small differences, all numerical rounding
# mean((colMeans(surv.boots[,-(1:2)]) - s_s$fit$par)^2)
# mean((colMeans(grow.boots[,-(1:2)]) - g_st.ty$fit$par)^2)
# mean((colMeans(phen.effect.boots[,-(1:2)]) - g_phen$fit$par[c(7, 14)])^2)

### Writing to csvs

write.csv(
  surv.boots, na = '', row.names = FALSE,
  '03_construct_kernels/bootstrapped_model_coefs/surv_boot_coefs.csv'
)

write.csv(
  grow.boots, na = '', row.names = FALSE,
  '03_construct_kernels/bootstrapped_model_coefs/grow_vegt_boot_coefs.csv'
)

write.csv(
  phen.effect.boots, na = '', row.names = FALSE,
  '03_construct_kernels/bootstrapped_model_coefs/grow_phen_boot_coefs.csv'
)

# surv.boots = read.csv('03_construct_kernels/bootstrapped_model_coefs/surv_boot_coefs.csv')
# grow.boots = read.csv('03_construct_kernels/bootstrapped_model_coefs/grow_vegt_boot_coefs.csv')
# phen.effect.boots = read.csv('03_construct_kernels/bootstrapped_model_coefs/grow_phen_boot_coefs.csv')

cat('Done.\n')

# --- Get bootstrapped kernels *for all phenology* -----------------------------

# Data frame for making predictions once every week over growing season
bootstrap.full.backbone = expand.grid(
  size.prev = (5:60)/10,
  size.cur = (5:60)/10,
  phen.c = (-2:2) * 7,
  trt = c('control', 'drought', 'irrigated')
)

# List for storing each bootstrapped sample in
boots.full.list = vector('list', length = n.straps)

cat('Fitting bootstrapped growth+survival subkernels...\n')

# Do kernel estimation on each bootstrapped set of parameters

# for (i in 1:n.straps) {
#   
#   boots.full.list[[i]] = bootstrap.full.backbone %>%
#     # Predicted survival
#     mutate(
#       pred.surv = predict(
#         newdata = .,
#         object = s_s, type = 'response',
#         # survival prediction made with ith bootstrap parameter set
#         newparams = surv.boots[i, -(1:2)],
#         re.form = ~ 0, allow.new.levels = TRUE
#       )
#     ) %>%
#     # Predicted growth
#     mutate(
#       # Growth without phenoloyg (applied to non-flowering plants)
#       pred.grow.mean = predict(
#         newdata = .,
#         object = g_st.ty, type = 'response',
#         # growth prediction made with ith bootstrap parameter set
#         newparams = grow.boots[i, -(1:2)],
#         re.form = ~ 0, allow.new.levels = TRUE
#       )
#     ) %>%
#     # Model with phenology
#     # Need to change name of year column to get annual predictions
#     mutate(phen.grow.mean = pred.grow.mean + phen.effect.boots$beta[i] * phen.c) %>%
#     # Predicted distribution of sizes in next time step
#     # OLD # note: 'betad' parameter is the log of the model's estimated residual variance term
#     # OLD mutate(p.grow.size = 0.1 * dnorm(size.cur, mean = pred.grow.mean, sd = sqrt(exp(grow.boots$betad[i])))) %>%
#     # OLD # Combine all together to get overall size distribution in next time step
#     # OLD # mutate(p.size.cur = pred.surv * p.grow.size) %>%
#     # Predicted distribution of sizes in next time step
#     mutate(
#       pv.grow.size = 0.1 * dnorm(size.cur, mean = pred.grow.mean, sd = sqrt(exp(grow.boots$betad[i]))),
#       pf.grow.size = 0.1 * dnorm(size.cur, mean = phen.grow.mean, sd = sqrt(exp(grow.boots$betad[i])))
#     ) %>%
#     # Re-center phenology
#     # mutate(phen = phen.c + phen.ctrl.mean) %>%
#     # Remove unnecessary columns
#     select(-c(phen.grow.mean, pred.grow.mean)) %>%
#     # Label bootstrap number
#     mutate(boot = paste0('b', i))
#   
#   print(i)
#   
# }

boots.full.list = map(
  1:n.straps,
  \(i) bootstrap.full.backbone %>%
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
      # Growth without phenoloyg (applied to non-flowering plants)
      pred.grow.mean = predict(
        newdata = .,
        object = g_st.ty, type = 'response',
        # growth prediction made with ith bootstrap parameter set
        newparams = grow.boots[i, -(1:2)],
        re.form = ~ 0, allow.new.levels = TRUE
      )
    ) %>%
    # Model with phenology
    # Need to change name of year column to get annual predictions
    mutate(phen.grow.mean = pred.grow.mean + phen.effect.boots$beta[i] * phen.c) %>%
    # Predicted distribution of sizes in next time step
    mutate(
      pv.grow.size = 0.1 * dnorm(size.cur, mean = pred.grow.mean, sd = exp(grow.boots$betadisp[i])),
      pf.grow.size = 0.1 * dnorm(size.cur, mean = phen.grow.mean, sd = exp(grow.boots$betadisp[i]))
    ) %>%
    # Re-center phenology
    # mutate(phen = phen.c + phen.ctrl.mean) %>%
    # Remove unnecessary columns
    select(-c(phen.grow.mean, pred.grow.mean)) %>%
    # Label bootstrap number
    mutate(boot = paste0('b', i)),
  .progress = FALSE
)

# Combine kernels, convert to wide form, and export

do.call(rbind, boots.full.list) %>%
  # pivot_wider(names_from = boot, values_from = p.size.cur) %>%
  write.csv(
    '03_construct_kernels/out/deterministic_growsurv_bootstrap_allphen.csv',
    row.names = FALSE, na = ''
  )

# --- Get bootstrapped kernels *for LTRE phenology* -----------------------------

# Generate backbone, reading in a script with bootstrapped treatment means
# (estimated in an upstream file)
boot.ltre.backbone = expand.grid(
  size.prev = (5:60)/10,
  size.cur = (5:60)/10,
  boot = 1:n.straps,
  # This column will be used for manipulating the phenology date and the vital
  # rate estimation
  trt.phen.idx = 1:7
) %>%
  # Read in LTRE design information
  merge(read.csv('03_construct_kernels/ltre_treatment_key.csv')) %>%
  # Now merge in to get the buddates for the trt.phen column
  merge(phen.boot.trt, by.x = c('trt.phen', 'boot'), by.y = c('trt', 'boot')) %>%
  rename(trt = trt.rate) %>%
  mutate(phen.c = mean.phen - phen.ctrl.mean) %>%
  select(-c(mean.phen, trt.phen))

# List for storing each bootstrapped sample in
boots.ltre.list = vector('list', length = n.straps)

# Do kernel estimation on each bootstrapped set of parameters

# for (i in 1:n.straps) {
#   
#   boots.ltre.list[[i]] = boot.ltre.backbone %>%
#     # Give us only bootstrap rep (including boot phenology) i
#     filter(boot %in% i) %>%
#     # Predicted survival
#     mutate(
#       pred.surv = predict(
#         newdata = .,
#         object = s_s, type = 'response',
#         # survival prediction made with ith bootstrap parameter set
#         newparams = surv.boots[i, -(1:2)],
#         re.form = ~ 0, allow.new.levels = TRUE
#       )
#     ) %>%
#     # Predicted growth
#     mutate(
#       # Growth without phenoloyg (applied to non-flowering plants)
#       pred.grow.mean = predict(
#         newdata = .,
#         object = g_st.ty, type = 'response',
#         # growth prediction made with ith bootstrap parameter set
#         newparams = grow.boots[i, -(1:2)],
#         re.form = ~ 0, allow.new.levels = TRUE
#       )
#     ) %>%
#     # Model with phenology
#     # Need to change name of year column to get annual predictions
#     mutate(phen.grow.mean = pred.grow.mean + phen.effect.boots$beta[i] * phen.c) %>%
#     # Predicted distribution of sizes in next time step
#     # OLD # note: 'betad' parameter is the log of the model's estimated residual variance term
#     # OLD mutate(p.grow.size = 0.1 * dnorm(size.cur, mean = pred.grow.mean, sd = sqrt(exp(grow.boots$betad[i])))) %>%
#     # OLD # Combine all together to get overall size distribution in next time step
#     # OLD # mutate(p.size.cur = pred.surv * p.grow.size) %>%
#     # Predicted distribution of sizes in next time step
#     mutate(
#       pv.grow.size = 0.1 * dnorm(size.cur, mean = pred.grow.mean, sd = sqrt(exp(grow.boots$betad[i]))),
#       pf.grow.size = 0.1 * dnorm(size.cur, mean = phen.grow.mean, sd = sqrt(exp(grow.boots$betad[i])))
#     ) %>%
#     # Re-center phenology
#     mutate(phen = phen.c + phen.ctrl.mean) %>%
#     # Remove unnecessary columns
#     select(-c(phen.c, phen.grow.mean, pred.grow.mean, trt, phen))
#   
#   print(i)
#   
# }

boots.ltre.list = map(
  1:n.straps,
  \(i) boot.ltre.backbone %>%
    # Give us only bootstrap rep (including boot phenology) i
    filter(boot %in% i) %>%
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
      # Growth without phenoloyg (applied to non-flowering plants)
      pred.grow.mean = predict(
        newdata = .,
        object = g_st.ty, type = 'response',
        # growth prediction made with ith bootstrap parameter set
        newparams = grow.boots[i, -(1:2)],
        re.form = ~ 0, allow.new.levels = TRUE
      )
    ) %>%
    # Model with phenology
    # Need to change name of year column to get annual predictions
    mutate(phen.grow.mean = pred.grow.mean + phen.effect.boots$beta[i] * phen.c) %>%
    # Predicted distribution of sizes in next time step
    # OLD # note: 'betad' parameter is the log of the model's estimated residual variance term
    # OLD mutate(p.grow.size = 0.1 * dnorm(size.cur, mean = pred.grow.mean, sd = sqrt(exp(grow.boots$betad[i])))) %>%
    # OLD # Combine all together to get overall size distribution in next time step
    # OLD # mutate(p.size.cur = pred.surv * p.grow.size) %>%
    # Predicted distribution of sizes in next time step
    mutate(
      pv.grow.size = 0.1 * dnorm(size.cur, mean = pred.grow.mean, sd = exp(grow.boots$betadisp[i])),
      pf.grow.size = 0.1 * dnorm(size.cur, mean = phen.grow.mean, sd = exp(grow.boots$betadisp[i]))
    ) %>%
    # Re-center phenology
    mutate(phen = phen.c + phen.ctrl.mean) %>%
    # Remove unnecessary columns
    select(-c(phen.c, phen.grow.mean, pred.grow.mean, trt, phen)),
  .progress = FALSE
)

# Export
do.call(rbind, boots.ltre.list) %>%
  # Don't do the pivoting yet...
  # pivot_wider(names_from = boot, values_from = c(pred.surv, pv.grow.size, pf.grow.size)) %>%
  write.csv(
    file = '03_construct_kernels/out/deterministic_growsurv_bootstrap_ltre.csv',
    na = '', row.names = FALSE
  )

cat('Done.\n')

# ------
# ------ Perturbation bootstrapping
# ------ (repeat above procedure but with perturbations at each vital rate)
# ------

# Parameters to perturb:
# - Survival: nothing (no treatment effects in any terms)
# - Growth:
#   - Intercept
#   - Slope
# - Growth-phenology:
#   - Intercept

cat('Perturbing bootstrapped growth+survival subkernels...\n')

# Perturbation amount
delta = 0.001

# Get a list for outputs
gs.pert.boot = vector('list', n.straps)

# Loop through each boot sample
for (i in 1:n.straps) {
  
  # Set up list for outputs
  this.boot = vector('list', 3)
  
  # 1: Growth model intercept
  
  this.boot[[1]] = boot.ltre.backbone %>%
    # Get relevant phenology dates only
    filter(boot %in% i) %>%
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
    # Model with phenology
    mutate(phen.grow.mean = pred.grow.mean + phen.effect.boots$beta[i] * phen.c) %>%
    # Take the average of the growth kernel across years
    # group_by(size.prev, size.cur, trt, trt.phen, phen.c, pred.surv, pred.grow.mean) %>%
    # summarise(phen.grow.mean = mean(phen.grow.mean)) %>%
    # ungroup() %>%
    # Predicted distribution of sizes in next time step
    mutate(
      pv.grow.size = 0.1 * dnorm(size.cur, mean = pred.grow.mean, sd = exp(grow.boots$betadisp[i])),
      pf.grow.size = 0.1 * dnorm(size.cur, mean = phen.grow.mean, sd = exp(grow.boots$betadisp[i]))
    ) %>%
    # Add in perturbation information
    mutate(perturb.param = 'grow.int') %>%
    # Remove unneeded columns
    select(-c(trt, phen.c, pred.grow.mean, phen.grow.mean))
  
  # 2: Growth model slope
  
  this.boot[[2]] = boot.ltre.backbone %>%
    # Get relevant phenology dates only
    filter(boot %in% i) %>%
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
    # Model with phenology
    mutate(phen.grow.mean = pred.grow.mean + phen.effect.boots$beta[i] * phen.c) %>%
    # Take the average of the growth kernel across years
    # group_by(size.prev, size.cur, trt, trt.phen, phen.c, pred.surv, pred.grow.mean) %>%
    # summarise(phen.grow.mean = mean(phen.grow.mean)) %>%
    # ungroup() %>%
    # Predicted distribution of sizes in next time step
    mutate(
      pv.grow.size = 0.1 * dnorm(size.cur, mean = pred.grow.mean, sd = exp(grow.boots$betadisp[i])),
      pf.grow.size = 0.1 * dnorm(size.cur, mean = phen.grow.mean, sd = exp(grow.boots$betadisp[i]))
    ) %>%
    # Add in perturbation information
    mutate(perturb.param = 'grow.slope') %>%
    # Remove unneeded columns
    select(-c(trt, phen.c, pred.grow.mean, phen.grow.mean))
  
  # 3: growth model phenology effect
  
  this.boot[[3]] = boot.ltre.backbone %>%
    # Get relevant phenology dates only
    filter(boot %in% i) %>%
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
        newparams = grow.boots[i,-(1:2)],
        re.form = ~ 0, allow.new.levels = TRUE
      )
    ) %>%
    # Model with phenology
    mutate(phen.grow.mean = pred.grow.mean + delta + (phen.effect.boots$beta[i] * phen.c)) %>%
    # Take the average of the growth kernel across years
    # group_by(size.prev, size.cur, trt, trt.phen, phen.c, pred.surv, pred.grow.mean) %>%
    # summarise(phen.grow.mean = mean(phen.grow.mean)) %>%
    # ungroup() %>%
    # Predicted distribution of sizes in next time step
    mutate(
      pv.grow.size = 0.1 * dnorm(size.cur, mean = pred.grow.mean, sd = exp(grow.boots$betadisp[i])),
      pf.grow.size = 0.1 * dnorm(size.cur, mean = phen.grow.mean, sd = exp(grow.boots$betadisp[i]))
    ) %>%
    # Add in perturbation information
    mutate(perturb.param = 'phen.grow') %>%
    # Remove unneeded columns
    select(-c(trt, phen.c, pred.grow.mean, phen.grow.mean))
  
  
  gs.pert.boot[[i]] = do.call(rbind, this.boot) %>% mutate(boot = i)
  
  # print(i)
  
}

# Export perturbed boostrapped kernels

do.call(rbind, gs.pert.boot) %>%
  # pivot_wider(names_from = boot, values_from = p.size.cur) %>%
  write.csv(
    '03_construct_kernels/out/deterministic_growsurv_perturb_bootstraps.csv',
    row.names = FALSE, na = ''
  )


# Export perturbed parameters (for vital rate differences in LTRE)

# Create a wide and centered version of the phenology bootstraps for indexing
phen.boot.trt.export = phen.boot.trt %>%
  mutate(mean.phen = mean.phen - phen.ctrl.mean) %>%
  pivot_wider(names_from = trt, values_from = mean.phen)

cbind(
  boot = 1:n.straps,
  grow.int_control = unlist(grow.boots[,-(1:2)][1]),
  grow.int_drought = unlist(grow.boots[,-(1:2)][1] + grow.boots[,-(1:2)][3]),
  grow.int_irrigated = unlist(grow.boots[,-(1:2)][1] + grow.boots[,-(1:2)][4]),
  grow.slope_control = unlist(grow.boots[,-(1:2)][2]),
  grow.slope_drought = unlist(grow.boots[,-(1:2)][2] + grow.boots[,-(1:2)][5]),
  grow.slope_irrigated = unlist(grow.boots[,-(1:2)][2] + grow.boots[,-(1:2)][6]),
  phen.grow_control = unlist(grow.boots[,-(1:2)][1] + phen.effect.boots$beta * phen.boot.trt.export$control),
  phen.grow_drought = unlist(grow.boots[,-(1:2)][1] + phen.effect.boots$beta * phen.boot.trt.export$drought),
  phen.grow_irrigated = unlist(grow.boots[,-(1:2)][1] + phen.effect.boots$beta * phen.boot.trt.export$irrigated)
) %>%
  write.csv(
    file = '03_construct_kernels/out/deterministic_growsurv_bootstrapped_perturbed_params.csv',
    row.names = FALSE, na = ''
  )

cat('Done.\n')
