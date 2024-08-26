# Script for generating bootstrapped survival-growth subkernels, including
# perturbations
# Outputs:
# - 100 bootstrapped reproduction sub-kernels
# - 100 bootstrapped perturbed reproduction sub-kernels at weekly-intervals
# - 100 parameter differences used to generate sub-kernels at observed treatment means
# NOTE: the bootstrapped predictions take some time to compute (because they are
# not parallelized and because they are done for multiple phenology levels) -
# each of these loops takes ~20 minutes on my computer

# --- Setup ---------------------------------------------------------

library(ggplot2)
library(dplyr)
library(tidyr)
library(parallel)
library(glmmTMB)

# Get rid of this super annoying feature
options(dplyr.summarise.inform = FALSE)

rm(list = ls())

# Run wrapper script to prepare demo + phen data data
source('03_construct_kernels/prepare_demo_data_repr.R')

# --- Fit original mods -----------------------------------------------
# # These are needed for either getting treatment phen (bud dates) and for
# centering bootstrap estimates

# === Phenology model ===
# Response: day of umbel budding (continuous julian date)
# Predictors: year (factor/categorical) and treatment (factor)
d_t = glmmTMB(
  phen.julian ~ trt + Year + (1 | Plot / plantid),
  data = phen
)

# === Flowering/umbel count model ===
u_s_s.ty = glmmTMB(
  No.umbels ~ size + (1 | Year) + (1 | Plot / plantid),
  family = 'truncated_poisson',
  ziformula = ~ size + trt + (1 | Year) + (1 | Year:trt) + (1 | Plot / plantid),
  data = demo.flow
)

# === Seed model ===
# Response: number of seeds (negative binomial distribution) of an umbel
# Predictors: treatment (categorical), year (factor) , mean budding date of
# plant (centered, continuous), plant size (continuous)
s_st.p_s.u.p2 = glmmTMB(
  no.seeds ~ trt * size + Year + phen.c + (1 | Plot / plantid),
  family = 'nbinom2',
  ziformula = ~ size + phen.umbels + Year + poly(phen.c, 2) + (1 | Plot / plantid),
  data = seed
)

# === Recruit size model ===
r_t.y = glmmTMB(size ~ trt + (1 | Year) + (1 | Plot), data = demo.recr)

# --- Boots ---------------------------------------------------------
# Generate bootstrapped parameter sets for each set of parameters
# Because the prob flower + umbel number, and prob. umbel success + seed set
# models are combined, these parameter estimates are estimated in tandem (i.e.,
# joint bootstrapped parameter estimates with both the zero-inflation terms and
# conditional model terms estimated in on the *same* bootstrapped sample)

set.seed(940820)

# Number of bootstraps
n.straps = 100

# --- Flowering and umbel production model bootstrap

flow.numb.boot = demo.flow %>%
  # Do bootstrapped resampling
  # copy the data frame (each row duplicated)
  uncount(weights = n.straps) %>%
  # label these entries with a `samp` (sample) column to delineate different
  # bootstrap samples
  group_by(plantid, Year) %>%
  mutate(samp = 1:n.straps) %>%
  # Perform the resampling, preserving plot and year structure
  group_by(Plot, Year, samp) %>%
  sample_n(size = n(), replace = TRUE) %>%
  ungroup() %>%
  # Split the dataset by each sample and re-fit the flowering/umbel count model
  split(.$samp) %>%
  mclapply(
    function(df) {
      glmmTMB(
        No.umbels ~ size + (1 | Year) + (1 | Plot / plantid),
        family = 'truncated_poisson',
        ziformula = ~ size + trt + (1 | Year) + (1 | Year:trt) + (1 | Plot / plantid),
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
  mutate(mod = 'flow.numb', i = 1:n.straps) %>%
  select(mod, i, everything())

# --- Umbel success and seed production

succ.seed.boot = seed %>%
  # Do bootstrapped resampling
  # copy the data frame (each row duplicated)
  # Doing an extra step in here - adding an obs.no column to do proper sample
  # labeling (this is because in the seed/umbel dataset there are multiple rows
  # per individual plant, but individual umbels are not uniquely identified)
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
        no.seeds ~ trt * size + Year + phen.c + (1 | Plot / plantid),
        family = 'nbinom2',
        ziformula = ~ size + phen.umbels + Year + poly(phen.c, 2) + (1 | Plot / plantid),
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
  mutate(mod = 'succ.seed', i = 1:n.straps) %>%
  select(mod, i, everything())

# --- Recruit size

recr.boot = demo.recr %>%
  # Do bootstrapped resampling
  # copy the data frame (each row duplicated)
  uncount(weights = n.straps) %>%
  # Label these entries with a `samp` (sample) column to delineate different
  # bootstrap samples
  group_by(plantid, Year) %>%
  mutate(samp = 1:n.straps) %>%
  # Perform the resampling, preserving plot and survival structure
  group_by(Plot, Year, samp) %>%
  sample_n(size = n(), replace = TRUE) %>%
  ungroup() %>%
  # Split the dataset by each sample and re-fit the recruit size
  split(.$samp) %>%
  mclapply(
    function(df) {
      glmmTMB(size ~ trt + (1 | Year) + (1 | Plot), data = df) %>%
        # extract model parameters
        (function(mod) mod$fit$par)
    },
    mc.cores = 6
  ) %>%
  # Combine into a single data frame
  do.call(rbind, .) %>%
  data.frame() %>%
  mutate(mod = 'recr', i = 1:n.straps) %>%
  select(mod, i, everything())

# Phenology bootstrapping

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

# # Un-biasing the bootstrapped estimates
# Manually shifting individual columns of the bootstrap to have the mean of
# boostrapped samples be identical to the parameter values of the true model

# Flower/umbel bootstraps
flow.numb.boot[,-(1:2)] = flow.numb.boot[,-(1:2)] + matrix(
  (u_s_s.ty$fit$par - colMeans(flow.numb.boot[,-(1:2)])), 
  nrow = n.straps, ncol = length(u_s_s.ty$fit$par), byrow = TRUE,
)

# Success/seed bootstraps
succ.seed.boot[,-(1:2)] = succ.seed.boot[,-(1:2)] + matrix(
  (s_st.p_s.u.p2$fit$par - colMeans(succ.seed.boot[,-(1:2)])), 
  nrow = n.straps, ncol = length(s_st.p_s.u.p2$fit$par), byrow = TRUE,
)

# Recruit size bootstraps
recr.boot[,-(1:2)] = recr.boot[-(1:2)] + matrix(
  (r_t.y$fit$par - colMeans(recr.boot[,-(1:2)])),
  nrow = n.straps, ncol = length(r_t.y$fit$par), byrow = TRUE
)

# Phenology estimate bootstraps
phen.boot[,-(1:2)] = phen.boot[-(1:2)] + matrix(
  (d_t$fit$par - colMeans(phen.boot[,-(1:2)])),
  nrow = n.straps, ncol = length(d_t$fit$par), byrow = TRUE
)


# very small differences, all numerical rounding
(colMeans(flow.numb.boot[,-(1:2)]) - u_s_s.ty$fit$par)
(colMeans(succ.seed.boot[,-(1:2)]) - s_st.p_s.u.p2$fit$par)
(colMeans(recr.boot[,-(1:2)]) - r_t.y$fit$par)
(colMeans(phen.boot[,-(1:2)]) - d_t$fit$par)

# --- Full-phenology kernel --------------------------------------------------

# # Doing 100 boot straps for each day is computationally way too expensive
# So instead doing one bootstrap per week
# # This will be used for getting uncertainty in the lambda-vs-phenology script

bootstrap.full.backbone = expand.grid(
  size = (5:60)/10,
  size.nex = (5:60)/10,
  trt = c('control', 'drought', 'irrigated'),
  phen.c = (-4:4) * 7,
  Year = 2021:2024
)

# Designate an output list 
boots.full.list = vector('list', length = n.straps)

for (i in 1:n.straps) {
  
  boots.full.list[[i]] = bootstrap.full.backbone %>%
    # Rename to not put the year random effect in these predictions
    rename(year = Year) %>%
    mutate(
      # Umbel count
      phen.umbels = predict(
        u_s_s.ty, newdata = ., 
        # use bootstrapped parameters
        newparams = flow.numb.boot[i, -(1:2)],
        allow.new.levels = TRUE, re.form = ~ 0, type = 'response'
      )
    ) %>%
    rename(Year = year) %>%
    mutate(
      # Model predictions at treatment means
      # (doing this on linear scale for each for easier averaging)
      seeds.zinf.linear =  predict(
        s_st.p_s.u.p2, newdata = ., allow.new.levels = TRUE, re.form = ~ 0,
        # use bootstrapped parameters
        newparams = succ.seed.boot[i,-(1:2)],
        type = 'zlink'
      ),
      seeds.seed.linear =  predict(
        s_st.p_s.u.p2, newdata = ., allow.new.levels = TRUE, re.form = ~ 0,
        # use bootstrapped parameters
        newparams = succ.seed.boot[i,-(1:2)],
        type = 'link'
      ),
    ) %>%
    # Taking out the size-zinf terms for 2021 - very extrapolatory, affects averages too much
    mutate(seeds.zinf.linear = ifelse(Year %in% 2021, NA, seeds.zinf.linear),) %>%
    group_by(size, size.nex, trt, phen.c, phen.umbels) %>%
    summarise(across(c(seeds.zinf.linear, seeds.seed.linear), ~ mean(.x, na.rm = TRUE))) %>%
    ungroup() %>%
    mutate(
      # Transform from linear scale to response scale
      seeds.per.umbel = (1 / (1 + exp(-seeds.zinf.linear))) * exp(seeds.seed.linear),
      # Total umbels per plant
      seeds.total = seeds.per.umbel * phen.umbels
    ) %>%
    select(-c(seeds.zinf.linear, seeds.seed.linear)) %>%
    mutate(
      # Mean recruit size
      recr.mean = predict(
        r_t.y, allow.new.levels = TRUE, re.form = ~ 0, newdata = .,
        # use bootstrapped parameters
        newparams = recr.boot[i, -(1:2)]
      ),
      # Get the number of seeds produced for each size grouping
      # (note: 'betad' parameter here is the log of the residual variance from the model fit)
      p.size.cur = 0.1 * seeds.total * dnorm(x = size.nex, mean = recr.mean, sd = sqrt(exp(recr.boot$betad[i])))
    ) %>%
    # Rename column
    rename(size.prev = size) %>%
    # Re-center phenology column
    mutate(mean.phen = phen.c + round(mean(seed$mean.phen))) %>%
    # Add a column for distinguishing bootstrapped samples
    mutate(boot = paste0('b', i)) %>%
    # Remove unnecessary columns (save space)
    select(-c(phen.umbels, seeds.per.umbel, seeds.total, recr.mean, phen.c))
  
    print(i)
  
}

# Combine into one data frame and export (pivot to wider so the file takes up
# less space)
do.call(rbind, boots.full.list) %>%
  pivot_wider(names_from = boot, values_from = p.size.cur) %>%
  write.csv(
    '03_construct_kernels/out/deterministic_reprod_bootstrap_allphen.csv',
    row.names = FALSE, na = ''
  )

# --- Treatment-phenology kernel -----------------------------------------------

# Here, performing bootstrapping only at the mean date for each treatment
# (not doing bootstrapping to incorporate the uncertainty in the phenology estimates, yet)
# These will be used for the LTRE analysis

# Get the mean bud dates for each treatment
trt.mean.buddates = expand.grid(trt = c('control', 'drought', 'irrigated'), Year = factor(2021:2024)) %>%
  mutate(
    mean.bud = predict(
      d_t, re.form = ~ 0, allow.new.levels = TRUE,
      newdata = expand.grid(trt = c('control', 'drought', 'irrigated'), Year = factor(2021:2024))
    )
  ) %>%
  group_by(trt.phen = trt) %>%
  summarise(mean.phen = mean(mean.bud))

# Boostrap out 

phen.backbone = expand.grid(
  trt = c('control', 'drought', 'irrigated'),
  Year = factor(2021:2024)
)

# There is probably a sleeker way to do this but I will leave it at this for now:

phen.list.out = vector('list', length = n.straps)

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
  
  # print(i)
}

phen.boot.trt = do.call(rbind, phen.list.out) %>% mutate(boot = as.numeric(boot))

boot.ltre.backbone = expand.grid(
  size = (5:60)/10,
  size.nex = (5:60)/10,
  Year = 2021:2024,
  boot = 1:n.straps,
  # This column will be used for manipulating the phenology date and the vital
  # rate estimation
  trt.phen.idx = 1:7
) %>%
  merge(
    data.frame(
      trt.phen.idx = 1:7,
      trt.phen = c('drought', 'control', 'irrigated', 'control', 'drought', 'irrigated', 'control'),
      trt.rate = c('drought', 'drought', 'irrigated', 'irrigated', 'control', 'control', 'control')
    )
  ) %>%
  # This column is superfluous now
  select(-trt.phen.idx) %>%
  # Now merge in to get the buddates for the trt.phen column
  merge(phen.boot.trt, by.x = c('trt.phen', 'boot'), by.y = c('trt', 'boot')) %>%
  rename(trt = trt.rate) %>%
  mutate(phen.c = mean.phen - round(mean(seed$mean.phen))) %>%
  select(-mean.phen)

# Use this to get a kernel backbone for each treatment-bud day combo
# bootstrap.ltre.backbone = expand.grid(
#   size = (5:60)/10,
#   size.nex = (5:60)/10,
#   trt = c('control', 'drought', 'irrigated'),
#   mean.phen = trt.mean.buddates$mean.bud,
#   Year = factor(2021:2024)
# ) %>%
#   # Going to get rid of the treatment-phen combos we did not observe
#   # (which will therefore not be included in the LTRE)
#   filter(!(trt %in% 'drought' & floor(mean.phen) == 126)) %>%
#   filter(!(trt %in% 'irrigated' & floor(mean.phen) == 122)) %>%

# bootstrap.ltre.backbone = expand.grid(
#   size = (5:60)/10,
#   size.nex = (5:60)/10,
#   Year = 2021:2024,
#   # This column will be used for manipulating the phenology date and the vital
#   # rate estimation
#   trt.phen.idx = 1:7
# ) %>%
#   merge(
#     data.frame(
#       trt.phen.idx = 1:7,
#       # Treatment associated with the buddates used for umbel success + seeds
#       trt.phen = c('drought', 'control', 'irrigated', 'control', 'drought', 'irrigated', 'control'),
#       # Treatment associated with direct treatment effects on vital rates
#       trt.rate = c('drought', 'drought', 'irrigated', 'irrigated', 'control', 'control', 'control')
#     )
#   ) %>%
#   merge(trt.mean.buddates) %>%
#   # center the phenology column and rename the `trt` column so it can be used in
#   # vital rate estimates
#   mutate(
#     phen.c = mean.phen - mean(round(seed$mean.phen))
#   )

# Start a list for storing the bootstrapped LTREs
boots.ltre.list = vector('list', length = n.straps)

for (i in 1:n.straps) {
  
  boots.ltre.list[[i]] = boot.ltre.backbone %>%
    filter(boot %in% i) %>%
    # Rename to not put the year random effect in these predictions
    rename(year = Year) %>%
    mutate(
      # Umbel count
      phen.umbels = predict(
        u_s_s.ty, newdata = ., 
        # use bootstrapped parameters
        newparams = flow.numb.boot[i, -(1:2)],
        allow.new.levels = TRUE, re.form = ~ 0, type = 'response'
      )
    ) %>%
    rename(Year = year) %>%
    mutate(
      # Model predictions at treatment means
      # (doing this on linear scale for each for easier averaging)
      seeds.zinf.linear =  predict(
        s_st.p_s.u.p2, newdata = ., allow.new.levels = TRUE, re.form = ~ 0,
        # use bootstrapped parameters
        newparams = succ.seed.boot[i,-(1:2)],
        type = 'zlink'
      ),
      seeds.seed.linear =  predict(
        s_st.p_s.u.p2, newdata = ., allow.new.levels = TRUE, re.form = ~ 0,
        # use bootstrapped parameters
        newparams = succ.seed.boot[i,-(1:2)],
        type = 'link'
      ),
    ) %>%
    # Taking out the size-zinf terms for 2021 - very extrapolatory, affects averages too much
    mutate(seeds.zinf.linear = ifelse(Year %in% 2021, NA, seeds.zinf.linear),) %>%
    group_by(size, size.nex, trt, trt.phen, phen.c, phen.umbels) %>%
    summarise(across(c(seeds.zinf.linear, seeds.seed.linear), ~ mean(.x, na.rm = TRUE))) %>%
    ungroup() %>%
    mutate(
      # Transform from linear scale to response scale
      seeds.per.umbel = (1 / (1 + exp(-seeds.zinf.linear))) * exp(seeds.seed.linear),
      # Total umbels per plant
      seeds.total = seeds.per.umbel * phen.umbels
    ) %>%
    select(-c(seeds.zinf.linear, seeds.seed.linear)) %>%
    mutate(
      # Mean recruit size
      recr.mean = predict(
        r_t.y, allow.new.levels = TRUE, re.form = ~ 0, newdata = .,
        # use bootstrapped parameters
        newparams = recr.boot[i, -(1:2)]
      ),
      # Get the number of seeds produced for each size grouping
      # (note: 'betad' parameter here is the log of the residual variance from the model fit)
      p.size.cur = 0.1 * seeds.total * dnorm(x = size.nex, mean = recr.mean, sd = sqrt(exp(recr.boot$betad[i])))
    ) %>%
    # Rename column
    rename(size.prev = size) %>%
    # Add a column for distinguishing bootstrapped samples
    mutate(boot = paste0('b', i)) %>%
    # Remove unnecessary columns
    select(-c(phen.umbels, seeds.per.umbel, seeds.total, recr.mean, phen.c))
    
  print(i)
  
}

# Combine into one data frame and export (pivot wider to save space)
do.call(rbind, boots.ltre.list) %>%
  pivot_wider(names_from = boot, values_from = p.size.cur) %>%
  write.csv(
    file = '03_construct_kernels/out/deterministic_reprod_bootstrap_ltre.csv',
    na = '', row.names = FALSE
  )

# ------ Perturbation bootstrapping --------------------------------------------

# We can use the same kernel backbone for generating these estimates

# Perturbation size
delta = 0.0001

# Make an output list (list of lists, will rbind after loop)
fr.pert.boot = vector('list', n.straps)

for (i in 1:n.straps) {
  
  # Set up a list for storing outputs
  # (will rbind this at the end of the loop to make one data frame per loop
  # iteration, then store that as the ith element of the list)
  this.boot = vector('list', 6)
  
  # Parameter 1: probability of flowering intercept

  this.boot[[1]] = boot.ltre.backbone %>%
    filter(boot %in% i) %>%
    rename(year = Year) %>%
    mutate(
      # Umbel count
      phen.umbels = predict(
        u_s_s.ty, newdata = .,
        newparams = flow.numb.boot[i,-(1:2)] %>%
          # apply the perturbation
          (function(x) {
            x[3] <- x[3] + delta
            return(x)
          }),
        allow.new.levels = TRUE, re.form = ~ 0, type = 'response'
      )
    ) %>%
    rename(Year = year) %>%
    mutate(
      # Model predictions for seed set on linear (link) scale for averaging
      seeds.zinf.linear =  predict(
        s_st.p_s.u.p2, newdata = ., allow.new.levels = TRUE, re.form = ~ 0,
        type = 'zlink',
        newparams = succ.seed.boot[i,-(1:2)]
      ),
      seeds.seed.linear =  predict(
        s_st.p_s.u.p2, newdata = ., allow.new.levels = TRUE, re.form = ~ 0,
        type = 'link',
        newparams = succ.seed.boot[i,-(1:2)]
      ),
    ) %>%
    # Taking out the size-zinf terms for 2021 - very extrapolatory, affects averages too much
    mutate(seeds.zinf.linear = ifelse(Year %in% 2021, NA, seeds.zinf.linear),) %>%
    group_by(size, size.nex, trt, trt.phen, phen.c, phen.umbels) %>%
    summarise(across(c(seeds.zinf.linear, seeds.seed.linear), ~ mean(.x, na.rm = TRUE))) %>%
    ungroup() %>%
    mutate(
      # Transform from linear scale to response scale
      seeds.per.umbel = (1 / (1 + exp(-seeds.zinf.linear))) * exp(seeds.seed.linear),
      # Total umbels per plant
      seeds.total = seeds.per.umbel * phen.umbels
    ) %>%
    select(-c(seeds.zinf.linear, seeds.seed.linear)) %>%
    mutate(
      # Mean recruit size
      recr.mean = predict(
        r_t.y, allow.new.levels = TRUE, re.form = ~ 0, newdata = .,
        newparams = recr.boot[i,-(1:2)]
      ),
      # Get the number of seeds produced for each size grouping
      p.size.cur = 0.1 * seeds.total * dnorm(x = size.nex, mean = recr.mean, sd = sqrt(exp(recr.boot$betad[i])))
    ) %>%
    # Rename column
    rename(size.prev = size) %>%
    select(-c(phen.umbels, phen.c, seeds.per.umbel, seeds.total, recr.mean)) %>%
    mutate(param = 'flow.int')
  
  # Seed set intercept (conditional model)
  
  this.boot[[2]] = boot.ltre.backbone %>%
    filter(boot %in% i) %>%
    rename(year = Year) %>%
    mutate(
      # Umbel count
      phen.umbels = predict(
        u_s_s.ty, newdata = .,
        allow.new.levels = TRUE, re.form = ~ 0, type = 'response',
        newparams = flow.numb.boot[i,-(1:2)]
      )
    ) %>%
    rename(Year = year) %>%
    mutate(
      # Model predictions for seed set on linear (link) scale for averaging
      seeds.zinf.linear =  predict(
        s_st.p_s.u.p2, newdata = ., allow.new.levels = TRUE, re.form = ~ 0,
        type = 'zlink',
        newparams = succ.seed.boot[i,-(1:2)]
      ),
      seeds.seed.linear =  predict(
        s_st.p_s.u.p2, newdata = ., allow.new.levels = TRUE, re.form = ~ 0,
        newparams = succ.seed.boot[i,-(1:2)] %>%
          # apply the perturbation
          (function(x) {
            x[1] <- x[1] + delta
            return(x)
          }),
        type = 'link'
      ),
    ) %>%
    # Taking out the size-zinf terms for 2021 - very extrapolatory, affects averages too much
    mutate(seeds.zinf.linear = ifelse(Year %in% 2021, NA, seeds.zinf.linear),) %>%
    group_by(size, size.nex, trt, trt.phen, phen.c, phen.umbels) %>%
    summarise(across(c(seeds.zinf.linear, seeds.seed.linear), ~ mean(.x, na.rm = TRUE))) %>%
    ungroup() %>%
    mutate(
      # Transform from linear scale to response scale
      seeds.per.umbel = (1 / (1 + exp(-seeds.zinf.linear))) * exp(seeds.seed.linear),
      # Total umbels per plant
      seeds.total = seeds.per.umbel * phen.umbels
    ) %>%
    select(-c(seeds.zinf.linear, seeds.seed.linear)) %>%
    mutate(
      # Mean recruit size
      recr.mean = predict(
        r_t.y, allow.new.levels = TRUE, re.form = ~ 0, newdata = .,
        newparams = recr.boot[i,-(1:2)]
      ),
      # Get the number of seeds produced for each size grouping
      p.size.cur = 0.1 * seeds.total * dnorm(x = size.nex, mean = recr.mean, sd = sqrt(exp(recr.boot$betad[i])))
    ) %>%
    # Rename column
    rename(size.prev = size) %>%
    select(-c(phen.umbels, phen.c, seeds.per.umbel, seeds.total, recr.mean)) %>%
    mutate(param = 'seed.int')
  
  # Seed set slope (conditional model) 
  
  this.boot[[3]] = boot.ltre.backbone %>%
    filter(boot %in% i) %>%
    rename(year = Year) %>%
    mutate(
      # Umbel count
      phen.umbels = predict(
        u_s_s.ty, newdata = .,
        allow.new.levels = TRUE, re.form = ~ 0, type = 'response',
        newparams = flow.numb.boot[i,-(1:2)]
      )
    ) %>%
    rename(Year = year) %>%
    mutate(
      # Model predictions for seed set on linear (link) scale for averaging
      seeds.zinf.linear =  predict(
        s_st.p_s.u.p2, newdata = ., allow.new.levels = TRUE, re.form = ~ 0,
        type = 'zlink',
        newparams = succ.seed.boot[i,-(1:2)]
      ),
      seeds.seed.linear =  predict(
        s_st.p_s.u.p2, newdata = ., allow.new.levels = TRUE, re.form = ~ 0,
        newparams = succ.seed.boot[i,-(1:2)] %>%
          # apply the perturbation
          (function(x) {
            x[4] <- x[4] + delta
            return(x)
          }),
        type = 'link'
      ),
    ) %>%
    # Taking out the size-zinf terms for 2021 - very extrapolatory, affects averages too much
    mutate(seeds.zinf.linear = ifelse(Year %in% 2021, NA, seeds.zinf.linear),) %>%
    group_by(size, size.nex, trt, trt.phen, phen.c, phen.umbels) %>%
    summarise(across(c(seeds.zinf.linear, seeds.seed.linear), ~ mean(.x, na.rm = TRUE))) %>%
    ungroup() %>%
    mutate(
      # Transform from linear scale to response scale
      seeds.per.umbel = (1 / (1 + exp(-seeds.zinf.linear))) * exp(seeds.seed.linear),
      # Total umbels per plant
      seeds.total = seeds.per.umbel * phen.umbels
    ) %>%
    select(-c(seeds.zinf.linear, seeds.seed.linear)) %>%
    mutate(
      # Mean recruit size
      recr.mean = predict(
        r_t.y, allow.new.levels = TRUE, re.form = ~ 0, newdata = .,
        newparams = recr.boot[i,-(1:2)]
      ),
      # Get the number of seeds produced for each size grouping
      p.size.cur = 0.1 * seeds.total * dnorm(x = size.nex, mean = recr.mean, sd = sqrt(exp(recr.boot$betad[i])))
    ) %>%
    # Rename column
    rename(size.prev = size) %>%
    select(-c(phen.umbels, phen.c, seeds.per.umbel, seeds.total, recr.mean)) %>%
    mutate(param = 'seed.slope')
  
  
  # 4: recruit size intercept
  
  this.boot[[4]] = boot.ltre.backbone %>%
    filter(boot %in% i) %>%
    rename(year = Year) %>%
    mutate(
      # Umbel count
      phen.umbels = predict(
        u_s_s.ty, newdata = .,
        allow.new.levels = TRUE, re.form = ~ 0, type = 'response',
        newparams = flow.numb.boot[i,-(1:2)]
      )
    ) %>%
    rename(Year = year) %>%
    mutate(
      # Model predictions for seed set on linear (link) scale for averaging
      seeds.zinf.linear =  predict(
        s_st.p_s.u.p2, newdata = ., allow.new.levels = TRUE, re.form = ~ 0,
        type = 'zlink',
        newparams = succ.seed.boot[i,-(1:2)]
      ),
      seeds.seed.linear =  predict(
        s_st.p_s.u.p2, newdata = ., allow.new.levels = TRUE, re.form = ~ 0,
        type = 'link',
        newparams = succ.seed.boot[i,-(1:2)]
      ),
    ) %>%
    # Taking out the size-zinf terms for 2021 - very extrapolatory, affects averages too much
    mutate(seeds.zinf.linear = ifelse(Year %in% 2021, NA, seeds.zinf.linear),) %>%
    group_by(size, size.nex, trt, trt.phen, phen.c, phen.umbels) %>%
    summarise(across(c(seeds.zinf.linear, seeds.seed.linear), ~ mean(.x, na.rm = TRUE))) %>%
    ungroup() %>%
    mutate(
      # Transform from linear scale to response scale
      seeds.per.umbel = (1 / (1 + exp(-seeds.zinf.linear))) * exp(seeds.seed.linear),
      # Total umbels per plant
      seeds.total = seeds.per.umbel * phen.umbels
    ) %>%
    select(-c(seeds.zinf.linear, seeds.seed.linear)) %>%
    mutate(
      # Mean recruit size
      recr.mean = predict(
        r_t.y, allow.new.levels = TRUE, re.form = ~ 0, newdata = .,
        newparams = recr.boot[i,-(1:2)] %>%
          # apply the perturbation
          (function(x) {
            x[1] <- x[1] + delta
            return(x)
          })
      ),
      # Get the number of seeds produced for each size grouping
      p.size.cur = 0.1 * seeds.total * dnorm(x = size.nex, mean = recr.mean, sd = sqrt(exp(recr.boot$betad[i])))
    ) %>%
    # Rename column
    rename(size.prev = size) %>%
    select(-c(phen.umbels, phen.c, seeds.per.umbel, seeds.total, recr.mean)) %>%
    mutate(param = 'recr.int')
  
  # 5: phenology effects on umbel success
  
  this.boot[[5]] = boot.ltre.backbone %>%
    filter(boot %in% i) %>%
    rename(year = Year) %>%
    mutate(
      # Umbel count
      phen.umbels = predict(
        u_s_s.ty, newdata = .,
        allow.new.levels = TRUE, re.form = ~ 0, type = 'response',
        newparams = flow.numb.boot[i,-(1:2)]
      )
    ) %>%
    rename(Year = year) %>%
    # Perturb phen variable (for only the zinf term)
    mutate(phen.c = phen.c + delta) %>%
    mutate(
      # Model predictions for seed set on linear (link) scale for averaging
      seeds.zinf.linear =  predict(
        s_st.p_s.u.p2, newdata = ., allow.new.levels = TRUE, re.form = ~ 0,
        type = 'zlink',
        newparams = succ.seed.boot[i,-(1:2)]
      )
    ) %>%
    # Reset the phen variable (so conditional is unaffected)
    mutate(phen.c = phen.c - delta) %>%
    mutate(
      seeds.seed.linear =  predict(
        s_st.p_s.u.p2, newdata = ., allow.new.levels = TRUE, re.form = ~ 0,
        type = 'link',
        newparams = succ.seed.boot[i,-(1:2)]
      ),
    ) %>%
    # Taking out the size-zinf terms for 2021 - very extrapolatory, affects averages too much
    mutate(seeds.zinf.linear = ifelse(Year %in% 2021, NA, seeds.zinf.linear),) %>%
    group_by(size, size.nex, trt, trt.phen, phen.c, phen.umbels) %>%
    summarise(across(c(seeds.zinf.linear, seeds.seed.linear), ~ mean(.x, na.rm = TRUE))) %>%
    ungroup() %>%
    mutate(
      # Transform from linear scale to response scale
      seeds.per.umbel = (1 / (1 + exp(-seeds.zinf.linear))) * exp(seeds.seed.linear),
      # Total umbels per plant
      seeds.total = seeds.per.umbel * phen.umbels
    ) %>%
    select(-c(seeds.zinf.linear, seeds.seed.linear)) %>%
    mutate(
      # Mean recruit size
      recr.mean = predict(
        r_t.y, allow.new.levels = TRUE, re.form = ~ 0, newdata = .,
        newparams = recr.boot[i,-(1:2)]
      ),
      # Get the number of seeds produced for each size grouping
      p.size.cur = 0.1 * seeds.total * dnorm(x = size.nex, mean = recr.mean, sd = sqrt(exp(recr.boot$betad[i])))
    ) %>%
    # Rename column
    rename(size.prev = size) %>%
    select(-c(phen.umbels, phen.c, seeds.per.umbel, seeds.total, recr.mean)) %>%
    mutate(param = 'phen.succ')
  
  # 6: phen effects on seed set
  
  this.boot[[6]] = boot.ltre.backbone %>%
    filter(boot %in% i) %>%
    rename(year = Year) %>%
    mutate(
      # Umbel count
      phen.umbels = predict(
        u_s_s.ty, newdata = .,
        allow.new.levels = TRUE, re.form = ~ 0, type = 'response',
        newparams = flow.numb.boot[i,-(1:2)]
      )
    ) %>%
    rename(Year = year) %>%
    mutate(
      # Model predictions for seed set on linear (link) scale for averaging
      seeds.zinf.linear =  predict(
        s_st.p_s.u.p2, newdata = ., allow.new.levels = TRUE, re.form = ~ 0,
        type = 'zlink',
        newparams = succ.seed.boot[i,-(1:2)]
      )
    ) %>%
    # Perturb phen variable (for only the cond term)
    mutate(phen.c = phen.c + delta) %>%
    mutate(
      seeds.seed.linear =  predict(
        s_st.p_s.u.p2, newdata = ., allow.new.levels = TRUE, re.form = ~ 0,
        type = 'link',
        newparams = succ.seed.boot[i,-(1:2)]
      ),
      # Reset the phen variable
      phen.c = phen.c - delta
    ) %>%
    # Taking out the size-zinf terms for 2021 - very extrapolatory, affects averages too much
    mutate(seeds.zinf.linear = ifelse(Year %in% 2021, NA, seeds.zinf.linear),) %>%
    group_by(size, size.nex, trt, trt.phen, phen.c, phen.umbels) %>%
    summarise(across(c(seeds.zinf.linear, seeds.seed.linear), ~ mean(.x, na.rm = TRUE))) %>%
    ungroup() %>%
    mutate(
      # Transform from linear scale to response scale
      seeds.per.umbel = (1 / (1 + exp(-seeds.zinf.linear))) * exp(seeds.seed.linear),
      # Total umbels per plant
      seeds.total = seeds.per.umbel * phen.umbels
    ) %>%
    select(-c(seeds.zinf.linear, seeds.seed.linear)) %>%
    mutate(
      # Mean recruit size
      recr.mean = predict(
        r_t.y, allow.new.levels = TRUE, re.form = ~ 0, newdata = .,
        newparams = recr.boot[i,-(1:2)]
      ),
      # Get the number of seeds produced for each size grouping
      p.size.cur = 0.1 * seeds.total * dnorm(x = size.nex, mean = recr.mean, sd = sqrt(exp(recr.boot$betad[i])))
    ) %>%
    # Rename column
    rename(size.prev = size) %>%
    select(-c(phen.umbels, phen.c, seeds.per.umbel, seeds.total, recr.mean)) %>%
    mutate(param = 'phen.seed')
  
  print(i)
  
  # combine into data frame and add a labeling column
  fr.pert.boot[[i]] = do.call(rbind, this.boot) %>% mutate(boot = paste0('b', i))
  
}

# Export data frame (wide-pivoted)
do.call(rbind, fr.pert.boot) %>%
  pivot_wider(names_from = boot, values_from = p.size.cur) %>%
  write.csv(
    file = '03_construct_kernels/out/deterministic_reprod_perturb_bootstraps.csv',
    row.names = FALSE, na = ''
  )

# # # # Code recycled above
# # # # Phenology bootstrapping
# 
# # Predict the mean buddate from each bootstrapped model
# 
# phen.backbone = expand.grid(
#   trt = c('control', 'drought', 'irrigated'),
#   Year = factor(2021:2024)
# )
# 
# # I'm sure there is a sleeker way to do this... but, feeling lazy now.
# 
# phen.list.out = vector('list', length = n.straps)
# 
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
#   # print(i)
# }
# 
# phen.boot.trt = do.call(rbind, phen.list.out) %>% mutate(boot = as.numeric(boot))
# 
# phen.boot.backbone = expand.grid(
#     size = (5:60)/10,
#     size.nex = (5:60)/10,
#     Year = 2021:2024,
#     boot = 1:n.straps,
#     # This column will be used for manipulating the phenology date and the vital
#     # rate estimation
#     trt.phen.idx = 1:7
# ) %>%
#   merge(
#     data.frame(
#       trt.phen.idx = 1:7,
#       trt.phen = c('drought', 'control', 'irrigated', 'control', 'drought', 'irrigated', 'control'),
#       trt.rate = c('drought', 'drought', 'irrigated', 'irrigated', 'control', 'control', 'control')
#     )
#   ) %>%
#   # This column is superfluous now
#   select(-trt.phen.idx) %>%
#   # Now merge in to get the buddates for the trt.phen column
#   merge(phen.boot.trt, by.x = c('trt.phen', 'boot'), by.y = c('trt', 'boot')) %>%
#   rename(trt = trt.rate) %>%
#   mutate(phen.c = mean.phen - round(mean(seed$mean.phen))) %>%
#   select(-mean.phen)
  
# List to store outputs
phen.boot.list = vector('list', length = n.straps)

# Bootstrap with observed values:

for (i in 1:n.straps) {
  
  phen.boot.list[[i]] = phen.boot.backbone %>%
    filter(boot %in% i) %>%
    # Rename to not put the year random effect in these predictions
    rename(year = Year) %>%
    mutate(
      # Umbel count
      phen.umbels = predict(
        u_s_s.ty, newdata = ., 
        # use bootstrapped parameters
        newparams = flow.numb.boot[i, -(1:2)],
        allow.new.levels = TRUE, re.form = ~ 0, type = 'response'
      )
    ) %>%
    rename(Year = year) %>%
    mutate(
      # Model predictions at treatment means
      # (doing this on linear scale for each for easier averaging)
      seeds.zinf.linear =  predict(
        s_st.p_s.u.p2, newdata = ., allow.new.levels = TRUE, re.form = ~ 0,
        # use bootstrapped parameters
        newparams = succ.seed.boot[i,-(1:2)],
        type = 'zlink'
      ),
      seeds.seed.linear =  predict(
        s_st.p_s.u.p2, newdata = ., allow.new.levels = TRUE, re.form = ~ 0,
        # use bootstrapped parameters
        newparams = succ.seed.boot[i,-(1:2)],
        type = 'link'
      ),
    ) %>%
    # Taking out the size-zinf terms for 2021 - very extrapolatory, affects averages too much
    mutate(seeds.zinf.linear = ifelse(Year %in% 2021, NA, seeds.zinf.linear),) %>%
    group_by(size, size.nex, trt, trt.phen, phen.c, phen.umbels) %>%
    summarise(across(c(seeds.zinf.linear, seeds.seed.linear), ~ mean(.x, na.rm = TRUE))) %>%
    ungroup() %>%
    mutate(
      # Transform from linear scale to response scale
      seeds.per.umbel = (1 / (1 + exp(-seeds.zinf.linear))) * exp(seeds.seed.linear),
      # Total umbels per plant
      seeds.total = seeds.per.umbel * phen.umbels
    ) %>%
    select(-c(seeds.zinf.linear, seeds.seed.linear)) %>%
    mutate(
      # Mean recruit size
      recr.mean = predict(
        r_t.y, allow.new.levels = TRUE, re.form = ~ 0, newdata = .,
        # use bootstrapped parameters
        newparams = recr.boot[i, -(1:2)]
      ),
      # Get the number of seeds produced for each size grouping
      # (note: 'betad' parameter here is the log of the residual variance from the model fit)
      p.size.cur = 0.1 * seeds.total * dnorm(x = size.nex, mean = recr.mean, sd = sqrt(exp(recr.boot$betad[i])))
    ) %>%
    # Rename column
    rename(size.prev = size) %>%
    # Re-center phenology column and add 'boot' column
    mutate(
      mean.phen = phen.c + round(mean(seed$mean.phen)),
      boot = paste0('b', i)
    ) %>%
    # Remove unnecessary columns
    select(-c(phen.umbels, seeds.per.umbel, seeds.total, recr.mean, phen.c))
  
  print(i)
  
}

do.call(rbind, phen.boot.list) %>%
  select(-mean.phen) %>%
  pivot_wider(names_from = boot, values_from = p.size.cur) %>%
  write.csv(
    '03_construct_kernels/out/deterministic_reprod_phen_bootstrap.csv',
    na = '', row.names = FALSE
  )

# Bootstrap perturbations

phen.boot.pert.list = vector('list', n.straps)

for (i in 1:n.straps) {
  
  # Set up a list for storing outputs
  # (will rbind this at the end of the loop to make one data frame per loop
  # iteration, then store that as the ith element of the list)
  this.boot = vector('list', 2)
  
  # Perturb phen for umbel success
  
  this.boot[[1]] = phen.boot.backbone %>%
    filter(boot %in% i) %>%
    rename(year = Year) %>%
    mutate(
      # Umbel count
      phen.umbels = predict(
        u_s_s.ty, newdata = .,
        allow.new.levels = TRUE, re.form = ~ 0, type = 'response',
        newparams = flow.numb.boot[i,-(1:2)]
      )
    ) %>%
    rename(Year = year) %>%
    # Perturb phen variable (for only the zinf term)
    mutate(phen.c = phen.c + delta) %>%
    mutate(
      # Model predictions for seed set on linear (link) scale for averaging
      seeds.zinf.linear =  predict(
        s_st.p_s.u.p2, newdata = ., allow.new.levels = TRUE, re.form = ~ 0,
        type = 'zlink',
        newparams = succ.seed.boot[i,-(1:2)]
      )
    ) %>%
    # Reset the phen variable (so conditional is unaffected)
    mutate(phen.c = phen.c - delta) %>%
    mutate(
      seeds.seed.linear =  predict(
        s_st.p_s.u.p2, newdata = ., allow.new.levels = TRUE, re.form = ~ 0,
        type = 'link',
        newparams = succ.seed.boot[i,-(1:2)]
      ),
    ) %>%
    # Taking out the size-zinf terms for 2021 - very extrapolatory, affects averages too much
    mutate(seeds.zinf.linear = ifelse(Year %in% 2021, NA, seeds.zinf.linear),) %>%
    group_by(size, size.nex, trt, phen.c, trt.phen, phen.umbels) %>%
    summarise(across(c(seeds.zinf.linear, seeds.seed.linear), ~ mean(.x, na.rm = TRUE))) %>%
    ungroup() %>%
    mutate(
      # Transform from linear scale to response scale
      seeds.per.umbel = (1 / (1 + exp(-seeds.zinf.linear))) * exp(seeds.seed.linear),
      # Total umbels per plant
      seeds.total = seeds.per.umbel * phen.umbels
    ) %>%
    select(-c(seeds.zinf.linear, seeds.seed.linear)) %>%
    mutate(
      # Mean recruit size
      recr.mean = predict(
        r_t.y, allow.new.levels = TRUE, re.form = ~ 0, newdata = .,
        newparams = recr.boot[i,-(1:2)]
      ),
      # Get the number of seeds produced for each size grouping
      p.size.cur = 0.1 * seeds.total * dnorm(x = size.nex, mean = recr.mean, sd = sqrt(exp(recr.boot$betad[i])))
    ) %>%
    # Rename column
    rename(size.prev = size) %>%
    select(-c(phen.umbels, seeds.per.umbel, seeds.total, recr.mean)) %>%
    mutate(param = 'phen.succ')
  
  # 6: phen effects on seed set
  
  this.boot[[2]] = phen.boot.backbone %>%
    filter(boot %in% i) %>%
    rename(year = Year) %>%
    mutate(
      # Umbel count
      phen.umbels = predict(
        u_s_s.ty, newdata = .,
        allow.new.levels = TRUE, re.form = ~ 0, type = 'response',
        newparams = flow.numb.boot[i,-(1:2)]
      )
    ) %>%
    rename(Year = year) %>%
    mutate(
      # Model predictions for seed set on linear (link) scale for averaging
      seeds.zinf.linear =  predict(
        s_st.p_s.u.p2, newdata = ., allow.new.levels = TRUE, re.form = ~ 0,
        type = 'zlink',
        newparams = succ.seed.boot[i,-(1:2)]
      )
    ) %>%
    # Perturb phen variable (for only the cond term)
    mutate(phen.c = phen.c + delta) %>%
    mutate(
      seeds.seed.linear =  predict(
        s_st.p_s.u.p2, newdata = ., allow.new.levels = TRUE, re.form = ~ 0,
        type = 'link',
        newparams = succ.seed.boot[i,-(1:2)]
      ),
      # Reset the phen variable
      phen.c = phen.c - delta
    ) %>%
    # Taking out the size-zinf terms for 2021 - very extrapolatory, affects averages too much
    mutate(seeds.zinf.linear = ifelse(Year %in% 2021, NA, seeds.zinf.linear),) %>%
    group_by(size, size.nex, trt, phen.c, trt.phen, phen.umbels) %>%
    summarise(across(c(seeds.zinf.linear, seeds.seed.linear), ~ mean(.x, na.rm = TRUE))) %>%
    ungroup() %>%
    mutate(
      # Transform from linear scale to response scale
      seeds.per.umbel = (1 / (1 + exp(-seeds.zinf.linear))) * exp(seeds.seed.linear),
      # Total umbels per plant
      seeds.total = seeds.per.umbel * phen.umbels
    ) %>%
    select(-c(seeds.zinf.linear, seeds.seed.linear)) %>%
    mutate(
      # Mean recruit size
      recr.mean = predict(
        r_t.y, allow.new.levels = TRUE, re.form = ~ 0, newdata = .,
        newparams = recr.boot[i,-(1:2)]
      ),
      # Get the number of seeds produced for each size grouping
      p.size.cur = 0.1 * seeds.total * dnorm(x = size.nex, mean = recr.mean, sd = sqrt(exp(recr.boot$betad[i])))
    ) %>%
    # Rename column
    rename(size.prev = size) %>%
    select(-c(phen.umbels, seeds.per.umbel, seeds.total, recr.mean)) %>%
    mutate(param = 'phen.seed')
  
  print(i)
  
  # combine into data frame and add a labeling column
  phen.boot.pert.list[[i]] = do.call(rbind, this.boot) %>% mutate(boot = paste0('b', i))
  
}

do.call(rbind, phen.boot.pert.list) %>%
  select(-phen.c) %>%
  pivot_wider(names_from = boot, values_from = p.size.cur) %>%
  write.csv(
    '03_construct_kernels/out/deterministic_reprod_phen_perturb_boostraps.csv',
    na = '', row.names = FALSE
  )

# Export the parameter estimates from each bootstrap

cbind(
  # Global id for bootstrap number
  boot = 1:n.straps,
  flow.int_control = unlist(flow.numb.boot[,-(1:2)][3]),
  flow.int_drought = unlist(flow.numb.boot[,-(1:2)][3] + flow.numb.boot[,-(1:2)][5]),
  flow.int_irrigated = unlist(flow.numb.boot[,-(1:2)][3] + flow.numb.boot[,-(1:2)][6]),
  seed.int_control = unlist(succ.seed.boot[,-(1:2)][1]),
  seed.int_drought = unlist(succ.seed.boot[,-(1:2)][1] + succ.seed.boot[,-(1:2)][2]),
  seed.int_irrigated = unlist(succ.seed.boot[,-(1:2)][1] + succ.seed.boot[,-(1:2)][3]),
  seed.slope_control = unlist(succ.seed.boot[,-(1:2)][4]),
  seed.slope_drought = unlist(succ.seed.boot[,-(1:2)][4] + succ.seed.boot[,-(1:2)][9]),
  seed.slope_irrigated = unlist(succ.seed.boot[,-(1:2)][4] + succ.seed.boot[,-(1:2)][10]),
  recr.int_control = unlist(recr.boot[,-(1:2)][1]),
  recr.int_drought = unlist(recr.boot[,-(1:2)][1] + recr.boot[,-(1:2)][2]),
  recr.int_irrigated = unlist(recr.boot[,-(1:2)][1] + recr.boot[,-(1:2)][3]),
  phen_control = phen.boot.trt %>% filter(trt %in% 'control') %>% pull(mean.phen),
  phen_drought = phen.boot.trt %>% filter(trt %in% 'drought') %>% pull(mean.phen),
  phen_irrigated = phen.boot.trt %>% filter(trt %in% 'irrigated') %>% pull(mean.phen)
) %>%
  write.csv(
    file = '03_construct_kernels/out/reprod_bootstrapped_perturbed_params.csv',
    row.names = FALSE, na = ''
  )
