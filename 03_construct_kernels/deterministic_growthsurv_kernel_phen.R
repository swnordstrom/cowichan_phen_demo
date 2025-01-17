# Script for re-running reproductive vital rate models and performing model
# selection and exporting a deterministic, survival e kernel for
# IPM analysis.
# Reads in processed demo/seed (including phenology) data, 2016-2024
# Reads in processed demo/seed (including phenology) data, 2016-2024
# based on the script `deterministic_growthsurv_kernel.R` which does not have
# phenology.
# (sn init jan 2025)
# --- Setup ---------------------------------------------------------

library(ggplot2)
library(dplyr)
library(tidyr)
library(glmmTMB)

rm(list = ls())

# Read in demo data and merge with treatment info
all.data = merge(
  x = read.csv('01_data_cleaning/out/demo_phen_seed_2016-2024_final.csv'),
  y = read.csv('00_raw_data/plot_treatments.csv'),
  by.x = 'Plot', by.y = 'plot'
)

nrow(all.data)
head(all.data)

# Read in phenology means
phen.treatment.means = read.csv('03_construct_kernels/phen_treatment_means.csv')
# Mean that will be used for centering
phen.ctrl.mean = phen.treatment.means$mean.phen[phen.treatment.means$trt %in% 'control']

# We're interested only in plants that are in demo for this analysis
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
) %>%
  # Change years to factors
  mutate(across(contains('year'), as.factor))

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
  filter(!(prev.year %in% 2016)) %>%
  # Add size columns
  mutate(size.prev = log(No.leaves.pre * Leaf.length.pre))

nrow(demo.surv.sizes)
table(demo.surv.sizes$surv, useNA = 'always')

# Subset for growth estimation
demo.grow = demo.surv.sizes %>% 
  filter(surv, !is.na(Leaf.length) & !is.na(No.leaves) & Leaf.length > 0 & No.leaves > 0) %>%
  mutate(size.cur = log(Leaf.length * No.leaves))

# Finally: subset surv dataset to not include 2023-2024 surv
demo.surv.sizes = demo.surv.sizes %>% filter(!(surv.year %in% 2024))

# Get phenology dataset
# this is for merging in with growth dataset, so relevant measure is by plant (not by umbel)
phen.by.plant.for.growth = all.data %>%
  # Give me plants that are in phenology
  filter(in.phen) %>%
  # Extract the bud dates on file
  separate_wider_delim(phen.julis, names = paste0('uu', 1:12), delim = ';', too_few = 'align_start') %>%
  pivot_longer(starts_with('uu'), names_to = 'umbel.number', values_to = 'phen.julian') %>%
  filter(!is.na(phen.julian)) %>%
  # convert to julian date
  mutate(phen.julian = as.numeric(gsub('\\s', '', phen.julian))) %>%
  # give me the mean bud date for each plant
  group_by(plantid, Year) %>%
  summarise(phen.mean = mean(phen.julian)) %>%
  ungroup() %>%
  # Change year to factor
  mutate(Year = as.factor(Year))

demo.grow = merge(
  demo.grow, phen.by.plant.for.growth, 
  by.x = c('prev.year', 'plantid'), by.y = c('Year', 'plantid'),
  all.x = TRUE, all.y = FALSE
) %>%
  # Center mean around control
  mutate(phen.c = phen.mean - phen.ctrl.mean)


# ------------------------------------------------
# Construct kernels

# --- Refit final models

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


# --- Extract parameters needed

# Residual variance in growth models
gv.sd = summary(g_st.ty)$sigma
gf.sd = summary(g_phen)$sigma

# Phenology effect per day
phen.effect = g_phen$fit$par[7]

# --- Construct and work with data frame

# Get a scaffold
grow.surv.kernel = expand.grid(
  size.prev = (5:60)/10,
  size.cur = (5:60)/10,
  trt = c('control', 'drought', 'irrigated'),
  phen.c = -28:28
)

grow.surv.kernel = grow.surv.kernel %>%
  # Predicted survival
  mutate(
    pred.surv = predict(
      newdata = .,
      object = s_s, type = 'response',
      re.form = ~ 0, allow.new.levels = TRUE
    )
  ) %>%
  # Predicted growth
  mutate(
    # Model without phenology
    pred.grow.mean = predict(
      newdata = .,
      object = g_st.ty, type = 'response',
      re.form = ~ 0, allow.new.levels = TRUE
    )
  ) %>%
  # Model with phenology
  # Need to change name of year column to get annual predictions
  mutate(phen.grow.mean = pred.grow.mean + phen.effect * phen.c) %>%
  # # Take the average of the growth kernel across years
  # group_by(size.prev, size.cur, trt, phen.c, pred.surv, pred.grow.mean) %>%
  # summarise(phen.grow.mean = mean(phen.grow.mean)) %>%
  # ungroup() %>%
  # Predicted distribution of sizes in next time step
  mutate(
    pv.grow.size = 0.1 * dnorm(size.cur, mean = pred.grow.mean, sd = gv.sd),
    pf.grow.size = 0.1 * dnorm(size.cur, mean = phen.grow.mean, sd = gf.sd)
  )

grow.surv.kernel %>%
  filter(phen.c %in% c(-28, 28)) %>%
  mutate(p.size.cur = pred.surv * pf.grow.size) %>%
  ggplot(aes(x = size.prev, y = size.cur)) +
  geom_tile(aes(fill = p.size.cur)) +
  scale_y_reverse() +
  scale_fill_viridis_c() +
  facet_wrap(trt ~ phen.c, nrow = 3)
# Not that different

grow.surv.kernel %>%
  filter(phen.c %in% 0) %>%
  pivot_longer(c(pv.grow.size, pf.grow.size), names_to = 'model', values_to = 'p.grow.size') %>%
  mutate(p.size.cur = pred.surv * p.grow.size) %>%
  ggplot(aes(x = size.prev, y = size.cur)) +
  geom_tile(aes(fill = p.size.cur)) +
  scale_fill_viridis_c() +
  scale_y_reverse() +
  facet_wrap(trt ~ model, nrow = 3)
# Hmm...

write.csv(
  grow.surv.kernel %>%
    mutate(phen = phen.c + phen.ctrl.mean) %>%
    select(-c(phen.c, phen.grow.mean, pred.grow.mean)),
  file = '03_construct_kernels/out/deterministic_growsurv_kernel_phen.csv',
  row.names = FALSE
)


# --- Lambda estimates for LTRE

# LTRE backbone (combinations of treatment used to estimate vital rate and
# treatment used for phenology)
ltre.backbone = expand.grid(
  size.prev = (5:60)/10,
  size.cur  = (5:60)/10,
  # year = 2021:2023,
  # This column will be used for manipulating the phenology date and the vital
  # rate estimation
  trt.phen.idx = 1:7
) %>%
  merge(
    data.frame(
      trt.phen.idx = 1:7,
      # Treatment associated with the buddates used for umbel success + seeds
      trt.phen = c('drought', 'control', 'irrigated', 'control', 'drought', 'irrigated', 'control'),
      # Treatment associated with direct treatment effects on vital rates
      trt.rate = c('drought', 'drought', 'irrigated', 'irrigated', 'control', 'control', 'control')
    )
  ) %>%
  # Merge with estimated mean bud date per treatment
  merge(phen.treatment.means, by.x = 'trt.phen', by.y = 'trt') %>%
  # Rename trt column so it is used in models
  rename(trt = trt.rate) %>% 
  # center the phenology column and rename the `trt` column so it can be used in
  # vital rate estimates
  mutate(phen.c = mean.phen - phen.ctrl.mean)

# Get kernel
ltre.kernel = ltre.backbone %>%
  # Predicted survival
  mutate(
   pred.surv = predict(
    newdata = .,
    object = s_s, type = 'response',
    re.form = ~ 0, allow.new.levels = TRUE
  )
) %>%
  # Predicted growth
  mutate(
    # Model without phenology
    pred.grow.mean = predict(
      newdata = .,
      object = g_st.ty, type = 'response',
      re.form = ~ 0, allow.new.levels = TRUE
    )
  ) %>%
  # Model with phenology
  mutate(phen.grow.mean = pred.grow.mean + phen.effect * phen.c) %>%
  # # Take the average of the growth kernel across years
  # group_by(size.prev, size.cur, trt, trt.phen, phen.c, pred.surv, pred.grow.mean) %>%
  # summarise(phen.grow.mean = mean(phen.grow.mean)) %>%
  # ungroup() %>%
  # Predicted distribution of sizes in next time step
  mutate(
    pv.grow.size = 0.1 * dnorm(size.cur, mean = pred.grow.mean, sd = gv.sd),
    pf.grow.size = 0.1 * dnorm(size.cur, mean = phen.grow.mean, sd = gf.sd)
  )

# Export
write.csv(
  ltre.kernel %>%
    mutate(phen = phen.c + phen.ctrl.mean) %>%
    select(-c(phen.c, mean.phen, phen.grow.mean, pred.grow.mean, trt.phen.idx)),
  file = '03_construct_kernels/out/deterministic_growsurv_kernel_phen_ltre.csv',
  row.names = FALSE
)


# --- Sensitivities

# Parameters of interest to us (vary by either/both of treatment and budding phenology)
# - Growth model (no phenology):
#   - Treatment intercept term
#   - Treatment slope (trt:size.prev term)
# - Growth model (with phenology)
#   - Phenology term
# (conundrum... are the treatment intercepts in the two models independent of
# each other?) (I suppose we can treat these as growth of a flowering plant and
# growth of a vegetative plant...)

# Perturbation amount
delta = 0.0001

# Get a list for outputs
outputs = vector('list', 3)

# Start the perturbations

# 1: Growth model (no phen) treatment intercept term

# intercept (control) is first parameter listed ([1]), 
# drought and control effects on intercept are [3] and [4]

outputs[[1]] = ltre.backbone %>%
  # Predicted survival
  mutate(
    pred.surv = predict(
      newdata = .,
      object = s_s, type = 'response',
      re.form = ~ 0, allow.new.levels = TRUE
    )
  ) %>%
  # Predicted growth
  mutate(
    pred.grow.mean = predict(
      newdata = .,
      object = g_st.ty, type = 'response',
      newparams = g_st.ty$fit$par %>%
        (function(x) {
          x[1] <- x[1] + delta
          return(x)
        }),
      re.form = ~ 0, allow.new.levels = TRUE
    )
  ) %>%
  # Model with phenology
  mutate(phen.grow.mean = pred.grow.mean + phen.effect * phen.c) %>%
  # Take the average of the growth kernel across years
  # group_by(size.prev, size.cur, trt, trt.phen, phen.c, pred.surv, pred.grow.mean) %>%
  # summarise(phen.grow.mean = mean(phen.grow.mean)) %>%
  # ungroup() %>%
  # Predicted distribution of sizes in next time step
  mutate(
    pv.grow.size = 0.1 * dnorm(size.cur, mean = pred.grow.mean, sd = gv.sd),
    pf.grow.size = 0.1 * dnorm(size.cur, mean = phen.grow.mean, sd = gf.sd)
  ) %>%
  # Add in perturbation information
  mutate(
    perturb.param = 'grow.int',
    orig.par.val = case_match(
      trt,
      'control' ~ g_st.ty$fit$par[1],
      'drought' ~ g_st.ty$fit$par[1] + g_st.ty$fit$par[3],
      'irrigated' ~ g_st.ty$fit$par[1] + g_st.ty$fit$par[4]
    ) # g_st.ty$fit$par[1]
  )

# 2: Growth model (no phen) slope; slope parameter (control) is [2] param listed,
# the drought and irrigated slope effect terms resp. are [5] and [6]

outputs[[2]] = ltre.backbone %>%
  # Predicted survival
  mutate(
    pred.surv = predict(
      newdata = .,
      object = s_s, type = 'response',
      re.form = ~ 0, allow.new.levels = TRUE
    )
  ) %>%
  # Predicted growth
  mutate(
    # Model with no phenology
    pred.grow.mean = predict(
      newdata = .,
      object = g_st.ty, type = 'response',
      newparams = g_st.ty$fit$par %>%
        (function(x) {
          x[2] <- x[2] + delta
          return(x)
        }),
      re.form = ~ 0, allow.new.levels = TRUE
    )
  ) %>%
  # Model with phenology
  mutate(phen.grow.mean = pred.grow.mean + phen.effect * phen.c) %>%
  # # Take the average of the growth kernel across years
  # group_by(size.prev, size.cur, trt, trt.phen, phen.c, pred.surv, pred.grow.mean) %>%
  # summarise(phen.grow.mean = mean(phen.grow.mean)) %>%
  # ungroup() %>%
  # Predicted distribution of sizes in next time step
  mutate(
    pv.grow.size = 0.1 * dnorm(size.cur, mean = pred.grow.mean, sd = gv.sd),
    pf.grow.size = 0.1 * dnorm(size.cur, mean = phen.grow.mean, sd = gf.sd)
  ) %>%
  # Add in perturbation information
  mutate(
    perturb.param = 'grow.slope',
    orig.par.val = case_match(
      trt,
      'control' ~ g_st.ty$fit$par[2],
      'drought' ~ g_st.ty$fit$par[2] + g_st.ty$fit$par[5],
      'irrigated' ~ g_st.ty$fit$par[2] + g_st.ty$fit$par[6]
    ) # g_st.ty$fit$par[2]
  )

# 3: Growth model phenology effect

outputs[[3]] = ltre.backbone %>%
  # Predicted survival
  mutate(
    pred.surv = predict(
      newdata = .,
      object = s_s, type = 'response',
      re.form = ~ 0, allow.new.levels = TRUE
    )
  ) %>%
  # Predicted growth
  mutate(
    # Model with no phenology
    pred.grow.mean = predict(
      newdata = .,
      object = g_st.ty, type = 'response',
      re.form = ~ 0, allow.new.levels = TRUE
    )
  ) %>%
  # Model with phenology
  mutate(phen.grow.mean = pred.grow.mean + (phen.c + delta) * phen.effect) %>%
  # # Take the average of the growth kernel across years
  # group_by(size.prev, size.cur, trt, trt.phen, phen.c, pred.surv, pred.grow.mean) %>%
  # summarise(phen.grow.mean = mean(phen.grow.mean)) %>%
  # ungroup() %>%
  # Predicted distribution of sizes in next time step
  mutate(
    pv.grow.size = 0.1 * dnorm(size.cur, mean = pred.grow.mean, sd = gv.sd),
    pf.grow.size = 0.1 * dnorm(size.cur, mean = phen.grow.mean, sd = gf.sd)
  ) %>%
  # Add in perturbation information
  mutate(
    perturb.param = 'phen.grow',
    orig.par.val = phen.c
  )


# Bind them all together
outputs.all = do.call(rbind, outputs) %>%
  # De-center phenology
  mutate(phen = phen.c + phen.ctrl.mean)

write.csv(
  outputs.all %>% 
    select(
      size.prev, size.cur, trt, trt.phen, phen, 
      pred.surv, pv.grow.size, pf.grow.size, perturb.param, orig.par.val
    ),
  file = '03_construct_kernels/out/deterministic_grow_coef_perturbation_phen.csv',
  row.names = FALSE
)
