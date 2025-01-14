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
  ungroup()

demo.grow = merge(
  demo.grow, phen.by.plant.for.growth, 
  by.x = c('prev.year', 'plantid'), by.y = c('Year', 'plantid'),
  all.x = TRUE, all.y = FALSE
) %>%
  mutate(phen.mean.c = phen.mean - round(mean(phen.mean, na.rm = TRUE)))

# Phenology for mean bud date by phenology
# i.e., for umbel
# This is used for getting dates for LTRE
# Umbel-level budding phenology data

phen.by.umbel.for.ltre = all.data %>% 
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
  size.cur ~ size.prev * prev.year + trt * prev.year + phen.mean.c + (1 | Plot / plantid),
  data = demo.grow %>% filter(!is.na(phen.mean.c))
)

# Phen (umbel, treatment effects) model
d_t = glmmTMB(
  phen.julian ~ trt + Year + (1 | Plot / plantid),
  data = phen.by.umbel.for.ltre
)

# --- Extract parameters needed

# Residual variance in growth models
gv.sd = summary(g_st.ty)$sigma
gf.sd = summary(g_phen)$sigma

# --- Construct and work with data frame

# Get a scaffold
grow.surv.kernel = expand.grid(
  size.prev = (5:60)/10,
  size.cur = (5:60)/10,
  trt = c('control', 'drought', 'irrigated'),
  year = 2021:2023,
  phen.mean.c = -28:28
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
  rename(prev.year = year) %>%
  mutate(
    phen.grow.mean = predict(
      newdata = .,
      object = g_phen, type = 'response',
      re.form = ~ 0, allow.new.levels = TRUE
    )
  ) %>%
  # Take the average of the growth kernel across years
  group_by(size.prev, size.cur, trt, phen.mean.c, pred.surv, pred.grow.mean) %>%
  summarise(phen.grow.mean = mean(phen.grow.mean)) %>%
  ungroup() %>%
  # Predicted distribution of sizes in next time step
  mutate(
    pv.grow.size = 0.1 * dnorm(size.cur, mean = pred.grow.mean, sd = gv.sd),
    pf.grow.size = 0.1 * dnorm(size.cur, mean = phen.grow.mean, sd = gf.sd)
  )

grow.surv.kernel %>%
  filter(phen.mean.c %in% c(-28, 28)) %>%
  mutate(p.size.cur = pred.surv * pf.grow.size) %>%
  ggplot(aes(x = size.prev, y = size.cur)) +
  geom_tile(aes(fill = p.size.cur)) +
  scale_y_reverse() +
  scale_fill_viridis_c() +
  facet_wrap(phen.mean.c ~ trt)
# Not that different

grow.surv.kernel %>%
  filter(phen.mean.c %in% 0) %>%
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
    mutate(phen = phen.mean.c + round(mean(phen.by.plant.for.growth$phen.mean))) %>%
    select(-c(phen.mean.c, phen.grow.mean, pred.grow.mean)),
  file = '03_construct_kernels/out/deterministic_growsurv_kernel_phen.csv',
  row.names = FALSE
)


# --- Lambda estimates for LTRE

# Mean bud date for each treatment
trt.mean.buddates = expand.grid(trt = c('control', 'drought', 'irrigated'), Year = factor(2021:2024)) %>%
  mutate(
    mean.bud = predict(
      d_t, re.form = ~ 0, allow.new.levels = TRUE,
      newdata = expand.grid(trt = c('control', 'drought', 'irrigated'), Year = factor(2021:2024))
    )
  ) %>%
  group_by(trt) %>%
  summarise(mean.phen = mean(mean.bud))

# LTRE backbone (combinations of treatment used to estimate vital rate and
# treatment used for phenology)
ltre.backbone = expand.grid(
  size.prev = (5:60)/10,
  size.cur  = (5:60)/10,
  year = 2021:2023,
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
  merge(trt.mean.buddates, by.x = 'trt.phen', by.y = 'trt') %>%
  # Rename trt column so it is used in models
  rename(trt = trt.rate) %>% 
  # center the phenology column and rename the `trt` column so it can be used in
  # vital rate estimates
  mutate(phen.mean.c = mean.phen - round(mean(phen.by.plant.for.growth$phen.mean)))

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
  # Need to change name of year column to get annual predictions
  rename(prev.year = year) %>%
  mutate(
    phen.grow.mean = predict(
      newdata = .,
      object = g_phen, type = 'response',
      re.form = ~ 0, allow.new.levels = TRUE
    )
  ) %>%
  # Take the average of the growth kernel across years
  group_by(size.prev, size.cur, trt, trt.phen, phen.mean.c, pred.surv, pred.grow.mean) %>%
  summarise(phen.grow.mean = mean(phen.grow.mean)) %>%
  ungroup() %>%
  # Predicted distribution of sizes in next time step
  mutate(
    pv.grow.size = 0.1 * dnorm(size.cur, mean = pred.grow.mean, sd = gv.sd),
    pf.grow.size = 0.1 * dnorm(size.cur, mean = phen.grow.mean, sd = gf.sd)
  )

# Export
write.csv(
  ltre.kernel %>%
    mutate(phen = phen.mean.c + round(mean(phen.by.plant.for.growth$phen.mean))) %>%
    select(-c(phen.mean.c, phen.grow.mean, pred.grow.mean)),
  file = '03_construct_kernels/out/deterministic_growsurv_kernel_phen_ltre.csv',
  row.names = FALSE
)


# --- Sensitivities

# Parameters of interest:
# - Survival model
#   - Intercept
#   - Slope (size-dependence)
# - Growth model:
#   - Intercept
#   - Slope (size-dependence)
# - 

# Perturbation amount
delta = 0.0001

# Get a list for outputs
outputs = vector('list', 5)

# Residual variance in growth models
grow.sd = summary(g_st.ty)$sigma

# Data frame to generate predictions for
grow.surv.kernel = expand.grid(
  size.prev = (5:60)/10,
  size.cur = (5:60)/10,
  trt = c('control', 'drought', 'irrigated')
)

# Start the perturbations

# 1: Survival model intercept
outputs[[1]] = grow.surv.kernel %>%
  # Predicted survival
  mutate(
    pred.surv = predict(
      newdata = .,
      object = s_s, type = 'response',
      newparams = s_s$fit$par %>%
        (function(x) {
          x[1] <- x[1] + delta
          return(x)
        }),
      re.form = ~ 0, allow.new.levels = TRUE
    )
  ) %>%
  # Predicted growth
  mutate(
    pred.grow.mean = predict(
      newdata = .,
      object = g_st.ty, type = 'response',
      re.form = ~ 0, allow.new.levels = TRUE
    )
  ) %>%
  mutate(
    p.grow.size = 0.1 * dnorm(size.cur, pred.grow.mean, grow.sd)
  ) %>%
  # Combine all together to get overall size distribution in next time step
  mutate(p.size.cur = pred.surv * p.grow.size) %>%
  mutate(
    perturb.param = 'surv_int',
    # no treatment effects on survival
    orig.par.val = s_s$fit$par[1]
  )

# 2: Survival model 
outputs[[2]] = grow.surv.kernel %>%
  # Predicted survival
  mutate(
    pred.surv = predict(
      newdata = .,
      object = s_s, type = 'response',
      newparams = s_s$fit$par %>%
        (function(x) {
          x[2] <- x[2] + delta
          return(x)
        }),
      re.form = ~ 0, allow.new.levels = TRUE
    )
  ) %>%
  # Predicted growth
  mutate(
    pred.grow.mean = predict(
      newdata = .,
      object = g_st.ty, type = 'response',
      re.form = ~ 0, allow.new.levels = TRUE
    )
  ) %>%
  mutate(
    p.grow.size = 0.1 * dnorm(size.cur, pred.grow.mean, grow.sd)
  ) %>%
  # Combine all together to get overall size distribution in next time step
  mutate(p.size.cur = pred.surv * p.grow.size) %>%
  mutate(
    perturb.param = 'surv_slope',
    # no treatment effects on survival size dependence
    orig.par.val = s_s$fit$par[2]
  )

# 3: Growth model intercept

outputs[[3]] = grow.surv.kernel %>%
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
  mutate(
    p.grow.size = 0.1 * dnorm(size.cur, pred.grow.mean, grow.sd)
  ) %>%
  # Combine all together to get overall size distribution in next time step
  mutate(p.size.cur = pred.surv * p.grow.size) %>%
  mutate(
    perturb.param = 'grow_int',
    orig.par.val = case_when(
      trt %in% 'control' ~ g_st.ty$fit$par[1],
      trt %in% 'drought' ~ g_st.ty$fit$par[1] + g_st.ty$fit$par[3],
      trt %in% 'irrigated' ~ g_st.ty$fit$par[1] + g_st.ty$fit$par[4]
    ) # g_st.ty$fit$par[1]
  )


# 4: Growth model intercept

outputs[[4]] = grow.surv.kernel %>%
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
          x[2] <- x[2] + delta
          return(x)
        }),
      re.form = ~ 0, allow.new.levels = TRUE
    )
  ) %>%
  mutate(
    p.grow.size = 0.1 * dnorm(size.cur, pred.grow.mean, grow.sd)
  ) %>%
  # Combine all together to get overall size distribution in next time step
  mutate(p.size.cur = pred.surv * p.grow.size) %>%
  mutate(
    perturb.param = 'grow_slope',
    orig.par.val = case_when(
      trt %in% 'control' ~ g_st.ty$fit$par[2],
      trt %in% 'drought' ~ g_st.ty$fit$par[2] + g_st.ty$fit$par[5],
      trt %in% 'irrigated' ~ g_st.ty$fit$par[2] + g_st.ty$fit$par[6]
    ) # g_st.ty$fit$par[2]
  )

# 5: Growth model standard deviation

outputs[[5]] = grow.surv.kernel %>%
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
      re.form = ~ 0, allow.new.levels = TRUE
    )
  ) %>%
  mutate(
    p.grow.size = 0.1 * dnorm(size.cur, pred.grow.mean, grow.sd + delta)
  ) %>%
  # Combine all together to get overall size distribution in next time step
  mutate(p.size.cur = pred.surv * p.grow.size) %>%
  mutate(
    perturb.param = 'grow_sigma',
    orig.par.val = grow.sd
  )

outputs.all = do.call(rbind, outputs)

write.csv(
  outputs.all %>% select(size.prev, size.cur, trt, p.size.cur, perturb.param, orig.par.val),
  file = '03_construct_kernels/out/deterministic_grow_coef_perturbation_no_phen.csv',
  row.names = FALSE
)
