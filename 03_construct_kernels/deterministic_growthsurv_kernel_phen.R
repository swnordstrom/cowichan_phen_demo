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

source('03_construct_kernels/prepare_demo_data_growsurv.R')

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
# gf.sd = summary(g_phen)$sigma

# Phenology effect per day
phen.effect = g_phen$fit$par[7]

# --- Construct and work with data frame

# Get a scaffold
grow.surv.kernel = expand.grid(
  size.prev = (5:60)/10,
  size.cur = (5:60)/10,
  trt = c('control', 'drought', 'irrigated'),
  phen.c = -14:14
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
    pf.grow.size = 0.1 * dnorm(size.cur, mean = phen.grow.mean, sd = gv.sd)
  )

grow.surv.kernel %>%
  filter(phen.c %in% c(-14, 14)) %>%
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
    # mutate(phen = phen.c + phen.ctrl.mean) %>%
    select(-c(phen.grow.mean, pred.grow.mean)),
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
  # Read in combinations of treatment and phenology for LTRE desigm
  merge(read.csv('03_construct_kernels/ltre_treatment_key.csv')) %>%
  # Merge with estimated mean bud date per treatment
  # 'trt.phen' column is the treatment for which mean phenology is used
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
    pf.grow.size = 0.1 * dnorm(size.cur, mean = phen.grow.mean, sd = gv.sd)
  )

# Export
write.csv(
  ltre.kernel %>%
    select(size.prev, size.cur, trt.phen.idx, pred.surv, pv.grow.size, pf.grow.size),
    # mutate(phen = phen.c + phen.ctrl.mean) %>%
    # select(-c(phen.c, mean.phen, phen.grow.mean, pred.grow.mean, trt.phen.idx)),
  file = '03_construct_kernels/out/deterministic_growsurv_kernel_phen_ltre.csv',
  row.names = FALSE
)

cat('Exported growth+survival subkernel\n')


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
delta = 0.001

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
    pf.grow.size = 0.1 * dnorm(size.cur, mean = phen.grow.mean, sd = gv.sd)
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
    pf.grow.size = 0.1 * dnorm(size.cur, mean = phen.grow.mean, sd = gv.sd)
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
  mutate(phen.grow.mean = pred.grow.mean + delta + (phen.c * phen.effect)) %>%
  # # Take the average of the growth kernel across years
  # group_by(size.prev, size.cur, trt, trt.phen, phen.c, pred.surv, pred.grow.mean) %>%
  # summarise(phen.grow.mean = mean(phen.grow.mean)) %>%
  # ungroup() %>%
  # Predicted distribution of sizes in next time step
  mutate(
    pv.grow.size = 0.1 * dnorm(size.cur, mean = pred.grow.mean, sd = gv.sd),
    pf.grow.size = 0.1 * dnorm(size.cur, mean = phen.grow.mean, sd = gv.sd)
  ) %>%
  # Add in perturbation information
  mutate(
    perturb.param = 'phen.grow',
    orig.par.val = (phen.c) * phen.effect
  )


# Bind them all together
outputs.all = do.call(rbind, outputs) # %>%
  # De-center phenology
  # mutate(phen = phen.c + phen.ctrl.mean)

write.csv(
  outputs.all %>% 
    select(
      size.prev, size.cur, trt.phen.idx, 
      pred.surv, pv.grow.size, pf.grow.size, perturb.param, orig.par.val
    ),
  file = '03_construct_kernels/out/deterministic_grow_coef_perturbation_phen.csv',
  row.names = FALSE
)

cat('Exported growth+survival perturbed subkernels\n')
