library(ggplot2)
library(dplyr)
library(tidyr)
library(glmmTMB)

# =====================================
#  ------------------------------------
# Data read in and processing/reshaping
#  ------------------------------------
# =====================================

rm(list = ls())

# Read in data with wrapper script
source('03_construct_kernels/prepare_demo_data_repr.R')

# =======================================================
# -------------------------------------------------------
# Fit models (model selection performed in other scripts)
# -------------------------------------------------------
# =======================================================

# === Phenology model ===
# Response: day of umbel budding (continuous julian date)
# Predictors: year (factor/categorical) and treatment (factor)
d_t = glmmTMB(
  phen.julian ~ trt + Year + (1 | Plot / plantid),
  data = phen
)

# summary(d_t)

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

# =========================
# -------------------------
# Generate kernel estimates
# -------------------------
# =========================

# Strategy here:
# For each combination of (1) current size, (2) subsequent size, (3) treatment,
# and (4) mean bud date, get an estimate of a transition probability
# Because this is the recruitment subkernel, (2) is recruit size and (4) is
# proportional to number of recruits in that size class
# Because the seed model has year as a fixed effect, we need to also generate
# yearly estimates (which will then be averaged together)

# Get SD for the new recruit size
sigma.recr = sigma.recr = summary(r_t.y)$sigma

# === A kernel for all bud dates === 

backbone = expand.grid(
  size = (5:60)/10,
  size.nex = (5:60)/10,
  trt = c('control', 'drought', 'irrigated'),
  phen.c = -28:28,
  year = 2021:2024
)

nrow(backbone)
# massive, massive data frame...

# Generate the kernel

all.phen.kernel = backbone %>%
  mutate(
    # Umbel count
    phen.umbels = predict(
      u_s_s.ty, newdata = ., 
      allow.new.levels = TRUE, re.form = ~ 0, type = 'response'
    )
  ) %>%
  mutate(Year = year) %>%
  mutate(
   seeds.zinf.linear =  predict(
      s_st.p_s.u.p2, newdata = ., allow.new.levels = TRUE, re.form = ~ 0,
      type = 'zlink'
    ),
    seeds.seed.linear =  predict(
      s_st.p_s.u.p2, newdata = ., allow.new.levels = TRUE, re.form = ~ 0,
      type = 'link'
    ),
  ) %>%
  # Taking out the size-zinf terms for 2021 - very extrapolatory, affects averages too much
  mutate(seeds.zinf.linear = ifelse(Year %in% 2021, NA, seeds.zinf.linear)) %>%
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
      r_t.y, allow.new.levels = TRUE, re.form = ~ 0, newdata = .
    ),
    # Get the number of seeds produced for each size grouping
    p.size.cur = 0.1 * seeds.total * dnorm(x = size.nex, mean = recr.mean, sd = sigma.recr)
  ) %>%
  # Rename column
  rename(size.prev = size)

# Export
write.csv(
  all.phen.kernel %>% 
    mutate(phen = phen.c + round(mean(seed$mean.phen))) %>% 
    select(-c(phen.umbels, recr.mean, seeds.per.umbel, seeds.total, phen.c)),
  file = '03_construct_kernels/out/deterministic_reprod_kernel_phen.csv',
  row.names = FALSE, na = ''
)


# === Kernel at observed phenolgoy dates (for LTRE) ===

# This is actually a better way to do the LTRE backbone I think:
trt.mean.buddates = expand.grid(trt = c('control', 'drought', 'irrigated'), Year = factor(2021:2024)) %>%
  mutate(
    mean.bud = predict(
      d_t, re.form = ~ 0, allow.new.levels = TRUE,
      newdata = expand.grid(trt = c('control', 'drought', 'irrigated'), Year = factor(2021:2024))
    )
  ) %>%
  group_by(trt) %>%
  summarise(mean.phen = mean(mean.bud))

ltre.backbone = expand.grid(
  size = (5:60)/10,
  size.nex = (5:60)/10,
  Year = 2021:2024,
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
  mutate(phen.c = mean.phen - mean(round(seed$mean.phen)))

ltre.kernel = ltre.backbone %>%
  # Rename to not put the year random effect in these predictions
  rename(year = Year) %>%
  mutate(
    # Umbel count
    phen.umbels = predict(
      u_s_s.ty, newdata = ., 
      allow.new.levels = TRUE, re.form = ~ 0, type = 'response'
    )
  ) %>%
  rename(Year = year) %>%
  mutate(
    # Model predictions at treatment means
    # (doing this on linear scale for each for easier averaging)
    seeds.zinf.linear =  predict(
      s_st.p_s.u.p2, newdata = ., allow.new.levels = TRUE, re.form = ~ 0,
      type = 'zlink'
    ),
    seeds.seed.linear =  predict(
      s_st.p_s.u.p2, newdata = ., allow.new.levels = TRUE, re.form = ~ 0,
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
      r_t.y, allow.new.levels = TRUE, re.form = ~ 0, newdata = .
    ),
    # Get the number of seeds produced for each size grouping
    p.size.cur = 0.1 * seeds.total * dnorm(x = size.nex, mean = recr.mean, sd = sigma.recr)
  ) %>%
  # Rename column
  rename(size.prev = size)

# # Do the same thing but for drought and irrigation but at teh control mean phenology
# ltre.control.backbone = ltre.backbone %>%
#   select(-phen.c) %>%
#   pivot_wider(names_from = trt, values_from = mean.phen) %>%
#   mutate(drought = control, irrigated = control) %>%
#   # (We don't need the control estimates)
#   select(-control) %>%
#   pivot_longer(c(drought, irrigated), names_to = 'trt', values_to = 'mean.phen') %>%
#   mutate(phen.c = mean.phen - round(mean(seed$mean.phen)))
# 
# ltre.control.kernel = ltre.control.backbone %>%
#   # Rename to not put the year random effect in these predictions
#   rename(year = Year) %>%
#   mutate(
#     # Umbel count
#     phen.umbels = predict(
#       u_s_s.ty, newdata = ., 
#       allow.new.levels = TRUE, re.form = ~ 0, type = 'response'
#     )
#   ) %>%
#   rename(Year = year) %>%
#   mutate(
#     # Model predictions at *control means*
#     seeds.zinf.linear =  predict(
#       s_st.p_s.u.p2, newdata = ., allow.new.levels = TRUE, re.form = ~ 0,
#       type = 'zlink'
#     ),
#     seeds.seed.linear =  predict(
#       s_st.p_s.u.p2, newdata = ., allow.new.levels = TRUE, re.form = ~ 0,
#       type = 'link'
#     ),
#   ) %>%
#   # Taking out the size-zinf terms for 2021 - very extrapolatory, affects averages too much
#   mutate(seeds.zinf.linear = ifelse(Year %in% 2021, NA, seeds.zinf.linear),) %>%
#   group_by(size, size.nex, trt, phen.c, phen.umbels) %>%
#   summarise(across(c(seeds.zinf.linear, seeds.seed.linear), ~ mean(.x, na.rm = TRUE))) %>%
#   ungroup() %>%
#   mutate(
#     # Transform from linear scale to response scale
#     seeds.per.umbel = (1 / (1 + exp(-seeds.zinf.linear))) * exp(seeds.seed.linear),
#     # Total umbels per plant
#     seeds.total = seeds.per.umbel * phen.umbels
#   ) %>%
#   select(-c(seeds.zinf.linear, seeds.seed.linear)) %>%
#   mutate(
#     # Mean recruit size
#     recr.mean = predict(
#       r_t.y, allow.new.levels = TRUE, re.form = ~ 0, newdata = .
#     ),
#     # Get the number of seeds produced for each size grouping
#     p.size.cur = 0.1 * seeds.total * dnorm(x = size.nex, mean = recr.mean, sd = sigma.recr)
#   ) %>%
#   # Rename column
#   rename(size.prev = size)

head(ltre.kernel)
head(ltre.control.kernel)

# Do some formatting and export
# rbind(ltre.kernel, ltre.control.kernel) %>%
ltre.kernel %>%
  mutate(phen = phen.c + round(mean(seed$mean.phen))) %>%
  select(-c(phen.umbels, recr.mean, seeds.per.umbel, seeds.total, phen.c)) %>%
  write.csv(
    file = '03_construct_kernels/out/determinstic_reprod_kernel_phen_ltre.csv',
    row.names = FALSE, na = ''
  )


# ==========================================
# ------------------------------------------
# Do perturbations for sensitivity estimates
# ------------------------------------------
# ==========================================

### Parameters to perturb
# - Probability of flowering (zero inflation part of umbel mod)
#   - intercept
#   - **slope does not vary by treatment, so no need for LTRE**
# - Umbels produced (conditional part of umbel mod)
#   - **does not vary according to treatment, so no need for LTRE**
# - Umbel success (zero inflation part of seed set model)
#   - **does not vary according to treatment, so no need for LTRE**
# - Seed set (conditional part of seed set model)
#   - intercept
#   - slope
# - Recruit size
#   - intercept
# - Phenology
#   - phen effect on umbel success
#   - phen effect on seed set

# Overall number of effects: six to test here

perturb.list = vector(length = 6, mode = 'list')

delta = 0.0001

sigma.recr = summary(r_t.y)$sigma

# Parameter 1: probability of flowering intercept

u_s_s.ty$fit$par
# zi intercept is [3], drought term is [5], irr is [6]

# Rename to not put the year random effect in these predictions
perturb.list[[1]] = ltre.backbone %>%
  rename(year = Year) %>%
  mutate(
    # Umbel count
    phen.umbels = predict(
      u_s_s.ty, newdata = .,
      newparams = u_s_s.ty$fit$par %>%
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
      type = 'zlink'
    ),
    seeds.seed.linear =  predict(
      s_st.p_s.u.p2, newdata = ., allow.new.levels = TRUE, re.form = ~ 0,
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
      r_t.y, allow.new.levels = TRUE, re.form = ~ 0, newdata = .
    ),
    # Get the number of seeds produced for each size grouping
    p.size.cur = 0.1 * seeds.total * dnorm(x = size.nex, mean = recr.mean, sd = sigma.recr)
  ) %>%
  # Rename column
  rename(size.prev = size) %>%
  mutate(
    param = 'flow.int',
    orig.parval = case_when(
      trt %in% 'control' ~ u_s_s.ty$fit$par[3],
      trt %in% 'drought' ~ u_s_s.ty$fit$par[3] + u_s_s.ty$fit$par[5],
      trt %in% 'irrigated' ~ u_s_s.ty$fit$par[3] + u_s_s.ty$fit$par[6]
    )
  )

# Seed set intercept (conditional model)

s_st.p_s.u.p2$fit$par
# Intercept is [1], drought is [2], irr is [3]

perturb.list[[2]] = ltre.backbone %>%
  rename(year = Year) %>%
  mutate(
    # Umbel count
    phen.umbels = predict(
      u_s_s.ty, newdata = .,
      allow.new.levels = TRUE, re.form = ~ 0, type = 'response'
    )
  ) %>%
  rename(Year = year) %>%
  mutate(
    # Model predictions for seed set on linear (link) scale for averaging
    seeds.zinf.linear =  predict(
      s_st.p_s.u.p2, newdata = ., allow.new.levels = TRUE, re.form = ~ 0,
      type = 'zlink'
    ),
    seeds.seed.linear =  predict(
      s_st.p_s.u.p2, newdata = ., allow.new.levels = TRUE, re.form = ~ 0,
      newparams = s_st.p_s.u.p2$fit$par %>%
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
      r_t.y, allow.new.levels = TRUE, re.form = ~ 0, newdata = .
    ),
    # Get the number of seeds produced for each size grouping
    p.size.cur = 0.1 * seeds.total * dnorm(x = size.nex, mean = recr.mean, sd = sigma.recr)
  ) %>%
  # Rename column
  rename(size.prev = size) %>%
  mutate(
    param = 'seed.int',
    orig.parval = case_when(
      trt %in% 'control' ~ s_st.p_s.u.p2$fit$par[1],
      trt %in% 'drought' ~ s_st.p_s.u.p2$fit$par[1] + s_st.p_s.u.p2$fit$par[2],
      trt %in% 'irrigated' ~ s_st.p_s.u.p2$fit$par[1] + s_st.p_s.u.p2$fit$par[3]
    )
  )

# Seed set slope (conditional model) 

s_st.p_s.u.p2$fit$par
# slope is 4, drought is 9, irr is 10

perturb.list[[3]] = ltre.backbone %>%
  rename(year = Year) %>%
  mutate(
    # Umbel count
    phen.umbels = predict(
      u_s_s.ty, newdata = .,
      allow.new.levels = TRUE, re.form = ~ 0, type = 'response'
    )
  ) %>%
  rename(Year = year) %>%
  mutate(
    # Model predictions for seed set on linear (link) scale for averaging
    seeds.zinf.linear =  predict(
      s_st.p_s.u.p2, newdata = ., allow.new.levels = TRUE, re.form = ~ 0,
      type = 'zlink'
    ),
    seeds.seed.linear =  predict(
      s_st.p_s.u.p2, newdata = ., allow.new.levels = TRUE, re.form = ~ 0,
      newparams = s_st.p_s.u.p2$fit$par %>%
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
      r_t.y, allow.new.levels = TRUE, re.form = ~ 0, newdata = .
    ),
    # Get the number of seeds produced for each size grouping
    p.size.cur = 0.1 * seeds.total * dnorm(x = size.nex, mean = recr.mean, sd = sigma.recr)
  ) %>%
  # Rename column
  rename(size.prev = size) %>%
  mutate(
    param = 'seed.slope',
    orig.parval = case_when(
      trt %in% 'control' ~ s_st.p_s.u.p2$fit$par[4],
      trt %in% 'drought' ~ s_st.p_s.u.p2$fit$par[4] + s_st.p_s.u.p2$fit$par[9],
      trt %in% 'irrigated' ~ s_st.p_s.u.p2$fit$par[4] + s_st.p_s.u.p2$fit$par[10]
    )
  )


# 4: recruit size intercept

r_t.y$fit$par
# control is 1, drought 2, irr 3

perturb.list[[4]] = ltre.backbone %>%
  rename(year = Year) %>%
  mutate(
    # Umbel count
    phen.umbels = predict(
      u_s_s.ty, newdata = .,
      allow.new.levels = TRUE, re.form = ~ 0, type = 'response'
    )
  ) %>%
  rename(Year = year) %>%
  mutate(
    # Model predictions for seed set on linear (link) scale for averaging
    seeds.zinf.linear =  predict(
      s_st.p_s.u.p2, newdata = ., allow.new.levels = TRUE, re.form = ~ 0,
      type = 'zlink'
    ),
    seeds.seed.linear =  predict(
      s_st.p_s.u.p2, newdata = ., allow.new.levels = TRUE, re.form = ~ 0,
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
      newparams = r_t.y$fit$par %>%
        (function(x) {
          x[1] <- x[1] + delta
          return(x)
        })
    ),
    # Get the number of seeds produced for each size grouping
    p.size.cur = 0.1 * seeds.total * dnorm(x = size.nex, mean = recr.mean, sd = sigma.recr)
  ) %>%
  # Rename column
  rename(size.prev = size) %>%
  mutate(
    param = 'recr.int',
    orig.parval = case_when(
      trt %in% 'control' ~ r_t.y$fit$par[1],
      trt %in% 'drought' ~ r_t.y$fit$par[1] + r_t.y$fit$par[2],
      trt %in% 'irrigated' ~ r_t.y$fit$par[1] + r_t.y$fit$par[3]
    )
  )

# 5: phenology effects on umbel success

perturb.list[[5]] = ltre.backbone %>%
  rename(year = Year) %>%
  mutate(
    # Umbel count
    phen.umbels = predict(
      u_s_s.ty, newdata = .,
      allow.new.levels = TRUE, re.form = ~ 0, type = 'response'
    )
  ) %>%
  rename(Year = year) %>%
  # Perturb phen variable (for only the zinf term)
  mutate(phen.c = phen.c + delta) %>%
  mutate(
      # Model predictions for seed set on linear (link) scale for averaging
    seeds.zinf.linear =  predict(
      s_st.p_s.u.p2, newdata = ., allow.new.levels = TRUE, re.form = ~ 0,
      type = 'zlink'
    )
    # Reset the phen variable (so conditional is unaffected)
  ) %>%
  mutate(phen.c = phen.c - delta) %>%
  mutate(
    seeds.seed.linear =  predict(
      s_st.p_s.u.p2, newdata = ., allow.new.levels = TRUE, re.form = ~ 0,
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
      r_t.y, allow.new.levels = TRUE, re.form = ~ 0, newdata = .
    ),
    # Get the number of seeds produced for each size grouping
    p.size.cur = 0.1 * seeds.total * dnorm(x = size.nex, mean = recr.mean, sd = sigma.recr)
  ) %>%
  # Rename column
  rename(size.prev = size) %>%
  mutate(
    param = 'phen.succ',
    orig.parval = phen.c
  )

# 6: phen effects on seed set

perturb.list[[6]] = ltre.backbone %>%
  rename(year = Year) %>%
  mutate(
    # Umbel count
    phen.umbels = predict(
      u_s_s.ty, newdata = .,
      allow.new.levels = TRUE, re.form = ~ 0, type = 'response'
    )
  ) %>%
  rename(Year = year) %>%
  mutate(
    # Model predictions for seed set on linear (link) scale for averaging
    seeds.zinf.linear =  predict(
      s_st.p_s.u.p2, newdata = ., allow.new.levels = TRUE, re.form = ~ 0,
      type = 'zlink'
    )
  ) %>%
  # Perturb phen variable (for only the cond term)
  mutate(phen.c = phen.c + delta) %>%
  mutate(
    seeds.seed.linear =  predict(
      s_st.p_s.u.p2, newdata = ., allow.new.levels = TRUE, re.form = ~ 0,
      type = 'link'
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
      r_t.y, allow.new.levels = TRUE, re.form = ~ 0, newdata = .
    ),
    # Get the number of seeds produced for each size grouping
    p.size.cur = 0.1 * seeds.total * dnorm(x = size.nex, mean = recr.mean, sd = sigma.recr)
  ) %>%
  # Rename column
  rename(size.prev = size) %>%
  mutate(
    param = 'phen.seed',
    orig.parval = phen.c
  )

perturb.df = do.call(rbind, perturb.list) %>%
  mutate(mean.phen = phen.c + round(mean(seed$mean.phen))) %>%
  select(-c(phen.c, phen.umbels, seeds.per.umbel, seeds.total, recr.mean))

write.csv(
  perturb.df,
  file = '03_construct_kernels/out/deterministic_repr_coef_perturbation_phen.csv',
  row.names = FALSE, na = ''
)
