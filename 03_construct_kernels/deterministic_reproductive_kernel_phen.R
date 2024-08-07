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

# Read in data
all.data = merge(
  x = read.csv('01_data_cleaning/out/demo_phen_seed_2016-2024_final.csv'),
  y = read.csv('00_raw_data/plot_treatments.csv'),
  by.x = 'Plot', by.y = 'plot'
)

# Dataset for flowering (additional processing below)
demo.flow = all.data %>%
  # Get only plants that are:
  # - known/estimated to be alive
  # - Has size measurements (from demography) and has non-zero leaf counts and lengths (for size measurements)
  # - Is not missing umbel counts (because we're subsetting this for umbel count analysis)
  filter(
    in.demo, Year > 2016, surv,
    !is.na(No.leaves) & !is.na(Leaf.length) & No.leaves > 0 & Leaf.length > 0,
    !is.na(No.umbels)
  ) %>%
  # Get a size column
  mutate(size = log(No.leaves * Leaf.length))

# Dataset for estimating recruit size distribution
demo.recr = all.data %>%
  arrange(Year) %>%
  # Want demo records (to get size distribution)
  filter(in.demo) %>%
  # Give me only the first sighting for each plant
  distinct(plantid, .keep_all = TRUE) %>%
  # Give me only records after 2017 (these are more likely to be seen for the
  # first time)
  filter(Year > 2017) %>%
  # Give me only records with leaf counts and lengths (for estimating size)
  filter(!is.na(No.leaves), !is.na(Leaf.length)) %>%
  # Give me plants that did not flower (plants unlikely to flower in first year)
  filter(!No.umbels & !is.na(No.umbels)) %>%
  # Okay... going to assume that only plants with one leaf are new recruits
  filter(No.leaves < 2) %>%
  mutate(size = log(Leaf.length))

# Finalize the flowering demo dataset 
demo.flow = demo.flow %>%
  # (currently there is one row per umbel for plants in the seed set dataset -
  # give me just one record per plant)
  distinct(Year, plantid, .keep_all = TRUE) %>%
  # Remove plot 2 in 2024 (missing a lot of vegetative plant sizes)
  filter(!(Plot %in% 2 & Year %in% 2024))

# Umbel-level budding phenology data
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

# All plants in phen, seed, and demo
seed = all.data %>% 
  filter(in.phen, in.seed, in.demo) %>%
  # Need size data
  filter(!is.na(No.leaves), !is.na(Leaf.length), No.leaves > 0, Leaf.length > 0) %>%
  # Hmm... okay, separae() and across() doesn't work, so I guess I can merge
  # this in with the data frame above, summarised by mean phen date per plant
  merge(
    y = phen %>% 
      group_by(plantid, Year) %>% 
      summarise(
        mean.phen = mean(phen.julian, na.rm = TRUE), 
        sept.phen = length(unique(phen.julian))
      ) %>%
      ungroup(),
    all.x = TRUE, all.y = FALSE
  )

seed = rbind(
  seed,
  seed %>%
    group_by(Year = factor(Year), plantid) %>%
    mutate(miss.umbel = ifelse(phen.umbels < n(), 0, phen.umbels - n())) %>%
    ungroup() %>%
    distinct(Year, plantid, .keep_all = TRUE) %>%
    uncount(miss.umbel) %>%
    mutate(no.seeds = 0)
) %>%
  mutate(phen.c = mean.phen - round(mean(mean.phen))) %>%
  mutate(size = log(No.leaves * Leaf.length))


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

summary(d_t)

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

# Get mean phen bud date for each year
# (because year is categorical in phenology model, gettting one estimate for each year and then averaging over treatment)
ltre.backbone = expand.grid(
  size = (5:60)/10,
  size.nex = (5:60)/10,
  trt = c('control', 'drought', 'irrigated'),
  Year = factor(2021:2024)
) %>%
  # Model predictions (from d_t model)
  mutate(mean.phen = predict(d_t, newdata = ., re.form = ~ 0, allow.new.levels = TRUE)) %>%
  # Get mean across treatment 
  group_by(trt) %>%
  # (using `mutate` instead of summarise because we'll still want the years for
  # the seed set model)
  mutate(mean.phen = mean(mean.phen)) %>%
  ungroup() %>%
  # Center (because predictor in phen-seed model is centered)
  mutate(phen.c = mean.phen - round(mean(seed$mean.phen)))


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

# Do the same thing but for drought and irrigation but at teh control mean phenology
ltre.control.backbone = ltre.backbone %>%
  select(-phen.c) %>%
  pivot_wider(names_from = trt, values_from = mean.phen) %>%
  mutate(drought = control, irrigated = control) %>%
  # (We don't need the control estimates)
  select(-control) %>%
  pivot_longer(c(drought, irrigated), names_to = 'trt', values_to = 'mean.phen') %>%
  mutate(phen.c = mean.phen - round(mean(seed$mean.phen)))

ltre.control.kernel = ltre.control.backbone %>%
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
    # Model predictions at *control means*
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

head(ltre.kernel)
head(ltre.control.kernel)

# Do some formatting and export
rbind(ltre.kernel, ltre.control.kernel) %>%
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

delta = 0.0001



