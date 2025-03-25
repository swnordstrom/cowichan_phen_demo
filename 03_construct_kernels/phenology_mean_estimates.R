# --- Setup ---------------------------------------------------------

library(ggplot2)
library(dplyr)
library(tidyr)
library(glmmTMB)

rm(list = ls())

### --- Read in all data ---

all.data = merge(
  x = read.csv('01_data_cleaning/out/demo_phen_seed_2016-2024_final.csv'),
  y = read.csv('00_raw_data/plot_treatments.csv'),
  by.x = 'Plot', by.y = 'plot'
)

### --- Extract umbel phenology from phen data ----------------------

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

### --- Means -------------------------------------------------------

# Model mean (from umbel-level observations) bud dates by each treatment
d_t = glmmTMB(
  phen.julian ~ trt + Year + (1 | Plot / plantid),
  data = phen.by.umbel
)

# Compare with null model (hypothesis test)
d_0 = glmmTMB(
  phen.julian ~ Year + (1 | Plot / plantid),
  data = phen.by.umbel
)

# Compare with treatment-year effect
d_ty = glmmTMB(
  phen.julian ~ trt * Year + (1 | Plot / plantid),
  data = phen.by.umbel
)

AIC(d_0, d_t, d_ty) %>%
  mutate(daic = round(AIC - min(AIC), 2))

# Excellent - evidence of treatment effect, interaction not supported

# Estimate means from model
trt.mean.buddates = expand.grid(trt = c('control', 'drought', 'irrigated'), Year = factor(2021:2024)) %>%
  # For each combo of treatment and year (above), estimate mean from model
  mutate(
    mean.bud = predict(
      d_t, re.form = ~ 0, allow.new.levels = TRUE,
      newdata = expand.grid(trt = c('control', 'drought', 'irrigated'), Year = factor(2021:2024))
    )
  ) %>%
  # Now average across years to get a mean
  group_by(trt) %>%
  summarise(mean.phen = mean(mean.bud))

trt.mean.buddates

# Export
write.csv(trt.mean.buddates, row.names = FALSE, '03_construct_kernels/phen_treatment_means.csv')


### Summary statistics

# Confidence intervals on effect sizes
summary(d_t)$coefficients$cond['trtdrought',1] + c(-1, 1) * 1.96 * summary(d_t)$coefficients$cond['trtdrought',2]
summary(d_t)$coefficients$cond['trtirrigated',1] + c(-1, 1) * 1.96 * summary(d_t)$coefficients$cond['trtirrigated',2]

# Variance explained by treatment (check indexing every time re-run!)
vars_t = exp(c(1, 2, 2) * d_t$fit$par[7:9])
vars_0 = exp(c(1, 2, 2) * d_0$fit$par[5:7])

# Pseudo-R^2
# (one source for this, https://web.pdx.edu/~newsomj/mlrclass/ho_r2.pdf)
# (see: Snijders and Bosker 1998)
1 - (sum(vars_t) / sum(vars_0))
# (is very small)

# As expected, minimal change to individual- and among-plant plots
# Reduction in spatial (plot-level) variance:
1 - (vars_t[3] / vars_0[3])

# Looking at year effects
as.Date(d_t$fit$par[1] + c(0, d_t$fit$par[4:6]))
