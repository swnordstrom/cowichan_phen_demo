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

# However, for our purposes we want to get the mean umbel emergence date by plant

phen.by.plant = phen.by.umbel %>%
  group_by(Plot, plantid, trt, Year) %>%
  summarise(phen.julian = mean(phen.julian), n.umbel = n()) %>%
  mutate(Year = factor(Year)) %>%
  ungroup()


head(phen.by.plant)


### --- Model selection ---------------------------------------------

# Model mean (from umbel-level observations) bud dates by each treatment
d_t = glmmTMB(
  phen.julian ~ trt + Year + (1 | Plot / plantid),
  data = phen.by.plant
)

# Compare with null model (hypothesis test)
d_0 = glmmTMB(
  phen.julian ~ Year + (1 | Plot / plantid),
  data = phen.by.plant
)

# Compare with treatment-year effect
d_ty = glmmTMB(
  phen.julian ~ trt * Year + (1 | Plot / plantid),
  data = phen.by.plant
)

AIC(d_0, d_t, d_ty) %>% mutate(daic = round(AIC - min(AIC), 2))

# Excellent - evidence of treatment effect, interaction not supported

# Compare with a random-effects model
d_t_r = glmmTMB(
  phen.julian ~ trt + (1 | Year) + (1 | Plot / plantid),
  data = phen.by.plant
)

AIC(d_t_r, d_t)
# Interestingly, random effects model does very poorly.

summary(d_t)
# Drought effect: 3.29 day advance
# Irrigation effect: 1.31 day delay

phen.by.plant %>%
  mutate(dt_resid = residuals(d_t)) %>%
  ggplot(aes(x = dt_resid)) +
  geom_histogram(aes(group = trt), alpha = 0.25, position = 'identity')
# Residuals look fine to me

# And the random effects
hist(ranef(d_t)$cond$`plantid:Plot`[,1])
# Maybe slightly skewed but otherwise fine.
hist(ranef(d_t)$cond$Plot[,1])

### Run model (selection code above)
# Full model:
d_t = glmmTMB(
  phen.julian ~ trt + Year + (1 | Plot / plantid),
  data = phen.by.plant
)

### --- Model means -------------------------------------------------

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

# Make folder for data output files if not already existing


# Export
write.csv(trt.mean.buddates, row.names = FALSE, '03_construct_kernels/out/phen_treatment_means.csv')

cat('Exported mean phenology treatments.\n')

### --- Annual means -------------------------------------------------

# Estimate means from model
trt.annual.buddates = expand.grid(trt = c('control', 'drought', 'irrigated'), Year = factor(2021:2024)) %>%
  # For each combo of treatment and year (above), estimate mean from model
  mutate(
    mean.bud = predict(
      d_t, re.form = ~ 0, allow.new.levels = TRUE,
      newdata = expand.grid(trt = c('control', 'drought', 'irrigated'), Year = factor(2021:2024))
    )
  )

# Export
write.csv(trt.annual.buddates, row.names = FALSE, '03_construct_kernels/out/phen_annual_treatment_means.csv')

### --- Summary statistics -------------------------------------------

# Confidence intervals on effect sizes
round(summary(d_t)$coefficients$cond['trtdrought',1] + c(-1, 0,  1) * 1.96 * summary(d_t)$coefficients$cond['trtdrought',2], 2)
round(summary(d_t)$coefficients$cond['trtirrigated',1] + c(-1, 0,  1) * 1.96 * summary(d_t)$coefficients$cond['trtirrigated',2], 2)

# Variance explained by treatment (check indexing every time re-run!)
vars_t = exp(2 * d_t$fit$par[7:9])
vars_0 = exp(2 * d_0$fit$par[5:7])

# Pseudo-R^2
# (one source for this, https://web.pdx.edu/~newsomj/mlrclass/ho_r2.pdf)
# (see: Snijders and Bosker 1999)
# 1 - (sum(vars_t) / sum(vars_0))
# 6.3%

x0 = model.matrix(~ Year, d_0$frame)
b0 = fixef(d_0)$cond
vars_0 = exp(2 * d_0$fit$par[5:7])

colnames(x0) == names(b0)

sigma2_f0 = var(x0 %*% as.matrix(b0))

rsq0 = sigma2_f0 / (sigma2_f0 + sum(vars_0))
# Year explains 27.4%

x1 = model.matrix(~ trt + Year, d_t$frame)
b1 = fixef(d_t)$cond
vars_t = exp(2 * d_t$fit$par[7:9])

colnames(x1) == names(b1)

sigma2_f1 = var(x1 %*% as.matrix(b1))

rsq1 = sigma2_f / (sigma2_f + sum(vars_t))
# Adding years gives gives 32.6%
# (but it includes the year effects...)
# (how to describe this...)

(rsq1 - rsq0) / (1 - rsq0)

# As expected, minimal change to individual- and among-plant plots
# Reduction in spatial (plot-level) variance:
# 1 - (vars_t[3] / vars_0[3])
# 67%

# Looking at year effects
(as.Date(d_t$fit$par[1] + c(0, d_t$fit$par[4:6]))) %>% range()

### Likelihood ratio tests on effect sizes
# anova(d_0, d_t) # p = 0.0036
# anova(d_t, d_ty) # p = 0.04338
# (pseudo-R^2 is 0.69)
