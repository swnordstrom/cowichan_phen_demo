# Script for creating growth/survival kernels for IPMs
# Reads in processed datasets for growth and survival,
# outputs one kernel for treatment
# No climate in these kernels so far
# - init SN 4 Mar 2024

# ------------------------------------------------
# Setup
# - Packages, namespace
# - Read in and prepare raw data

# Load packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(glmmTMB) # won't work with automation, lmao
library(lme4)

# Load namespace
rm(list = ls())

# Load in processed demo data
# demo = read.csv('01_data_cleaning/out/demo_imputed_survival.csv')
demo = read.csv('01_data_cleaning/out/Lomatium_demo_2016-2024.csv')

# Create demo dataset for survival
# we'll say that survival in year $t$ means surviving from $t$ to $t+1$
demo.for.surv = merge(
  # x - status in year $t$ - only interested in plants that survived in t
  x = demo %>% 
    filter(surv) %>% 
    select(Year, No.leaves, Leaf.length, No.umbels, demo.note, proc.note, plantid, Plot, trt),
  # y - whether plant was observed alive in year $t+1$
  y = demo %>%
    mutate(p.year = Year - 1) %>%
    select(p.year, plantid, surv, No.leaves, Leaf.length),
  by.x = c('Year', 'plantid'), by.y = c('p.year', 'plantid'),
  # suffixes are t and t+1
  suffixes = c('.t', '.tp1')
)

# Further subset demo with *only plants with usable sizes* prior to the transition
# (I imagine this is what will be used in the final analysis)
# This also means cutting out the 2016 plants because their size measures are not good
demo.for.surv.sizes = demo.for.surv %>%
  # Remove plants iwth no measurements and plants in 2016
  filter(!is.na(No.leaves.t) & !is.na(Leaf.length.t) & No.leaves.t > 0, Year > 2016) %>%
  # add a size column
  mutate(size.t = log(No.leaves.t * Leaf.length.t))

# Get a dataset of demo for growth - subsetting the data frame above
# to get this, we need plants that survived *and* have measurements
demo.for.growth = demo.for.surv.sizes %>%
  filter(!is.na(No.leaves.tp1) & !is.na(Leaf.length.tp1) & No.leaves.tp1 > 0) %>%
  # add a size column
  mutate(size.tp1 = log(No.leaves.tp1 * Leaf.length.tp1))

# Remove 2023 from surv because of NAs (plants seen in 2023 but not 2024 -
# unsure if alive or just missed)
demo.for.surv = demo.for.surv %>% 
  filter(!Year %in% 2023) %>%
  # Change year to factor
  mutate(Year = factor(Year))

demo.for.surv.sizes = demo.for.surv.sizes %>% 
  filter(!Year %in% 2023) %>%
  # Change year to factor
  mutate(Year = factor(Year))

# Remove plot 2 in 2024 from growth (because of censoring
demo.for.growth = demo.for.growth %>% filter(!(Year %in% 2023 & Plot %in% 2)) %>%   # Change year to factor
  mutate(Year = factor(Year))

# ------------------------------------------------
# Fit survival models

# First question: do we need size and/or year
# (I suspect we will)

s_0 = glmmTMB(
  formula = surv ~ (1 | Plot),
  family = 'binomial',
  data = demo.for.surv.sizes
)

s_y = glmmTMB(
  formula = surv ~ (1 | Year) + (1 | Plot),
  family = 'binomial',
  data = demo.for.surv.sizes
)

s_s = glmmTMB(
  formula = surv ~ size.t + (1 | Plot),
  family = 'binomial',
  data = demo.for.surv.sizes
)

s_sy = glmmTMB(
  formula = surv ~ size.t + (1 | Year) + (1 | Plot),
  family = 'binomial',
  data = demo.for.surv.sizes
)

AIC(s_0, s_s, s_y, s_sy)
# Unsurprisingly, size effect present.
# BUT, amazingly, the year effect is negligible (!)

# # Also check for random effects vs. fixed effects?
# 
# s_s_fix = glmmTMB(
#   formula = surv ~ size.t + Year + (1 | Plot),
#   family = 'binomial',
#   data = demo.for.surv.sizes
# )
# 
# AIC(s_s, s_s_fix)
# # Hmm... when removing individual-level random effects in survival model,
# # random effect-year model out-performs the fixed effect year model
# # (good!)

# Models to test here:
# - treatment effect
# - size-treatment interaction
# - year-varying size effect
# - (model will not converge for treatment effect varying by year)

surv.mod.forms = c(
  'surv ~ size.t',
  'surv ~ size.t + trt',
  'surv ~ size.t * trt'
  # old code for models with year effects (no longer supported)
  # 'surv ~ size.t + (1 | Year)',
  # 'surv ~ size.t + trt + (1 | Year)',
  # 'surv ~ size.t * trt + (1 | Year)',
  # 'surv ~ (size.t | Year)',
  # 'surv ~ trt + (size.t | Year)' ,
  # 'surv ~ size.t + trt + (1 | Year) + (1 | Year:trt)'
) %>%
  paste('+ (1 | Plot)')

surv.mod.list = vector('list', length = length(surv.mod.forms))

for (j in 1:length(surv.mod.forms)) {
  print(j)
  surv.mod.list[[j]] = glmer(
    formula = surv.mod.forms[j],
    family = 'binomial',
    data = demo.for.surv.sizes,
    control = glmerControl(optimizer = 'bobyqa')
  )
}
# Convergence issues... (come on man how does this happen?)

data.frame(
  AIC = sapply(surv.mod.list, AIC) %>% unlist(),
  form = surv.mod.forms
) %>%
  mutate(daic = round(AIC - min(AIC), 2))
# Best supported model has no treatment effects present


# Plot model outputs:
# expand.grid(size.t = (5:60)/10, Year = factor(2017:2022)) %>%
expand.grid(size.t = (5:60)/10) %>%
  mutate(
    # pred_y = predict(
    #   s_s, newdata = ., type = 'response',
    #   re.form = ~ (1 | Year), allow.new.levels = TRUE
    # ),
    pred = predict(
      s_s, newdata = ., type = 'response',
      re.form = ~ 0, allow.new.levels = TRUE
    )
  ) %>%
  # pivot_longer(c(pred_y, pred), names_to = 'model', values_to = 'pred') %>%
  ggplot(aes(x = size.t, y = pred)) +
  geom_point(
    data = demo.for.surv.sizes,
    aes(x = size.t, y = as.numeric(surv)),
    size = 3, alpha = 0.1
  ) +
  # geom_line(aes(linetype = model)) +
  geom_line() +
  labs(x = 'size', y = 'survival') +
  facet_wrap(~ Year)
# Looks good/plausible to me

# ------------------------------------------------
# Fit growth models

# Size effect of course will be here.

grow.mod.forms = c(
  'size.tp1 ~ size.t + (1 | Year)',
  'size.tp1 ~ size.t + trt + (1 | Year)',
  'size.tp1 ~ size.t * trt + (1 | Year)',
  'size.tp1 ~ (size.t | Year)',
  'size.tp1 ~ trt + (size.t | Year)',
  'size.tp1 ~ size.t + trt + (1 | Year) + (1 | Year:trt)',
  'size.tp1 ~ trt + (size.t | Year) + (1 | Year:trt)',
  'size.tp1 ~ trt + (size.t | Year) + (size.t | Year:trt)',
  'size.tp1 ~ size.t * trt + (1 | Year) + (1 | Year:trt)'
) %>%
  paste('+ (1 | Plot / plantid)')

grow.mod.list = vector('list', length = length(grow.mod.forms))

for (j in 1:length(grow.mod.forms)) {
  print(j)
  grow.mod.list[[j]] = lmer(
    formula = grow.mod.forms[j],
    data = demo.for.growth,
    control = lmerControl(optimizer = 'bobyqa')
  )
}
# Quick! Thank you Gaussian link!

data.frame(
  AIC = sapply(grow.mod.list, AIC) %>% unlist(),
  form = grow.mod.forms
) %>%
  mutate(daic = round(AIC - min(AIC), 2)) %>%
  arrange(daic)

# Wow. Tear-varying treatment effect dominates here.
# Only vital rate I have seen like this!
# Furthermore definitely no support for year-varying size effects

g_st_ty = lmer(
  size.tp1 ~ size.t * trt + (1 | Year) + (1 | Year:trt) + (1 | Plot / plantid),
  data = demo.for.growth,
  control = lmerControl(optimizer = 'bobyqa')
)

summary(g_st_ty)
# Treatment effects are incredibly small, but...
# there is a slight benefit to both treatments... lol

# Residuals:
demo.for.growth %>%
  mutate(resid = residuals(g_st_ty)) %>%
  ggplot(aes(x = size.t, y = resid)) +
  geom_point() +
  facet_wrap( ~ Year)
# Some negative outliers (negative tail looking bleh), and 2018 is weird, 2021 a
# little strange, maybe a handful of outliers in 2022
# otherwise okay
# Not worried about heteroskedasticity here

# The year-treatment random effects
ranef(g_st_ty)$`Year:trt` %>%
  mutate(year.trt = row.names(.)) %>%
  separate(year.trt, into = c('year', 'trt'), sep = ':') %>%
  rename(intercept = `(Intercept)`) %>%
  ggplot(aes(x = year, y = intercept, colour = trt)) +
  geom_point() +
  geom_line(aes(group = trt)) +
  scale_colour_manual(values = c('black', 'red', 'blue'))
# Not seeing clear patterns here
# Maybe some negative temporal autocorrelation?
# Controls are pretty noisy though

# Model predictions, year-varying first
expand.grid(
  size.t = (9:60)/10, 
  Year = factor(2017:2023), 
  trt = c('control', 'drought', 'irrigated')
) %>%
  mutate(
    pred_y = predict(
      g_st_ty,
      newdata = .,
      re.form = ~ (1 | Year),
      type = 'response',
      allow.new.levels = TRUE
    ),
    pred_mean = predict(
      g_st_ty,
      newdata = .,
      re.form = ~ 0,
      type = 'response',
      allow.new.levels = TRUE
    )
  ) %>%
  pivot_longer(c(pred_mean, pred_y), names_to = 'model', values_to = 'pred') %>%
  ggplot(aes(x = size.t, y = pred)) +
  geom_point(
    data = demo.for.growth,
    aes(x = size.t, y = size.tp1, colour = trt),
    inherit.aes = FALSE,
    size = 3, alpha = 0.1
  ) +
  geom_segment(aes(x = .9, xend = 6, y = .9, yend = 6), linetype = 2, colour = 'gray77') +
  geom_line(aes(colour = trt, linetype = model)) +
  scale_colour_manual(values = c('black', 'red', 'blue')) +
  facet_wrap(~ Year)
# That 2020 fit looks bad... maybe it's a plot-variance effect?
# Yeesh.

# Ignoring year effects
expand.grid(
  size.t = (9:60)/10, 
  Year = factor(2017:2023), 
  trt = c('control', 'drought', 'irrigated')
) %>%
  mutate(
    pred = predict(
      g_st_ty,
      newdata = .,
      re.form = ~ 0,
      type = 'response',
      allow.new.levels = TRUE
    )
  ) %>%
  # pivot_longer(c(pred_mean, pred_y), names_to = 'model', values_to = 'pred') %>%
  ggplot(aes(x = size.t, y = pred)) +
  geom_point(
    data = demo.for.growth,
    aes(x = size.t, y = size.tp1, colour = trt),
    inherit.aes = FALSE,
    size = 3, alpha = 0.1
  ) +
  geom_segment(aes(x = .9, xend = 6, y = .9, yend = 6), linetype = 2, colour = 'gray77') +
  geom_line(aes(colour = trt)) +
  scale_colour_manual(values = c('black', 'red', 'blue'))
# Still looks from this like the overall effect of drought is... positive
# Larger plants will shrink less with abundant water...

expand.grid(
  size.t = (9:60)/10, 
  Year = factor(2017:2023), 
  trt = c('control', 'drought', 'irrigated')
) %>%
  mutate(
    pred_y = predict(
      g_st_ty,
      newdata = .,
      re.form = ~ (1 | Year),
      type = 'response',
      allow.new.levels = TRUE
    )
  ) %>%
  ggplot(aes(x = size.t, y = pred_y)) +
  geom_segment(aes(x = .9, xend = 6, y = .9, yend = 6), linetype = 2, colour = 'gray77') +
  geom_line(aes(colour = trt, group = interaction(Year, trt))) +
  scale_colour_manual(values = c('black', 'red', 'blue'))
# Hmm... oay not a super cool plot...

# ------------------------------------------------
# Construct kernels

# Residual variance in growth models
grow.sd = summary(g_s_ty)$sigma

# Get a scaffold
grow.surv.kernel = expand.grid(
  size.t = (5:60)/10,
  size.tp1 = (5:60)/10,
  trt = c('control', 'drought', 'irrigated')
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
    pred.grow.mean = predict(
      newdata = .,
      object = g_s_ty, type = 'response',
      re.form = ~ 0, allow.new.levels = TRUE
    )
  ) %>%
  # Predicted distribution of sizes in next time step
  mutate(p.grow.size.tp1 = 0.1 * dnorm(size.tp1, mean = pred.grow.mean, sd = grow.sd)) %>%
  # Combine all together to get overall size distribution in next time step
  mutate(p.size.tp1 = pred.surv * p.grow.size.tp1)

ggplot(grow.surv.kernel, aes(x = size.t, y = size.tp1)) +
  geom_tile(aes(fill = p.size.tp1)) +
  scale_y_reverse() +
  facet_wrap(~ trt)
# Neat

# write.csv(
#   grow.surv.kernel,
#   file = '03_construct_kernels/out/test_survgrow_kernel.csv',
#   row.names = FALSE
# )
