#-------------------------------------------------------------------
# Model selection script for somatic (survival, growth) vital rate models
# Involves following models
# - Probability of survival
# - Growth function (growth, conditioned on survival)
#   - Growth for all surviving plant (vegetative)
#   - Growth for flowering plants, conditioned on phenology (flowering)
# - Size of recruits
#-------------------------------------------------------------------

#-------------------------------------------------------------------
# Setup

##### Load packages

library(glmmTMB)
library(ggplot2)
library(lme4)
library(dplyr)
library(tidyr)
# library(cowplot)

##### Clear namespace
rm(list = ls())

##### Run script for loading/processing data

source('03_construct_kernels/prepare_demo_data_growsurv.R')

#-------------------------------------------------------------------
# Fit models for probability of survival per plant

# Datasource is `demo.surv.sizes` which has observations of plants with recorded
# size and their fate in the following year (including plants with imputed
# survival)

# Testing to see if survival is dependent on size and for best random effects terms:

s_0 = glmmTMB(
  formula = surv ~ (1 | Plot),
  family = 'binomial',
  data = demo.surv.sizes
)

s_y = glmmTMB(
  formula = surv ~ (1 | surv.year) + (1 | Plot),
  family = 'binomial',
  data = demo.surv.sizes
)

s_s = glmmTMB(
  formula = surv ~ size.prev + (1 | Plot),
  family = 'binomial',
  data = demo.surv.sizes
)

s_s_y = glmmTMB(
  formula = surv ~ size.prev + (1 | surv.year) + (1 | Plot),
  family = 'binomial',
  data = demo.surv.sizes
)

AIC(s_0, s_y, s_s, s_s_y) %>% mutate(daic = round(AIC - min(AIC), 2))

# Surprisingly, no year effects!

summary(s_s)
# And unsurprisingly the effect of previous size is positive.

# Models to test:
# - Treatment effect
# - Size-treatment interaction

s_t = glmmTMB(
  formula = surv ~ size.prev + trt + (1 | Plot),
  family = 'binomial',
  data = demo.surv.sizes
)

s_ts = glmmTMB(
  formula = surv ~ size.prev * trt + (1 | Plot),
  family = 'binomial',
  data = demo.surv.sizes
)

AIC(s_s, s_t, s_ts) %>% mutate(daic = round(AIC - min(AIC), 2))
# No treatment effects.

# Plot of this model:
surv.preds = data.frame(size.prev = (5:60) / 10) %>%
  mutate(pred = predict(s_s, re.form = ~ 0, newdata = ., allow.new.levels = TRUE, type = 'response'))

surv.preds %>%
  ggplot(aes(x = size.prev, y = pred)) +
  geom_point(
    data = demo.surv.sizes,
    aes(y = as.numeric(surv)),
    position = position_jitter(height = 0.05),
    size = 2, alpha = 0.1
  ) +
  geom_line()

# Plausible.

# Table with AIC for all:
AIC(s_0, s_s, s_y, s_s_y, s_t, s_ts) %>%
  mutate(daic = AIC - min(AIC)) %>%
  mutate(across(c(AIC, daic), ~ round(., 2)))


#-------------------------------------------------------------------
# Fit models for growth conditioned on survival

# 
g_0 = glmmTMB(
  size.cur ~ size.prev + (1 | Plot / plantid),
  data = demo.grow
)

g_y = glmmTMB(
  size.cur ~ size.prev + (1 | prev.year) + (1 | Plot / plantid),
  data = demo.grow
)

AIC(g_y, g_0)
# Yep. Support for a year effect.

# Treatment effect

g_t = glmmTMB(
  size.cur ~ size.prev + trt + (1 | prev.year) + (1 | Plot / plantid),
  data = demo.grow
)

AIC(g_t, g_y)
# Wow... no treatment effect.

# There are a couple of other things to try:
# - treatment-year
# - treatment-size
# - size-year

g_ty = glmmTMB(
  size.cur ~ size.prev + trt + (1 | prev.year) + (1 | prev.year:trt) + (1 | Plot / plantid),
  data = demo.grow
)

g_st.ty = glmmTMB(
  size.cur ~ size.prev + size.prev * trt + (1 | prev.year) + (1 | prev.year:trt) + (1 | Plot / plantid),
  data = demo.grow
)

g_st = glmmTMB(
  size.cur ~ size.prev * trt + (1 | prev.year) + (1 | Plot / plantid),
  data = demo.grow
)

g_sy = glmmTMB(
  size.cur ~ size.prev + (size.prev | prev.year) + (1 | Plot / plantid),
  data = demo.grow
)


# Final AIC table
AIC(g_0, g_y, g_t, g_ty, g_st.ty, g_st, g_sy) %>% 
  mutate(daic = AIC - min(AIC)) %>%
  mutate(across(c(AIC, daic), ~ round(., 2)))

# next:
# - visualization of growth function
# - phen/growth
# - flesh out appendix section on growth/surv models
# - recruit size models + tables in this script
