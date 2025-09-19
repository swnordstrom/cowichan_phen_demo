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

#----------------------------------------------
# Now, the survival-growth model

# Check to make sure there are no year or treatment effects in this subset
sp_0 = glmmTMB(
  surv ~ size.prev + (1 | Plot),
  family = 'binomial',
  data = demo.surv.sizes %>% filter(!is.na(phen.c))
)

sp_y = glmmTMB(
  surv ~ size.prev + surv.year + (1 | Plot),
  family = 'binomial',
  data = demo.surv.sizes %>% filter(!is.na(phen.c))
)

sp_t = glmmTMB(
  surv ~ size.prev + trt + (1 | Plot),
  family = 'binomial',
  data = demo.surv.sizes %>% filter(!is.na(phen.c))
)

AIC(sp_0, sp_y, sp_t)
# Good - no evidence of either (not testing for interactions)

# Now, testing for phenology effects

sp_m = glmmTMB(
  surv ~ size.prev + phen.c + (1 | Plot),
  family = 'binomial',
  data = demo.surv.sizes %>% filter(!is.na(phen.c))
)

sp_m2 = glmmTMB(
  surv ~ size.prev + poly(phen.c, 2) + (1 | Plot),
  family = 'binomial',
  data = demo.surv.sizes %>% filter(!is.na(phen.c))
)


AIC(sp_0, sp_m, sp_m2) %>% mutate(daic = round(AIC - min(AIC), 2))
# No evidence of a phenology effect on growth.

# All models in one table:
AIC(sp_0, sp_y, sp_t, sp_m, sp_m2) %>%
  mutate(daic = AIC - min(AIC)) %>%
  mutate(across(c(AIC, daic), ~ round(., 2)))

#-------------------------------------------------------------------
# Fit models for growth conditioned on survival

#----------------------------------------------
# First, non-phenology model

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
# No treatment effect.

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

### Visualizing predictions from final model:
growth.viz = expand.grid(
  size.prev = (5:60)/10, 
  trt = c('drought', 'irrigated', 'control'), 
  prev.year = unique(demo.grow$prev.year)
) %>%
  mutate(Plot = 0, plantid = 0) %>%
  mutate(
    pred.mean = predict(g_st.ty, newdata = ., re.form = ~ 0, allow.new.levels = TRUE),
    pred.year = predict(g_st.ty, newdata = ., re.form = NULL, allow.new.levels = TRUE)
  )

growth.viz %>%
  ggplot(aes(x = size.prev)) +
  geom_point(data = demo.grow, aes(y = size.cur, colour = trt), size = 2, alpha = 0.1) +
  geom_line(aes(y = pred.mean, group = trt, colour = trt), linewidth = 1.2) +
  geom_line(
    aes(y = pred.year, group = interaction(trt, prev.year), colour = trt), 
    linewidth = 1.2, linetype = 2
  ) +
  scale_colour_manual(values = c('black', 'goldenrod', 'dodgerblue')) +
  facet_wrap(~ prev.year)
# Dotted lines: annual means (so there is certainly some annual variation...)
# Solid lines: overall means


#----------------------------------------------
# Now, the phenology-growth model

gp_0a = glmmTMB(
  size.cur ~ size.prev + prev.year + (1 | Plot / plantid),
  data = demo.grow %>% filter(!is.na(phen.c))
)

gp_0b = glmmTMB(
  size.cur ~ size.prev * prev.year + (1 | Plot / plantid),
  data = demo.grow %>% filter(!is.na(phen.c))
)

AIC(gp_0a, gp_0b)
# Okay - null model will have size-year interaction (different slopes for each year)

# Treatment effects if they are here
gp_t = glmmTMB(
  size.cur ~ size.prev * prev.year + trt + (1 | Plot / plantid),
  data = demo.grow %>% filter(!is.na(phen.c))
)

gp_ty = glmmTMB(
  size.cur ~ size.prev * prev.year + trt * prev.year + (1 | Plot / plantid),
  data = demo.grow %>% filter(!is.na(phen.c))
)

gp_t.sy = glmmTMB(
  size.cur ~ size.prev * prev.year + trt * size.prev + (1 | Plot / plantid),
  data = demo.grow %>% filter(!is.na(phen.c))
)
  
AIC(gp_0b, gp_t, gp_ty, gp_t.sy) %>% mutate(daic = round(AIC - min(AIC), 2))
# Same model structure as before (good).

gp_m = glmmTMB(
  size.cur ~ size.prev * prev.year + trt * prev.year + phen.c + (1 | Plot / plantid),
  data = demo.grow %>% filter(!is.na(phen.c))
)

gp_m2 = glmmTMB(
  size.cur ~ size.prev * prev.year + trt * prev.year + poly(phen.c, 2) + (1 | Plot / plantid),
  data = demo.grow %>% filter(!is.na(phen.c))
)

AIC(gp_ty, gp_m, gp_m2) %>% mutate(daic = round(AIC - min(AIC), 2))

# Just in case, test for a year-varying phenology effect
gp_my = glmmTMB(
  size.cur ~ size.prev * prev.year + trt * prev.year + phen.c * prev.year + (1 | Plot / plantid),
  data = demo.grow %>% filter(!is.na(phen.c))
)

AIC(gp_m, gp_my)
# Year-effects are n.s.

# Table with all AICs:
AIC(
  gp_0a, gp_0b,
  gp_t, gp_ty, gp_t.sy,
  gp_m, gp_m2, gp_my
) %>%
  mutate(daic = AIC - min(AIC)) %>%
  mutate(across(c(AIC, daic), ~ round(., 2)))

summary(gp_m)


#-------------------------------------------------------------------
# Fit models for size of recruits

# (first need to run a script to prepare recruit data)
# (running this script clears the namespace)
source('03_construct_kernels/prepare_demo_data_repr.R')

# Set year to factor (probably not needed but still good to be safe)
demo.recr$Year = factor(demo.recr$Year)

# Testing for best random effects
r_0 = glmmTMB(size ~ 1, data = demo.recr)
r_p = glmmTMB(size ~ (1 | Plot), data = demo.recr)

AIC(r_0, r_p)
# Yes, plot-level random effects

# Year-level random effects
r_y = glmmTMB(size ~ (1 | Plot) + (1 | Year), data = demo.recr)

AIC(r_y, r_p)
# Yes, year-level random effects

# Test for effects of treatment
r_t.y = glmmTMB(size ~ trt + (1 | Year) + (1 | Plot), data = demo.recr)

AIC(r_t.y, r_y)
# Yes (although not incredibly strong)

# Test for year-varying treatment effects
r_ty = glmmTMB(size ~ trt + (1 | trt:Year) + (1 | Year) + (1 | Plot), data = demo.recr)

AIC(r_ty, r_t.y)
# No evidence of year-varying treatment effects

# Whole table:
AIC(r_0, r_p, r_y, r_t.y, r_ty) %>%
  mutate(daic = AIC - min(AIC)) %>%
  mutate(across(c(AIC, daic), ~ round(., 2)))

summary(r_t.y)
