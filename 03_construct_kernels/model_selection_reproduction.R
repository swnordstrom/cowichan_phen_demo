#-------------------------------------------------------------------
# Model selection script for reproductive vital rate models
# Involves multiple steps:
# - Probability of flowering
# - Number of umbels produced
# - Probability of umbels surviving
# - Seeds per umbel
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

source('03_construct_kernels/prepare_demo_data_repr.R')

#-------------------------------------------------------------------
# Fit models for probability of flowering + umbels per plant
# Strategy here: use a hurdle model
# - Response is number of umbels per plant
# - Hurdle portion (through zero-inflation) is the probability of (not) flowering
# - Umbel count is modeled with a zero-inflated Poisson distribution

### Predictors:
# - size: size of plant in year of observation
# - plot: plot containing plant
# - plantid: unique identifier for each plant
# - year: year of observation
# - trt: drought, irrigated, control (categorical); control is reference


### Test first for size-effects on both parts of the model
# Then test for treatment effects
# Then test for potential interactions


u_0 = glmmTMB(
  No.umbels ~ (1 | Year) + (1 | Plot / plantid),
  family = 'truncated_poisson',
  ziformula = ~ (1 | Year) + (1 | Plot / plantid),
  data = demo.flow
)

u_s_0 = glmmTMB(
  No.umbels ~ size + (1 | Year) + (1 | Plot / plantid),
  family = 'truncated_poisson',
  ziformula = ~ (1 | Year) + (1 | Plot / plantid),
  data = demo.flow
)

u_0_s = glmmTMB(
  No.umbels ~ (1 | Year) + (1 | Plot / plantid),
  family = 'truncated_poisson',
  ziformula = ~ size + (1 | Year) + (1 | Plot / plantid),
  data = demo.flow
)

u_s_s = glmmTMB(
  No.umbels ~ size + (1 | Year) + (1 | Plot / plantid),
  family = 'truncated_poisson',
  ziformula = ~ size + (1 | Year) + (1 | Plot / plantid),
  data = demo.flow
)

AIC(u_0, u_0_s, u_s_0, u_s_s) %>% mutate(daic = round(AIC - min(AIC), 2))
# Keep size in the models

### Test for treatment effects on either rate

u_s.t_s = glmmTMB(
  No.umbels ~ trt + size + (1 | Year) + (1 | Plot / plantid),
  family = 'truncated_poisson',
  ziformula = ~ size + (1 | Year) + (1 | Plot / plantid),
  data = demo.flow
)

u_s_s.t = glmmTMB(
  No.umbels ~ size + (1 | Year) + (1 | Plot / plantid),
  family = 'truncated_poisson',
  ziformula = ~ trt + size + (1 | Year) + (1 | Plot / plantid),
  data = demo.flow
)

# Size treatment interactions
u_st_s = glmmTMB(
  No.umbels ~ size * trt + (1 | Year) + (1 | Plot / plantid),
  family = 'truncated_poisson',
  ziformula = ~ size + (1 | Year) + (1 | Plot / plantid),
  data = demo.flow
)

u_s_st = glmmTMB(
  No.umbels ~ size + (1 | Year) + (1 | Plot / plantid),
  family = 'truncated_poisson',
  ziformula = ~ trt * size + (1 | Year) + (1 | Plot / plantid),
  data = demo.flow
)

# Year-varying treatment effects
u_s.ty_s = glmmTMB(
  No.umbels ~ size + trt + (1 | Year) + (1 | Year:trt) + (1 | Plot / plantid),
  family = 'truncated_poisson',
  ziformula = ~ size + (1 | Year) + (1 | Plot / plantid),
  data = demo.flow
)

u_s_s.ty = glmmTMB(
  No.umbels ~ size + (1 | Year) + (1 | Plot / plantid),
  family = 'truncated_poisson',
  ziformula = ~ size + trt + (1 | Year) + (1 | Year:trt) + (1 | Plot / plantid),
  data = demo.flow
)

u_s.ty_s.ty = glmmTMB(
  No.umbels ~ size + trt + (1 | Year) + (1 | Year:trt) + (1 | Plot / plantid),
  family = 'truncated_poisson',
  ziformula = ~ size + trt + (1 | Year) + (1 | Year:trt) + (1 | Plot / plantid),
  data = demo.flow
)

AIC(u_s.t_s, u_s_s.t, u_st_s, u_s_st, u_s_s, u_s.ty_s, u_s_s.ty, u_s.ty_s.ty) %>%
  mutate(daic = round(AIC - min(AIC), 2))
# Best performing model has treatment-year effect for probability of flowering
# But how large is the effect size?

summary(u_s_s.ty)

# Testing for size-treatment effects with the year-varying treatment effects in
# the zero-inflation (prob. flowering) model
u_st_s.ty = glmmTMB(
  No.umbels ~ size * trt + (1 | Year) + (1 | Plot / plantid),
  family = 'truncated_poisson',
  ziformula = ~ size + trt + (1 | Year) + (1 | Year:trt) + (1 | Plot / plantid),
  data = demo.flow
)

u_s.t_st.ty = glmmTMB(
  No.umbels ~ size + trt + (1 | Year) + (1 | Plot / plantid),
  family = 'truncated_poisson',
  ziformula = ~ size * trt + (1 | Year) + (1 | Year:trt) + (1 | Plot / plantid),
  data = demo.flow
)

u_s_st.ty = glmmTMB(
  No.umbels ~ size + (1 | Year) + (1 | Plot / plantid),
  family = 'truncated_poisson',
  ziformula = ~ size * trt + (1 | Year) + (1 | Year:trt) + (1 | Plot / plantid),
  data = demo.flow
)

AIC(u_st_s.ty, u_s.t_st.ty, u_s_st.ty, u_s_s.ty) %>% mutate(daic = round(AIC - min(AIC), 2))
# No size-treatment effects

# Final AIC table
AIC(
  u_0, u_0_s, u_s_0, u_s_s, u_s_s.t, u_s.t_s, u_s_st, u_st_s, u_s_s.ty, 
  u_s.ty_s, u_s.ty_s.ty, u_s.t_st.ty, u_st_s.ty, u_s_st.ty
) %>%
  mutate(dAIC = round(AIC - min(AIC), 2))

# Summary of best performing model
summary(u_s_s.ty)

# Plot visuals:
expand.grid(size = (5:60)/10, trt = c('control', 'drought', 'irrigated')) %>%
  mutate(
    pred.umbl = predict(
      u_s_s.ty, newdata = ., allow.new.levels = TRUE, re.form = ~ 0, type = 'response'
    )
  ) %>%
  ggplot(aes(x = size, y = pred.umbl, colour = trt)) +
  geom_point(
    data = demo.flow,
    aes(x = size, y = No.umbels),
    size = 2, position = position_jitter(height = 0.25), alpha = 0.1
  ) +
  geom_line(linewidth = 1.2) + # scale_y_log10()
  labs(x = 'Size', y = 'Number of umbels') +
  scale_colour_manual(values = c('black', 'goldenrod', 'dodgerblue'))

# Looking at just the zero-inflation model (probability of flowering) - what
# does the zero-inflation model look like?

expand.grid(
  size = (5:60)/10, 
  trt = c('control', 'drought', 'irrigated'),
  Year = 2017:2024,
  Plot = 0, plantid = 0
) %>%
  mutate(
    p.zero.year = predict(u_s_s.ty, newdata = ., allow.new.levels = TRUE, re.form = NA, type = 'zprob')# ,
    # p.zero.mean = predict(u_s_s.ty, newdata = ., allow.new.levels = TRUE, re.form = NULL, type = 'zprob')
  ) %>%
  ggplot(aes(x = size, colour = trt)) +
  geom_line(aes(y = 1 - p.zero.year), linetype = 3) +
  # geom_line(aes(y = 1 - p.zero.mean), linewidth = 1.2) +
  labs(x = 'Size', y = 'Number of umbels') +
  scale_colour_manual(values = c('black', 'goldenrod', 'dodgerblue')) +
  facet_wrap(~ Year)

#-------------------------------------------------------------------
# # Here, we'll also use a zero-inflated model
# # - Probability of producing zero seed estimated with zero inflation term
# # - Distribution of seeds after the zero inflation is modeled by a Negative Binomial
# #   - doing this because it accurately captures overdispersion and because the 
# #     data-generating process can also account for umbels that were simply not 
# #     sufficiently pollinated (i.e., including zeros)

### Testing for size effects (expecting they will be sign. in both models)

s_0 = glmmTMB(
  no.seeds ~ Year + (1 | Plot / plantid),
  family = 'nbinom2',
  ziformula = ~ Year + (1 | Plot / plantid),
  data = seed
)

s_s_0 = glmmTMB(
  no.seeds ~ size + Year + (1 | Plot / plantid),
  family = 'nbinom2',
  ziformula = ~ Year + (1 | Plot / plantid),
  data = seed
)

s_0_s = glmmTMB(
  no.seeds ~ Year + (1 | Plot / plantid),
  family = 'nbinom2',
  ziformula = ~ size + Year + (1 | Plot / plantid),
  data = seed
)

s_s_s = glmmTMB(
  no.seeds ~ size + Year + (1 | Plot / plantid),
  family = 'nbinom2',
  ziformula = ~ size + Year + (1 | Plot / plantid),
  data = seed
)

AIC(s_0, s_s_0, s_0_s, s_s_s) %>% arrange(AIC) %>% mutate(daic = round(AIC - min(AIC), 2))
# As expected, size in both components (zero-inflation and seed count)

### Testing for treatmemnt effects (including interactions with size),
# on both/either model components

s_s.t_s = glmmTMB(
  no.seeds ~ trt + size + Year + (1 | Plot / plantid),
  family = 'nbinom2',
  ziformula = ~ size + Year + (1 | Plot / plantid),
  data = seed
)

s_s_s.t = glmmTMB(
  no.seeds ~ size + Year + (1 | Plot / plantid),
  family = 'nbinom2',
  ziformula = ~ trt + size + Year + (1 | Plot / plantid),
  data = seed
)

s_s.t_s.t = glmmTMB(
  no.seeds ~ trt + size + Year + (1 | Plot / plantid),
  family = 'nbinom2',
  ziformula = ~ trt + size + Year + (1 | Plot / plantid),
  data = seed
)

s_st_s = glmmTMB(
  no.seeds ~ trt * size + Year + (1 | Plot / plantid),
  family = 'nbinom2',
  ziformula = ~ size + Year + (1 | Plot / plantid),
  data = seed
)

s_s_st = glmmTMB(
  no.seeds ~ size + Year + (1 | Plot / plantid),
  family = 'nbinom2',
  ziformula = ~ trt * size + Year + (1 | Plot / plantid),
  data = seed
)

s_st_st = glmmTMB(
  no.seeds ~ trt * size + Year + (1 | Plot / plantid),
  family = 'nbinom2',
  ziformula = ~ trt * size + Year + (1 | Plot / plantid),
  data = seed
)

AIC(s_s_s, s_s.t_s, s_st_s, s_s_s.t, s_s_st, s_s.t_s.t, s_st_st) %>%
  mutate(daic = round(AIC - min(AIC), 2)) %>%
  arrange(daic)
# Best of these models has size-treatment interaction for seeds
# Check to make sure there's no linear treatment effect on the ZI term

s_st_s.t = glmmTMB(
  no.seeds ~ trt * size + Year + (1 | Plot / plantid),
  family = 'nbinom2',
  ziformula = ~ trt + size + Year + (1 | Plot / plantid),
  data = seed
)

AIC(s_st_s, s_st_s.t)
# Yes - no treatment effects on the zero-inflation term

### Test for treatment-year effects

s_st.ty_s = glmmTMB(
  no.seeds ~ trt * size + trt * Year + (1 | Plot / plantid),
  family = 'nbinom2',
  ziformula = ~ size + Year + (1 | Plot / plantid),
  data = seed
)

s_st_s.ty = glmmTMB(
  no.seeds ~ trt * size + Year + (1 | Plot / plantid),
  family = 'nbinom2',
  ziformula = ~ size + trt * Year + (1 | Plot / plantid),
  data = seed
)

s_st.ty_s.ty = glmmTMB(
  no.seeds ~ trt * size + trt * Year + (1 | Plot / plantid),
  family = 'nbinom2',
  ziformula = ~ size + trt * Year + (1 | Plot / plantid),
  data = seed
)

AIC(s_st_s, s_st_s.t, s_st.ty_s, s_st_s.ty, s_st.ty_s.ty) %>%
  mutate(daic = round(AIC - min(AIC), 2))
# so no treatment-year effects, very cool

### Look for effects of number of umbels produced

s_st.u_s = glmmTMB(
  no.seeds ~ trt * size + phen.umbels + Year + (1 | Plot / plantid),
  family = 'nbinom2',
  ziformula = ~ size + Year + (1 | Plot / plantid),
  data = seed
)

s_st_s.u = glmmTMB(
  no.seeds ~ trt * size + Year + (1 | Plot / plantid),
  family = 'nbinom2',
  ziformula = ~ size + phen.umbels + Year + (1 | Plot / plantid),
  data = seed
)

s_st.u_s.u = glmmTMB(
  no.seeds ~ trt * size + phen.umbels + Year + (1 | Plot / plantid),
  family = 'nbinom2',
  ziformula = ~ size + phen.umbels + Year + (1 | Plot / plantid),
  data = seed
)

AIC(s_st_s, s_st.u_s, s_st_s.u, s_st.u_s.u) %>%
  mutate(daic = round(AIC - min(AIC), 2))
# Umbel count in zero inflation term (this makes sense)
# However, support for the umbel count in the conditional term is very weak

### Tests for effects of phenology
# Start by looking at effects of phen on the count component (seed set)

s_st.p_s.u = glmmTMB(
  no.seeds ~ trt * size + Year + phen.c + (1 | Plot / plantid),
  family = 'nbinom2',
  ziformula = ~ size + phen.umbels + Year + (1 | Plot / plantid),
  data = seed
)

s_st.p2_s.u = glmmTMB(
  no.seeds ~ trt * size + Year + poly(phen.c, 2) + (1 | Plot / plantid),
  family = 'nbinom2',
  ziformula = ~ size + phen.umbels + Year + (1 | Plot / plantid),
  data = seed
)

s_st.tp_s.u = glmmTMB(
  no.seeds ~ trt * size + Year + trt * phen.c + (1 | Plot / plantid),
  family = 'nbinom2',
  ziformula = ~ size + phen.umbels + Year + (1 | Plot / plantid),
  data = seed
)

s_st.tp2_s.u = glmmTMB(
  no.seeds ~ trt * size + Year + trt * poly(phen.c, 2) + (1 | Plot / plantid),
  family = 'nbinom2',
  ziformula = ~ size + phen.umbels + Year + (1 | Plot / plantid),
  data = seed
)

s_st.py_s.u = glmmTMB(
  no.seeds ~ trt * size + Year * phen.c + (1 | Plot / plantid),
  family = 'nbinom2',
  ziformula = ~ size + phen.umbels + Year + (1 | Plot / plantid),
  data = seed
)

s_st.p2y_s.u = glmmTMB(
  no.seeds ~ trt * size + Year * poly(phen.c, 2) + (1 | Plot / plantid),
  family = 'nbinom2',
  ziformula = ~ size + phen.umbels + Year + (1 | Plot / plantid),
  data = seed
)

AIC(s_st_s.u, s_st.p_s.u, s_st.p2_s.u, s_st.tp_s.u, s_st.py_s.u, s_st.tp2_s.u, s_st.p2y_s.u) %>%
  mutate(daic = round(AIC - min(AIC), 2)) %>%
  arrange(daic)
# Evidence of linear phenology effect on seed count, not failure
# Very weak evidence of interactions with treatment or year, 
# parsimonious model assumes effect is linear

# Look for phenology effects on the zero-inflation term

s_st.p_s.u.p = glmmTMB(
  no.seeds ~ trt * size + Year + phen.c + (1 | Plot / plantid),
  family = 'nbinom2',
  ziformula = ~ size + phen.umbels + Year + phen.c + (1 | Plot / plantid),
  data = seed
)

s_st.p_s.u.tp = glmmTMB(
  no.seeds ~ trt * size + Year + phen.c + (1 | Plot / plantid),
  family = 'nbinom2',
  ziformula = ~ size + phen.umbels + Year + trt * phen.c + (1 | Plot / plantid),
  data = seed
)

s_st.p_s.u.py = glmmTMB(
  no.seeds ~ trt * size + Year + phen.c + (1 | Plot / plantid),
  family = 'nbinom2',
  ziformula = ~ size + phen.umbels + Year * phen.c + (1 | Plot / plantid),
  data = seed
)

s_st.p_s.u.p2 = glmmTMB(
  no.seeds ~ trt * size + Year + phen.c + (1 | Plot / plantid),
  family = 'nbinom2',
  ziformula = ~ size + phen.umbels + Year + poly(phen.c, 2) + (1 | Plot / plantid),
  data = seed
)

s_st.p_s.u.tp2 = glmmTMB(
  no.seeds ~ trt * size + Year + phen.c + (1 | Plot / plantid),
  family = 'nbinom2',
  ziformula = ~ size + phen.umbels + Year + trt * poly(phen.c, 2) + (1 | Plot / plantid),
  data = seed
)

s_st.p_s.u.p2y = glmmTMB(
  no.seeds ~ trt * size + Year + phen.c + (1 | Plot / plantid),
  family = 'nbinom2',
  ziformula = ~ size + phen.umbels + Year * poly(phen.c, 2) + (1 | Plot / plantid),
  data = seed
)

AIC(s_st.p_s.u, s_st.p_s.u.p, s_st.p_s.u.p2, s_st.p_s.u.tp, s_st.p_s.u.py, s_st.p_s.u.tp2, s_st.p_s.u.p2y) %>% 
  mutate(daic = round(AIC - min(AIC), 2)) %>%
  arrange(daic)
# There is definitely an effect
# Model with lowest AIC has polynomial term
# But, delta AIC against more parsimonious model with linear term is 1.27
# So, go with the linear model.

summary(s_st.p_s.u.p)

# Best parsimonious model is s_st.y_s.y

AIC(
  s_0, s_s_0, s_0_s, s_s_s,
  s_s.t_s, s_s_s.t, s_s.t_s.t, s_st_s, s_s_st, s_st_st, s_st_s.t,
  s_st.ty_s, s_st_s.ty, s_st.ty_s.ty,
  s_st.u_s, s_st_s.u, s_st.u_s.u,
  s_st.p_s.u, s_st.p2_s.u, s_st.tp_s.u, s_st.tp2_s.u, s_st.py_s.u, s_st.p2y_s.u,
  s_st.p_s.u.p, s_st.p_s.u.tp, s_st.p_s.u.py, s_st.p_s.u.p2, s_st.p_s.u.tp2, s_st.p_s.u.p2y
) %>%
  mutate(
    daic = round(AIC - min(AIC), 2),
    AIC = round(AIC, 2)
  )


#-------------------------------------------------------------------
# Getting some quantities reported in the manuscript.

### Estimate of phenology effects (one week acceleration of flowering)
expand.grid(phen.c = 0:-7, Year = 2021:2024) %>%
  mutate(size = 3.9, phen.umbels = 1, trt = 'control') %>%
  mutate(
    lin.zinf = predict(s_st.p_s.u.p, newdata = ., allow.new.levels = TRUE, re.form = ~ 0, type = 'zlink'),
    lin.cond = predict(s_st.p_s.u.p, newdata = ., allow.new.levels = TRUE, re.form = ~ 0, type = 'link')
  ) %>%
  group_by(phen.c, trt) %>%
  summarise(across(c(lin.zinf, lin.cond), mean)) %>%
  mutate(
    # Transform from linear scale to response scale
    seeds.per.umbel = (1 / (1 + exp(lin.zinf))) * exp(lin.cond),
  )


# Model comparisons: treatment-phenology effects on reproduction
anova(s_st.p_s.u.p, s_st.tp_s.u.p)
anova(s_st.p_s.u.p, s_st.p_s.u.tp)

# Model comparisons: quadratic terms for phenology
# Tests of these above were did not include fits of full model + quadratic terms
# Fitting those here for comparison:

s_st.p2_s.u.p = glmmTMB(
  no.seeds ~ trt * size + Year + poly(phen.c, 2) + (1 | Plot / plantid),
  family = 'nbinom2',
  ziformula = ~ size + phen.umbels + Year + phen.c + (1 | Plot / plantid),
  data = seed
)

s_st.p_s.u.p2 = glmmTMB(
  no.seeds ~ trt * size + Year + phen.c + (1 | Plot / plantid),
  family = 'nbinom2',
  ziformula = ~ size + phen.umbels + Year + poly(phen.c, 2) + (1 | Plot / plantid),
  data = seed
)

anova(s_st.p_s.u.p, s_st.p2_s.u.p)
anova(s_st.p_s.u.p, s_st.p_s.u.p2)
