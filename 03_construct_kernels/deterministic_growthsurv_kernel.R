# Script for re-running reproductive vital rate models and performing model
# selection and exporting a deterministic, survival e kernel for
# IPM analysis.
# Reads in processed demo/seed (including phenology) data, 2016-2024
# (sn init july 2024)

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

### ----- MODEL SELECTION -----

# # ------------------------------------------------
# # Fit survival models
# 
# # First question: do we need size and/or year
# # (I suspect we will)
# 
# s_0 = glmmTMB(
#   formula = surv ~ (1 | Plot),
#   family = 'binomial',
#   data = demo.surv.sizes
# )
# 
# s_y = glmmTMB(
#   formula = surv ~ (1 | surv.year) + (1 | Plot),
#   family = 'binomial',
#   data = demo.surv.sizes
# )
# 
# s_s = glmmTMB(
#   formula = surv ~ size.prev + (1 | Plot),
#   family = 'binomial',
#   data = demo.surv.sizes
# )
# 
# s_s_y = glmmTMB(
#   formula = surv ~ size.prev + (1 | surv.year) + (1 | Plot),
#   family = 'binomial',
#   data = demo.surv.sizes
# )
# 
# AIC(s_0, s_y, s_s, s_s_y) %>% mutate(daic = round(AIC - min(AIC), 2))
# 
# # Interesting! No support for year effects.
# 
# summary(s_s)
# # Increasing size increases odds of survival, unsurprisingly!
# 
# # Models to test:
# # - Treatment effect
# # - Size-treatment interaciton
# 
# s_t = glmmTMB(
#   formula = surv ~ size.prev + trt + (1 | Plot),
#   family = 'binomial',
#   data = demo.surv.sizes
# )
# 
# s_ts = glmmTMB(
#   formula = surv ~ size.prev * trt + (1 | Plot),
#   family = 'binomial',
#   data = demo.surv.sizes
# )
# 
# AIC(s_s, s_t, s_ts) %>% mutate(daic = round(AIC - min(AIC), 2))
# # No treatment effects.
# 
# # Plot of this model:
# surv.preds = data.frame(size.prev = (5:60) / 10) %>%
#   mutate(pred = predict(s_s, re.form = ~ 0, newdata = ., allow.new.levels = TRUE, type = 'response'))
# 
# head(surv.preds)
# 
# demo.surv.sizes %>%
#   mutate(surv = as.numeric(surv)) %>%
#   ggplot(aes(x = size.prev, y = surv)) +
#   geom_point(position = position_jitter(height = 0.2), alpha = 0.1, size = 3) +
#   geom_line(data = surv.preds, aes(y = pred))
# 
# # # Random effects:
# 
# # Plot level effects:
# data.frame(ranef(s_s_y)$cond$Plot) %>%
#   mutate(Plot = row.names(.)) %>%
#   rename(incp = `X.Intercept.`) %>%
#   merge(distinct(demo.surv.sizes, Plot, trt)) %>%
#   mutate(Plot = factor(Plot)) %>%
#   ggplot(aes(x = Plot, y = incp, colour = trt)) +
#   geom_point(size = 4) +
#   scale_colour_manual(values = c('black', 'red', 'blue'))
# # Definitely not seeing a clear pattern here
# 
# # Maybe slightly higher intercepts for irrigated but I trust the models
# # suggesting no treatment effect
# 
# 
# # ------------------------------------------------
# # Fit growth models
# 
# # Test for year effects
# 
# g_0 = glmmTMB(
#   size.cur ~ size.prev + (1 | Plot / plantid),
#   data = demo.grow
# )
# 
# g_y = glmmTMB(
#   size.cur ~ size.prev + (1 | prev.year) + (1 | Plot / plantid),
#   data = demo.grow
# )
# 
# AIC(g_y, g_0)
# # Yep. Support for a year effect.
# 
# # Treatment effect
# 
# g_t = glmmTMB(
#   size.cur ~ size.prev + trt + (1 | prev.year) + (1 | Plot / plantid),
#   data = demo.grow
# )
# 
# AIC(g_t, g_y)
# # Wow... no treatment effect.
# 
# # There are a couple of other things to try:
# # - treatment-year
# # - treatment-size
# # - size-year
# 
# g_ty = glmmTMB(
#   size.cur ~ size.prev + trt + (1 | prev.year) + (1 | prev.year:trt) + (1 | Plot / plantid),
#   data = demo.grow
# )
# 
# g_st.ty = glmmTMB(
#   size.cur ~ size.prev + size.prev * trt + (1 | prev.year) + (1 | prev.year:trt) + (1 | Plot / plantid),
#   data = demo.grow
# )
# 
# g_st = glmmTMB(
#   size.cur ~ size.prev * trt + (1 | prev.year) + (1 | Plot / plantid),
#   data = demo.grow
# )
# 
# g_sy = glmmTMB(
#   size.cur ~ size.prev + (size.prev | prev.year) + (1 | Plot / plantid),
#   data = demo.grow
# )
# 
# AIC(g_y, g_t, g_ty, g_st.ty, g_st, g_sy) %>% mutate(daic = round(AIC - min(AIC), 2))
# # Yep - as before, overwhelming support for the year-treatment interaction
# 
# summary(g_st.ty)
# summary(g_ty)
# # These feel like *very* weak effects on average...
# # (We also have slight advantages to both drought and irrigation...)
# 
# # glmmTMB won't let me generate predictions with half-specified random effects
# data.frame(ranef(g_st.ty)$cond$`prev.year:trt`) %>%
#   rename(intcp = `X.Intercept.`) %>%
#   mutate(year.trt = row.names(.)) %>%
#   separate(year.trt, into = c('year', 'trt'), sep = ':') %>%
#   ggplot(aes(x = year, y = intcp, colour = trt, group = trt)) +
#   geom_line() +
#   scale_colour_manual(values = c('black', 'red', 'blue'))
# # Hmm... yep pretty noisy  
# # Although if the mean effect of drought is +.12... that's actually kinda big
# # compared to the magnitude of fluctuations?
# 
# # What would results look like if I made predictions for the median plant+plot?
# abs(ranef(g_st.ty)$cond$Plot) %>% arrange(`(Intercept)`) # median plot is 13, eff is +.0222
# ranef(g_st.ty)$cond$`plantid:Plot` %>% 
#   mutate(pla = row.names(.)) %>% 
#   filter(grepl('\\_13\\:', pla)) %>% 
#   mutate(eff = `(Intercept)` + .0222) %>% 
#   arrange(abs(eff))
# # seems like median plant is our old friend 3989_13:13
# 
# expand.grid(
#   size.prev = (5:60)/10,
#   prev.year = factor(2017:2023),
#   trt = c('control', 'drought', 'irrigated')
# ) %>%
#   mutate(Plot = 13, plantid = '3989_13') %>%
#   mutate(pp = predict(g_st.ty, newdata = ., re.form = NULL)) %>%
#   ggplot(aes(x = size.prev, y = pp, colour = trt, group = trt)) +
#   geom_segment(aes(x = 0.5, y = 0.5, xend = 6, yend = 6), linetype = 2, colour = 'gray77') +
#   geom_line() +
#   scale_colour_manual(values = c('black', 'red', 'blue')) +
#   facet_wrap(~ prev.year)
# # Usually, but not always, control performing below drought...
# # (It would be cool to iterate individuals through the years...)
# 
# # Compare predictions from our two different models
# expand.grid(
#   size.prev = (5:60)/10,
#   trt = c('control', 'drought', 'irrigated')
# ) %>%
#   mutate(
#     pred.yrf = predict(g_st.ty, newdata = ., re.form = ~ 0, allow.new.levels = TRUE),
#     pred.nrf = predict(g_st, newdata = ., re.form = ~ 0, allow.new.levels = TRUE)
#   ) %>%
#   ggplot(aes(x = pred.nrf, y = pred.yrf, colour = trt)) +
#   geom_segment(aes(x = 2, y = 2, xend = 4.5, yend = 4.5), linetype = 2, colour = 'gray77') +
#   geom_point(size = 3, alpha = 0.25) +
#   labs(x = 'non-year varying', y = 'year-varying') +
#   scale_colour_manual(values = c('black', 'red', 'blue'))
# # above diagonal: greater predicted growth in year-varying model
# # below diagonal: greater predicted growth in non-year varying model
# # there's the *least* change in the control treatment...
# 
# expand.grid(
#   size.prev = (5:60)/10,
#   trt = c('control', 'drought', 'irrigated')
# ) %>%
#   mutate(
#     pred.yrf = predict(g_st.ty, newdata = ., re.form = ~ 0, allow.new.levels = TRUE),
#     pred.nrf = predict(g_st, newdata = ., re.form = ~ 0, allow.new.levels = TRUE)
#   ) %>%
#   pivot_longer(cols = c(pred.yrf, pred.nrf), names_to = 'model', values_to = 'pred') %>%
#   ggplot(aes(x = size.prev, y = pred, group = interaction(trt, model), colour = trt)) +
#   geom_segment(aes(x = .5, y = .5, xend = 6, yend = 6), linetype = 2, colour = 'gray77') +
#   geom_line(aes(linetype = model)) +
#   labs(x = 'non-year varying', y = 'year-varying') +
#   scale_colour_manual(values = c('black', 'red', 'blue'))
# # The model with the year-varying effects has less shrinkage for large plants, less growth for small plants
# 
# # I just don't trust these random effects
# AIC(g_y, g_st, g_sy) %>% mutate(daic = round(AIC - min(AIC), 2))
# 
# # Out of these, g_st appears to be best supported
# summary(g_st.ty)
# 
# expand.grid(trt = c('control', 'drought', 'irrigated'), size.prev = (5:60)/10) %>%
#   mutate(size.cur = predict(g_st.ty, re.form = ~ 0, newdata = ., allow.new.levels = TRUE)) %>%
#   ggplot(aes(x = size.prev, y = size.cur, colour = trt, group = trt)) +
#   geom_segment(aes(x = .5, y = .5, xend = 6, yend = 6), linetype = 2, linewidth = 0.25) +
#   geom_point(data = demo.grow, size = 3, alpha = 0.1) +
#   geom_line() +
#   scale_colour_manual(values = c('black', 'red', 'blue')) +
#   facet_wrap(~ trt)
# # Hmm...
# 
# # Look at residuals
# demo.grow %>%
#   mutate(resids = residuals(g_st.ty)) %>%
#   ggplot(aes(x = size.prev, y = resids, colour = trt)) +
#   geom_segment(aes(x = .5, xend = 6, y = 0, yend = 0), colour = 'gray', linetype = 2, linewidth = 0.5) +
#   scale_colour_manual(values = c('black' ,'red', 'blue')) +
#   geom_point(alpha = 0.1, size = 3)
# # Eh... a few outliers for medium-large plants...
# 
# demo.grow %>%
#   mutate(resids = residuals(g_st)) %>%
#   ggplot(aes(x = prev.year, y = resids)) +
#   geom_segment(aes(x = 2017, xend = 2023, y = 0, yend = 0), linetype = 2, linewidth = 0.5) +
#   geom_point(position = position_jitter(width = 0.4), alpha = 0.05, size = 3)
# # A few negative outliers in 2022, 2023 but nothing awful.

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

# --- Extract parameters needed

# Residual variance in growth models
grow.sd = summary(g_st.ty)$sigma

# --- Construct and work with data frame

# Get a scaffold
grow.surv.kernel = expand.grid(
  size.prev = (5:60)/10,
  size.cur = (5:60)/10,
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
      object = g_st.ty, type = 'response',
      re.form = ~ 0, allow.new.levels = TRUE
    )
  ) %>%
  # Predicted distribution of sizes in next time step
  mutate(p.grow.size = 0.1 * dnorm(size.cur, mean = pred.grow.mean, sd = grow.sd)) %>%
  # Combine all together to get overall size distribution in next time step
  mutate(p.size.cur = pred.surv * p.grow.size)

ggplot(grow.surv.kernel, aes(x = size.prev, y = size.cur)) +
  geom_tile(aes(fill = p.size.cur)) +
  scale_y_reverse() +
  facet_wrap(~ trt)

# Neat

write.csv(
  grow.surv.kernel,
  file = '03_construct_kernels/out/deterministic_growsurv_kernel.csv',
  row.names = FALSE
)
