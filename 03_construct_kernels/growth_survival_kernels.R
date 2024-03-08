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
demo = read.csv('01_data_cleaning/out/demo_imputed_survival.csv')

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
  mutate(size.t = log(No.leaves.t * Leaf.length.t)) %>%
  # Change year to factor
  mutate(Year = factor(Year))

# Get a dataset of demo for growth - subsetting the data frame above
# to get this, we need plants that survived *and* have measurements
demo.for.growth = demo.for.surv.sizes %>%
  filter(!is.na(No.leaves.tp1) & !is.na(Leaf.length.tp1) & No.leaves.tp1 > 0) %>%
  # add a size column
  mutate(size.tp1 = log(No.leaves.tp1 * Leaf.length.tp1))

# ------------------------------------------------
# Fit survival models

# First question: do we need size
# (I suspect we will)

s_0 = glmmTMB(
  formula = surv ~ (1 | Year) + (1 | Plot / plantid),
  family = 'binomial',
  data = demo.for.surv.sizes
)

s_s = glmmTMB(
  formula = surv ~ size.t + (1 | Year) + (1 | Plot / plantid),
  family = 'binomial',
  data = demo.for.surv.sizes
)

AIC(s_0, s_s)
# Unsurprisingly, size effect present.

# Also check for random effects vs. fixed effects?

s_s_fix = glmmTMB(
  formula = surv ~ size.t + Year + (1 | Plot / plantid),
  family = 'binomial',
  data = demo.for.surv.sizes
)

AIC(s_s, s_s_fix)
# Fixed effects model performs slightly better
# BUT in looking for inter-annual variation it just means too many parameters
# So continue with fixed effects

# Models to test here:
# - treatment effect
# - size-treatment interaction
# - year-varying size effect
# - (model will not converge for treatment effect varying by year)

surv.mod.forms = c(
  'surv ~ size.t + (1 | Year)',
  'surv ~ size.t + trt + (1 | Year)',
  'surv ~ size.t * trt + (1 | Year)',
  'surv ~ (size.t | Year)',
  'surv ~ trt + (size.t | Year)'# ,
  # 'surv ~ size.t + (trt | Year)'
) %>%
  paste('+ (1 | Plot / plantid)')

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
# No convergence issues here

data.frame(
  AIC = sapply(surv.mod.list, AIC) %>% unlist(),
  form = surv.mod.forms
) %>%
  mutate(daic = round(AIC - min(AIC), 2))
# Very funny! AIC ordering is same as the order in which I entered models

# Okay, but now the fixed effects model is still the best-performing one.

# Inspect residuals
demo.for.surv.sizes %>%
  mutate(
    r.ran = residuals(s_s),
    r.fix = residuals(s_s_fix)
) %>%
  pivot_longer(cols = c(r.ran, r.fix), names_to = 'mod', values_to = 'resid') %>%
  ggplot(aes(x = size.t, y = resid)) +
  geom_point() +
  facet_wrap(Year ~ mod)
# Ah right... residuals are not normally distributed here
# But these look very similar otherwise...

# Comparison of magnitude of annual variation:
data.frame(
  fix.y.var = fixef(s_s_fix)$cond %>%
    (function(x) x[grepl('Year', names(x))]) %>%
    (function(x) c(0, x)) %>%
    (function(x) x - mean(x)),
  ran.y.var = ranef(s_s)$cond$Year %>% unlist()
) %>%
  ggplot(aes(x = fix.y.var, y = ran.y.var)) +
  geom_segment(aes(x = -0.3, xend = 0.3, y = -0.3, yend = 0.3), linetype = 2) +
  geom_point(size = 3)
# Fixed effects have a lot more variation in magnitude

# Going to go with... fixed effects here
# Will estimate linear predictors and add random noise to it

# Plot model outputs:
expand.grid(size.t = (5:60)/10, Year = factor(2017:2022)) %>%
  mutate(
    pred = predict(
      s_s_fix, newdata = ., type = 'response',
      re.form = ~ 0, allow.new.levels = TRUE
    )
  ) %>%
  ggplot(aes(x = size.t, y = pred)) +
  geom_point(
    data = demo.for.surv.sizes,
    aes(x = size.t, y = as.numeric(surv)),
    size = 3, alpha = 0.1
  ) +
  geom_line() +
  labs(x = 'size', y = 'survival') +
  facet_wrap(~ Year)
# Ah right... a lot of behavior here being driven by plants in the tail.
# Okay, well, that happens sometimes.

# And now one for the "average" year:
# it'll be easier if I make an object here with the year random effects

s_s_fix_year_terms = fixef(s_s_fix)$cond %>%
  (function(x) x[grepl('Year', names(x))]) %>%
  (function(x) c(0, x))

expand.grid(size.t = (5:60)/10, Year = factor(2017)) %>%
  mutate(
    lin.pred = predict(
      s_s_fix, newdata = ., type = 'link',
      re.form = ~ 0, allow.new.levels = TRUE
    ),
    lin.pred = lin.pred + mean(s_s_fix_year_terms),
    pred = 1 / (1 + exp(-lin.pred))
  ) %>%
  ggplot(aes(x = size.t, y = pred, group = Year)) +
  geom_point(
    data = demo.for.surv.sizes,
    aes(x = size.t, y = as.numeric(surv)),
    size = 3, alpha = 0.1
  ) +
  geom_line() +
  labs(x = 'size', y = 'survival')

# Okay.

# ------------------------------------------------
# Fit growth models

# Size effect of course will be here.
# First test is whether annual variation should be fixed or random:

g_ran = lmer(
  formula = size.tp1 ~ size.t + (1 | Year) + (1 | Plot / plantid),
  data = demo.for.growth
)

g_fix = lmer(
  formula = size.tp1 ~ size.t + Year + (1 | Plot / plantid),
  data = demo.for.growth
)

AIC(g_ran, g_fix)
# Slight advantage to the random effects model here.
# So go with the random effects. Good!

grow.mod.forms = c(
  'size.tp1 ~ size.t + (1 | Year)',
  'size.tp1 ~ size.t + trt + (1 | Year)',
  'size.tp1 ~ size.t * trt + (1 | Year)',
  'size.tp1 ~ (size.t | Year)',
  'size.tp1 ~ trt + (size.t | Year)'
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
# Quick!

data.frame(
  AIC = sapply(grow.mod.list, AIC) %>% unlist(),
  form = grow.mod.forms
) %>%
  mutate(daic = round(AIC - min(AIC), 2)) %>%
  arrange(daic)

# Once again, no treatment effect, no good evidence of year-varying effects of
# size

# Residuals:
demo.for.growth %>%
  mutate(resid = residuals(g_ran)) %>%
  ggplot(aes(x = size.t, y = resid)) +
  geom_point() +
  facet_wrap( ~ Year)
# Some negative outliers (negative tail looking bleh), and 2018 is weird, but
# otherwise okay
# No visual evidence of heteroskedasticity

# Model predictions, year-varying first
expand.grid(size.t = (9:60)/10, Year = factor(2017:2022)) %>%
  mutate(
    pred = predict(
      g_ran,
      newdata = .,
      re.form = ~ (1 | Year),
      type = 'response',
      allow.new.levels = TRUE
    )
  ) %>%
  ggplot(aes(x = size.t, y = pred)) +
  geom_point(
    data = demo.for.growth,
    aes(x = size.t, y = size.tp1),
    inherit.aes = FALSE,
    size = 3, alpha = 0.1
  ) +
  geom_segment(aes(x = .9, xend = 6, y = .9, yend = 6), linetype = 2, colour = 'gray77') +
  geom_line() +
  facet_wrap(~ Year)
# Alright, okey dokey, etc.

# Now plot grand mean:
expand.grid(size.t = (9:60)/10) %>%
  mutate(
    pred = predict(
      g_ran,
      newdata = .,
      re.form = ~ 0,
      type = 'response',
      allow.new.levels = TRUE
    )
  ) %>%
  ggplot(aes(x = size.t, y = pred)) +
  geom_point(
    data = demo.for.growth,
    aes(x = size.t, y = size.tp1),
    inherit.aes = FALSE,
    size = 3, alpha = 0.1
  ) +
  geom_segment(aes(x = .9, xend = 6, y = .9, yend = 6), linetype = 2, colour = 'gray77') +
  geom_line()
# Agh. It does look like it over-predicts for low sizes. Ugh.
# Oh well.


# ------------------------------------------------
# Construct kernels

### Some things I'll need 

# Distribution of year-fixed effects in year model
s_s_fix_year_terms = fixef(s_s_fix)$cond %>%
  (function(x) x[grepl('Year', names(x))]) %>%
  (function(x) c(0, x))

# Residual variance in growth models
grow.sd = summary(g_ran)$sigma

# Get a scaffold
kernel.scaffold = expand.grid(
  size.t = (5:60)/10,
  size.tp1 = (5:60)/10
)

nrow(kernel.scaffold)

kernel.surv = kernel.scaffold %>%
  mutate(Year = factor(2017)) %>%
  mutate(
    linpred = predict(
      newdata = .,
      object = s_s_fix, type = 'link',
      re.form = ~ 0, allow.new.levels = TRUE
    ),
    linpred = linpred + mean(s_s_fix_year_terms),
    pred.surv = 1 / (1+exp(-linpred))
  ) %>%
  select(-c(linpred, Year))

head(kernel.surv)

kernel.grow = kernel.scaffold %>%
  mutate(
    pred.mean = predict(
      newdata = .,
      object = g_ran, type = 'response',
      re.form = ~ 0, allow.new.levels = TRUE
    ),
    p.grow.tp1 = 0.1 * dnorm(size.tp1, mean = pred.mean, sd = grow.sd)
  )

head(kernel.grow)

kernel.sg = merge(x = kernel.surv, y = kernel.grow, by = c('size.t', 'size.tp1')) %>%
  mutate(p.size.tp1 = pred.surv * p.grow.tp1)

ggplot(kernel.sg, aes(x = size.t, y = size.tp1)) +
  geom_tile(aes(fill = p.size.tp1)) +
  scale_y_reverse()

# write.csv(
#   kernel.sg,
#   file = '03_construct_kernels/out/test_survgrow_kernel.csv',
#   row.names = FALSE
# )
