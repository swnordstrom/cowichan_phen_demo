########################################################
# Script with preliminary no-climate analysis of survival and growth
# Hopefully also an attempt to build growth kernels?
# SN - init 8 feb 2024
########################################################

library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyr)
library(glmmTMB)

rm(list = ls())


########################################################
# Read in data

# Demo data
demo = read.csv('01_data_cleaning/out/demo_imputed_survival.csv')


########################################################
# Process data as needed

# Demo dataset, subsetted for survival
# we'll say that survival in year $y$ means surviving from $y$ to $y+1$
demo.for.surv = merge(
  # x - status in year $y$ - only interested in plants that survived in y
  x = demo %>% 
    filter(surv) %>% 
    select(Year, No.leaves, Leaf.length, No.umbels, demo.note, proc.note, plantid, Plot, trt),
  # y - whether plant was observed alive in year $y+1$
  y = demo %>%
    mutate(p.year = Year - 1) %>%
    select(p.year, plantid, surv, No.leaves, Leaf.length),
  by.x = c('Year', 'plantid'), by.y = c('p.year', 'plantid'),
  suffixes = c('.t', '.tp1')
)

head(demo.for.surv)
nrow(demo.for.surv)
length(unique(demo.for.surv$plantid))
table(demo.for.surv$Year)

# Further subset demo with *only plants with usable sizes* prior to the transition
# (I imagine this is what will be used in the final analysis)
# This also means cutting out the 2016 plants because their size measures are not good
demo.for.surv.sizes = demo.for.surv %>%
  filter(!is.na(No.leaves.t) & !is.na(Leaf.length.t) & No.leaves.t > 0, Year > 2016) %>%
  # add a size column
  mutate(size.t = log(No.leaves.t * Leaf.length.t))

nrow(demo.for.surv.sizes)
length(unique(demo.for.surv.sizes$plantid)) # still have 940 plants

# Further subsetting the above dataset to get growth (conditioned on survival)
# to get this, we need plants that survived *and* have measurements
demo.for.growth = demo.for.surv.sizes %>%
  filter(!is.na(No.leaves.tp1) & !is.na(Leaf.length.tp1) & No.leaves.tp1 > 0) %>%
  # add a size column
  mutate(size.tp1 = log(No.leaves.tp1 * Leaf.length.tp1))


nrow(demo.for.growth)
length(unique(demo.for.growth$plantid))


########################################################
# Visualize raw data

# First, survival data
demo.for.surv.sizes %>%
  mutate(surv = as.numeric(surv) - 0.125 * as.numeric(trt %in% 'drought') + 0.125 * as.numeric(trt %in% 'irrigated')) %>%
  ggplot(aes(x = size.t, y = surv, colour = trt)) +
  geom_point(size = 3, alpha = 0.125, position = position_jitter(height = (1/24))) +
  scale_colour_manual(values = c('black', 'red', 'blue')) +
  scale_y_continuous(breaks = (0:4)/4, labels = (0:4)/4) +
  facet_wrap(~ Year) +
  theme(
    panel.background = element_blank(),
    legend.position = 'none'
  )

# Growth 
demo.for.growth %>%
  ggplot(aes(x = size.t, y = size.tp1, colour = trt)) +
  geom_segment(aes(x = 1, xend = 5, y = 1, yend = 5), linetype = 2) +
  geom_point(size = 3, alpha = 0.125) +
  scale_colour_manual(values = c('black', 'red', 'blue')) +
  facet_wrap(~ Year) +
  theme(
    panel.background = element_blank(),
    legend.position = 'none'
  )
# Interesting. Def some year-to-year variation!

demo.for.growth %>%
  ggplot(aes(x = size.t, y = size.tp1 / size.t, colour = trt)) +
  geom_segment(aes(x = 1, xend = 5, y = 1, yend = 1), linetype = 2) +
  geom_point(size = 3, alpha = 0.125) +
  scale_colour_manual(values = c('black', 'red', 'blue')) +
  facet_wrap(~ Year) +
  theme(
    panel.background = element_blank(),
    legend.position = 'none'
  )
# Interesting non-linearity here... what does it mean...?

########################################################
# Fit some models

# First, probability of survival #######################

# First thing's first,
s_0 = glmmTMB(
  surv ~ (1 | Plot / plantid) + (1 | Year),
  data = demo.for.surv.sizes,
  family = 'binomial'
)

s_s = glmmTMB(
  surv ~ size.t + (1 | Plot / plantid) + (1 | Year),
  data = demo.for.surv.sizes,
  family = 'binomial'
)

AIC(s_0, s_s)
# Unsurprisingly, effect of size

s_sy = glmmTMB(
  surv ~ (1 | Plot / plantid) + (size.t | Year),
  data = demo.for.surv.sizes,
  family = 'binomial'
)

AIC(s_sy, s_s)
# Interesting. No evidence of year effects!

# Polynomial?
s_s2 = glmmTMB(
  surv ~ poly(size.t, 2) + (1 | Plot / plantid) + (1 | Year),
  data = demo.for.surv.sizes,
  family = 'binomial'
)

AIC(s_s2, s_s)
# no polynomial

# Look for treatment effects

s_s_t = glmmTMB(
  surv ~ size.t + trt + (1 | Plot / plantid) + (1 | Year),
  data = demo.for.surv.sizes,
  family = 'binomial'
)

s_st = glmmTMB(
  surv ~ size.t * trt + (1 | Plot / plantid) + (1 | Year),
  data = demo.for.surv.sizes,
  family = 'binomial'
)

AIC(s_s_t, s_st, s_s)
# no treatment effects, once again...

summary(s_s)
# one unit (lol) increase in size means survival odds go up by exp(0.647) =
# 1.909 or ~91%!


# Next, growth kernel (conditioned on survival) #######

g_0 = glmmTMB(
  size.tp1 ~ (1 | Plot / plantid) + (1 | Year),
  data = demo.for.growth
)

g_s = glmmTMB(
  size.tp1 ~ size.t + (1 | Plot / plantid) + (1 | Year),
  data = demo.for.growth
)

AIC(g_s, g_0)
# very unsurprising, there is an effect of current size on next year's size...

# yearly-varying effect of previous size?
g_sy = glmmTMB(
  size.tp1 ~ (1 | Plot / plantid) + (size.t | Year),
  data = demo.for.growth
)

AIC(g_sy, g_s)
# yearly-varying effects... not in this model!

# 
g_s2 = glmmTMB(
  size.tp1 ~ poly(size.t, 2) + (1 | Plot / plantid) + (1 | Year),
  data = demo.for.growth
)

AIC(g_s2, g_s)
# okay... a polynomial effect
summary(g_s2) # but it does look kinda weak rel. to uncertainty

g_s_t = glmmTMB(
  size.tp1 ~ size.t + trt + (1 | Plot / plantid) + (1 | Year),
  data = demo.for.growth
)

g_st = glmmTMB(
  size.tp1 ~ size.t * trt + (1 | Plot / plantid) + (1 | Year),
  data = demo.for.growth
)

g_s2_t = glmmTMB(
  size.tp1 ~ poly(size.t, 2) + trt + (1 | Plot / plantid) + (1 | Year),
  data = demo.for.growth
)

g_s2t = glmmTMB(
  size.tp1 ~ poly(size.t, 2) * trt + (1 | Plot / plantid) + (1 | Year),
  data = demo.for.growth
)

# g_s2y_t = glmmTMB(
#   size.tp1 ~ trt + (1 | Plot / plantid) + (poly(size.t, 2) | Year),
#   data = demo.for.growth
# ) # convergence issue, what a bummer

# g_s2y = glmmTMB(
#   size.tp1 ~ (1 | Plot / plantid) + (poly(size.t, 2) | Year),
#   data = demo.for.growth
# )

AIC(g_s, g_s2, g_s_t, g_st, g_s2_t, g_s2t) %>%
  mutate(daic = round(AIC - min(AIC), 2)) %>%
  arrange(daic)

# ah very cool... there is a treatment effect!

summary(g_s2t)

# Crude plot of this...
expand.grid(
  size.t = (23:43)/10,
  trt = c('control', 'drought', 'irrigated')
) %>%
  mutate(
    pred.size.tp1 = predict(
      object = g_s2t,
      newdata = .,
      re.form = ~ 0,
      allow.new.levels = TRUE
    )
  ) %>%
  ggplot(aes(x = size.t, y = pred.size.tp1)) +
  geom_segment(aes(x = 2.3, xend = 4.3, y = 2.3, yend = 4.3), linetype = 2) +
  geom_line(aes(group = trt, colour = trt)) +
  scale_colour_manual(values = c('black', 'red', 'blue'))
# LMAO more growth in drought than in control...
# still kinda surprised by this

# Random effects for the growth model
ranef(g_s2t)$cond$Year %>% unlist() %>% plot()
# interesting... does look like there's an increase over time
ranef(g_s2t)$cond$Plot %>% unlist() %>% hist()
# plausibly normally distributed
ranef(g_s2t)$cond$`plantid:Plot` %>% unlist() %>% hist()
# hey cool that looks super normally distributed! hell yeah bro
ranef(g_s2t)$cond$`plantid:Plot` %>% unlist() %>% qqnorm()
ranef(g_s2t)$cond$`plantid:Plot` %>% unlist() %>% qqline()
# hmm upper tail of these is off but otherwise looks good
residuals(g_s2t) %>% qqnorm()
residuals(g_s2t) %>% qqline()
# right tail also looks off but otherwise good

demo.for.growth %>%
  mutate(resid = residuals(g_s2t)) %>%
  ggplot(aes(x = size.t, y = resid)) +
  geom_point(aes(colour = trt)) +
  geom_segment(aes(x = 1, xend = 5, y = 0, yend = 0), linetype = 2) +
  scale_colour_manual(values = c('black', 'red', 'blue')) +
  facet_wrap(~ Year)
# Honestly this looks fine to me. Looks like residual variance is similar across years.
# Also doesn't look like the variance changes with size.


########################################################
# Plot model predictions

##### Plot survival

surv.preds = data.frame(size.t = (10:50)/10) %>%
  mutate(
    pred.surv = predict(
      object = s_s,
      newdata = .,
      allow.new.levels = TRUE,
      re.form = ~ 0,
      type = 'response'
    )
  )

head(surv.preds)

surv.data.base = demo.for.surv.sizes %>%
  mutate(surv = as.numeric(surv) - 0.125 * as.numeric(trt %in% 'drought') + 0.125 * as.numeric(trt %in% 'irrigated')) %>%
  ggplot(aes(x = size.t, y = surv, colour = trt)) +
  geom_point(size = 3, alpha = 0.125, position = position_jitter(height = (1/24))) +
  labs(x = 'Probability of survival', y = 'Size after transition') +
  scale_colour_manual(values = c('black', 'red', 'blue')) +
  scale_y_continuous(breaks = (0:4)/4, labels = (0:4)/4) +
  facet_wrap(~ Year) +
  theme(
    panel.background = element_blank(),
    legend.position = 'none'
  )

surv.data.base +
  geom_line(
    data = surv.preds,
    inherit.aes = FALSE,
    aes(x = size.t, y = pred.surv)
  )
# Interesting.


##### Plot growth

growth.preds = expand.grid(
  size.t = (10:50)/10,
  trt = c('control', 'drought', 'irrigated')
) %>%
  mutate(
    pred.size.tp1 = predict(
      object = g_s2t,
      newdata = .,
      re.form = ~ 0,
      allow.new.levels = TRUE
    )
  ) 

growth.data.base = demo.for.growth %>%
  ggplot(aes(x = size.t, y = size.tp1, colour = trt)) +
  # geom_segment(aes(x = 1, xend = 5, y = 1, yend = 5), linetype = 2) +
  geom_point(size = 3, alpha = 0.125) +
  labs(x = 'Size before transition', y = 'Size after transition') +
  scale_colour_manual(values = c('black', 'red', 'blue')) +
  facet_wrap(~ Year) +
  theme(
    panel.background = element_blank(),
    legend.position = 'none'
  )

growth.data.base +
  geom_line(
    data = growth.preds,
    inherit.aes = FALSE,
    aes(x = size.t, y = pred.size.tp1, group = trt, colour = trt)
  )
