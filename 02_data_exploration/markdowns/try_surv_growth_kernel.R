#' ---
#' title: 'First stab at a survival+growth kernel'
#' date: '9 Feb 2024'
#' output:
#'   pdf_document:
#'     keep_tex: true
#' ---
#' 
#' Here is a first attempt at estimating a survival + growth kernel for each treatment. 

#+ setup, echo = FALSE, warning = FALSE, message = FALSE

# Load packages
library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyr)
library(glmmTMB)

# Clear namespace
rm(list = ls())

# Set wd for markdown
# knitr::opts_knit$set(root.dir = rprojroot::find_rstudio_root_file()) # why is this not working...?
knitr::opts_knit$set(root.dir = '~/repos/cowichan_phen_demo/') # I hate hard-coding this in but whatever


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

# Further subset demo with *only plants with usable sizes* prior to the transition
# (I imagine this is what will be used in the final analysis)
# This also means cutting out the 2016 plants because their size measures are not good
demo.for.surv.sizes = demo.for.surv %>%
  filter(!is.na(No.leaves.t) & !is.na(Leaf.length.t) & No.leaves.t > 0, Year > 2016) %>%
  # add a size column
  mutate(size.t = log(No.leaves.t * Leaf.length.t))

# Further subsetting the above dataset to get growth (conditioned on survival)
# to get this, we need plants that survived *and* have measurements
demo.for.growth = demo.for.surv.sizes %>%
  filter(!is.na(No.leaves.tp1) & !is.na(Leaf.length.tp1) & No.leaves.tp1 > 0) %>%
  # add a size column
  mutate(size.tp1 = log(No.leaves.tp1 * Leaf.length.tp1))


#' ### Survival models
#' 
#' I fit a bunch of survival models looking for fixed effects of size and treament.
#' 
#' Here is an AIC table comparing the models:
#' 
#+ fit-surv-models, echo = FALSE, message = FALSE

########################################################
# Fit some models

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

s_sy = glmmTMB(
  surv ~ (1 | Plot / plantid) + (size.t | Year),
  data = demo.for.surv.sizes,
  family = 'binomial'
)

s_s2 = glmmTMB(
  surv ~ poly(size.t, 2) + (1 | Plot / plantid) + (1 | Year),
  data = demo.for.surv.sizes,
  family = 'binomial'
)

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

AIC(s_0, s_s, s_s_t, s_st, s_sy, s_s2) %>%
  select(AIC) %>%
  mutate(delta.AIC = round(AIC - min(AIC), 2)) %>%
  mutate(form = sapply(X = list(s_0, s_s, s_s_t, s_st, s_sy, s_s2), FUN = formula)) %>%
  arrange(delta.AIC)

#' The best performing model to predict survival has only a fixed effect of size
#' (before the transition). Treatment does not have a discernible effect. There
#' is no evidence of a quadratic or year-varying effect of size on survival
#' odds.
#' 
#' Here is the model summary for the best-performing model:

#+ surv-mod-summary, echo = FALSE

summary(s_s)

#' An increase in size of one unit corresponds to an increase of exp(0.647) =
#' 1.90, or 90%, in the odds of survival. The standard deviation for sizes
#' across the whole dataset is approx 0.8, in which case an increase of one
#' standard deviation in size corresponds to a 67% increase in the odds of
#' survival.
#' 
#' ### Growth models
#' 
#' I fit several models to estimate growth (conditioned on survival).
#' Here is an AIC table comparing the models:

#+ fit-growth-mods, echo = FALSE

g_0 = glmmTMB(
  size.tp1 ~ (1 | Plot / plantid) + (1 | Year),
  data = demo.for.growth
)

g_s = glmmTMB(
  size.tp1 ~ size.t + (1 | Plot / plantid) + (1 | Year),
  data = demo.for.growth
)

g_sy = glmmTMB(
  size.tp1 ~ (1 | Plot / plantid) + (size.t | Year),
  data = demo.for.growth
)

g_s2 = glmmTMB(
  size.tp1 ~ poly(size.t, 2) + (1 | Plot / plantid) + (1 | Year),
  data = demo.for.growth
)

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


AIC(g_0, g_s, g_sy, g_s2, g_s_t, g_st, g_s2_t, g_s2t) %>%
  select(AIC) %>%
  mutate(delta.aic = round(AIC - min(AIC), 2)) %>%
  mutate(formula = sapply(list(g_0, g_s, g_sy, g_s2, g_s_t, g_st, g_s2_t, g_s2t), formula)) %>%
  arrange(delta.aic)

#' The best performing margin, by a considerable margin, contains a quadratic
#' effect of size (before transition) that varies by treatment. I am not sure
#' what the comma in the year-level random effect is doing - all of these models
#' have a simple year-varying intercept except for the model `g.sy` (which is has
#' a $\Delta$ AIC of ~25).
#' 
#' Here is the model output for the optimal model:

#+ summarise-growth, echo = FALSE

summary(g_s2t)

#' Here is a visualization of model predictions for each of the three treatments
#' (a subsequent figure will include raw data)

#+ vis-growth-prediction, echo = FALSE

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
  geom_segment(aes(x = 2.3, xend = 4.3, y = 2.3, yend = 4.3), colour = 'white', linetype = 2) +
  geom_line(aes(group = trt, colour = trt)) +
  scale_colour_manual(values = c('black', 'red', 'blue')) +
  labs(x = 'Size before transition', y = 'Size after transition')

#' The dashed line is the 1-1 line. This plot does suggest that large plants
#' will, on average, shrink (i.e., for large sizes pre-transition, the
#' expectation is below the 1-1 line). 
#' 
#' I'll plot the model output against raw data below. Before then, I'll plot the
#' residuals of the best-performing growth model. This is a test of
#' heteroskedasticity across size, years, or treatments.

#+ residual-check, echo = FALSE

demo.for.growth %>%
  mutate(resid = residuals(g_s2t)) %>%
  ggplot(aes(x = size.t, y = resid)) +
  geom_point(aes(colour = trt)) +
  geom_segment(aes(x = 1, xend = 5, y = 0, yend = 0), linetype = 2) +
  scale_colour_manual(values = c('black', 'red', 'blue')) +
  facet_wrap(~ Year)

#' There is some funny stuff happening in 2018 but otherwise this looks fine to
#' me.
#' 
#' ### Plotting model predictions
#' 
#' Now I'll plot the model predictions against data.
#' 
#' First, here is the survival data. The curves here are the same for each year
#' - the predicted probability for the "average" year. I fit these with a
#' different package and I'm having trouble getting the year-specific
#' predictions working.

#+ plot-survival-against-data, echo = FALSE

########################################################
# Plotting model output

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

#' Treatment was not significant in the model, but I dodged them on the y-axis
#' here just to confirm this visually. It looks correct to me - there doesn't
#' seem to be any difference among the treatments here.
#' 
#' Very small plants have a lower hance of survival - this may be a data
#' sparseness issue (we have few observations of small plants). Otherwise though
#' this looks okay to me.
#' 
#' Here is the growth model output:

#+ plot-growth-against-data, echo = FALSE

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

#' That quadratic fit for the drought treatment might be due to outliers from a
#' handful of plants. Also note again that the year-to-year variation is not
#' reflected in the predictions.
#' 
#' 
#' ### Putting together a growth/survival kernel
#' 
#' Here is a visualization of a growth-survival kernel:

#+ kernel, echo = FALSE

########################################################
# Trying to construct a growth kernel

growth.resid.sd = summary(g_s2t)$sigma

kernels = expand.grid(
  size.t = (45:595)/100,
  size.tp1 = (45:595)/100,
  trt = c('control', 'drought', 'irrigated')
)

kernels = kernels %>%
  mutate(
    p.survive = predict(
      object = s_s,
      newdata = .,
      allow.new.levels = TRUE,
      re.form = ~ 0,
      type = 'response'
    ),
    exp.size.tp1 = predict(
      object = g_s2t,
      newdata = .,
      allow.new.levels = TRUE,
      re.form = ~ 0,
      type = 'response'
    )
  )

kernels = kernels %>%
  mutate(
    p.grow.size.tp1 = dnorm(x = size.tp1, mean = exp.size.tp1, sd = growth.resid.sd),
    p.size.tp1 = p.grow.size.tp1 * p.survive
  )

ggkernel = ggplot(kernels, aes(x = size.t, y = size.tp1, fill = p.size.tp1))

ggkernel +
  geom_tile() +
  geom_segment(aes(x = .5, xend = 5.5, y = .5, yend = 5.5), colour = 'white', linetype = 2) +
  scale_y_reverse() +
  scale_fill_viridis_c('probability') +
  labs(x = 'Size in t', y = 'Size in t+1') +
  coord_equal() +
  facet_wrap(~ trt) +
  theme(
    legend.position = 'bottom',
    panel.background = element_blank()
  )

#' The white line is the 1-1 line (i.e., the matrix diagonal), so density above
#' the white line is shrinkage and density below the line is growth. That
#' quadratic effect in the drought treatment is quite strong. We'll see if it
#' appears in subsequent models.

# /* 
# knitr::spin('02_data_exploration/markdowns/try_surv_growth_kernel.R')
# rmarkdown::render('02_data_exploration/markdowns/try_surv_growth_kernel.R') # <- this one is better
# */