library(ggplot2)
library(dplyr)
library(tidyr)
library(lme4)

rm(list = ls())

phen = merge(
  x = read.csv('01_data_cleaning/out/phenology_all_cleaned.csv'),
  y = read.csv('00_raw_data/plot_treatments.csv')
)

head(phen)

# Putting aside issues with tags between the multiple demo datasets...

#### Visualizations

phen %>%
  select(-c(init.wk, fina.wk)) %>%
  # I think there are multiple rows in here just due to umbel counts...
  # I don't think these rows actually reflect umbel tracking...
  distinct(plantid, year, .keep_all = TRUE) %>%
  pivot_longer(
    cols = contains('doy'),
    names_to = 'marker',
    values_to = 'period'
  ) %>%
  ggplot(aes(x = plot, y = period, colour = trt)) +
  geom_point(aes(shape = marker), position = position_jitter(width = 0.25, height = 2)) +
  scale_shape_manual(values = c(21, 19)) +
  scale_colour_manual(values = c('black', 'red', 'blue')) +
  facet_wrap(~ year)
# not super helpful...

phen %>%
  # I think there are multiple rows in here just due to umbel counts...
  # I don't think these rows actually reflect umbel tracking...
  distinct(plantid, year, .keep_all = TRUE) %>%
  ggplot(aes(x = plot, y = init.doy, colour = trt)) +
  geom_point(alpha = 0.5, position = position_jitter(width = 0.25, height = 2)) +
  scale_colour_manual(values = c('black', 'red', 'blue')) +
  facet_wrap(~ year)
# also doesn't look super great...

phen.mod = phen %>% 
  mutate(across(contains('doy'), ~ . - round(mean(.x, na.rm = TRUE))))

### Fit some models

# Null model
d_0 = lmer(
  init.doy ~ (1 | plot) + (1 | year),
  data = phen.mod
)

c(2.5569^2, 4.5077^2, 9.4504^2) %>% (function(x) x/sum(x))
# okay, so absent any other effects estimated, only ~24% of variance is plot and year-level
# (seems... not super convincing...)
# (although - our estimates of year-to-year variance is going to be super uncertain)

summary(d_0)

# Treatment effect
d_t = lmer(
  init.doy ~ trt + (1 | plot) + (1 | year),
  data = phen.mod
)

summary(d_t)
# ate away at some of the (already small) plot-to-plot variance
# does appear that drought means earlier flowering though!

# Try a likelihood ratio test (chi-squared test statistic)
anova(d_t, d_0)

# Okay! Does indeed look like there are effects of treatment!

# Could now see if treatment effects vary by year (which they surely do...)

d_yt = lmer(
  init.doy ~ (1 | plot) + (trt | year),
  data = phen.mod
)
# argh... singularity
# probably just not enough years to do this well.


# Look at residuals of the treatment model.
# (are residuals plot/year centered?)
predict(d_t) + residuals(d_t)
predict(d_t, re.form = NA) + residuals(d_t) 
predict(d_t, re.form = NULL) + residuals(d_t) 
# default is NULL, so residuals are plot/year centered

# Group = plot
phen.mod %>%
  mutate(resid = residuals(d_t)) %>%
  ggplot(aes(x = resid, group = plot, fill = trt)) +
  geom_histogram(position = position_identity(), alpha = 0.5, binwidth = 4) +
  scale_fill_manual(values = c('black', 'red', 'blue')) +
  facet_wrap(~ year, ncol = 1)
# oh my... lots of variance in 2023

# Group = trt
phen.mod %>%
  mutate(resid = residuals(d_t)) %>%
  ggplot(aes(x = resid, group = trt, fill = trt)) +
  geom_histogram(position = position_identity(), alpha = 0.5, binwidth = 4) +
  scale_fill_manual(values = c('black', 'red', 'blue')) +
  facet_wrap(~ year, ncol = 1)

# So inference here... not ideal.
# (this may explain why the treatment-by-year model didn't behave)

phen.mod %>%
  mutate(resid = residuals(d_t)) %>%
  ggplot(aes(x = plot, y = resid, colour = trt)) +
  geom_segment(aes(x = 0.5, xend = 15.5, y = 0, yend = 0), linetype = 2, colour = 'gray55') +
  geom_point(position = position_jitter(width = 0.25)) +
  scale_colour_manual(values = c('black', 'red', 'blue')) +
  facet_wrap(~ year, ncol = 1)
# doesn't look so bad, although certainly more variance in 2022

# Look at effects
ranef(d_t)$plot %>%
  mutate(plot = row.names(.)) %>%
  merge(distinct(phen.mod, plot, trt)) %>%
  ggplot(aes(x = plot, y = `(Intercept)`, colour = trt)) +
  geom_point(size = 3) +
  scale_colour_manual(values = c('black', 'red', 'blue'))
# Not very much variation in the control plots, but considerable variation in
# treatment plots
# Also nothing super consistent within treatment groups.

summary(d_t)

coef(d_t)
coefficients(d_t)

summary(d_t)$coefficients %>%
  as.data.frame() %>%
  mutate(param = row.names(.)) %>%
  ggplot(aes(y = param, colour = param)) +
  geom_point(aes(x = Estimate), size = 4) +
  geom_segment(aes(yend = param, x = Estimate - `Std. Error`, xend = Estimate + `Std. Error`)) +
  scale_colour_manual(values = c('black', 'red', 'blue'))

# but maybe I could look just at start dates for individual umbels?

###################################
###################################
###################################

rm(list = ls())

phen = merge(
  x = read.csv('01_data_cleaning/out/phenology_all_ind_cleaned.csv'),
  y = read.csv('00_raw_data/plot_treatments.csv')
)

head(phen)

phen.mod = phen %>%
  mutate(tagplot = gsub('\\_\\d{1,2}[A-Z]', '', plantid)) %>%
  mutate(init.doy = init.doy - round(mean(init.doy)))

head(phen.mod)
nrow(phen.mod)

d_0 = lmer(
  formula = init.doy ~ (1 | plot / tagplot) + (1 | year),
  data = phen.mod
)

summary(d_0)

d_t = lmer(
  formula = init.doy ~ trt + (1 | plot / tagplot) + (1 | year),
  data = phen.mod
)

summary(d_t)


anova(d_t, d_0)
# still significant, which is cool...
# here though mean effect of irrigation is positive... although this is 

phen.preds = expand.grid(
  year = 2021:2023,
  trt  = c('control', 'drought', 'irrigated')
) %>%
  mutate(
    pred.0 = predict(d_t, newdata = ., re.form = ~ 0),
    pred.y = predict(d_t, newdata = ., re.form = ~ (1 | year))
  )

phen.preds %>%
  mutate(year = factor(year)) %>%
  ggplot(aes(x = trt, colour = trt)) +
  geom_point(aes(y = pred.0), size = 5, shape = 21) +
  geom_point(
    aes(y = pred.y, shape = year),
    size = 3,
    position = position_dodge(width = 0.5)
  ) +
  scale_colour_manual(values = c('black', 'red', 'blue'))

d_pi.0 = bootMer(
  d_t, nsim = 1000,
  FUN = function(x) predict(x, newdata = phen.preds, re.form = ~ 0),
  seed = 405
)
# extremely slow for some reason...

dim(d_ci.0$t)

d_ci.0 = bootMer(
  d_t, nsim = 1000, predict, re.form = NULL,
  seed = 22901
)
# but this has different dimension than the new data frame... ugh

phen.preds %>%
  cbind(. , t(apply(d_pi.0$t, 2, function(x) quantile(x, c(0.025, 0.975))))) %>%
  mutate(year = factor(year)) %>%
  ggplot(aes(x = trt, colour = trt)) +
  geom_segment(aes(xend = trt, y = `2.5%`, yend = `97.5%`)) +
  geom_point(aes(y = pred.0), size = 5) +
  scale_colour_manual(values = c('black', 'red', 'blue'))
# oh... lmao, damn wtf

# hmm...

# ### some housekeeping
# 
# # The demo.seed is upstream of the imputed survival data frame I've been using
# # It's possible that some plantids changed when doing the survival imputing
# # Let's do a quick comparison here
# 
# # Read in both data frames
# demo.impu = read.csv('01_data_cleaning/out/demo_imputed_survival.csv')
# demo.seed = read.csv('01_data_cleaning/out/demo_seed_v1.csv')
# 
# head(demo.seed)
# table(demo.seed$year)
# # the seed data frame has two different plantid columns
# demo.seed %>% group_by(is.same = plantid.demo == plantid.seed) %>% summarise(n = n())
# # in most cases  they're different...
# demo.seed %>% filter(plantid.demo != plantid.seed) %>% head()
# # looks like the differences are mostly in the coordinates (lol)
# 
# # Get a data frame to compare them
# demo.compare = merge(
#   # add an `in.demo` column
#   x = demo.impu %>% filter(Year > 2020) %>% select(Year, plantid) %>% mutate(in.demo = TRUE),
#   y = demo.seed %>% distinct(year, plantid.demo),
#   by.x = c('Year', 'plantid'), by.y = c('year', 'plantid.demo'),
#   # get all plants from the seed dataset
#   # if there's a mis-alignment, then we'll get NAs in the `in.demo` column
#   all.y = TRUE
# )
# 
# demo.compare %>% filter(is.na(in.demo)) %>% nrow()
# # 38 plants
# 
# demo.compare %>% filter(is.na(in.demo))
# # why are there NAs...
# 
