# Quick-and-dirty script for looking at leave-one-year-out on survival data
# To do:
# - add year-random effects to formula ist?
# - add new predictors to formula list
#   (remember to re-scale precipitation)

library(ggplot2)
library(dplyr)
library(tidyr)
library(lme4)

rm(list = ls())

# Read in demography data
demo = read.csv('01_data_cleaning/out/demo_imputed_survival.csv')

head(demo)
nrow(demo)

# Add a "flowering" column
# Check to make sure this works
demo %>%
  mutate(
    flowering = case_when(
      No.umbels > 0 ~ TRUE,
      !No.umbels ~ FALSE,
      is.na(No.umbels) ~ FALSE,
      .default = NA
  )
  ) %>%
  distinct(No.umbels, flowering, surv)
# Looks good.

# Add column to data frame
demo = demo %>%
  mutate(
    flowering = case_when(
      No.umbels > 0 ~ TRUE,
      !No.umbels ~ FALSE,
      is.na(No.umbels) ~ FALSE,
      .default = NA
    )
  )

### Get this year's survival and last year's demographic data
demo.surv = merge(
  # x is response
  # a plant surviving to year t; we'll use convention that it survived in year t-1
  # so need to subtract 1 from year
  x = demo %>%
    mutate(Year = Year - 1) %>%
    select(Year, plantid, Plot, trt, surv, demo.note, proc.note, edited),
  # y is previous year's data
  # give *only* the plants that have size measurements
  # NOTE: this will cut out some plants from analysis if we don't have their
  # size in the previous year
  y = demo %>%
    # Filter out plants only that were alive in previous year (and have measurements)
    filter(!is.na(No.leaves) & !is.na(Leaf.length) & No.leaves > 0) %>%
    # Filter out only plants known to be alive
    filter(surv) %>%
    # Add a previous size column
    mutate(size.prev = log(No.leaves * Leaf.length)) %>%
    select(Year, plantid, size.prev, flowering, demo.note, proc.note, edited),
  all = FALSE, by = c("Year", "plantid"), suffixes = c("", ".prev")
) %>%
  mutate(Plot = factor(Plot))

nrow(demo.surv)

head(demo.surv)

# # Remove 2016 plants # (actually I think this is handled later...)
# demo.surv = demo.surv %>% filter(Year > 2016)

length(unique(demo.surv$plantid))
table(demo.surv$surv)
with(demo.surv, table(Year, surv))

# Should only be one dead record per plant
demo.surv %>%
  group_by(plantid) %>%
  filter(sum(!surv) > 1)
# good

### Read in climate data
clim = read.csv('01_data_cleaning/out/climate_summary.csv') %>%
  # fix freeze date column (change to numeric)
  mutate(last.freeze = as.numeric(as.Date(last.freeze)) - as.numeric(as.Date('2019-12-31'))) %>%
  # center temperature and precip variables at mean
  group_by(pd.label) %>%
  mutate(across(contains('Temp'), ~ (. - mean(., na.rm = TRUE)))) %>%
  mutate(across(contains('prec.sum'), ~ (. - mean(., na.rm = TRUE)) / 100 )) %>%
  ungroup() %>%
  # center last freeze date at some neutral-ish date
  mutate(last.freeze = last.freeze - 90)


head(clim)
# oh... should think carefully about the merge...
# survival from 2016 to 2017 should depend on growing season 2016, winter 2016,
# early season 2017...
# I think the pgrow from the flowering script should work

clim.merge = rbind(
  clim,
  clim %>% filter(pd.label %in% 'grow') %>% mutate(Year = Year + 1, pd.label = 'pgrow')
)

clim.merge = clim.merge %>%
  # the window start dates aren't necessary
  # also remove the last.freeze column (for now) because it makes things difficult
  select(-c(period.start, last.freeze)) %>%
  pivot_longer(c(prec.sum, contains("Temp")), names_to = "stat", values_to = "statVal") %>%
  pivot_wider(names_from = c(pd.label, stat), values_from = statVal) %>%
  merge(y = clim %>% filter(!pd.label %in% 'pgrow') %>% distinct(Year, last.freeze)) %>%
  select(
    Year, last.freeze, grow_meanTemp, early_prec.sum, early_minnTemp,
    winter_prec.sum, winter_meanTemp, summer_meanTemp, pgrow_meanTemp
  )

head(clim.merge)

demo.surv.cl = merge(x = demo.surv, y = clim.merge, by = 'Year') %>%
  filter(Year > 2016)

head(demo.surv.cl)
nrow(demo.surv.cl)

### Start fitting models

s_0 = glmer(
  formula = surv ~ (1 | Year) + (1 | Plot),
  data = demo.surv.cl,
  family = 'binomial'
)

s_s = glmer(
  formula = surv ~ size.prev + (1 | Year) + (1 | Plot),
  data = demo.surv.cl,
  family = 'binomial'
)
  
s_s.f = glmer(
  formula = surv ~ size.prev + flowering + (1 | Year) + (1 | Plot),
  data = demo.surv.cl,
  family = 'binomial'
)


AIC(s_0, s_s, s_s.f)
# flowering... not super helpful!

s_sf = glmer(
  formula = surv ~ size.prev * flowering + (1 | Year) + (1 | Plot),
  data = demo.surv.cl,
  family = 'binomial'
)

AIC(s_0, s_s, s_s.f, s_sf)
# hmm... 

s_s.f.t = glmer(
  formula = surv ~ size.prev + flowering * trt + (1 | Year) + (1 | Plot),
  data = demo.surv.cl,
  family = 'binomial'
)

s_st.f = glmer(
  formula = surv ~ size.prev * trt + flowering + (1 | Year) + (1 | Plot),
  data = demo.surv.cl,
  family = 'binomial'
)

s_s.ft = glmer(
  formula = surv ~ size.prev + flowering * trt + (1 | Year) + (1 | Plot),
  data = demo.surv.cl,
  family = 'binomial'
)

s_st.ft = glmer(
  formula = surv ~ size.prev * trt + flowering * trt + (1 | Year) + (1 | Plot),
  data = demo.surv.cl,
  family = 'binomial'
)

AIC(s_0, s_s, s_s.f, s_sf, s_s.f.t, s_st.f, s_s.ft, s_st.ft) %>%
  mutate(daic = round(AIC - min(AIC)))
# No evidence of treatment effects mattering here!
# but which of the models with and without flowering to include...
# try both I suppose

##### Wrapper functions for model fitting

# Inverse logit wrapper
ilogit = function(x) (1+exp(-x))^(-1)

data.splitter = function(df) {
  
  yrs = unique(df$Year)
  lis = vector('list', length = length(yrs))
  
  for (y in 1:length(yrs)) {
    lis[[y]]$train = df %>% filter(!(Year %in% yrs[y]))
    lis[[y]]$test  = df %>% filter(Year %in% yrs[y])
  }
  
  return(lis)
  
}

demo.sv.split = data.splitter(demo.surv.cl)

str(demo.sv.split)
length(demo.sv.split)

fit.split.model = function(dfs, mod.formula) {
  
  try( 
    mod <- glmer(
      formula = mod.formula,
      family = 'binomial',
      data = dfs$train,
      control = glmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 1e6))
    )
  )
  
  # NOTE: here, we're doing the predictions including the plot-level random effects.
  if (exists('mod')) { 
    mod.pred  = predict(mod, newdata = dfs$test, re.form = ~ (1 | Plot)) %>% ilogit()
    mod.score = mean((dfs$test$surv - mod.pred))^2
    return(mod.score)
  } else { return(NA) }
}

sapply(
  demo.sv.split, 
  fit.split.model, 
  mod.formula = 'surv ~ size.prev + (1 | Year) + (1 | Plot)'
)
# cool

##### Try some leave-one-year-out validation...

forms = c(
  'surv ~ size.prev + (1 | Year) + (1 | Plot)',
  'surv ~ size.prev + last.freeze + (1 | Plot)',
  'surv ~ size.prev + early_prec.sum + (1 | Plot)',
  'surv ~ size.prev + early_minnTemp + (1 | Plot)',
  'surv ~ size.prev + grow_meanTemp + (1 | Plot)',
  'surv ~ size.prev + summer_meanTemp + (1 | Plot)',
  'surv ~ size.prev + winter_meanTemp + (1 | Plot)',
  'surv ~ size.prev + winter_prec.sum + (1 | Plot)',
  'surv ~ size.prev + pgrow_meanTemp + (1 | Plot)',
  'surv ~ size.prev + flowering + (1 | Year) + (1 | Plot)',
  'surv ~ size.prev + flowering + last.freeze + (1 | Plot)',
  'surv ~ size.prev + flowering + early_prec.sum + (1 | Plot)',
  'surv ~ size.prev + flowering + early_minnTemp + (1 | Plot)',
  'surv ~ size.prev + flowering + grow_meanTemp + (1 | Plot)',
  'surv ~ size.prev + flowering + summer_meanTemp + (1 | Plot)',
  'surv ~ size.prev + flowering + winter_meanTemp + (1 | Plot)',
  'surv ~ size.prev + flowering + winter_prec.sum + (1 | Plot)',
  'surv ~ size.prev + flowering + pgrow_meanTemp + (1 | Plot)',
  'surv ~ size.prev + last.freeze * trt + (1 | Plot)',
  'surv ~ size.prev + early_prec.sum * trt + (1 | Plot)',
  'surv ~ size.prev + early_minnTemp * trt + (1 | Plot)',
  'surv ~ size.prev + grow_meanTemp * trt + (1 | Plot)',
  'surv ~ size.prev + summer_meanTemp * trt + (1 | Plot)',
  'surv ~ size.prev + winter_meanTemp * trt + (1 | Plot)',
  'surv ~ size.prev + winter_prec.sum * trt + (1 | Plot)',
  'surv ~ size.prev + pgrow_meanTemp * trt + (1 | Plot)',
  'surv ~ size.prev + flowering + last.freeze * trt + (1 | Plot)',
  'surv ~ size.prev + flowering + early_prec.sum * trt + (1 | Plot)',
  'surv ~ size.prev + flowering + early_minnTemp * trt + (1 | Plot)',
  'surv ~ size.prev + flowering + grow_meanTemp * trt + (1 | Plot)',
  'surv ~ size.prev + flowering + summer_meanTemp * trt + (1 | Plot)',
  'surv ~ size.prev + flowering + winter_meanTemp * trt + (1 | Plot)',
  'surv ~ size.prev + flowering + winter_prec.sum * trt + (1 | Plot)',
  'surv ~ size.prev + flowering + pgrow_meanTemp * trt + (1 | Plot)',
  'surv ~ size.prev * last.freeze + (1 | Plot)',
  'surv ~ size.prev * early_prec.sum + (1 | Plot)',
  'surv ~ size.prev * early_minnTemp + (1 | Plot)',
  'surv ~ size.prev * grow_meanTemp + (1 | Plot)',
  'surv ~ size.prev * summer_meanTemp + (1 | Plot)',
  'surv ~ size.prev * winter_meanTemp + (1 | Plot)',
  'surv ~ size.prev * winter_prec.sum + (1 | Plot)',
  'surv ~ size.prev * pgrow_meanTemp + (1 | Plot)'
)

# form.scores = vector(mode = 'double', length = length(forms))
form.scores = matrix(ncol = length(demo.sv.split), nrow = length(forms))

for (form in 1:length(forms)) {
  form.scores[form,] = sapply(
    X = demo.sv.split,
    FUN = fit.split.model,
    mod.formula = forms[form]
  )
  print(form)
}
# got singularity issues on all models (except for [1]) when year random effect was included
# I'm not sure why - this worked for the probability of flowering model
# convergence warning on final formula...

form.scores
# but... still got scores...

data.frame(form.scores) %>%
  mutate(form = 1:nrow(.)) %>%
  pivot_longer(-form, names_to = 'holdout', values_to = 'score') %>%
  ggplot(aes(x = form, y = score, group = holdout)) +
  geom_line(aes(colour = holdout))

apply(form.scores, 1, mean) %>% (function(x) x - min(x))
# these are pretty similar
# although there are two years that get biffed a few times...

forms[order(rowSums(form.scores))]
# interesting...

demo.surv.cl %>%
  ggplot(aes(x = Year, y = size.prev)) +
  geom_point(
    aes(colour = surv, shape = flowering, group = interaction(flowering, surv)),
    position = position_jitterdodge(jitter.width = .1, dodge.width = .5)
  ) +
  facet_wrap(~ trt, ncol = 1)
# hmm...

demo.surv.cl %>%
  ggplot(aes(x = size.prev, y = winter_meanTemp)) +
  geom_point(
    aes(colour = surv), 
    position = position_jitter(height = 0.02), 
    alpha = 0.75, size = 2
  ) +
  facet_wrap(~ trt)
# doesn't look all that strong...

demo.surv.cl %>%
  ggplot(aes(x = size.prev, y = winter_meanTemp)) +
  geom_point(
    aes(colour = surv), 
    position = position_jitter(height = 0.01), 
    alpha = 0.5, size = 3
  )

demo.surv.cl %>%
  ggplot(aes(x = size.prev, y = grow_meanTemp)) +
  geom_point(
    aes(colour = surv, alpha = surv), 
    position = position_jitter(height = 0.04), 
    size = 2
  ) +
  scale_alpha_manual(values = c(1, 0.5)) +
  facet_wrap(~ trt)
# maybe?

demo.surv.cl %>%
  ggplot(aes(x = size.prev, y = last.freeze)) +
  geom_point(
    aes(colour = surv, alpha = surv), 
    position = position_jitter(height = 0.25), 
    size = 2
  ) +
  scale_alpha_manual(values = c(1, 0.5)) +
  facet_wrap(~ trt)
# eh... maybe...

demo.surv.cl %>%
  ggplot(aes(x = size.prev, y = summer_meanTemp)) +
  geom_point(
    aes(colour = surv, alpha = surv), 
    position = position_jitter(height = 0.05), 
    size = 2
  ) +
  scale_alpha_manual(values = c(1, 0.5)) +
  facet_wrap(~ trt)
# actually... maybe?

demo.surv.cl %>%
  ggplot(aes(x = size.prev, y = winter_prec.sum)) +
  geom_point(
    aes(colour = surv, alpha = surv), 
    position = position_jitter(height = 0.05), 
    size = 2
  ) +
  scale_alpha_manual(values = c(1, 0.5)) +
  facet_wrap(~ trt)
# hmm... having trouble seeing an effect here...

### Make a raw data plot

surv.year.plot = demo.surv.cl %>%
  mutate(surv = as.numeric(surv)) %>%
  ggplot(aes(x = size.prev)) +
  geom_point(aes(y = surv, colour = trt), position = position_jitter(height = 0.1), size = 3) +
  scale_colour_manual(values = c('black', 'red', 'blue')) +
  facet_wrap(~ paste(Year, trt), ncol = 3)

# there's some effect of size here... but it's really small...

surv.year.fl.plot = demo.surv.cl %>%
  mutate(surv = as.numeric(surv)) %>%
  ggplot(aes(x = size.prev)) +
  geom_point(aes(y = surv, colour = trt), position = position_jitter(height = 0.1), size = 3) +
  scale_colour_manual(values = c('black', 'red', 'blue')) +
  facet_wrap(~ paste(Year, trt, flowering), ncol = 6)
# imbalance in here

### Get model predictions

s_s.f.wp = glmer(
  formula = surv ~ size.prev + flowering + winter_prec.sum + (1 | Plot),
  data = demo.surv.cl,
  family = 'binomial'
)

s_s.wp = glmer(
  formula = surv ~ size.prev + winter_prec.sum + (1 | Plot),
  data = demo.surv.cl,
  family = 'binomial'
)

s_s.wpt = glmer(
  formula = surv ~ size.prev + winter_prec.sum * trt + (1 | Plot),
  data = demo.surv.cl,
  family = 'binomial'
)

AIC(s_s.f.wp, s_s.wp, s_s.wpt) # hmm..
# flowering + winter precip seems to perform best

summary(s_s.f.wint)
summary(s_s.wint)
summary(s_swint)

winter.preds = expand.grid(
  Year = 2017:2022,
  trt = c("control", "drought", "irrigated"),
  size.prev = (5:50)/10,
  flowering = c(TRUE, FALSE)
) %>%
  merge(y = demo.surv.cl %>% distinct(Year, winter_prec.sum)) %>%
  mutate(
    pred.s.f.wp = predict(s_s.f.wp, newdata = ., re.form = NA) %>% ilogit(),
    pred.s.wp = predict(s_s.wp, newdata = ., re.form = NA) %>% ilogit(),
    pred.s.wpt  = predict(s_s.wpt,  newdata = ., re.form = NA) %>% ilogit()
  ) %>%
  pivot_longer(contains('pred.'), names_to = 'model', values_to = 'prediction') %>%
  mutate(model = gsub('pred\\.', '', model))

surv.year.plot +
  geom_line(
    data = winter.preds %>% filter(!grepl('f', model)),
    aes(y = prediction, group = model, linetype = model)
  )
# Slight differences in the left tail (where we have evidence of only a few plants)

surv.year.fl.plot +
  geom_line(
    data = winter.preds,
    aes(y = prediction, group = model, linetype = model)
  )
# not loving the flowering-separate models; few observations to support differences...

summary(s_s.f.wp)
summary(s_s.wp)
summary(s_s.wpt)
# whole lotta n.s. terms in the treatment-wp model
# (slight possibility that wp only has effets on the drought and irrigation models...)

winter.extrap.preds = expand.grid(
  size.prev = (5:50)/10,
  winter_prec.sum = -1:2
) %>%
  mutate(pred = predict(s_s.wp, newdata = ., re.form = NA) %>% ilogit())

demo.surv.cl %>%
  mutate(surv = as.numeric(surv)) %>%
  ggplot(aes(x = size.prev)) +
  geom_point(
    aes(y = surv),
    position = position_jitter(height = 0.1), size = 2, alpha = 0.5
  ) +
  geom_line(
    data = winter.extrap.preds,
    aes(y = pred, group = winter_prec.sum, colour = winter_prec.sum)
  )


### Try growing season models and their fits...

s_s.f.g = glmer(
  formula = surv ~ size.prev + flowering + grow_meanTemp + (1 | Plot),
  data = demo.surv.cl,
  family = 'binomial'
)

s_sg = glmer(
  formula = surv ~ size.prev * grow_meanTemp + (1 | Plot),
  data = demo.surv.cl,
  family = 'binomial'
)

s_s.g = glmer(
  formula = surv ~ size.prev + grow_meanTemp + (1 | Plot),
  data = demo.surv.cl,
  family = 'binomial'
)

s_s.f.gt = glmer(
  formula = surv ~ size.prev + flowering + trt * grow_meanTemp + (1 | Plot),
  data = demo.surv.cl,
  family = 'binomial'
)

growtemp.preds = expand.grid(
  Year = 2017:2022,
  trt = c("control", "drought", "irrigated"),
  size.prev = (5:50)/10,
  flowering = c(TRUE, FALSE)
) %>%
  merge(y = demo.surv.cl %>% distinct(Year, grow_meanTemp)) %>%
  mutate(
    pred.s.f.g  = predict(s_s.f.g, newdata = ., re.form = NA) %>% ilogit(),
    pred.sg     = predict(s_sg, newdata = ., re.form = NA) %>% ilogit(),
    pred.s.g    = predict(s_s.g, newdata = ., re.form = NA) %>% ilogit(),
    pred.s.f.gt = predict(s_s.f.gt, newdata = ., re.form = NA) %>% ilogit()
  ) %>%
  pivot_longer(contains('pred'), names_to = 'model', values_to = 'pred')

surv.year.plot +
  geom_line(
    data = growtemp.preds %>% filter(!grepl('f', model)),
    aes(x = size.prev, y = pred, group = model, linetype = model),
    colour = 'gray44'
  )

surv.year.fl.plot +
  geom_line(
    data = growtemp.preds,
    aes(x = size.prev, y = pred, group = model, linetype = model),
    colour = 'gray44'
  )
