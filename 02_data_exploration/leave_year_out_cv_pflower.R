# Trying a "leave one year out" cross validation approach
# for fitting models of probability of flowering

library(ggplot2)
library(dplyr)
library(tidyr)
library(lme4)

rm(list = ls())

# Inverse logit wrapper
ilogit = function(x) (1+exp(-x))^(-1)

# Read in demography data
demo = read.csv('01_data_cleaning/out/demo_imputed_survival.csv')

head(demo)
nrow(demo)

# Remove plants that weren't surveyed in 2020
demo = demo %>%
  separate(plantid, into = c('tag', 'plot', 'coord'), remove = FALSE, sep = '_') %>%
  filter(
    !( (Year %in% 2020 & plot %in% c(1:2, 13, 15) & grepl('\\d{2}', coord)) |
       (Year %in% 2020 & plot %in% 3 & (grepl('\\d{2}', coord) | grepl('[7890]', coord))))
  ) %>%
  select(-c(tag, plot, coord))

# Add a "flowering" column
# Check to make sure this works
demo %>%
  mutate(flowering = case_when(
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

# Attach this year's flowering to prev year's demographic data
demo.flower = merge(
  # x is year of flowering
  # here: need *surviving* plants only
  # also want flowering column
  x = demo %>%
    filter(surv) %>%
    select(Year, plantid, Plot, trt, flowering, demo.note, proc.note, edited),
  # y is previous year's data
  # do this by adding 1 to year (will allow matching)
  # give *only* the plants that have size measurements
  # NOTE: this will cut out some plants from analysis if we don't have their
  # size in the previous year
  y = demo %>%
    mutate(Year = Year + 1) %>%
    filter(!is.na(No.leaves) & !is.na(Leaf.length) & No.leaves > 0) %>%
    mutate(size.prev = log(No.leaves * Leaf.length)) %>%
    select(Year, plantid, size.prev, flowering, demo.note, proc.note, edited),
  all = FALSE, by = c("Year", "plantid"), suffixes = c("", ".prev")
) %>%
  mutate(Plot = factor(Plot))

# Look at this subcase to make sure we got this right
demo %>% filter(grepl('1411', plantid))
demo.flower %>% filter(grepl('1411', plantid))

head(demo.flower)
nrow(demo.flower)

# Read in climate data
clim = read.csv('01_data_cleaning/out/climate_summary.csv') %>%
  # fix freeze date column (change to numeric)
  mutate(last.freeze = as.numeric(as.Date(last.freeze)) - as.numeric(as.Date('2019-12-31'))) %>%
  # center variables at mean
  group_by(pd.label) %>%
  mutate(across(contains('Temp'), ~ (. - mean(., na.rm = TRUE)))) %>%
  mutate(prec.sum = prec.sum - mean(prec.sum, na.rm = TRUE)) %>%
  ungroup() %>%
  # center last freeze date at some neutral date... like 90
  mutate(last.freeze = last.freeze - 90)
  
head(clim)
str(clim)

# We would expect flowering probabilities to probably be influenced by the 12
# months previous to ~April...
# E.g., in year 2020, we would expect growing season 2020 (?), early season
# 2020, winter 2019/2020, and summer 2019 to influence flowering probabilities
# Maybe also 2019 growing season too? might be a good idea to add this as
# 'pgrow' for 'previous growing season'...

rbind(
  clim,
  clim %>% filter(pd.label %in% 'grow') %>% mutate(Year = Year + 1, pd.label = 'pgrow')
) %>%
  arrange(Year, period.start)
# looks good!

clim.merge = rbind(
  clim,
  clim %>% filter(pd.label %in% 'grow') %>% mutate(Year = Year + 1, pd.label = 'pgrow')
)

# For merging, we want a wide-version of data frame
# can then go through and subset relevant columns

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

clim.merge
# looks okay to me!

demo.fl.clim = merge(x = demo.flower, y = clim.merge, by = 'Year')

head(demo.fl.clim)
nrow(demo.fl.clim)
# Good

demo.fl.clim %>% filter(Year %in% 2022) %>% head()
# Looks right to me, mostly

table(demo.fl.clim$flowering)
table(demo.fl.clim$flowering.prev)
with(demo.fl.clim, table(flowering, flowering.prev))
# wow... that's a lot on the diagonal
# lots of plants that are repeat-flowerers!
table(demo.fl.clim$Year)
# remember - 2016 size data is suspect
# in fact... maybe just take it out? has potential to influence year effects

demo.fl.clim = demo.fl.clim %>% filter(Year > 2017)

##### Want to try a cross-val approach

data.splitter = function(df) {

  yrs = unique(df$Year)
  lis = vector('list', length = length(yrs))
  
  for (y in 1:length(yrs)) {
    lis[[y]]$train = df %>% filter(!(Year %in% yrs[y]))
    lis[[y]]$test  = df %>% filter(Year %in% yrs[y])
  }
  
  return(lis)
  
}

demo.fl.split = data.splitter(demo.fl.clim)

str(demo.fl.split)
length(demo.fl.split)

mod.fit.test = function(mod.formula, data.split) {
  
 try(
   mod <- glmer(
    formula = mod.formula,
    data = data.split$train,
    family = 'binomial'
    )
  )
  
  if (exists('mod')) {
    pred = predict(mod, newdata = data.split$test, re.form = NA)
  
    mod.score = mean(abs(data.split$test$flowering - ilogit(pred)))
  
    if (!length(summary(mod)$fitMsgs)) { return(mod.score)
    } else return(NA)
  } else return(NA)

}

mod.fit.test('flowering ~ size.prev + (1 | Plot / plantid)', demo.fl.split[[1]])
sapply(demo.fl.split, mod.fit.test, mod.formula = 'flowering ~ size.prev + (1 | Plot / plantid)') %>%
  mean()
sapply(demo.fl.split, mod.fit.test, mod.formula = 'flowering ~ size.prev + flowering.prev + (1 | Plot / plantid)') %>%
  mean()
# welp... singularities
# okay well probably get rid of plant-level random effects

sapply(demo.fl.split, mod.fit.test, mod.formula = 'flowering ~ size.prev + (1 | Plot)') %>%
  mean()
sapply(demo.fl.split, mod.fit.test, mod.formula = 'flowering ~ size.prev + flowering.prev + (1 | Plot)') %>%
  mean()
# cool

### Fit some base models

f_0 = glmer(
  formula = flowering ~ (1 | Year) + (1 | Plot),
  data = demo.fl.clim,
  family = 'binomial'
)

f_s = glmer(
  formula = flowering ~ size.prev + (1 | Year) + (1 | Plot),
  data = demo.fl.clim,
  family = 'binomial'
)

f_s.f = glmer(
  formula = flowering ~ size.prev + flowering.prev + (1 | Year) + (1 | Plot),
  data = demo.fl.clim,
  family = 'binomial'
)

f_sf = glmer(
  formula = flowering ~ size.prev * flowering.prev + (1 | Year) + (1 | Plot),
  data = demo.fl.clim,
  family = 'binomial'
)

AIC(f_sf, f_s.f, f_s, f_0) %>% mutate(daic = AIC - min(AIC))

f_s.f.t = glmer(
  formula = flowering ~ size.prev + flowering.prev + trt + (1 | Year) + (1 | Plot),
  data = demo.fl.clim,
  family = 'binomial'
)

f_s.ft = glmer(
  formula = flowering ~ size.prev + flowering.prev * trt + (1 | Year) + (1 | Plot),
  data = demo.fl.clim,
  family = 'binomial'
)

f_st.f = glmer(
  formula = flowering ~ size.prev * trt + flowering.prev + (1 | Year) + (1 | Plot),
  data = demo.fl.clim,
  family = 'binomial'
)

f_st.ft = glmer(
  formula = flowering ~ size.prev * trt + flowering.prev * trt + (1 | Year) + (1 | Plot),
  data = demo.fl.clim,
  family = 'binomial'
)

AIC(f_s.f.t, f_st.f, f_s.ft, f_st.ft, f_s.f) %>% mutate(daic = AIC - min(AIC))
# interesting...
# effects of flowering very by treatment?

summary(f_s.ft)
# so experimental treatments have lower flowering probabilities if they did not
# flower last year
# but it looks like flowering previously negates that?

##### Anyway let's now try fitting some models...

forms = c(
  'flowering ~ size.prev + flowering.prev * trt + (1 | Year) + (1 | Plot)',
  'flowering ~ size.prev + flowering.prev * trt + last.freeze + (1 | Plot)',
  'flowering ~ size.prev + flowering.prev * trt + early_prec.sum + (1 | Plot)',
  'flowering ~ size.prev + flowering.prev * trt + early_minnTemp + (1 | Plot)',
  'flowering ~ size.prev + flowering.prev * trt + grow_meanTemp + (1 | Plot)',
  'flowering ~ size.prev + flowering.prev * trt + summer_meanTemp + (1 | Plot)',
  'flowering ~ size.prev + flowering.prev * trt + winter_prec.sum + (1 | Plot)',
  'flowering ~ size.prev + flowering.prev * trt + winter_meanTemp + (1 | Plot)',
  'flowering ~ size.prev + flowering.prev * trt + pgrow_meanTemp + (1 | Plot)',
  'flowering ~ size.prev + flowering.prev * trt + (trt | Year) + (1 | Plot)',
  'flowering ~ size.prev + flowering.prev * trt + last.freeze * trt + (1 | Plot)',
  'flowering ~ size.prev + flowering.prev * trt + early_prec.sum * trt + (1 | Plot)',
  'flowering ~ size.prev + flowering.prev * trt + early_minnTemp * trt + (1 | Plot)',
  'flowering ~ size.prev + flowering.prev * trt + grow_meanTemp * trt + (1 | Plot)',
  'flowering ~ size.prev + flowering.prev * trt + summer_meanTemp * trt + (1 | Plot)',
  'flowering ~ size.prev + flowering.prev * trt + winter_prec.sum * trt + (1 | Plot)',
  'flowering ~ size.prev + flowering.prev * trt + winter_meanTemp * trt + (1 | Plot)',
  'flowering ~ size.prev + flowering.prev * trt + pgrow_meanTemp * trt + (1 | Plot)',
  'flowering ~ flowering.prev * trt + (size.prev | Year) + (1 | Plot)',
  'flowering ~ size.prev * last.freeze + flowering.prev * trt + (1 | Plot)',
  'flowering ~ size.prev * early_prec.sum + flowering.prev * trt + (1 | Plot)',
  'flowering ~ size.prev * early_minnTemp + flowering.prev * trt + (1 | Plot)',
  'flowering ~ size.prev * grow_meanTemp + flowering.prev * trt + (1 | Plot)',
  'flowering ~ size.prev * summer_meanTemp + flowering.prev * trt + (1 | Plot)',
  'flowering ~ size.prev * winter_prec.sum + flowering.prev * trt + (1 | Plot)',
  'flowering ~ size.prev * winter_meanTemp + flowering.prev * trt + (1 | Plot)',
  'flowering ~ size.prev * pgrow_meanTemp + flowering.prev * trt + (1 | Plot)',
  'flowering ~ size.prev * trt + flowering.prev * trt + (1 | Year) + (1 | Plot)',
  'flowering ~ size.prev * last.freeze * trt + flowering.prev * trt + (1 | Plot)',
  'flowering ~ size.prev * early_prec.sum * trt + flowering.prev * trt + (1 | Plot)',
  'flowering ~ size.prev * early_minnTemp * trt + flowering.prev * trt + early_minnTemp + (1 | Plot)',
  'flowering ~ size.prev * grow_meanTemp * trt + flowering.prev * trt + (1 | Plot)',
  'flowering ~ size.prev * summer_meanTemp * trt + flowering.prev * trt + (1 | Plot)',
  'flowering ~ size.prev * winter_meanTemp * trt + flowering.prev * trt + (1 | Plot)',
  'flowering ~ size.prev * winter_prec.sum * trt + flowering.prev * trt + (1 | Plot)',
  'flowering ~ size.prev * pgrow_meanTemp * trt + flowering.prev * trt + (1 | Plot)'
)

# form.scores = vector(mode = 'double', length = length(forms))
form.scores = matrix(ncol = length(demo.fl.split), nrow = length(forms))

for (form in 1:length(forms)) {
  form.scores[form,] = sapply(
    X = demo.fl.split,
    FUN = mod.fit.test,
    mod.formula = forms[form])
}
# many convergence warnings
form.scores
# uh... okay well I guess all of these ran?

data.frame(form.scores) %>%
  mutate(form = 1:nrow(.)) %>%
  pivot_longer(-form, names_to = 'holdout', values_to = 'score') %>%
  ggplot(aes(x = form, y = score, group = holdout)) +
  geom_line(aes(colour = holdout))
# whoa!

# data.frame(score = form.scores) %>% mutate(score.delt = score - min(score))
# interesting...
order(rowSums(form.scores))
forms[order(rowSums(form.scores))]
# it looks like summer temps are the most predictive
# followed by early precipitation

best.mod.q = sapply(
  X = demo.fl.split,
  FUN = mod.fit.test,
  mod.formula = forms[which.min(rowSums(form.scores))]
)
# one convergence issue
best.mod.q
# still considerable variation...

f_s.ft.summert = glmer(
  formula = flowering ~ size.prev + flowering.prev * trt + summer_meanTemp * trt + (1 | Plot),
  family = binomial,
  data = demo.fl.clim,
  control = glmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 1e6))
)

summary(f_s.ft.summert)

f_sst.ft = glmer(
  formula = flowering ~ size.prev * summer_meanTemp * trt + flowering.prev * trt + (1 | Plot),
  family = binomial,
  data = demo.fl.clim,
  control = glmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 1e6))
)

f_s.ft.summer = glmer(
  formula = flowering ~ size.prev + flowering.prev * trt + summer_meanTemp + (1 | Plot),
  family = binomial,
  data = demo.fl.clim,
  control = glmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 1e6))
)

f_st.ft.summer = glmer(
  formula = flowering ~ size.prev * summer_meanTemp + flowering.prev * trt + (1 | Plot),
  family = binomial,
  data = demo.fl.clim,
  control = glmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 1e6))
)

AIC(f_s.ft.summer, f_s.ft.summert, f_sst.ft, f_st.ft.summer) %>%
  mutate(daic = AIC - min(AIC))
# okay... from this, it looks like the full model (all interactions)
# is best

summary(f_sst.ft)
# yeah that's a lot of stuff, huh

predict.backbone = expand.grid(
  size.prev = c(1:50)/10,
  trt = c('control', 'drought', 'irrigated'),
  flowering.prev = c(TRUE, FALSE),
  summer_meanTemp = (-2:2)
) %>%
  mutate(pred = predict(f_sst.ft, newdata = ., re.form = NA) %>% ilogit())

head(predict.backbone)

predict.backbone %>%
  ggplot(aes(x = size.prev, y = pred)) +
  geom_line(
    aes(
      colour = trt, 
      linetype = flowering.prev,
      group = interaction(flowering.prev, trt)
    )
  ) +
  scale_colour_manual(values = c('black', 'red', 'blue')) +
  facet_wrap(~ summer_meanTemp)
# warm summers increase probability of flowering...
# control plants especially unlikely to flower? I don't buy it.
# hot summers - drought plants more likely to flower than irrigated?
# don't buy it.

# let's try now with the second-best supported model...

predict.backbone = expand.grid(
  size.prev = c(1:50)/10,
  trt = c('control', 'drought', 'irrigated'),
  flowering.prev = c(TRUE, FALSE),
  summer_meanTemp = .5 + (-2:1)
) %>%
  mutate(pred = predict(f_s.ft.summert, newdata = ., re.form = NA) %>% ilogit())

head(predict.backbone)

predict.backbone %>%
  ggplot(aes(x = size.prev, y = pred)) +
  geom_line(
    aes(
      colour = trt, 
      linetype = flowering.prev,
      group = interaction(flowering.prev, trt)
    )
  ) +
  scale_colour_manual(values = c('black', 'red', 'blue')) +
  facet_wrap(~ summer_meanTemp)
# jesus
# I could believe that warmer years mean more flowering
# but here, drought is favored in warm years? simply can't be true!

# Try some sort of raw data plot...
# bin by size?

demo.fl.clim %>%
  mutate(size.round = round(size.prev, 1)) %>%
  ggplot(aes(x = size.round)) +
  geom_histogram(binwidth = 0.1) +
  facet_grid(rows = vars(trt), cols = vars(Year))

demo.fl.clim %>%
  mutate(size.round = round(size.prev, 1)) %>%
  group_by(Year, summer_meanTemp, size.round,  trt, flowering.prev) %>%
  summarise(
    n = n(),
    p.flower = mean(flowering)
  ) %>%
  ggplot(aes(x = size.round, y = p.flower)) +
  geom_line(aes(colour = trt, group = interaction(flowering.prev, trt)), linewidth = 0.25) +
  geom_point(aes(colour = trt, shape = flowering.prev), size = 3, alpha = 0.5) +
  scale_colour_manual(values = c('black', 'red', 'blue')) +
  facet_wrap(~ Year)

# Try a plot with summer temp on the y-axis?
# color can tell us flowering status
demo.fl.clim %>%
  ggplot(aes(x = size.prev, y = summer_meanTemp)) +
  geom_point(
    aes(colour = flowering), 
    size = 3, alpha = 0.5, position = position_jitter(height = 0.05)
  ) +
  facet_wrap(flowering.prev ~ trt)
# Yeah... looks like the warmer years are indeed more conducive to flowering
# especially for previously-flowering plants (interesting!)
# the two cooler years seem to have the least flowering, although it looks like
# one of these years is much cooler than the other

demo.fl.clim %>%
  ggplot(aes(x = size.prev, y = early_prec.sum)) +
  geom_point(
    aes(colour = flowering), 
    size = 3, alpha = 0.5, position = position_jitter(height = 2)
  ) +
  facet_wrap(flowering.prev ~ trt)
# Not such a clear relationship here with early precipitation

# Previous year's growing season
demo.fl.clim %>%
  ggplot(aes(x = size.prev, y = pgrow_meanTemp)) +
  geom_point(
    aes(colour = flowering), 
    size = 3, alpha = 0.5, position = position_jitter(height = 0.05)
  ) +
  facet_wrap(flowering.prev ~ trt)
# also does not look very strong

# Kinda want to see if the ranges of previous flowering and size overlap?
# maybe the flowering previous is just an effect of size?
# (but then again... previously-flowering plants would grow more???)
demo.fl.clim %>%
  ggplot(aes(y = size.prev, x = Year, group = flowering.prev)) +
  geom_point(
    aes(colour = flowering), 
    size = 3, alpha = 0.5, position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.5)
  ) +
  facet_wrap(~ trt, ncol = 3)

# what if we removed that size effect...

demo.fl.clim %>%
  # group_by(Year, trt, flowering.prev) %>%
  # mutate(size.prev = size.prev - mean(size.prev)) %>%
  # ungroup() %>%
  ggplot(aes(y = size.prev, x = Year, group = interaction(flowering, flowering.prev))) +
  geom_point(
    aes(colour = flowering, shape = flowering.prev), 
    size = 3, alpha = 0.5, position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75)
  ) +
  facet_wrap(~ trt, nrow = 3)

##### Plot of two different model fits against raw data

# Model with size-trt-weather interaction
flower.backbone.year = expand.grid(
  size.prev = (5:50)/10,
  Year = 2017:2023,
  trt = c('control', 'drought', 'irrigated'),
  flowering.prev = c(TRUE, FALSE)
) %>%
  merge(y = demo.fl.clim %>% distinct(Year, summer_meanTemp)) %>%
  mutate(pred = predict(f_sst.ft, newdata = ., re.form = NA) %>% ilogit())

# Model without size-trt-weather interaction
flower.backbone.year2 = expand.grid(
  size.prev = (5:50)/10,
  Year = 2017:2023,
  trt = c('control', 'drought', 'irrigated'),
  flowering.prev = c(TRUE, FALSE)
) %>%
  merge(y = demo.fl.clim %>% distinct(Year, summer_meanTemp)) %>%
  mutate(pred = predict(f_s.ft.summert, newdata = ., re.form = NA) %>% ilogit())

# Plot
demo.fl.clim %>%
  mutate(flowering = as.numeric(flowering)) %>%
  ggplot(aes(x = size.prev)) +
  geom_point(
    aes(y = flowering, shape = flowering.prev, colour = trt),
    position = position_jitter(height = 0.1), alpha = 0.5, size = 3
  ) +
  geom_line(
    data = flower.backbone.year,
    aes(y = pred, colour = trt)
  ) +
  geom_line(
    data = flower.backbone.year2,
    aes(y = pred, colour = trt),
    linetype = 2
  ) +
  scale_colour_manual(values = c('black', 'red', 'blue')) +
  facet_wrap(~ paste(Year, trt, flowering.prev))
  # facet_wrap(~ Year + trt + flowering.prev)

# Well... 
# The differences in these two models looks to all be in the control treatment
# and the differences here seem to be in two years (2020, 2021), which kinda have sparse data
# and some of the discrepancy comes in a zone of extrapolation...

############################################################
##### Jennifer's suggestion of pulling out new recruits...

# Get demo that has all years (incl. 2016)
demo.fl.dupe = merge(x = demo.flower, y = clim.merge, by = 'Year')

nrow(demo.fl.dupe)

demo.fl.dupe = demo.fl.dupe %>%
  arrange(plantid, Year) %>%
  mutate(first.obs = !duplicated(plantid))

table(demo.fl.dupe$first.obs)

demo.fl.dupe %>%
  ggplot(aes(x = Year, y = size.prev, group = first.obs)) +
  geom_point(
    aes(colour = flowering, shape = first.obs), 
    size = 2, alpha = 0.75,
    position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.5)
  )
# hmm...
# new recruits are definitely smaller, although in some years they are smaller than others

# Let's try fitting models now with only the already-observed plants

demo.fl.dupe = demo.fl.dupe %>%
  filter(!first.obs) %>%
  select(-first.obs)

nrow(demo.fl.dupe)
head(demo.fl.dupe)
table(demo.fl.dupe$Year) # we lost a year
table(demo.fl.dupe$flowering)

d_0 = glmer(
  formula = flowering ~ (1 | Year) + (1 | Plot),
  data = demo.fl.dupe,
  family = 'binomial'
)

d_s = glmer(
  formula = flowering ~ size.prev + (1 | Year) + (1 | Plot),
  data = demo.fl.dupe,
  family = 'binomial'
)

d_s.f = glmer(
  formula = flowering ~ size.prev + flowering.prev + (1 | Year) + (1 | Plot),
  data = demo.fl.dupe,
  family = 'binomial'
)

d_sf = glmer(
  formula = flowering ~ size.prev * flowering.prev + (1 | Year) + (1 | Plot),
  data = demo.fl.dupe,
  family = 'binomial'
)

AIC(d_0, d_s, d_s.f, d_sf)
# qualitatively similar to above

summary(d_s.f)
# still a positive effect
# coefficients are very similar
