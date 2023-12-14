library(ggplot2)
library(dplyr)
library(tidyr)
library(lme4)

rm(list = ls())

# Logit wrapper function
ilogit = function(x) (1+exp(-x))^-1

# Read in demography data
demo = read.csv('01_data_cleaning/out/demo_imputed_survival.csv') %>%
  mutate(
    flowering = case_when(
      No.umbels > 0 ~ TRUE,
      !No.umbels ~ FALSE,
      is.na(No.umbels) ~ FALSE,
      .default = NA
    )
  )

head(demo)
nrow(demo)

### Get this year's growth and last year's demographic data
demo.grow = merge(
  # x is response - size in year t; we'll use convention that it grew in year t-1
  # so need to subtract 1 from year
  # also want to get only plants with leaf counts + measures,
  # and only living plants (those with non-zero leaf counts)
  x = demo %>%
    filter(!is.na(Leaf.length) & !is.na(No.leaves)) %>%
    filter(No.leaves > 0, surv) %>%
    mutate(
      Year = Year - 1,
      size = log(Leaf.length * No.leaves)
    ) %>%
    select(Year, plantid, Plot, trt, size, flowering, demo.note, proc.note, edited),
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

nrow(demo.grow)
head(demo.grow)

# # Remove 2016 plants # (actually I think this is handled later...)

length(unique(demo.grow$plantid))
table(demo.grow$Year)


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
# think carefully about the merge...
# growth from 2016 to 2017 should depend on growing season 2016, winter 2016,
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

demo.grow.cl = merge(x = demo.grow, y = clim.merge, by = 'Year') %>%
  filter(Year > 2016)

head(demo.grow.cl)
nrow(demo.grow.cl)

##### Some plots

ggplot(demo.grow.cl) +
  geom_point(aes(x = size.prev, y = size, colour = trt), alpha = 0.5, size = 3) +
  geom_segment(aes(x = 1, xend = 5, y = 1, yend = 5), linetype = 2) +
  facet_wrap(~ Year)

# Growth differences in flowering *this year*
demo.grow.cl %>%
  mutate(gr.lin = size - size.prev) %>%
  ggplot() +
  geom_point(
    aes(x = flowering, y = gr.lin, colour = trt), 
    alpha = 0.5, size = 3,
    position = position_jitter(width = 0.25)
  ) +
  geom_segment(aes(x = 1, xend = 2, y = 0, yend = 0), linetype = 2) +
  facet_wrap(~ Year)

# Growth differences in flowering *last year*
demo.grow.cl %>%
  mutate(gr.lin = size - size.prev) %>%
  ggplot() +
  geom_point(
    aes(x = flowering.prev, y = gr.lin, colour = trt), 
    alpha = 0.5, size = 3,
    position = position_jitter(width = 0.25)
  ) +
  geom_segment(aes(x = 1, xend = 2, y = 0, yend = 0), linetype = 2) +
  facet_wrap(~ Year)

##### Try fitting base models

g_0 = lmer(
  formula = size ~ size.prev + (1 | Year) + (1 | Plot),
  data = demo.grow.cl
)

summary(g_0) # good
ranef(g_0)$Year %>% unlist() %>% plot(type = 'l')
# hmm... possibly autocorrelated...

g_f0 = lmer(
  formula = size ~ size.prev + flowering.prev + (1 | Year) + (1 | Plot),
  data = demo.grow.cl
)

AIC(g_0, g_f0)
# keep flowering in prior year

g_f0.f1 = lmer(
  formula = size ~ size.prev + flowering + flowering.prev + (1 | Year) + (1 | Plot),
  data = demo.grow.cl
)

AIC(g_0, g_f0, g_f0.f1)
# flowering *this year* seems to have a huge impact

g_f0.f1.t = lmer(
  formula = size ~ size.prev + flowering + flowering.prev + trt + (1 | Year) + (1 | Plot),
  data = demo.grow.cl
)

AIC(g_0, g_f0, g_f0.f1, g_f0.f1.t)
# no straightforward treatment effect...

g_f0t.f1 = lmer(
  formula = size ~ size.prev + flowering + flowering.prev * trt + (1 | Year) + (1 | Plot),
  data = demo.grow.cl
)

g_f0.f1t  = lmer(
  formula = size ~ size.prev + flowering * trt + flowering.prev + (1 | Year) + (1 | Plot),
  data = demo.grow.cl
)

g_st.f0.f1 = lmer(
  formula = size ~ size.prev * trt + flowering + flowering.prev + (1 | Year) + (1 | Plot),
  data = demo.grow.cl
)

g_f0t.f1t = lmer(
  formula = size ~ size.prev + flowering * trt + flowering.prev * trt + (1 | Year) + (1 | Plot),
  data = demo.grow.cl
)

AIC(g_f0.f1, g_st.f0.f1, g_f0t.f1, g_f0.f1t, g_f0t.f1t)
# no treatment effects

# Size-flowering interactions

g_f0.sf1 = lmer(
  formula = size ~ size.prev * flowering  + flowering.prev  + (1 | Year) + (1 | Plot),
  data = demo.grow.cl
)

g_sf0.f1 = lmer(
  formula = size ~ flowering  + size.prev * flowering.prev + (1 | Year) + (1 | Plot),
  data = demo.grow.cl
)

AIC(g_f0.f1, g_sf0.f1, g_f0.sf1)
# size interactions with both terms plausible...

g_sf0.sf1 = lmer(
  formula = size ~ size.prev * flowering  + size.prev * flowering.prev  + (1 | Year) + (1 | Plot),
  data = demo.grow.cl
)

AIC(g_f0.f1, g_sf0.sf1, g_sf0.f1, g_f0.sf1)
# well... no clear best, but all are better than interaction-less model

g_sf0.sf1t = lmer(
  formula = size ~ size.prev * flowering * trt + size.prev * flowering.prev  + (1 | Year) + (1 | Plot),
  data = demo.grow.cl
)

g_sf0t.sf1 = lmer(
  formula = size ~ size.prev * flowering + size.prev * flowering.prev * trt  + (1 | Year) + (1 | Plot),
  data = demo.grow.cl
)

g_sf0t.sf1t = lmer(
  formula = size ~ size.prev * flowering * trt + size.prev * flowering.prev * trt  + (1 | Year) + (1 | Plot),
  data = demo.grow.cl
)

AIC(g_f0.f1, g_sf0.sf1, g_sf0.f1, g_f0.sf1, g_sf0t.sf1, g_sf0.sf1t, g_sf0t.sf1t) %>%
  arrange(AIC)
# still no treatment effects!

# Okay... which of these to keep?

summary(g_f0.sf1)
ranef(g_f0.sf1)$Plot %>% unlist() %>% hist()
ranef(g_f0.sf1)$Year %>% unlist() %>% hist()
ranef(g_f0.sf1)$Year %>% unlist() %>% plot(type = 'l')

summary(g_sf0.f1)
ranef(g_sf0.f1)$Plot %>% unlist() %>% hist()
ranef(g_sf0.f1)$Year %>% unlist() %>% hist()
ranef(g_sf0.f1)$Year %>% unlist() %>% plot(type = 'l')
# random effects etc. look near-identical

# Previous size and previous flowering make more sense to me.
# Maybe do some plotting...

noyear.preds = expand.grid(
  Year = 2017:2022,
  size.prev = (5:50)/10,
  flowering = c(TRUE, FALSE),
  flowering.prev = c(TRUE, FALSE)
) %>%
  mutate(
    pred.sf0.f1 = predict(g_sf0.f1, newdata = ., re.form = NA),
    pred.f0.sf1 = predict(g_f0.sf1, newdata = ., re.form = NA),
  ) %>%
  pivot_longer(contains('pred'), names_to = 'model', values_to = 'prediction') %>%
  mutate(model = gsub('pred\\.', '', model))

demo.grow.cl %>%
  ggplot(aes(x = size.prev)) +
  geom_point(aes(y = size), colour = 'orange') +
  geom_line(
    data = noyear.preds,
    aes(y = prediction, group = model, linetype = model)
  ) +
  facet_wrap(~ paste(Year, 'p', flowering.prev, flowering), ncol = 4)

# Ugh... idk man.
# I guess try fitting with each.
# Lotsa models.

### Split dataset by year

demo.years = sort(unique(demo.grow.cl$Year))
demo.grow.split = vector(mode = 'list', length = length(demo.years))

for (i in 1:length(demo.years)) {
  demo.grow.split[[i]]$train = demo.grow.cl %>% filter(!(Year %in% demo.years[i]))
  demo.grow.split[[i]]$test  = demo.grow.cl %>% filter(Year %in% demo.years[i])
  demo.grow.split[[i]]$holdout = demo.years[i] 
}

### Wrapper function for model fitting
fit.split.model = function(dfs, mod.formula) {
  
  try( 
    mod <- lmer(
      formula = mod.formula,
      data = dfs$train #,
      # control = lmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 1e6))
    )
  )
  
  # NOTE: here, we're doing the predictions including the plot-level random effects.
  if (exists('mod')) { 
    mod.pred  = predict(mod, newdata = dfs$test, re.form = ~ (1 | Plot)) %>% ilogit()
    mod.score = mean((dfs$test$size - mod.pred))^2
    return(mod.score)
  } else { return(NA) }
}

fit.split.model(demo.grow.split[[1]], mod.formula = 'size ~ size.prev + (1 | Year) + (1 | Plot)')
# seems alright

### Define formulas

forms = c(
  'size ~ size.prev * flowering.prev + flowering + (size.prev | Year) + (1 | Plot)',
  'size ~ size.prev * flowering.prev + flowering + last.freeze + (1 | Year) + (1 | Plot)',
  'size ~ size.prev * flowering.prev + flowering + grow_meanTemp + (1 | Year) + (1 | Plot)',
  'size ~ size.prev * flowering.prev + flowering + early_prec.sum + (1 | Year) + (1 | Plot)',
  'size ~ size.prev * flowering.prev + flowering + early_minnTemp + (1 | Year) + (1 | Plot)',
  'size ~ size.prev * flowering.prev + flowering + winter_prec.sum + (1 | Year) + (1 | Plot)',
  'size ~ size.prev * flowering.prev + flowering + winter_meanTemp + (1 | Year) + (1 | Plot)',
  'size ~ size.prev * flowering.prev + flowering + summer_meanTemp + (1 | Year) + (1 | Plot)',
  'size ~ size.prev * flowering.prev + flowering + pgrow_meanTemp + (1 | Year) + (1 | Plot)',
  'size ~ size.prev * flowering + flowering.prev + (1 | Year) + (1 | Plot)',
  'size ~ size.prev * flowering + flowering.prev + last.freeze + (1 | Year) + (1 | Plot)',
  'size ~ size.prev * flowering + flowering.prev + grow_meanTemp + (1 | Year) + (1 | Plot)',
  'size ~ size.prev * flowering + flowering.prev + early_prec.sum + (1 | Year) + (1 | Plot)',
  'size ~ size.prev * flowering + flowering.prev + early_minnTemp + (1 | Year) + (1 | Plot)',
  'size ~ size.prev * flowering + flowering.prev + winter_prec.sum + (1 | Year) + (1 | Plot)',
  'size ~ size.prev * flowering + flowering.prev + winter_meanTemp + (1 | Year) + (1 | Plot)',
  'size ~ size.prev * flowering + flowering.prev + summer_meanTemp + (1 | Year) + (1 | Plot)',
  'size ~ size.prev * flowering + flowering.prev + pgrow_meanTemp + (1 | Year) + (1 | Plot)',
  'size ~ size.prev * flowering.prev * last.freeze + flowering + (1 | Year) + (1 | Plot)',
  'size ~ size.prev * flowering.prev * grow_meanTemp + flowering + (1 | Year) + (1 | Plot)',
  'size ~ size.prev * flowering.prev * early_prec.sum + flowering + (1 | Year) + (1 | Plot)',
  'size ~ size.prev * flowering.prev * early_minnTemp + flowering + (1 | Year) + (1 | Plot)',
  'size ~ size.prev * flowering.prev * winter_prec.sum + flowering + (1 | Year) + (1 | Plot)',
  'size ~ size.prev * flowering.prev * winter_meanTemp + flowering + (1 | Year) + (1 | Plot)',
  'size ~ size.prev * flowering.prev * summer_meanTemp + flowering + (1 | Year) + (1 | Plot)',
  'size ~ size.prev * flowering.prev * pgrow_meanTemp + flowering + (1 | Year) + (1 | Plot)',
  'size ~ size.prev * flowering * last.freeze + flowering.prev + (1 | Year) + (1 | Plot)',
  'size ~ size.prev * flowering * grow_meanTemp + flowering.prev + (1 | Year) + (1 | Plot)',
  'size ~ size.prev * flowering * early_prec.sum + flowering.prev + (1 | Year) + (1 | Plot)',
  'size ~ size.prev * flowering * early_minnTemp + flowering.prev + (1 | Year) + (1 | Plot)',
  'size ~ size.prev * flowering * winter_prec.sum + flowering.prev + (1 | Year) + (1 | Plot)',
  'size ~ size.prev * flowering * winter_meanTemp + flowering.prev + (1 | Year) + (1 | Plot)',
  'size ~ size.prev * flowering * summer_meanTemp + flowering.prev + (1 | Year) + (1 | Plot)',
  'size ~ size.prev * flowering * pgrow_meanTemp + flowering.prev + (1 | Year) + (1 | Plot)',
  'size ~ size.prev * flowering.prev + flowering * last.freeze + (1 | Year) + (1 | Plot)',
  'size ~ size.prev * flowering.prev + flowering * grow_meanTemp + (1 | Year) + (1 | Plot)',
  'size ~ size.prev * flowering.prev + flowering * early_prec.sum + (1 | Year) + (1 | Plot)',
  'size ~ size.prev * flowering.prev + flowering * early_minnTemp + (1 | Year) + (1 | Plot)',
  'size ~ size.prev * flowering.prev + flowering * winter_prec.sum + (1 | Year) + (1 | Plot)',
  'size ~ size.prev * flowering.prev + flowering * winter_meanTemp + (1 | Year) + (1 | Plot)',
  'size ~ size.prev * flowering.prev + flowering * summer_meanTemp + (1 | Year) + (1 | Plot)',
  'size ~ size.prev * flowering.prev + flowering * pgrow_meanTemp + (1 | Year) + (1 | Plot)',
  'size ~ size.prev * flowering + flowering.prev * last.freeze + (1 | Year) + (1 | Plot)',
  'size ~ size.prev * flowering + flowering.prev * grow_meanTemp + (1 | Year) + (1 | Plot)',
  'size ~ size.prev * flowering + flowering.prev * early_prec.sum + (1 | Year) + (1 | Plot)',
  'size ~ size.prev * flowering + flowering.prev * early_minnTemp + (1 | Year) + (1 | Plot)',
  'size ~ size.prev * flowering + flowering.prev * winter_prec.sum + (1 | Year) + (1 | Plot)',
  'size ~ size.prev * flowering + flowering.prev * winter_meanTemp + (1 | Year) + (1 | Plot)',
  'size ~ size.prev * flowering + flowering.prev * summer_meanTemp + (1 | Year) + (1 | Plot)',
  'size ~ size.prev * flowering + flowering.prev * pgrow_meanTemp + (1 | Year) + (1 | Plot)'
)

form.scores = matrix(NA, nrow = length(forms), ncol = length(demo.years))

for (i in 1:length(forms)) {
  form.scores[i,] = sapply(
    X = demo.grow.split,
    FUN = fit.split.model,
    mod.formula = forms[i]
  )
  print(i)
}
# why the fuck am I getting singularities
# fuck man
# I mean maybe this is good... it means no variance in random effects after incorporating slope...?

form.scores

form.scores %>%
  as.data.frame() %>%
  mutate(form.no = 1:nrow(.)) %>%
  pivot_longer(-form.no, names_to = 'holdout', values_to = 'score') %>%
  ggplot(aes(x = form.no, y = score, group = holdout, colour = holdout)) +
  geom_line()
# okay... really minimal differences among these...
# maybe there just isn't very much annual variation...
# even adding in these shitty motherfucking shits none of these fucking results change
# they're all identical to the year-only model
# what the fuck is happening?

data.frame(
  form = forms,
  mean.score = rowMeans(form.scores) %>% round(2)
) %>%
  arrange(mean.score)
# these are all fucking identical! what the absolute fuck? fuck this man.

g_f0.f1 = lmer(
  formula = size ~ flowering  + flowering.prev  + (size.prev | Year) + (1 | Plot),
  data = demo.grow.cl
)

g_f0.sf1 = lmer(
  formula = size ~ size.prev * flowering  + flowering.prev  + (size.prev | Year) + (1 | Plot),
  data = demo.grow.cl
)

g_sf0.f1 = lmer(
  formula = size ~ flowering + size.prev * flowering.prev + (size.prev | Year) + (1 | Plot),
  data = demo.grow.cl
)

AIC(g_f0.f1, g_f0.sf1, g_sf0.f1)
# okay... definitely some interaction here

summary(g_sf0.f1)
ranef(g_sf0.f1)
# really not that much variation in slope!

### Fit one or two of these models to see what the coefficients are (maybe they are shit!)

g_sf0.f1.emin = lmer(
  formula = size ~ size.prev * flowering.prev + flowering + early_minnTemp + (1 | Year) + (1 | Plot),
  data = demo.grow.cl
)

summary(g_f0.f1.emin)
# effect is very weak relative to uncertainty

g_sf0.f1.emin = lmer(
  formula = size ~ size.prev * flowering.prev + flowering + early_minnTemp + (size.prev | Year) + (1 | Plot),
  data = demo.grow.cl
)

summary(g_f0.f1.emin)
# variation in slope - still very small
# although now temperature effects are larger!

# plot against the least shitty climate variables

demo.grow.cl %>%
  mutate(g.lin = size - size.prev) %>%
  ggplot(aes(x = g.lin, fill = interaction(flowering.prev, flowering))) +
  geom_histogram(position = position_identity(), binwidth = 0.25, alpha = 0.5) +
  facet_wrap(~ Year) +
  theme(legend.position = 'bottom')

demo.grow.cl %>%
  mutate(g.lin = size - size.prev) %>%
  ggplot(aes(x = early_minnTemp)) +
  geom_segment(aes(x = -5, xend = 2, y = 0, yend = 0), linetype = 2) +
  geom_point(
    aes(y = g.lin, colour = interaction(flowering.prev, flowering)),
    alpha = 0.5,
    position = position_jitterdodge(dodge.width = 0.25),
  )

demo.grow.cl %>%
  mutate(g.lin = size - size.prev) %>%
  ggplot(aes(x = grow_meanTemp)) +
  geom_segment(aes(x = -3, xend = 1, y = 0, yend = 0), linetype = 2) +
  geom_point(
    aes(y = g.lin, colour = interaction(flowering.prev, flowering)),
    alpha = 0.5,
    position = position_jitterdodge(dodge.width = 0.25),
  )

demo.grow.cl %>%
  mutate(g.lin = size - size.prev) %>%
  ggplot(aes(x = pgrow_meanTemp)) +
  geom_segment(aes(x = -1, xend = 2, y = 0, yend = 0), linetype = 2) +
  geom_point(
    aes(y = g.lin, colour = interaction(flowering.prev, flowering)),
    alpha = 0.5,
    position = position_jitterdodge(dodge.width = 0.25),
  )

demo.grow.cl %>%
  ggplot(aes(x = size.prev, y = size, colour = grow_meanTemp)) +
  geom_point() +
  geom_segment(aes(x = 1, xend = 5, y = 1, yend = 5), linetype = 2, colour = 'white')
# lmao, fuck


