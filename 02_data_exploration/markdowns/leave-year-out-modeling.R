#' ---
#' title: 'A leave-year-out approach to model selection'
#' output: github_document
#' ---

#' We're trying to build models that link climate data to demographic data to
#' estimate vital rates and putative drivers. This is a non-trivial task. I will
#' set aside which specific variables will be included in our "candidate
#' variable" list for later. For now, I'll use a preliminary set of variables
#' that I estimated from Jenna's climate dataset. Ideally, if a different set of
#' candidate variables is desired, this workflow can be adapted to work with
#' those variables pretty easily.
#' 
#' In our dataset, we have fewer than ten years of data for estimating each
#' vital rate. I think with so few distinct "climates" observed in our dataset,
#' we should be conservative and only include one climate variable in our model.
#' The question, then, is which variable from our "candidate list" to include?
#' 
#' A useful paper for considering this is Tredennick et al. (2021, Ecology).
#' They identified three different modes or goals of data analysis: exploratory,
#' inferential, and predictive. They then analyzed the same dataset in three
#' different ways and highlighted how they might provide different results. In
#' adapting this to our project, I will first note that I think this project is
#' not inferential at all. In fact, the papers states that if you have more than
#' a handful of candidate models to evaluate, then you are not performing an
#' inferential analysis. This leaves us with exploratory analysis. I would argue
#' that because the project includes an experimental design to measure drought
#' conditions, that the project is suited to predict conditions under changing
#' precipitation regimes. In this case, it may make more sense to keep an eye
#' towards predictive analysis.
#' 
#' Another issue identified in a climate-demography working group is
#' "flickering" of predictors, i.e., that adding additional years or data to a
#' dataset may cause certain variables to come into and out of best models. This
#' problem seems difficult to avoid without an extensive climate dataset, which
#' we certainly do not have. A best approach, then, may be to see how well
#' different climate variables do at predicting out-of-sample years. We can try
#' to do this using a cross-validation approach.
#' 
#' Here, I implement a "leave one year out" cross validation approach. This
#' involves splitting the dataset into training and testing sets; the testing
#' dataset can be all observations in one year, whereas the training dataset is
#' the data from all remaining years. Note that one shortcoming of this is that
#' temporal autocorrelation among years is thus ignored by this model (Roberts
#' et al. 2017, Ecography, makes this argument). We then can measure in some way
#' the performance of a model in predicting the response in the held-out year.
#' This procedure can be repeated for all of the years in the dataset and
#' averaged across held-out years to get a sense of which predictors, on
#' average, are best at predicting held-out years in the dataset. Tredennick et
#' al. (2021) use a similar approach on 20 years of data (i.e., the model is fit
#' 20 times, each time with only 19 years of data, and with each one of those 20
#' models, they evalutate how well it predicts the held-out year).
#' 
#' I will implement this approach below. I'll do this on the probability of
#' flowering. In the markdown file I won't include the code to set up the data,
#' although this code will be under the hood in the .R file. I'll use some
#' tidyverse packages (`dplyr`, `tidyr`, and `ggplot2`) and will fit models in
#' `lme4`.

#+ r setup, include = FALSE, message = FALSE, warning = FALSE, echo = FALSE

library(ggplot2)
library(dplyr)
library(tidyr)
library(lme4)

# Clear namespace
rm(list = ls())

# Load in demographic data
# (survival is imputed here)
demo = read.csv('01_data_cleaning/out/demo_imputed_survival.csv')

# Remove plants that weren't surveyed in 2020
demo = demo %>%
  separate(plantid, into = c('tag', 'plot', 'coord'), remove = FALSE, sep = '_') %>%
  filter(
    !( (Year %in% 2020 & plot %in% c(1:2, 13, 15) & grepl('\\d{2}', coord)) |
       (Year %in% 2020 & plot %in% 3 & (grepl('\\d{2}', coord) | grepl('[7890]', coord))))
  ) %>%
  select(-c(tag, plot, coord))

# Add a column for whether the plant flowered
demo = demo %>%
  mutate(
    flowering = case_when(
      No.umbels > 0 ~ TRUE,
      !No.umbels ~ FALSE,
      is.na(No.umbels) ~ FALSE,
      .default = NA
    )
  )

# Combine observation year's data with previous year's demographic data
demo.fl = merge(
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

# Read in a climate dataset
# (csv is has one row for each "season" of each year, with columns for
# aggregated temperature and precipitation measurements in those years)
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

# For each year, get the previous year's growing season as well (not explicitly
# in current data frame)
clim = rbind(
  clim,
  clim %>% filter(pd.label %in% 'grow') %>% mutate(Year = Year + 1, pd.label = 'pgrow')
)

# Convert data frame to wide form (one row per year) and select some
# pre-screened variables
clim.wide = clim %>%
  # the window start dates aren't necessary
  # also remove the last.freeze column (for now) because it makes things difficult
  select(-c(period.start, last.freeze)) %>%
  pivot_longer(c(contains('prec.sum'), contains("Temp")), names_to = "stat", values_to = "statVal") %>%
  pivot_wider(names_from = c(pd.label, stat), values_from = statVal) %>%
  merge(y = clim %>% filter(!pd.label %in% 'pgrow') %>% distinct(Year, last.freeze)) %>%
  select(
    Year, last.freeze, grow_meanTemp, early_prec.sum, early_minnTemp,
    winter_prec.sum, winter_meanTemp, summer_meanTemp, pgrow_meanTemp
  )

# Merge together to get one data frame
demo.fl.clim = merge(x = demo.fl, y = clim.wide, by = 'Year') %>%
  # Remove 2017 data (because 2016 size measurements are unreliable)
  filter(Year > 2017)

#' The data frame is assigned to the variable `demo.fl.clim`.

#+ r show data frame, echo = TRUE, warning = FALSE, message = FALSE

# Columns of data frame
names(demo.fl.clim)
# Climate variables are `last.freeze` through `pgrow_meanTemp`.
# `size.prev` is size in the previous year
# `flowering.prev` is whether the plant flowered in the previous year
# `flowering` is our response (boolean)
# `plantid` is a unique identifier for each unique plant

demo.fl.clim[1:5, c(1:5, 15:18)]

# Total number of observations
nrow(demo.fl.clim)

# Number of unique plants
length(unique(demo.fl.clim$plantid))

#' First, I'll fit a model without climate-associated variables, simply
#' absorbing year-to-year variation by including year as a random effect. I'll
#' use this to find a model structure for non-climate related terms to include
#' in the models looking for climate effects. I'll use AIC to evaluate these
#' models. *Note: I'm not sure if this is the best way to evaluate these models
#' against each other.* This step seems exploratory. In Tredennick et al.
#' (2021), they perform exploratory analysis by doing null hypothesis
#' significance testing and then correcting for multiple comparisons. That
#' approach may be valid as well.
#' 
#' I'll do this step on the full dataset because we're not trying to predict
#' data for missing years in this step, just to get additional variables that
#' need to be controlled for in models of flowering probability.

#+ r fit base models, echo = TRUE, warning = FALSE, message = FALSE

# Null model (no fixed effects)
f_0 = glmer(
  formula = flowering ~ (1 | Year) + (1 | Plot),
  family = 'binomial',
  data = demo.fl.clim
)

# Model with size in previous year
f_s = glmer(
  formula = flowering ~ size.prev + (1 | Year) + (1 | Plot),
  family = 'binomial',
  data = demo.fl.clim
)

AIC(f_0, f_s)
# Unsurprisingly, considerably strong effect of including previous year's size.

# Try fitting models with probability of flowering in the previous year.

# Without interaction
f_s.f = glmer(
  formula = flowering ~ size.prev + flowering.prev + (1 | Year) + (1 | Plot),
  family = 'binomial',
  data = demo.fl.clim
)

# With interaction
f_sf = glmer(
  formula = flowering ~ size.prev * flowering.prev + (1 | Year) + (1 | Plot),
  family = 'binomial',
  data = demo.fl.clim
)

AIC(f_s, f_s.f, f_sf)
# From a multi-model perspective, either of these models (with or without
# flowering time) are plausible. Because we'll be fitting a lot of models later,
# I'll adopt the approach of simply choosing the most parsimonious of plausible
# models (i.e., simplest model with delta AIC within 2 of the best performing
# model), in which case I'll proceed without the interaction.

# Look for effects of treatment, and possible interactions

f_s.f.t = glmer(
  formula = flowering ~ size.prev + flowering.prev + trt + (1 | Year) + (1 | Plot),
  family = 'binomial',
  data = demo.fl.clim
)

f_st.f = glmer(
  formula = flowering ~ size.prev * trt + flowering.prev + (1 | Year) + (1 | Plot),
  family = 'binomial',
  data = demo.fl.clim
)

f_s.ft = glmer(
  formula = flowering ~ size.prev + flowering.prev * trt + (1 | Year) + (1 | Plot),
  family = 'binomial',
  data = demo.fl.clim
)

f_st.ft = glmer(
  formula = flowering ~ size.prev * trt + flowering.prev * trt + (1 | Year) + (1 | Plot),
  family = 'binomial',
  data = demo.fl.clim
)

AIC(f_s.f, f_s.f.t, f_st.f, f_s.ft, f_st.ft) %>% mutate(delta.aic = AIC - min(AIC))
# It looks like the most parsimonious, simplest model has a treatment effect
# interacting with flowering time, but not with size.

#' Based on this procedure, we'll fit all of our climate-variable models with a
#' fixed effect of size and interactions between fixed effects of flowering and
#' treatment.
#' 
#' Now, for the leave-one-out step
#' 
#' I'll make six different splits of the dataset, one with each year held out.
#' I'll store these in a list where each element of hte list contains the
#' testing and training dataset.

#+ r split data, warning = FALSE, message = FALSE

demo.years = sort(unique(demo.fl.clim$Year))
demo.fl.split = vector(mode = 'list', length = length(demo.years))

for (i in 1:length(demo.years)) {
  
  demo.fl.split[[i]]$train = demo.fl.clim %>% filter(!(Year %in% demo.years[i]))
  demo.fl.split[[i]]$test  = demo.fl.clim %>% filter(Year %in% demo.years[i])
  demo.fl.split[[i]]$holdout = demo.years[i]
  
}

length(demo.fl.split)

#' Before fitting models, it's worthwhile to think about how to evaluate models
#' at predicting our out-of-sample years. There are multiple ways to do this,
#' and I haven't thought incredibly hard about the advantages of each.
#' 
#' Two sensible options are to use the "mean absolute deviation" (MAD), i.e.,
#' `abs(y.true - y.pred)` averaged over all data points. Another is to use what
#' is used in OLS for linear models, the mean squared deviation, i.e., `(y.true
#' - y.pred)^2` averaged over all data points. Because we're dealing with
#' bernoulli data, `y.true` will be 0 or 1; we'll assume `y.pred` is a predicted
#' probability.
#' 
#' Here is a plot of each of those. Assume that `y.true` is zero (i.e., the
#' plant did not flower). The x-axis is model predictions, `y.pred`, and the
#' y-axis is the error associated with the given `y.pred`.

#+ r show errors, echo = FALSE, warning = FALSE, message = FALSE

data.frame(y.pred = (1:99) / 100) %>%
  mutate(abs.error = abs(0-y.pred), sqd.error = (0-y.pred)^2) %>%
  pivot_longer(-y.pred, names_to = 'error.function', values_to = 'error') %>%
  ggplot(aes(x = y.pred, y = error, group = error.function)) +
  geom_line(aes(linetype = error.function)) +
  theme(legend.position = 'bottom')

#' If `y.true` is 1, then the plot is identical but with the x-axis reversed. A
#' `y.pred` of 0.5 is predicting that either outcome is equally likely; in this
#' case where `y.true` is zero, then values less than 0.5 are predicting that
#' not flowering is predicted to be more likely, and values greater than 0.5
#' predict that flowering is more likely than not. The squared error function
#' gives a lower penalty for predicting 0.5, and gives a harsher penalty for
#' being on the wrong side of 0.5. Maybe for this reason it's a good idea to use
#' the squared error function.
#' 
#' I'll use `sapply` to run each element of the list through a wrapper function
#' that fits a the model on the training set, then makes model predictions on
#' the testing set, and finally computes a mean error on those predictions.

#+ r model fitting wrappers, echo = TRUE, warning = FALSE

# Wrapper function for inverse logit
ilogit = function(x) 1 / (1 + exp(-x))

# Wrapper function for fitting a model on the split dataset
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
    mod.score = mean((dfs$test$flowering - mod.pred))^2
    return(mod.score)
  } else { return(NA) }
}

# Example, run on the first split (holding out year 2017)
fit.split.model(
  dfs = demo.fl.split[[1]], 
  mod.formula = 'flowering ~ size.prev + flowering.prev * trt + (1 | Year) + (1 | Plot)'
)

#' Now, we need to evaluate a set of candidate models to try out. This list will
#' be quite long! i'll try fitting models with flowering as a response
#' (obvoiusly), but in each case, with different climate variables included.
#' We'll try also looking at interactions with treatment, interactions with
#' size, and interactions with flowering in the previous year. Following
#' Tredennick et al. (2021), I'll also fit the model with none of these terms as
#' a null model for comparison - if the null model outperforms all
#' climate-related models, then we have all bogus climate variables.
#' 
#' I'm making the choice to leave the year-level random effects in. I'm not
#' totally sure if this is the best idea, but my reasoning here is that even
#' though much of the year-to-year variation is accounted for in the climate
#' term, there will still be correlations in residuals within a year. A drawback
#' to this is that this may throw a bunch of additional convergence errors, in
#' which case we will not be totally confident that our training models are in
#' fact fit with the maximum likelihood estimators.
#' 

forms = c(
  # Null model
  'flowering ~ size.prev + flowering.prev * trt + (1 | Year) + (1 | Plot)', 
  # Models with fixed effects of climate variables
  'flowering ~ size.prev + flowering.prev * trt + (1 | Year) + (1 | Plot) + last.freeze',
  'flowering ~ size.prev + flowering.prev * trt + (1 | Year) + (1 | Plot) + grow_meanTemp',
  'flowering ~ size.prev + flowering.prev * trt + (1 | Year) + (1 | Plot) + early_prec.sum',
  'flowering ~ size.prev + flowering.prev * trt + (1 | Year) + (1 | Plot) + early_minnTemp',
  'flowering ~ size.prev + flowering.prev * trt + (1 | Year) + (1 | Plot) + winter_prec.sum',
  'flowering ~ size.prev + flowering.prev * trt + (1 | Year) + (1 | Plot) + winter_meanTemp',
  'flowering ~ size.prev + flowering.prev * trt + (1 | Year) + (1 | Plot) + summer_meanTemp',
  'flowering ~ size.prev + flowering.prev * trt + (1 | Year) + (1 | Plot) + pgrow_meanTemp',
  # Models with size-by-climate interactions
  # (Null model will have size.prev slope varying by year)
  'flowering ~ flowering.prev * trt + (size.prev | Year) + (1 | Plot)', 
  'flowering ~ size.prev * last.freeze + flowering.prev * trt + (1 | Year) + (1 | Plot)',
  'flowering ~ size.prev * grow_meanTemp + flowering.prev * trt + (1 | Year) + (1 | Plot)',
  'flowering ~ size.prev * early_prec.sum + flowering.prev * trt + (1 | Year) + (1 | Plot)',
  'flowering ~ size.prev * early_minnTemp + flowering.prev * trt + (1 | Year) + (1 | Plot)',
  'flowering ~ size.prev * winter_prec.sum + flowering.prev * trt + (1 | Year) + (1 | Plot)',
  'flowering ~ size.prev * winter_meanTemp + flowering.prev * trt + (1 | Year) + (1 | Plot)',
  'flowering ~ size.prev * summer_meanTemp + flowering.prev * trt + (1 | Year) + (1 | Plot)',
  'flowering ~ size.prev * pgrow_meanTemp + flowering.prev * trt + (1 | Year) + (1 | Plot)',
  # Models with treatment-by-climate interactions
  # # (Null model will have trt slope varying by year)
  # 'flowering ~ size.prev + flowering.prev * trt + (trt | Year) + (1 | Plot)', 
  'flowering ~ size.prev + flowering.prev * trt + trt * last.freeze + (1 | Year) + (1 | Plot)',
  'flowering ~ size.prev + flowering.prev * trt + trt * grow_meanTemp + (1 | Year) + (1 | Plot)',
  'flowering ~ size.prev + flowering.prev * trt + trt * early_prec.sum + (1 | Year) + (1 | Plot)',
  'flowering ~ size.prev + flowering.prev * trt + trt * early_minnTemp +  (1 | Year) + (1 | Plot)',
  'flowering ~ size.prev + flowering.prev * trt + trt * winter_prec.sum + (1 | Year) + (1 | Plot)',
  'flowering ~ size.prev + flowering.prev * trt + trt * winter_meanTemp + (1 | Year) + (1 | Plot)',
  'flowering ~ size.prev + flowering.prev * trt + trt * summer_meanTemp + (1 | Year) + (1 | Plot)',
  'flowering ~ size.prev + flowering.prev * trt + trt * pgrow_meanTemp + (1 | Year) + (1 | Plot)',
  # Models with flowering-by-treatment-by-climate interactions
  # # (Null model will have trt*flowering slope varying by year)
  # 'flowering ~ size.prev + (flowering.prev * trt | Year) + (1 | Plot)', 
  'flowering ~ size.prev + flowering.prev * trt * last.freeze + (1 | Year) + (1 | Plot)',
  'flowering ~ size.prev + flowering.prev * trt * grow_meanTemp + (1 | Year) + (1 | Plot)',
  'flowering ~ size.prev + flowering.prev * trt * early_prec.sum + (1 | Year) + (1 | Plot)',
  'flowering ~ size.prev + flowering.prev * trt * early_minnTemp +  (1 | Year) + (1 | Plot)',
  'flowering ~ size.prev + flowering.prev * trt * winter_prec.sum + (1 | Year) + (1 | Plot)',
  'flowering ~ size.prev + flowering.prev * trt * winter_meanTemp + (1 | Year) + (1 | Plot)',
  'flowering ~ size.prev + flowering.prev * trt * summer_meanTemp + (1 | Year) + (1 | Plot)',
  'flowering ~ size.prev + flowering.prev * trt * pgrow_meanTemp + (1 | Year) + (1 | Plot)'
)

# Initialize an object to store scores in
# rows correspond to models, cols correspond to hold-out years
form.scores = matrix(NA, nrow = length(forms), ncol = length(demo.years))

# Note: this loop takes a while because it fits ~200 models
for (i in 1:length(forms)) {
  form.scores[i,] = sapply(
    X = demo.fl.split, 
    FUN = fit.split.model,
    mod.formula = forms[i]
  )
}

# Distribution of the model scores
hist(rowMeans(form.scores))

# Look at the best performing models
data.frame(
  form  = forms,
  score = round(rowMeans(form.scores), 4)
) %>%
  arrange(score) %>%
  head(8)

#' From this, it appears that all of the best performing models all have mean
#' summer temperature (from the prior year) as a predictor. The differences
#' among models are simply due to interaction terms. The next-best set of
#' performing models all feature growing season temperature in the previous
#' year. Note that this near arrangement where the best performing models
#' feature the same sets of predictors is *not* guaranteed to happen. But, when
#' it does happen, it suggests that there may be predictive power in the
#' predictor.
#' 
#' To get coefficients, let's fit the models with summer temperature on the
#' whole dataset. We'll fit all of the mdoels for assessment.
#'

f_s.fts = glmer(
  formula = flowering ~ size.prev + flowering.prev * trt * summer_meanTemp + 
    (1 | Year) + (1 | Plot),
  family = binomial,
  data = demo.fl.clim
)

f_s.ft.ts = glmer(
  formula = flowering ~ size.prev + flowering.prev * trt + trt * summer_meanTemp + 
    (1 | Year) + (1 | Plot),
  family = binomial,
  data = demo.fl.clim
)

f_ss.ft = glmer(
  formula = flowering ~ size.prev * summer_meanTemp + flowering.prev * trt + 
    (1 | Year) + (1 | Plot),
  family = binomial,
  data = demo.fl.clim
)

f_s.ft.s = glmer(
  formula = flowering ~ size.prev + summer_meanTemp + flowering.prev * trt + 
    (1 | Year) + (1 | Plot),
  family = binomial,
  data = demo.fl.clim
)

# Look at AIC just to see how these models compare to each other
AIC(f_s.fts, f_s.ft.ts, f_ss.ft, f_s.ft.s) %>%
  mutate(delta.aic = AIC - min(AIC))
# The ordering of these models is roughly the same in as our procedure!

# Examine the coefficients in the "best" performing model
summary(f_s.fts)
# From these, it appears that:
# - Positive but n.s. effect of prior flowering in the controls (as in other models)
# - Positive but n.s. effect of summer temperature in previous year
# - For plants that flowered previously, the treatments have *positive* 
#   effects on flowering
# - Increasing summer temperature also has extra-strong effects on drought 
#   plants that did not flower in the prior year
# - For flowering plants, there may be a negative effect in irrigated plots for plants 
#   that flowered in the previous year (compared to those that did not)

#' Let's plot some of these model fits against the raw data to see how they differ from each other.
#' 

#+ r, plot model fits, echo = FALSE

summ.temp.predictions = expand.grid(
  size.prev = c(5:50)/10,
  flowering.prev = c(TRUE, FALSE),
  Year = 2017:2023,
  trt = c('control', 'drought', 'irrigated')
) %>%
  merge(y = demo.fl.clim %>% distinct(Year, summer_meanTemp)) %>%
  mutate(
    pred.s.fts = predict(f_s.fts, newdata = ., re.form = NA) %>% ilogit(),
    pred.s.ft.ts = predict(f_s.ft.ts, newdata = ., re.form = NA) %>% ilogit(),
    pred.s.ss.ft = predict(f_ss.ft, newdata = ., re.form = NA) %>% ilogit(),
    pred.s.ft.s = predict(f_s.ft.s, newdata = ., re.form = NA) %>% ilogit()
  ) %>%
  pivot_longer(contains('pred.'), names_to = 'model', values_to = 'prediction') %>%
  mutate(model = gsub('pred\\.', '', model))

demo.fl.clim %>%
  mutate(
    flowering = as.numeric(flowering),
    flowering.prev = ifelse(flowering.prev, 'veg.prior', 'flw.prior')
  ) %>%
  ggplot(aes(x = size.prev, colour = trt)) +
  geom_point(
    aes(y = flowering),
    size = 3, alpha = 0.5,
    position = position_jitter(height = 0.1)
  ) +
  geom_line(
    data = summ.temp.predictions %>%
      mutate(flowering.prev = ifelse(flowering.prev, 'veg.prior', 'flw.prior')),
    aes(y = prediction, group = model, linetype = model)
  ) +
  scale_colour_manual(values = c('black', 'red', 'blue')) +
  labs(x = 'Size in prev. year', y = 'Flowering probability') +
  facet_wrap(~ paste(Year, trt, flowering.prev)) +
  theme(legend.position = 'bottom')

#' Perhaps unsurprisingly (given how close all of the scores are), these models
#' appear to differ very little in their predictions!
#' 
#' It does look, though, like the differences mostly arise in areas where there
#' is no data!
