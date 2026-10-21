# Eyeballing germination rates from data recorded for Jenna's 2021-2022 seed
# addition study
# These data are probably not useful for formal vital rate estimation or
# projections (esp. considering they are for only one year) but they will be
# useful for justifying whatever value we may choose
# ---
library(ggplot2)
library(dplyr)
library(tidyr)
library(lme4)

# Clear namespace
rm(list = ls())

# Read in seedling (s-ling, or sling) data
sling = read.csv('00_raw_data/2022_seed_addition_counts.csv')

# Examine
head(sling)
nrow(sling)
table(sling$plot)
table(sling$quadrant)
table(sling$no.lomatium)
# Fully factorial - nice

# Are there NAs?
sling %>% filter(is.na(no.lomatium))
# hmm... yes, all in plot 7

# Remove these observations (no mo factorial)
sling = sling %>% filter(!is.na(no.lomatium))

# Merge in treatment data
sling = merge(sling, read.csv('00_raw_data/plot_treatments.csv'))

head(sling)

# Add treatment variables
sling = sling %>%
  mutate(loma.add = quadrant %in% c(3, 7)) %>%
  mutate(cage = quadrant %in% 5:8)

# First: check for overdispersion:

# Standard poisson model
s_0 = glmer(
  no.lomatium ~ (1 | plot),
  data = sling,
  family = 'poisson',
  control = glmerControl(optimizer = 'bobyqa')
)

s_n = glmer.nb(
  no.lomatium ~ (1 | plot),
  data = sling,
  control = glmerControl(optimizer = 'bobyqa')
)

AIC(s_0, s_n)
# Yes - much more support for the negative binomial model

summary(s_n)

# Maybe start first with cage/uncage
s_c = glmer.nb(
  no.lomatium ~ cage + (1 | plot),
  data = sling,
  control = glmerControl(optimizer = 'bobyqa')
)

AIC(s_n, s_c)
# No support for cage effects here

# Test for effects of addition
s_a = glmer.nb(
  no.lomatium ~ loma.add + (1 | plot),
  data = sling,
  control = glmerControl(optimizer = 'bobyqa')
)

AIC(s_a, s_n)
# Yes - there is a clear effect of addition (good)

summary(s_a)
# Addition increases the number of Lomatium in a plot by >20-fold
# (no shock here)

# Test for an interaction with the cage, just in case:
s_ac = glmer.nb(
  no.lomatium ~ loma.add * cage + (1 | plot),
  data = sling,
  control = glmerControl(optimizer = 'bobyqa')
)

AIC(s_a, s_ac, s_n)
# No evidence of a cage effect here


# Look for a treatment effect (takes a few seconds to run)
s_a_t = glmer.nb(
  no.lomatium ~ loma.add + trt + (1 | plot),
  data = sling,
  control = glmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 1e7), tolPwrss = 1e-6)
)

AIC(s_a, s_a_t)
# Great - no evidence of a treatment effect either

summary(s_a_t)

# AIC(s_t, s_n)
# anova(s_t, s_n) # models not nested
# Well, no support for treatment effect (although hard to fully rule out with conv. issue)


fits.errs = predict(
  s_a, 
  newdata = expand.grid(loma.add = c(FALSE, TRUE)), 
  re.form = ~ 0, se.fit = TRUE
) %>%
  do.call(what = 'cbind')

# (warning here is just to tell us that the SE is approximated)
fits.errs = data.frame(loma.add = c(FALSE, TRUE)) %>%
  cbind(fits.errs) %>%
  mutate(
    fit.mean = exp(fit),
    fit.cilo = exp(fit - 2*se.fit),
    fit.cihi = exp(fit + 2*se.fit)
  )

fits.errs %>%
  select(starts_with('fit.')) %>%
  summarise(across(everything(), ~ diff(.) / 30))
# Yeah. Lots of support for 0.01

# Some envelope math:
# the effect size for addition is exp(3.049) =  approx 21
# so there's a 21x increase in seedlings when 30 seeds are added
# let x be ambient seed rain (in expectation), y be number of seedlings, p be germination prob.
# y = px
# 21y = p(30+x)
# combine these to get .45y = 30p
# but also LHS is .45px = 30p
# so x = approx 30/.35 or x between 66 and 67
# y we have estimated from model outputs, and p = y/x

s_a_preds[-1,c("fit.cilo", "fit.mean", "fit.cihi")] / (66)
# (first row removed - this is from the addition+cage treatment where germ was higher)
# Hmm okay
# these numbers range from (roughlY) between .0005 and .005
# .0001 is in the middle of that range
# (.01 is way too high, .005 looks like it's pushing it)

