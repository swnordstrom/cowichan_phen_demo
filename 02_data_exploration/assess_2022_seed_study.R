# Eyeballing germination rates from data recorded for Jenna's 2021-2022 seed
# addition study
# These data are probably not useful for formal vital rate estimation or
# projections (esp. considering they are for only one year) but they will be
# useful for justifying whatever value we may choose
# ---

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

# Remove these observations
sling = sling %>% filter(!is.na(no.lomatium))

# Merge in treatment data
sling = merge(sling, read.csv('00_raw_data/plot_treatments.csv'))

head(sling)

# Let's try this.
sling = sling %>%
  mutate(loma.add = quadrant %in% c(3, 7)) %>%
  mutate(loma.cage = quadrant %in% 7) %>%
  mutate(loma.nocage = quadrant %in% 3) %>%
  mutate(cage = quadrant %in% (1:4)*2)

s_0 = glmer(
  no.lomatium ~ (1 | plot),
  data = sling,
  family = 'poisson'
)

summary(s_0)

s_t = glmer(
  no.lomatium ~ trt + (1 | plot),
  data = sling,
  family = 'poisson'
)

AIC(s_t, s_0)
anova(s_t, s_0)
# okay sweet, no solid treatment effect here

# Look at addition effects

s_a = glmer(
  no.lomatium ~ loma.add + (1 | plot),
  data = sling,
  family = 'poisson'
)

AIC(s_a, s_0)
anova(s_a, s_0)
# ah yes, there is an addition effect...
summary(s_a)
predict(s_a, data.frame(loma.add = c(TRUE, FALSE)), re.form = ~ 0, type = 'response')
# so ~3x higher in the addition plots than not...

s_c = glmer(
  no.lomatium ~ loma.add + cage + (1 | plot),
  data = sling,
  family = 'poisson'
)
# huh
anova(s_c, s_a)

s_c = glmer(
  no.lomatium ~ loma.add * cage + (1 | plot),
  data = sling,
  family = 'poisson'
)
anova(s_c, s_a)
AIC(s_c, s_a)
# oh... huh

s_cnc = glmer(
  no.lomatium ~ loma.cage + loma.nocage + (1 | plot),
  data = sling,
  family = 'poisson'
)

AIC(s_cnc, s_c)
# okay it's better to have a cage * addition interaction

summary(s_c)

# expand.grid(loma.add = c(TRUE, FALSE), cage = c(TRUE, FALSE)) %>%
#   mutate(pp = predict(s_c, newdata = ., re.form = ~ 0, type = 'response', se.fit = TRUE))
# above throws error


s_c_preds = expand.grid(loma.add = c(TRUE, FALSE), cage = c(TRUE, FALSE)) %>%
  cbind(
    predict(
      s_c, re.form = ~ 0, se.fit = TRUE,
      newdata = expand.grid(loma.add = c(TRUE, FALSE), cage = c(TRUE, FALSE))
    ) %>%
  sapply(cbind)
) %>%
  mutate(
    fit.mean = exp(fit),
    fit.cihi = exp(fit + se.fit),
    fit.cilo = exp(fit - se.fit)
  )

s_c_preds

s_c_preds %>%
  mutate(
    label = case_when(
      loma.add & cage ~ 'added with cage',
      loma.add & !cage ~ 'added no cage',
      !loma.add & cage ~ 'not added with cage',
      !loma.add & !cage ~ 'not added no cage'
    )
  ) %>%
  ggplot(aes(x = label)) +
  geom_segment(
    aes(xend = label, y = fit.cilo, yend = fit.cihi), colour = 'gray'
  ) +
  geom_point(aes(y = fit.mean), size = 4) +
  labs(x = 'treeatment', y = 'mean number of seeedlings')

# Hmm standard error estimates here are quite large.
# I suppose those are due in part to the standard effects?

# Cage seems to have a large effect! Possibly larger than the addition.
# But this effect is really only felt in the addition... i.e., interaction term

# Some envelope math:
# the effect size for addition+cage is ~ (-.51 - .61 + 1.5) approx .3725
# exp(.3725) = approx 1.45
# so there's a 1.45x increase in seedlings when 30 seeds are added (with frugifore exclosure)
# let x be ambient seed rain (in expectation), y be number of seedlings, p be germination prob.
# y = px
# 1.45y = p(30+x)
# combine these to get .45y = 30p
# but also LHS is .45px = 30p
# so x = approx 30/.35 or x between 66 and 67
# y we have estimated from model outputs, and p = y/x

s_c_preds[-1,c("fit.cilo", "fit.mean", "fit.cihi")] / (66)
# (first row removed - this is from the addition+cage treatment where germ was higher)
# Hmm okay
# these numbers range from (roughlY) between .0005 and .005
# .0001 is in the middle of that range
# (.01 is way too high, .005 looks like it's pushing it)

