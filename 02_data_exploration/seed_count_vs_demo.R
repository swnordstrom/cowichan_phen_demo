
##### Setup

# Loading packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(lme4)
library(glmmTMB)

# This message is super annoying and I hate it
options(dplyr.summarise.inform = FALSE)

# Clearing namespace
rm(list = ls())

# Load in raw data
demo.seed = read.csv('01_data_cleaning/out/demo_seed_v1.csv')

head(demo.seed)

# Load in plot info
plot.info = read.csv('00_raw_data/plot_treatments.csv')

# Load in an invlogit wrapper
ilogit = function(x) (1+exp(-x))^(-1)

##### Do cleaning as needed

### Main concern is plantids (for getting individual plants...)

demo.seed %>% 
  distinct(plantid.seed, plantid.demo) %>% 
  group_by(ids.match = plantid.seed == plantid.demo) %>%
  summarise(n = n())
# yes... matching is unbelievably poor

# Actually this is probably not a problem... I trust the demoIDs
# (demoIDs should be identical across years)

### Other thing is what to do about plants that were listed as not flowering...
demo.seed %>%
  filter(!flowering.in.demo)
# ah good - it's a small number
# okay - most of these are dead, eaten, etc. umbels
# I say ignore them

demo.seed = demo.seed %>% filter(flowering.in.demo)
nrow(demo.seed)

### Guess I should also check for NAs

apply(demo.seed, 2, function(x) sum(is.na(x)))
# This is good actually! Number of leaves and leaf length numbers are small!
# Can exclude these maybe on a case-by-case basis for now.

# While here let's also merge in treatment data
demo.seed = merge(demo.seed, plot.info)

##### Assessments

### How often to plants appear in multiple years?

demo.seed %>%
  group_by(plantid.demo) %>%
  summarise(n.years = length(unique(year))) %>%
  group_by(n.years) %>%
  summarise(nn = n())
# okay - some plants appearing multiple times
# (would be good to look at variation in these - are they iid?)

##### Some visualizations, first

### Number of plants per plot per year

demo.seed %>%
  distinct(year, plot, trt, plantid.demo) %>%
  group_by(year, plot, trt) %>%
  summarise(n = n()) %>%
  ggplot(aes(x = plot, y = n, colour = trt)) +
  geom_point() +
  facet_wrap(~ year)

# Ah right - sampling effort was low in 2021

### Number of umbels counted per plant, per plot?

# We probably only want this for 2022-2023
demo.seed %>%
  filter(year > 2021) %>%
  group_by(year = factor(year), plot, plantid.demo, trt) %>%
  summarise(n = n()) %>%
  ggplot(aes(x = plot, y = n, colour = trt)) +
  geom_point(position = position_jitter(width = 0.25, height = 0.25), alpha = 0.5) +
  facet_wrap(~ year, nrow = 1)

# Seems like we have multiple umbels for the majority of plants in here (good)

### Okay... time to look at the distribution of seed counts
# (there will be tons of zeros due to death or florivory)

demo.seed %>%
  mutate(year = factor(year)) %>%
  ggplot(aes(x = plot, y = no.seeds, colour = trt)) +
  geom_point(position = position_jitter(width = 0.25), alpha = 0.25) +
  facet_wrap(~ year)

# Hmm... perhaps not as many zeroes as expected...

# Proportion of zeroes per plot per year
zeros.per.year = demo.seed %>%
  mutate(year = factor(year)) %>%
  group_by(year, plot, trt) %>%
  summarise(p.zero = mean(!no.seeds, na.rm = TRUE)) 
# (wait... why are there NAs in the no.seeds column...?)

zeros.per.year %>%
  ggplot(aes(x = plot, y = p.zero)) +
  geom_point(aes(colour = year))

# lmao this is kind of a silly plot

zeros.per.year %>%
  ggplot(aes(x = year, y = p.zero, colour = trt)) +
  geom_point() +
  geom_line(aes(group = plot)) +
  scale_y_continuous(limits = 0:1)
# actually not seeing very much consistency in treatments
# (treatment x year does look like there's something happening though)
# there is some clear variation across years
# (and it's quite high in 2022-2023)

# *although* it's possible that the sampling in 2021 excluded plants with dead or eaten umbels...
# (phen data could answer this, perhaps?)

### Distribution of zeros within plants
# Does the number of failed umbels *on a plant* vary by treatment? Or number of umbels produced?

zeros.per.plant = demo.seed %>%
  group_by(
    year = factor(year), 
    plot, 
    trt,
    plantid = plantid.demo
  ) %>%
  summarise(p.zero = mean(!no.seeds, na.rm = TRUE), n = n()) 

head(zeros.per.plant)

zeros.per.plant %>%
  ggplot(aes(x = plot, y = p.zero)) +
  geom_point(
    aes(size = n, colour = trt),
    position = position_jitter(width = 0.25),
    alpha = 0.5
  ) +
  facet_wrap(~ year)
# hmm... not seeing a ton of variation across years I guess
# but it does look like most plants have at least one failed umbel!

zeros.per.plant %>%
  ggplot(aes(x = n, y = p.zero)) +
  geom_point(
    aes(colour = trt),
    position = position_jitter(width = 0.25),
    alpha = 0.5
  ) +
  facet_wrap(trt ~ year)
# cool... no clear pattern is visible...

### Leaf characters vs. seed count...
# idk, probably want umble count in here as a mediator
# (could run models that control for this... seeds per umbel?)

##### Try some models?

# Okay lol I suppose it's time to remove the NAs...
demo.seed %>% filter(is.na(no.seeds))
demo.seed = demo.seed %>% filter(!is.na(no.seeds))

### Okay... first I guess is probability of getting zeros (because that is
### important! will want a zero-inflated model)

p.zero.0 = glmer(
  formula = is.zero ~ (1 | plot / plantid),
  data = demo.seed %>% mutate(is.zero = !no.seeds, isnt.zero = no.seeds > 0, plantid = plantid.demo),
  family = 'binomial'
)
# singularity... rats
# standard deviation of zero for plant within plot are all zero
# mis-specification somehow but I'm not sure how?
# (oh and the plot-level random effects sure do seem to show a treatment effect so far!)

p.zero.0 = glmer(
  formula = cbind(is.zero, isnt.zero) ~ (1 | plot),
  data = demo.seed %>% 
    group_by(plantid.demo, plot, year) %>%
    summarise(is.zero = sum(!no.seeds), isnt.zero = sum(no.seeds > 0)),
  family = 'binomial'
)

p.zero.0
ranef(p.zero.0)
# cool
# (but there should be a random effect of plantid due to repeated sampling...
# grrr...)
# (trying the above with a random effect of plantid though gives us the
# singularity again though lmao)

p.zero.y = glmer(
  formula = cbind(is.zero, isnt.zero) ~ year + (1 | plot),
  data = demo.seed %>% 
    group_by(plantid.demo, plot, year = factor(year)) %>%
    summarise(is.zero = sum(!no.seeds), isnt.zero = sum(no.seeds > 0)),
  family = 'binomial'
)

p.zero.y
summary(p.zero.y)
# (I forget if z-values are appropriate here for significance testing though lmao)

AIC(p.zero.y, p.zero.0)
# yes - year effects model is better!

# Okay let's get spicier and add in treatment effects!

p.zero.y.t = glmer(
  formula = cbind(is.zero, isnt.zero) ~ year + trt + (1 | plot),
  data = demo.seed %>% 
    group_by(plantid.demo, plot, year = factor(year), trt) %>%
    summarise(is.zero = sum(!no.seeds), isnt.zero = sum(no.seeds > 0)),
  family = 'binomial'
)

summary(p.zero.y.t)
# wow... no treatment effect?

AIC(p.zero.y.t, p.zero.y)
# stick with the year model!

# but is there a year x trt effect? hopefully not... check just in case

p.zero.yt = glmer(
  formula = cbind(is.zero, isnt.zero) ~ year*trt + (1 | plot),
  data = demo.seed %>% 
    group_by(plantid.demo, plot, year = factor(year), trt) %>%
    summarise(is.zero = sum(!no.seeds), isnt.zero = sum(no.seeds > 0)),
  family = 'binomial'
) # convergence issues...

summary(p.zero.yt) # definitely getting increased s.e. of coefs
AIC(p.zero.yt, p.zero.y.t)
# 

### Okay... let's see if the non-zero plants are Poisson?

non.zeros = demo.seed %>% 
  filter(no.seeds > 0) %>%
  mutate(year = factor(year))

nrow(non.zeros)

nz.seedcount.0 = glmer(
  formula = no.seeds ~ (1 | plot / plantid.demo),
  data = non.zeros,
  family = 'poisson'
)

nz.seedcount.0 # considerable within-plant variance!
# definitely looks like overdispersion though... ugh

nz.seedcount.0 = glmer.nb(
  formula = no.seeds ~ (1 | plot / plantid.demo),
  data = non.zeros
)
# ahhhhhh fuck
summary(nz.seedcount.0)
# the within-plant variance is suddenly super small...
# come on man! shit.

nz.seedcount.0 = glmer.nb(
  formula = no.seeds ~ (1 | plot),
  data = non.zeros
)
# no convergence issues here though
summary(nz.seedcount.0)
# hmm...
ranef(nz.seedcount.0)

nz.seedcount.y = glmer.nb(
  formula = no.seeds ~ year + (1 | plot),
  data = non.zeros
)

# cool
summary(nz.seedcount.y)
# cool cool
ranef(nz.seedcount.y)
# treatment effects look less certain here...

AIC(nz.seedcount.y, nz.seedcount.0)
# yep - year effects are definitely needed here

non.zeros = non.zeros %>%
  mutate(pred.y = exp(predict(nz.seedcount.y, newdata = ., re.form = NA)))

# Visualizing this model (for now - absent treatment effects!)
non.zeros %>%
  ggplot(aes(x = year)) +
  geom_point(aes(y = no.seeds), position = position_jitter(width = 0.25), alpha = 0.5) +
  geom_point(aes(y = pred.y), shape = 'square', size = 3, colour = 'blue')
# cool

non.zeros %>%
  ggplot(aes(x = year)) +
  geom_point(aes(y = no.seeds, colour = trt), position = position_jitter(width = 0.25), alpha = 0.5) +
  geom_point(aes(y = pred.y), shape = 'square', size = 3, colour = 'blue')
# treatment effects not super obvious here... but this is due to plot construction

# Model with treatment effect
nz.seedcount.y.t = glmer.nb(
  formula = no.seeds ~ year + trt + (1 | plot),
  data = non.zeros
)

summary(nz.seedcount.y.t)
# also no strong treatment effects here...
ranef(nz.seedcount.y.t)

# Add model predictions
non.zeros = non.zeros %>%
  mutate(pred.y.t = exp(predict(nz.seedcount.y.t, newdata = ., re.form = NA)))

head(non.zeros)

AIC(nz.seedcount.y.t, nz.seedcount.y)
# no effect of treatment here either!

# What about a model with interaction of year and treatment? just to make sure...

nz.seedcount.yt = glmer.nb(
  formula = no.seeds ~ year * trt + (1 | plot),
  data = non.zeros
)

AIC(nz.seedcount.yt, nz.seedcount.y.t, nz.seedcount.y)
# so again, only year effects!

### Try a zero-inflated model now.
# (should also try a hurdle... could see either being better though...)

# Let's just add factors to the dataset
demo.seed = demo.seed %>% mutate(year = factor(year))

# (Hopefully my glmmTMB version is cool? hmm...)

seedcount.0.zip = glmmTMB(
  formula = no.seeds ~ (1 | plot / plantid.demo),
  family = 'nbinom2',
  ziformula = ~ 1,
  data = demo.seed
) # so fast!

seedcount.0.zip

seedcount.y.zip = glmmTMB(
  formula = no.seeds ~ year + (1 | plot / plantid.demo),
  family = 'nbinom2',
  ziformula = ~ 1,
  data = demo.seed
)

AIC(seedcount.y.zip, seedcount.0.zip)
# seems like year is good for seed counts
# also check to see if year should be put into the zero-inflation...

seedcount.0.zip.0 = glmmTMB(
  formula = no.seeds ~ (1 | plot / plantid.demo),
  family = 'nbinom2',
  ziformula = ~ .,
  data = demo.seed
)

seedcount.y.zip.y = glmmTMB(
  formula = no.seeds ~ year + (1 | plot / plantid.demo),
  family = 'nbinom2',
  ziformula = ~ .,
  data = demo.seed
)

AIC(seedcount.y.zip.y, seedcount.y.zip, seedcount.0.zip.0)
# Good to include random effect structure in both portions
# (let's just plan on doing this throughout)

seedcount.y.t.zip.y.t = glmmTMB(
  formula = no.seeds ~ year + trt + (1 | plot / plantid.demo),
  family = 'nbinom2',
  ziformula = ~ .,
  data = demo.seed
)

AIC(seedcount.y.t.zip.y.t, seedcount.y.zip.y) # no improvement!

# If removing trt from the seed count, but keeping it in zero-inflation...

seedcount.y.zip.y.t = glmmTMB(
  formula = no.seeds ~ year + (1 | plot / plantid.demo),
  family = 'nbinom2',
  ziformula = ~ year + trt + (1 | plot / plantid.demo),
  data = demo.seed
)

AIC(seedcount.y.t.zip.y.t, seedcount.y.zip.y, seedcount.y.zip.y.t)
# hmm... seems slightly better!

seedcount.y.t.zip.y = glmmTMB(
  formula = no.seeds ~ year + trt + (1 | plot / plantid.demo),
  family = 'nbinom2',
  ziformula = ~ year + (1 | plot / plantid.demo),
  data = demo.seed
)
# what will happen...
AIC(seedcount.y.t.zip.y.t, seedcount.y.zip.y, 
    seedcount.y.zip.y.t, seedcount.y.t.zip.y)
# cool - keep treatment in the zero-inflation!

# Okay... try some interactions now
seedcount.y.zip.yt = glmmTMB(
  formula = no.seeds ~ year + (1 | plot / plantid.demo),
  family = 'nbinom2',
  ziformula = ~ year * trt + (1 | plot / plantid.demo),
  data = demo.seed
)

AIC(seedcount.y.zip.y, seedcount.y.zip.yt)
# no treatment effects anywhere!

# interactions?
seedcount.yt.zip.yt = glmmTMB(
  formula = no.seeds ~ year * trt + (1 | plot / plantid.demo),
  family = 'nbinom2',
  ziformula = ~ year * trt + (1 | plot / plantid.demo),
  data = demo.seed
)

AIC(seedcount.y.zip.y, seedcount.y.zip.yt, seedcount.yt.zip.yt)
# leave treatment out of seed set altogether

# Also do a check of hurdle model
seedcount.y.hdl.y = glmmTMB(
  formula = no.seeds ~ year + (1 | plot / plantid.demo),
  family = 'truncated_nbinom2',
  ziformula = ~ .,
  data = demo.seed
)

AIC(seedcount.y.hdl.y, seedcount.y.zip.y) # about the same... let's go with the hurdlews

# Look at random effects of the best model!
ranef(seedcount.y.zip.y)
# man these are tiny effects of individual for zero-inflation

### Evaluate model predictions

demo.seed = demo.seed %>%
  mutate(pred.y_y = exp(predict(seedcount.y.zip.y, newdata = ., re.form = NA)))

# To assess: seed counts ~ year, different lines for each treatment?

demo.seed %>%
  ggplot(aes(x = year)) +
  geom_point(
    aes(y = no.seeds, colour = trt),
    position = position_jitter(width = 0.25),
    alpha = 1/3
  ) +
  geom_point(aes(y = pred.y_y), shape = 0, size = 4) +
  geom_line(aes(y = pred.y_y, group = trt)) +
  labs(y = 'seeds/umbel') +
  facet_wrap(~ trt)
# interesting!

ggsave('02_data_exploration/figs/seeds_per_umbel.png', width = 8, height = 5)

