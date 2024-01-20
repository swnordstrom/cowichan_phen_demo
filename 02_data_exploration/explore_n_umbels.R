# Init a script for looking at the number of umbels produced per plant

library(ggplot2)
library(dplyr)
library(tidyr)
library(glmmTMB)

rm(list = ls())

# Read in demo data
# also - add treatment data and clip the coordinates out of the plant ID
demo = read.csv('01_data_cleaning/out/demo_imputed_survival.csv') %>%
  merge(read.csv('00_raw_data/plot_treatments.csv') %>% rename(Plot = plot)) %>%
  separate(plantid, into = c('tag', 'plot', 'coord'), sep = '_', remove = FALSE) %>%
  select(-coord) %>%
  unite(tagplot, c(tag, plot))

head(demo)
length(unique(demo$plantid))
length(unique(demo$tagplot)) # hmm... okay maybe tagplot is not necessary

demo %>% 
  group_by(tagplot) %>% 
  filter(any(duplicated(Year))) %>% 
  arrange(tagplot, Year) #%>% View()
# that's it... maybe run this analysis just using plantid

demo.w.prev = merge(
  # This year's flowering info
  x = demo %>% select(-c(Stalk_Height, tagplot, obs.alive)),
  y = demo %>%
    # Align years
    mutate(Year = Year + 1) %>%
    select(Year, plantid, No.leaves, Leaf.length, No.umbels, demo.note, proc.note),
  by = c('plantid', 'Year'), suffixes = c('', '.prev')
)

head(demo.w.prev)

demo.w.prev = demo.w.prev %>%
  # Give me only the living plants with size measurements
  filter(!is.na(No.leaves), !is.na(Leaf.length), No.leaves > 0) %>%
  # Also give me only the plants where we have size measurements in the prior year
  filter(!is.na(No.leaves.prev), !is.na(Leaf.length.prev), No.leaves.prev > 0) %>%
  # Add previous size
  mutate(Size.prev = log(No.leaves.prev * Leaf.length.prev)) %>%
  # Change NA umbel count to zero
  mutate(No.umbels = ifelse(is.na(No.umbels), 0, No.umbels))

head(demo.w.prev)

# demo.w.prev %>% filter(is.na(No.umbels) | is.na(Size.prev))

demo.w.prev %>%
  # Forgot that 2016 sizes are bunk
  filter(Year > 2017) %>%
  ggplot(aes(x = Size.prev, y = No.umbels, colour = trt)) +
  geom_point(position = position_jitter(height = 0.25), size = 3) +
  scale_colour_manual(values = c('black', 'red', 'blue')) +
  facet_wrap(~ Year) +
  theme(legend.position = 'bottom')
# Sure looks lto me like there are more umbels on larger plants...

# Get a dataset with just flowering plants
flwr.w.prev = demo.w.prev %>% filter(No.umbels > 0, Year > 2017)

flwr.w.prev %>%
  ggplot(aes(x = Size.prev, y = No.umbels, colour = trt)) +
  geom_point(position = position_jitter(height = 0.25), size = 3) +
  scale_colour_manual(values = c('black', 'red', 'blue')) +
  scale_y_log10() +
  facet_wrap(~ Year) +
  theme(legend.position = 'bottom')

# Look at variance

flwr.w.prev %>%
  group_by(No.umbels) %>%
  summarise(n = n()) %>%
  mutate(p = n / sum(n)) %>%
  mutate(mean = sum(No.umbels * p)) %>%
  # wolfram alpha tells me this gives lambda approx 1.51
  mutate(
    lambda = 1.51,
    expected = lambda^No.umbels / ((exp(lambda) - 1) * factorial(No.umbels))
  ) %>%
  ggplot(aes(x = No.umbels)) +
  geom_point(aes(y = p), colour = 'blue') +
  geom_point(aes(y = expected), colour = 'red')
# not that far off...

var(flwr.w.prev$No.umbels) # observed variance ~ 1.472
with(l <- 1.51, ((l + l^2) / (1 - exp(-l))) - (l^2 / (1-exp(-l))^2)) # expected variance is ~1.10
# actually this is fine for a poisson I think

##### Run models

u_0 = glmmTMB(
  No.umbels ~ (1 | Plot / plantid) + (1 | Year),
  data = flwr.w.prev,
  ziformula = ~ 0,
  family = truncated_poisson()
)

summary(u_0)
# diagnose(u_0)

u_s = glmmTMB(
  No.umbels ~ Size.prev + (1 | Plot / plantid) + (1 | Year),
  data = flwr.w.prev,
  family = truncated_poisson()
)

summary(u_s) # significant effect
# diagnose(u_s)

u_s_t = glmmTMB(
  No.umbels ~ Size.prev + trt + (1 | Plot / plantid) + (1 | Year),
  data = flwr.w.prev,
  family = truncated_poisson()
)

summary(u_s_t)
# hmm...

anova(u_s_t, u_s, u_0)
# so... no treatment effects here, once again, lmao
# I also tried an interaction - also n.s.

# unless we look at treatment effects varying by year...
u_s_ty = glmmTMB(
  No.umbels ~ Size.prev + (1 | Plot / plantid) + (trt | Year),
  data = flwr.w.prev,
  family = truncated_poisson()
)
# convergence issue, unsurprisingly
rm(u_s_ty)

# Size relationship varying by year

u_sy = glmmTMB(
  No.umbels ~ (1 | Plot / plantid) + (Size.prev | Year),
  data = flwr.w.prev,
  family = truncated_poisson()
)

summary(u_sy)
anova(u_sy, u_s)
# lmao damn

summary(u_s)

hist(unlist(ranef(u_s)$cond$Year))
plot(unlist(ranef(u_s)$cond$Year))
# it is increasing with time...

hist(unlist(ranef(u_s)$cond$Plot))
ranef(u_s)$cond$Plot %>%
  cbind(plot = row.names(.)) %>%
  merge(y = read.csv('00_raw_data/plot_treatments.csv')) %>%
  ggplot(aes(x = trt, y = `(Intercept)`, colour = trt)) +
  geom_point(size = 4) +
  scale_colour_manual(values = c('black', 'red', 'blue'))
# interesting... two of the control plots have the lowest intercepts...

hist(unlist(ranef(u_s)$cond$`plantid:Plot`))
# ooof okay this is not normally distributed...
ranef(u_s)$cond$`plantid:Plot` %>% filter(`(Intercept)` > 0.3) %>% row.names()
# what's up with these ones?

demo %>% 
  filter(plantid %in% c('3123_13_5H', '3354_15_11B', '3610_1_0C', '3360_4_2G'))
# damn... these ones are just prolific I guess
# 11 umbels... damn
