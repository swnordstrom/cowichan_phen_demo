# Taking a first stab at looking at individual trade-offs in budding phenology
# i.e., is flowering (budding) earlier associated with a change in survival
# and/or growth?
# sn - init 2 apr 2024

library(ggplot2)
library(tidyr)
library(dplyr)
library(lme4)

rm(list = ls())

phen = read.csv('01_data_cleaning/out/seed_phen_demo_combined.csv')

head(phen)

demo = read.csv('01_data_cleaning/out/demo_imputed_survival.csv') %>%
  filter(Year > 2021) %>%
  mutate(Prev.year = Year - 1)

head(demo)

# Merging these together... is plantid column from demo? How often do they differ?

phen %>% filter(plantid != finalid.demo)
# 25 rows here - how many unique plants...?
phen %>% filter(plantid != finalid.demo) %>% distinct(Year, plantid, .keep_all = TRUE)
# 14 plants
# going to bet plantid is the demo key


phen.surv = merge(
  x = phen %>% 
    filter(Year < 2023) %>%
    group_by(Year, tagplot, plantid, n.phen.umbel, n.lost.umbel, n.demo.umbels, cur.size, mean.bud.day) %>%
    summarise(total.seeds = sum(no.seeds, na.rm = TRUE)) %>%
    ungroup(),
  y = demo %>%
    mutate(next.size = log(No.leaves * Leaf.length)) %>%
    select(Prev.year, Plot, trt, surv, next.size, plantid) %>%
    mutate(plantid = gsub('b', '', plantid)),
  by.x = c('Year', 'plantid'), by.y = c('Prev.year', 'plantid'),
  all.x = TRUE, all.y = FALSE
)

head(phen.surv)
nrow(phen.surv)

# Look for plants that are alive but missing a measurement
# (and remove them)

# There are 11
phen.surv %>% filter(surv & is.na(next.size))

# Remove them
phen.surv = phen.surv %>% filter(!(surv & is.na(next.size)))

#### Plots?

# Timing of budding vs. survival to next year
phen.surv %>%
  ggplot(aes(x = mean.bud.day, y = surv, colour = trt)) +
  geom_point(position = position_jitter(height = 1/8, width = 1), size = 3) +
  scale_colour_manual(values = c('black', 'red', 'blue'))
# eh - nothing visually obvious

# Plot of growth
phen.surv %>%
  # get rid of dead plants
  filter(surv) %>%
  # plot
  ggplot(aes(x = cur.size, y = next.size, colour = mean.bud.day)) +
  geom_segment(aes(x = 1, xend = 5.5, y = 1, yend = 5.5), linetype = 2) +
  geom_point(size = 3) +
  scale_colour_gradient(low = 'orange', high = 'royalblue') +
  facet_grid(Year ~ trt)
# Not super obvious to me here
# Growth looks pretty good 2021-2022 though...

# Another plot of growth
phen.surv %>%
  filter(surv, !is.na(cur.size)) %>%
  mutate(size.change = next.size - cur.size) %>%
  ggplot(aes(x = mean.bud.day, y = size.change)) +
  geom_point(aes(colour = trt), size = 3) +
  scale_colour_manual(values = c('black', 'red', 'blue'))
# huh.  

##### Models?

### Survival model

phen.surv = phen.surv %>% mutate(Year = factor(Year), Plot = factor(Plot))

s_0 = glm(
  surv ~ Year,
  family = 'binomial',
  data = phen.surv
)
# using a random effect here gives a singularity (probably some plot with no mortalities)

s_p = glm(
  surv ~ Year + mean.bud.day,
  family = 'binomial',
  data = phen.surv
)

anova(s_p, s_0, test = 'Chisq')
# oh... interesting
# but... conservativism due to lack of random effects? or does that not apply here?

summary(s_p)

# So a possible (marginally significant?) effect of phen on survival - 
# budding one day later associated with a reduction in survival by ~4%

### Growth model

phen.grow = phen.surv %>%
  filter(surv, !is.na(cur.size)) %>%
  mutate(
    size.change = next.size - cur.size,
    Year = factor(Year),
    Plot = factor(Plot)
  )

m_0 = lmer(
  size.change ~ Year + (1 | Plot),
  data = phen.grow
)

m_t = lmer(
  size.change ~ Year + trt + (1 | Plot),
  data = phen.grow
)

anova(m_t, m_0)
# interesting... no treatment effect...?

m_t_p = lmer(
  size.change ~ Year + trt + mean.bud.day + (1 | Plot),
  data = phen.grow
)

anova(m_t_p, m_t)
# So yes, does look like there is an effect even controlling for year and
# treatment

summary(m_t_p)

# and if we took treatment out?
m_p = lmer(
  size.change ~ Year + mean.bud.day + (1 | Plot),
  data = phen.grow
)

anova(m_t_p, m_p)
# Cool. Looks like there isn't a treatment effect at all.

# Year-treatment interactions?
m_py = lmer(
  size.change ~ Year * mean.bud.day + (1 | Plot),
  data = phen.grow
)

anova(m_py, m_p)
# Yes - very strong evidence.

summary(m_py)
# Oh so it's just an effect in one year?

phen.grow %>%
  mutate(pred.grow = predict(object = m_py, newdata = ., allow.new.levels = TRUE, re.form = ~ 0)) %>%
  ggplot(aes(x = mean.bud.day, colour = Year)) +
  geom_segment(aes(x = 105, xend = 160, y = 0, yend = 0), linetype = 2, colour = 'gray44') +
  geom_line(aes(y = pred.grow), linewidth = 1.2) +
  geom_point(aes(y = size.change))
# fits look fine
# flowering earlier (absolute date) means more growth in 2021?
# barely discernible effect in 2022

ranef(m_py)$Plot %>% unlist() %>% qqnorm()
ranef(m_py)$Plot %>% unlist() %>% qqline()
# woof... small sample size but that left tail is quite big

residuals(m_py) %>% hist()
# eh kinda big left tail but not super bad
qqnorm(residuals(m_py))
qqline(residuals(m_py))
# left tail deviates a bit but not awful
# (I guess we're not getting p-values though?)

phen.grow %>%
  mutate(resid = residuals(m_py)) %>%
  ggplot(aes(x = Year, y = resid)) +
  geom_point(position = position_jitter(width = 0.15))

phen.grow %>%
  mutate(resid = residuals(m_py)) %>%
  ggplot(aes(x = mean.bud.day, y = resid)) +
  geom_point()
# Not super concerned about heteroskedasticity here

phen.grow %>%
  mutate(resid = residuals(m_py)) %>%
  ggplot(aes(x = n.phen.umbel, y = resid)) +
  geom_point(position = position_jitter(width = 0.15))
# Looks like there may be a slight association with number of umbels?

phen.grow %>%
  mutate(resid = residuals(m_py)) %>%
  ggplot(aes(x = n.lost.umbel, y = resid)) +
  geom_point(position = position_jitter(width = 0.15))
# Meh

# Seems like a model 
