# Script for analysis of flowering phenology and flowering phenology vs. seed
# set.
# 2021-2024

library(ggplot2)
library(dplyr)
library(tidyr)
library(glmmTMB)
library(cowplot)

rm(list = ls())

# Read in data
all.data = merge(
  x = read.csv('01_data_cleaning/out/demo_phen_seed_2016-2024_final.csv'),
  y = read.csv('00_raw_data/plot_treatments.csv'),
  by.x = 'Plot', by.y = 'plot'
) %>%
  mutate(Year = factor(Year))

# Phen only
phen.each.umbel = all.data %>% 
  filter(in.phen) %>%
  # Split out the bud dates for bud date models; the most umbels seen in a
  # plant is 12, so use separate() to kick these out and then pivot_long to get
  # one row per umbel
  # (first - need to get one row per plant - do a distinct())
  distinct(Year, plantid, .keep_all = TRUE) %>%
  separate_wider_delim(phen.julis, names = paste0('uu', 1:12), delim = ';', too_few = 'align_start') %>%
  pivot_longer(starts_with('uu'), names_to = 'umbel.number', values_to = 'phen.julian') %>%
  filter(!is.na(phen.julian)) %>%
  mutate(phen.julian = as.numeric(gsub('\\s', '', phen.julian)))

head(phen.each.umbel)
nrow(phen.each.umbel)
sum(phen.each.umbel$phen.umbels)
hist(phen.each.umbel$phen.julian)

# All plants in phen, seed, and demo
seed = all.data %>% 
  filter(in.phen, in.seed, in.demo) %>%
  # Need size data
  filter(!is.na(No.leaves), !is.na(Leaf.length), No.leaves > 0, Leaf.length > 0) %>%
  # Hmm... okay, separae() and across() doesn't work, so I guess I can merge
  # this in with the data frame above, summarised by mean phen date per plant
  merge(
    y = phen.each.umbel %>% 
      group_by(plantid, Year) %>% 
      summarise(
        mean.phen = mean(phen.julian, na.rm = TRUE), 
        sept.phen = length(unique(phen.julian))
      ) %>%
      ungroup(),
    all.x = TRUE, all.y = FALSE
  )

# Now, need to impute in zeros for failed umbels that did not make it to seed counting

seed = rbind(
  seed,
  seed %>%
    group_by(Year, plantid) %>%
    mutate(miss.umbel = ifelse(phen.umbels < n(), 0, phen.umbels - n())) %>%
    ungroup() %>%
    distinct(Year, plantid, .keep_all = TRUE) %>%
    uncount(miss.umbel) %>%
    mutate(no.seeds = 0)
) %>%
  mutate(phen.c = mean.phen - round(mean(mean.phen))) %>%
  mutate(size = log(No.leaves * Leaf.length))

head(seed)
nrow(seed)
hist(seed$mean.phen)


# -----------------------------------------
# Analysis -------------
# -----------------------------------------

# ----------- Just phenology --------------

# Response: date of bud initiation per plant
# Gaussian probably works best

d_0 = glmmTMB(
  phen.julian ~ Year + (1 | Plot / plantid),
  data = phen.each.umbel
)

summary(d_0)
# Huge within-plant variation

d_t = glmmTMB(
  phen.julian ~ trt + Year + (1 | Plot / plantid),
  data = phen.each.umbel
)

summary(d_t)

AIC(d_t, d_0) # delta aic of 7
anova(d_t, d_0)

# Test for year-varying treatment effects

d_ty = d_0 = glmmTMB(
  phen.julian ~ Year * trt + (1 | Plot / plantid),
  data = phen.each.umbel
)

summary(d_ty)

AIC(d_ty, d_t)
anova(d_ty, d_t)
# marginally significant ANOVA, higher AIC for interaction model

# Best phen model has trt effect, same across years

# Plot model predictions:
# (I guess averaged across years...)
# (how to get an effect size plot here...)

# This is crude but it will get it done:
summary(d_t)$coefficients$cond %>%
  as.data.frame() %>%
  mutate(coef = row.names(.)) %>%
  filter(grepl('trt', coef)) %>%
  mutate(trt = gsub('trt', '', coef)) %>%
  ggplot(aes(y = trt)) +
  geom_segment(
    aes(
      x = Estimate - 1.96 * `Std. Error`,
      xend = Estimate + 1.96 * `Std. Error`, 
      yend = trt
    )
  ) +
  geom_point(aes(x = Estimate, colour = trt), size = 4) +
  scale_colour_manual(values = c('red', 'blue')) +
  labs(x = 'Effect size (days)', y = 'Treatment') +
  theme(legend.position = 'none')

# ----------- Seeds per umbel -----------

s_0 = glmmTMB(
  no.seeds ~ Year + (1 | Plot / plantid),
  family = 'nbinom2',
  ziformula = ~ Year + (1 | Plot / plantid),
  data = seed
)

s_s_0 = glmmTMB(
  no.seeds ~ size + Year + (1 | Plot / plantid),
  family = 'nbinom2',
  ziformula = ~ Year + (1 | Plot / plantid),
  data = seed
)

s_0_s = glmmTMB(
  no.seeds ~ Year + (1 | Plot / plantid),
  family = 'nbinom2',
  ziformula = ~ size + Year + (1 | Plot / plantid),
  data = seed
)

s_s_s = glmmTMB(
  no.seeds ~ size + Year + (1 | Plot / plantid),
  family = 'nbinom2',
  ziformula = ~ size + Year + (1 | Plot / plantid),
  data = seed
)

AIC(s_0, s_s_0, s_0_s, s_s_s) %>% arrange(AIC) %>% mutate(daic = round(AIC - min(AIC), 2))

# Size in all models

s_s.t_s = glmmTMB(
  no.seeds ~ trt + size + Year + (1 | Plot / plantid),
  family = 'nbinom2',
  ziformula = ~ size + Year + (1 | Plot / plantid),
  data = seed
)

s_s_s.t = glmmTMB(
  no.seeds ~ size + Year + (1 | Plot / plantid),
  family = 'nbinom2',
  ziformula = ~ trt + size + Year + (1 | Plot / plantid),
  data = seed
)

s_s.t_s.t = glmmTMB(
  no.seeds ~ trt + size + Year + (1 | Plot / plantid),
  family = 'nbinom2',
  ziformula = ~ trt + size + Year + (1 | Plot / plantid),
  data = seed
)

s_st_s = glmmTMB(
  no.seeds ~ trt * size + Year + (1 | Plot / plantid),
  family = 'nbinom2',
  ziformula = ~ size + Year + (1 | Plot / plantid),
  data = seed
)

s_s_st = glmmTMB(
  no.seeds ~ size + Year + (1 | Plot / plantid),
  family = 'nbinom2',
  ziformula = ~ trt * size + Year + (1 | Plot / plantid),
  data = seed
)

s_st_st = glmmTMB(
  no.seeds ~ trt * size + Year + (1 | Plot / plantid),
  family = 'nbinom2',
  ziformula = ~ trt * size + Year + (1 | Plot / plantid),
  data = seed
)

AIC(s_s_s, s_s.t_s, s_st_s, s_s_s.t, s_s_st, s_s.t_s.t, s_st_st) %>%
  mutate(daic = round(AIC - min(AIC), 2)) %>%
  arrange(daic)

# Best of these models has size-treatment interaction for seeds
# Check to make sure there's no linear treatment effect on the ZI term

s_st_s.t = glmmTMB(
  no.seeds ~ trt * size + Year + (1 | Plot / plantid),
  family = 'nbinom2',
  ziformula = ~ trt + size + Year + (1 | Plot / plantid),
  data = seed
)

# Test for treatment-year effects
s_st.ty_s = glmmTMB(
  no.seeds ~ trt * size + trt * Year + (1 | Plot / plantid),
  family = 'nbinom2',
  ziformula = ~ size + Year + (1 | Plot / plantid),
  data = seed
)

s_st_s.ty = glmmTMB(
  no.seeds ~ trt * size + Year + (1 | Plot / plantid),
  family = 'nbinom2',
  ziformula = ~ size + trt * Year + (1 | Plot / plantid),
  data = seed
)

s_st.ty_s.ty = glmmTMB(
  no.seeds ~ trt * size + Year + (1 | Plot / plantid),
  family = 'nbinom2',
  ziformula = ~ size + trt * Year + (1 | Plot / plantid),
  data = seed
)

AIC(s_st_s, s_st_s.t, s_st.ty_s, s_st_s.ty, s_st.ty_s.ty) %>%
  mutate(daic = round(AIC - min(AIC), 2))
# so no treatment-year effects, very cool

# Look for umbel count effects

s_st.u_s = glmmTMB(
  no.seeds ~ trt * size + phen.umbels + Year + (1 | Plot / plantid),
  family = 'nbinom2',
  ziformula = ~ size + Year + (1 | Plot / plantid),
  data = seed
)

s_st_s.u = glmmTMB(
  no.seeds ~ trt * size + Year + (1 | Plot / plantid),
  family = 'nbinom2',
  ziformula = ~ size + phen.umbels + Year + (1 | Plot / plantid),
  data = seed
)

s_st.u_s.u = glmmTMB(
  no.seeds ~ trt * size + phen.umbels + Year + (1 | Plot / plantid),
  family = 'nbinom2',
  ziformula = ~ size + phen.umbels + Year + (1 | Plot / plantid),
  data = seed
)

AIC(s_st_s, s_st.u_s, s_st_s.u, s_st.u_s.u) %>%
  mutate(daic = round(AIC - min(AIC), 2))

# Umbel count in zero inflation term (makes sense)

# NOW, tests for effects of phenology

s_st.p_s.u = glmmTMB(
  no.seeds ~ trt * size + Year + phen.c + (1 | Plot / plantid),
  family = 'nbinom2',
  ziformula = ~ size + phen.umbels + Year + (1 | Plot / plantid),
  data = seed
)

s_st.p2_s.u = glmmTMB(
  no.seeds ~ trt * size + Year + poly(phen.c, 2) + (1 | Plot / plantid),
  family = 'nbinom2',
  ziformula = ~ size + phen.umbels + Year + (1 | Plot / plantid),
  data = seed
)

s_st.tp_s.u = glmmTMB(
  no.seeds ~ trt * size + Year + trt * phen.c + (1 | Plot / plantid),
  family = 'nbinom2',
  ziformula = ~ size + phen.umbels + Year + (1 | Plot / plantid),
  data = seed
)

s_st.tp2_s.u = glmmTMB(
  no.seeds ~ trt * size + Year + trt * poly(phen.c, 2) + (1 | Plot / plantid),
  family = 'nbinom2',
  ziformula = ~ size + phen.umbels + Year + (1 | Plot / plantid),
  data = seed
)

s_st.py_s.u = glmmTMB(
  no.seeds ~ trt * size + Year * phen.c + (1 | Plot / plantid),
  family = 'nbinom2',
  ziformula = ~ size + phen.umbels + Year + (1 | Plot / plantid),
  data = seed
)

s_st.p2y_s.u = glmmTMB(
  no.seeds ~ trt * size + Year * phen.c + (1 | Plot / plantid),
  family = 'nbinom2',
  ziformula = ~ size + phen.umbels + Year * poly(phen.c, 2) + (1 | Plot / plantid),
  data = seed
)

AIC(s_st_s.u, s_st.p_s.u, s_st.p2_s.u, s_st.tp_s.u, s_st.py_s.u, s_st.tp2_s.u, s_st.p2y_s.u) %>%
  mutate(daic = round(AIC - min(AIC), 2))
# Evidence of linear phenology effect on seed count, not failure
# No evidence of interactions with treatment or year, or evidence of a quadratic relationship

# Look for phenology effects on the zero-inflation term

s_st.p_s.u.p = glmmTMB(
  no.seeds ~ trt * size + Year + phen.c + (1 | Plot / plantid),
  family = 'nbinom2',
  ziformula = ~ size + phen.umbels + Year + phen.c + (1 | Plot / plantid),
  data = seed
)

s_st.p_s.u.tp = glmmTMB(
  no.seeds ~ trt * size + Year + phen.c + (1 | Plot / plantid),
  family = 'nbinom2',
  ziformula = ~ size + phen.umbels + Year + trt * phen.c + (1 | Plot / plantid),
  data = seed
)

s_st.p_s.u.py = glmmTMB(
  no.seeds ~ trt * size + Year + phen.c + (1 | Plot / plantid),
  family = 'nbinom2',
  ziformula = ~ size + phen.umbels + Year * phen.c + (1 | Plot / plantid),
  data = seed
)

s_st.p_s.u.p2 = glmmTMB(
  no.seeds ~ trt * size + Year + phen.c + (1 | Plot / plantid),
  family = 'nbinom2',
  ziformula = ~ size + phen.umbels + Year + poly(phen.c, 2) + (1 | Plot / plantid),
  data = seed
)

s_st.p_s.u.tp2 = glmmTMB(
  no.seeds ~ trt * size + Year + phen.c + (1 | Plot / plantid),
  family = 'nbinom2',
  ziformula = ~ size + phen.umbels + Year + trt * poly(phen.c, 2) + (1 | Plot / plantid),
  data = seed
)

s_st.p_s.u.p2y = glmmTMB(
  no.seeds ~ trt * size + Year + phen.c + (1 | Plot / plantid),
  family = 'nbinom2',
  ziformula = ~ size + phen.umbels + Year * poly(phen.c, 2) + (1 | Plot / plantid),
  data = seed
)

AIC(s_st.p_s.u, s_st.p_s.u.p, s_st.p_s.u.p2, s_st.p_s.u.tp, s_st.p_s.u.py, s_st.p_s.u.tp2, s_st.p_s.u.p2y) %>% 
  mutate(daic = round(AIC - min(AIC), 2))
# Huh... polynomial phenology effect on zero-inflation term?

summary(s_st.p_s.u.p2)

# Examine model predictions

seed.all.preds = expand.grid(
  Year = factor(2021:2024),
  size = (5:60)/10,
  phen.umbels = c(1, 2, 5),
  trt = c('control', 'drought', 'irrigated'),
  phen.c = (-4:4)*7
) %>%
  mutate(
    zero = predict(s_st.p_s.u.p2, newdata = ., allow.new.levels = TRUE, re.form = ~ 0, type = 'zprob'),
    zlnk = predict(s_st.p_s.u.p2, newdata = ., allow.new.levels = TRUE, re.form = ~ 0, type = 'zlink'),
    cond = predict(s_st.p_s.u.p2, newdata = ., allow.new.levels = TRUE, re.form = ~ 0, type = 'conditional'),
    resp = predict(s_st.p_s.u.p2, newdata = ., allow.new.levels = TRUE, re.form = ~ 0, type = 'response'),
    link = predict(s_st.p_s.u.p2, newdata = ., allow.new.levels = TRUE, re.form = ~ 0, type = 'link'),
  )

# Okay... lots here... how to plot it all lmao

seed.all.preds %>%
  filter(phen.umbels < 2) %>%
  # select(-c(link, zlnk, zero, resp)) %>%
  # pivot_longer(c(zero, cond, resp), names_to = 'pred.type', values_to = 'pred') %>%
  ggplot(aes(x = size, y = resp, group = interaction(trt, phen.c))) +
  geom_line(aes(colour = phen.c)) +
  scale_colour_viridis_c() +
  scale_y_log10() +
  facet_wrap(trt ~ Year)
# So overall - more seeds in control, more seeds when *early* (not late), more
# seeds for larger plants but not by a ton

# Let's look at that zero inflation
seed.all.preds %>%
  filter(phen.umbels < 2) %>%
  # select(-c(link, zlnk, zero, resp)) %>%
  # pivot_longer(c(zero, cond, resp), names_to = 'pred.type', values_to = 'pred') %>%
  ggplot(aes(x = size, y = zero, group = interaction(trt, phen.c))) +
  geom_line(aes(colour = phen.c)) +
  scale_colour_viridis_c() +
  facet_wrap(trt ~ Year)
# Smaller plants more likely to have umbel failure (not surprising), umbel
# failure is more likely for plants with earlier mean-budding dates
# Earlier budding means more likely umbel failure

seed.all.preds %>%
  # Evaluated for a plant with one umbel at the IQR
  filter(phen.umbels < 2, size %in% c(3.5, 3.9, 4.3)) %>%
  # select(-c(link, zlnk, zero, resp)) %>%
  # pivot_longer(c(zero, cond, resp), names_to = 'pred.type', values_to = 'pred') %>%
  ggplot(aes(x = phen.c, y = resp, group = trt)) +
  geom_line(aes(colour = trt)) +
  scale_colour_manual(values = c('black', 'red', 'blue')) +
  scale_y_log10() +
  facet_wrap(size ~ Year)

# Overall: Most seeds for plants that bud early
# Fewest seeds for plants in drought 
# Smaller-sized plants: irrigation has slight advantage over control, larger
# plants: more seeds in controls

seed.all.preds %>%
  # Evaluated for a plant with one umbel at the IQR
  filter(phen.umbels < 2, size %in% c(3.5, 3.9, 4.3)) %>%
  # select(-c(link, zlnk, zero, resp)) %>%
  # pivot_longer(c(zero, cond, resp), names_to = 'pred.type', values_to = 'pred') %>%
  ggplot(aes(x = phen.c, y = zero)) +
  geom_line() +
  # scale_colour_manual(values = c('black', 'red', 'blue')) +
  facet_wrap(size ~ Year)

# Okay - let's try to get averages across years
seed.linear.means = seed.all.preds %>%
  select(-c(zero, resp, cond)) %>%
  group_by(phen.umbels, trt, size, phen.c) %>%
  summarise(across(c(zlnk, link), mean)) %>%
  ungroup()

seed.linear.means %>%
  filter(size %in% c(3.5, 3.9, 4.3)) %>%
  mutate(p.failure = 1 / (1 + exp(-zlnk))) %>%
  ggplot(aes(x = phen.c, y = p.failure, colour = trt)) +
  geom_line() +
  scale_colour_manual(values = c('black', 'red', 'blue')) +
  facet_wrap(size ~ phen.umbels)  

seed.linear.means %>%
  filter(size %in% c(3.5, 3.9, 4.3)) %>%
  mutate(mean.seeds = exp(link) / (1 + exp(zlnk))) %>%
  mutate(phen.umbels = factor(phen.umbels)) %>%
  ggplot(aes(x = phen.c, y = mean.seeds, colour = trt, linetype = phen.umbels)) +
  geom_line() +
  scale_colour_manual(values = c('black', 'red', 'blue')) +
  facet_wrap(~ size)  

summary(s_st.p_s.u.p2)

# okay - make a three-row plot (one per prediction type) to send to jenn

seed.resp.means = seed.linear.means %>%
  # filter(size %in% c(3.5, 3.9, 4.3), phen.umbels < 2) %>%
  mutate(
    prob.succ = 1 / (1 + exp(zlnk)),
    seed.cond = exp(link),
    mean.seed = prob.succ * seed.cond
  ) %>%
  mutate(phen = as.Date(phen.c + round(mean(seed$mean.phen)), format = '%b-%d')) # %>%
  # pivot_longer(c(prob.succ, seed.cond, mean.seed), names_to = 'pred.type', values_to = 'value')


p.succ.plot = seed.resp.means %>%
  filter(size %in% c(3.5, 3.9, 4.3), phen.umbels < 2) %>%
  ggplot(aes(x = phen, y = prob.succ, colour = trt, group = interaction(size, trt))) +
  geom_line() +
  scale_colour_manual(values = c('black', 'red', 'blue'), 'treatment') +
  labs(x = '', y = 'Prob. succesful umbel') +
  facet_wrap( ~ size) +
  theme(legend.position = 'none', axis.text.x = element_blank())

p.cond.plot = seed.resp.means %>%
  filter(size %in% c(3.5, 3.9, 4.3), phen.umbels < 2) %>%
  ggplot(aes(x = phen, y = seed.cond, colour = trt, group = interaction(size, trt))) +
  geom_line() +
  scale_colour_manual(values = c('black', 'red', 'blue')) +
  scale_y_log10() +
  labs(x = '', y = 'Seeds per succ. umble') +
  facet_wrap( ~ size) +
  theme(legend.position = 'none', axis.text.x = element_blank())

p.resp.plot = seed.resp.means %>%
  filter(size %in% c(3.5, 3.9, 4.3), phen.umbels < 2) %>%
  ggplot(aes(x = phen, y = mean.seed, colour = trt, group = interaction(size, trt))) +
  geom_line() +
  scale_colour_manual(values = c('black', 'red', 'blue')) +
  scale_y_log10() +
  labs(x = '', y = 'Overall seed set') +
  facet_wrap( ~ size) +
  theme(legend.position = 'none', axis.text.x = element_text(angle = 90))

p.plot.lgnd = get_legend(p.succ.plot + theme(legend.position = 'top'))

plot_grid(p.plot.lgnd, p.succ.plot, p.cond.plot, p.resp.plot, rel_heights = c(.1, 1, 1, 1), nrow = 4) %>%
  save_plot(filename = '02_data_exploration/figs/phen_seed_2021-2024.png', base_height = 8, base_width = 8)

