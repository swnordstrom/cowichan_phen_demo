# Script to get annual vital rate estimates from reproductive models
# These values will be used to get make year-specific kernels (if possible)

library(ggplot2)
library(tidyr)
library(dplyr)

rm(list = ls())

# --- Run script to prepare data

source('03_construct_kernels/prepare_demo_data_repr.R')

# --- Run models

u_s_s.ty = glmmTMB(
  No.umbels ~ size + (1 | Year) + (1 | Plot / plantid),
  family = 'truncated_poisson',
  ziformula = ~ size + trt + (1 | Year) + (1 | Year:trt) + (1 | Plot / plantid),
  data = demo.flow
)

s_st.y_sy.u = glmmTMB(
  no.seeds ~ trt * size + Year + (1 | Plot / plantid),
  ziformula = ~ size * Year + phen.umbels + (1 | Plot / plantid),
  family = 'nbinom2',
  data = demo.seed
)

r_t.y = glmmTMB(size ~ trt + (1 | Year) + (1 | Plot), data = demo.recr)

# --- Get predictions

# Oh I need to get annual random effects for this...
uss_ranef = data.frame(
  ranef(u_s_s.ty)$cond$Year %>% rename(numb = `(Intercept)`),
  ranef(u_s_s.ty)$zi$Year %>% rename(pflw = `(Intercept)`)
) %>%
  mutate(Year = row.names(.))

rty_ranef = data.frame(
  ranef(r_t.y)$cond$Year %>% rename(rcts = `(Intercept)`)
) %>%
  mutate(Year = row.names(.))

all_ranef = rbind(
  uss_ranef %>% pivot_longer(-Year, names_to = 'ranef', values_to = 'val'),
  rty_ranef %>% pivot_longer(-Year, names_to = 'ranef', values_to = 'val')
)

head(all_ranef)

# ### Find median plants in each model/model component
# WAIT - this probably won't work for multi-step models...
# 
# # # Flowering probability (zero inflation on umbel count model)
# ranef(u_s_s.ty)$zi$Plot %>% arrange(abs(`(Intercept)`))
# # Plot 11 (there's only one plot in plant 11 and its intercept is still quite large
# # next is 4
# ranef(u_s_s.ty)$zi$`plantid:Plot` %>% 
#   mutate(id = row.names(.)) %>% 
#   filter(grepl('\\:4$', id)) %>%
#   arrange(abs(`(Intercept)` + 0.13353997))
# # For flowering probability, median plant is 3571_4
# 
# # # Umbel count (conditional model on umbel count)
# ranef(u_s_s.ty)$cond$Plot %>% arrange(abs(`(Intercept)`))
# # plot 1
# ranef(u_s_s.ty)$zi$`plantid:Plot` %>%
#   mutate(id = row.names(.)) %>%
#   filter(grepl('\\:1$', id)) %>%
#   arrange(abs(`(Intercept)` - 0.02419042))
# # 3347_1 for umbel count
# 
# # # Recruit size
# ranef(r_t.y)$cond$Plot %>% arrange(abs(`(Intercept)`))
# # plot 4

# Initialize a data frame for generating predictions

backbone = expand.grid(
  size     = (5:60)/10,
  size.cur = (5:60)/10,
  trt = c('control', 'drought', 'irrigated'),
  Year = factor(2017:2024)
)

# --- Get probability of flowering and number of umbels

flow.numb = backbone %>%
  # Get mean zero inflation term (mean across all years) on the link scale
  mutate(
    flow.linear = predict(
      object = u_s_s.ty, newdata = ., re.form = ~ 0, allow.new.levels = TRUE,
      type = 'zlink'
    )
  ) %>%
  # Merge in to get the year random effects and add them in
  merge(all_ranef %>% filter(ranef %in% 'pflw')) %>%
  mutate(flow.linear.yr = flow.linear + val) %>%
  # estimate the PROBABILITY OF FLOWERING, NOT the probability of zero umbels
  mutate(p.flower.yr = 1 / (1 + exp(flow.linear.yr))) %>%
  # Remove unnecesary columns
  select(-c(ranef, val, flow.linear, flow.linear.yr)) %>%
  # Get mean number of umbels per plant (on link scale, in mean year)
  mutate(
    numb.linear = predict(
      object = u_s_s.ty, newdata = ., re.form = ~ 0, allow.new.levels = TRUE,
      type = 'link'
    )
  ) %>%
  # Merge in to get the random effect, add to linear predictor, get on natural scale
  merge(all_ranef %>% filter(ranef %in% 'numb')) %>%
  mutate(numb.linear.yr = numb.linear + val) %>%
  mutate(n.umbels.yr.cond = exp(numb.linear.yr)) %>%
  mutate(n.umbels.yr = p.flower.yr * n.umbels.yr.cond) %>%
  # Remove unnecessary columns (again)
  select(-c(ranef, val, numb.linear, numb.linear.yr))

head(flow.numb)

flow.numb %>%
  filter(size.cur %in% 0.5) %>%
  ggplot(aes(x = size, y = n.umbels.yr, group = trt)) +
  geom_line(aes(colour = trt)) +
  scale_colour_manual(values = c('black', 'red', 'blue')) +
  scale_y_log10() +
  facet_wrap(~ Year)

# Huh, you can see the treatment effects shifting around year-by-year

# Now get annual estimates of seed set
succ.seed = backbone %>%
  filter(as.numeric(as.character(Year)) > 2020) %>%
  # Merge in data frame above to get number of umbels per individual
  merge(
    flow.numb %>% select(-c(p.flower.yr, n.umbels.yr.cond)), 
    by = c('size', 'size.cur', 'trt', 'Year')
  ) %>%
  rename(phen.umbels = n.umbels.yr) %>%
  # Estimate umbel success and number of seeds
  mutate(
    p.umbel.failure.yr = predict(
      s_st.y_sy.u, newdata = ., allow.new.levels = TRUE, re.form = ~ 0,
      type = 'zprob'
    ),
    p.umbel.succ.yr = 1 - p.umbel.failure.yr,
    seeds.per.umbel.yr.cond = predict(
      s_st.y_sy.u, newdata = ., allow.new.levels = TRUE, re.form = ~ 0,
      type = 'cond'
    ),
    seeds.per.umbel.yr = predict(
      s_st.y_sy.u, newdata = ., allow.new.levels = TRUE, re.form = ~ 0,
      type = 'response'
    ),
    seeds.total.yr = phen.umbels * seeds.per.umbel.yr
  ) %>%
  # Remove unnecessary columns
  select(-c(phen.umbels, p.umbel.failure.yr))

succ.seed %>%
  filter(size.cur %in% 0.5) %>%
  ggplot(aes(x = size, y = seeds.total.yr, group = trt, colour = trt)) +
  geom_line() +
  scale_colour_manual(values = c('black', 'red', 'blue')) +
  scale_y_log10() +
  facet_wrap(~ Year)

# Huh. Cool.

# Next - get recruit size over time.
