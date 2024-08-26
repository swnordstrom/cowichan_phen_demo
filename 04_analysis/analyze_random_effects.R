library(ggplot2)
library(dplyr)
library(tidyr)
library(glmmTMB)
library(cowplot)

rm(list = ls())

# Run wrapper script to prepare demo + phen data data
source('03_construct_kernels/prepare_demo_data_repr.R')

# Run wrapper script to prepare growth + surv demo data
source('03_construct_kernels/prepare_demo_data_growsurv.R')

# --- Run models ----------------------------------------

# === Survival model ===
s_s = glmmTMB(
  formula = surv ~ size.prev + (1 | Plot),
  family = 'binomial',
  data = demo.surv.sizes
)

# === Growth model ===
g_st.ty = glmmTMB(
  size.cur ~ size.prev + size.prev * trt + (1 | prev.year) + (1 | prev.year:trt) + (1 | Plot / plantid),
  data = demo.grow
)

# === Phenology model ===
# Response: day of umbel budding (continuous julian date)
# Predictors: year (factor/categorical) and treatment (factor)
d_t = glmmTMB(
  phen.julian ~ trt + Year + (1 | Plot / plantid),
  data = phen
)

# === Flowering/umbel count model ===
u_s_s.ty = glmmTMB(
  No.umbels ~ size + (1 | Year) + (1 | Plot / plantid),
  family = 'truncated_poisson',
  ziformula = ~ size + trt + (1 | Year) + (1 | Year:trt) + (1 | Plot / plantid),
  data = demo.flow
)

# === Seed model ===
# Response: number of seeds (negative binomial distribution) of an umbel
# Predictors: treatment (categorical), year (factor) , mean budding date of
# plant (centered, continuous), plant size (continuous)
s_st.p_s.u.p2 = glmmTMB(
  no.seeds ~ trt * size + Year + phen.c + (1 | Plot / plantid),
  family = 'nbinom2',
  ziformula = ~ size + phen.umbels + Year + poly(phen.c, 2) + (1 | Plot / plantid),
  data = seed
)

# === Recruit size model ===
r_t.y = glmmTMB(size ~ trt + (1 | Year) + (1 | Plot), data = demo.recr)


# --- Extract effects ----------------------------------------

recr.year = ranef(r_t.y)$cond$Year %>%
  as.data.frame() %>%
  mutate(label = row.names(.)) %>%
  rename(int = `(Intercept)`) %>%
  mutate(type = 'ranef', mod = 'recr', var = 'Year')

recr.plot = ranef(r_t.y)$cond$Plot %>%
  as.data.frame() %>%
  mutate(label = row.names(.)) %>%
  rename(int = `(Intercept)`) %>%
  mutate(type = 'ranef', mod = 'recr', var = 'Plot')

seed.plot = ranef(s_st.p_s.u.p2)$cond$Plot %>%
  as.data.frame() %>%
  mutate(label = row.names(.)) %>%
  rename(int = `(Intercept)`) %>%
  mutate(type = 'ranef', mod = 'seed', var = 'Plot')

succ.plot = ranef(s_st.p_s.u.p2)$zi$Plot %>%
  as.data.frame() %>%
  mutate(label = row.names(.)) %>%
  rename(int = `(Intercept)`) %>%
  mutate(type = 'ranef', mod = 'succ', var = 'Plot')

numb.year = ranef(u_s_s.ty)$cond$Year %>%
  as.data.frame() %>%
  mutate(label = row.names(.)) %>%
  rename(int = `(Intercept)`) %>%
  mutate(type = 'ranef', mod = 'numb', var = 'Year')

numb.plot = ranef(u_s_s.ty)$cond$Plot %>%
  as.data.frame() %>%
  mutate(label = row.names(.)) %>%
  rename(int = `(Intercept)`) %>%
  mutate(type = 'ranef', mod = 'numb', var = 'Plot')

flow.plot = ranef(u_s_s.ty)$zi$Plot %>%
  as.data.frame() %>%
  mutate(label = row.names(.)) %>%
  rename(int = `(Intercept)`) %>%
  mutate(type = 'ranef', mod = 'flow', var = 'Plot')

flow.year = ranef(u_s_s.ty)$zi$Year %>%
  as.data.frame() %>%
  mutate(label = row.names(.)) %>%
  rename(int = `(Intercept)`) %>%
  mutate(type = 'ranef', mod = 'flow', var = 'Year')

grow.year = ranef(g_st.ty)$cond$prev.year %>%
  as.data.frame() %>%
  mutate(label = row.names(.)) %>%
  rename(int = `(Intercept)`) %>%
  mutate(type = 'ranef', mod = 'grow', var = 'Year')

grow.plot = ranef(g_st.ty)$cond$Plot %>%
  as.data.frame() %>%
  mutate(label = row.names(.)) %>%
  rename(int = `(Intercept)`) %>%
  mutate(type = 'ranef', mod = 'grow', var = 'Plot')

surv.plot = ranef(s_s)$cond$Plot %>%
  as.data.frame() %>%
  mutate(label = row.names(.)) %>%
  rename(int = `(Intercept)`) %>%
  mutate(type = 'ranef', mod = 'surv', var = 'Plot')

# seed.year = fixef(s_st.p_s.u.p2)$cond[c("(Intercept)", grep('^Year', names(fixef(s_st.p_s.u.p2)$cond), value = TRUE))] %>%
#   (function(x) c(x[1], x[1] + x[2], x[1] + x[3], x[1] + x[4])) %>%
#   data.frame(int = ., label = 2021:2024, type = 'fixef', mod = 'seed', var = 'Year')

seed.year = fixef(s_st.p_s.u.p2)$cond %>%
  data.frame(int = .) %>%
  mutate(label = row.names(.)) %>%
  filter(label %in% '(Intercept)' | grepl('^Year', label)) %>%
  mutate(int = int + ifelse(label %in% '(Intercept)', 0, int[label %in% '(Intercept)'])) %>%
  mutate(label = as.numeric(ifelse(label %in% '(Intercept)', 2021, gsub('Year', '', label)))) %>%
  mutate(type = 'fixef', mod = 'seed', var = 'Year')

succ.year = fixef(s_st.p_s.u.p2)$zi %>%
  data.frame(int = .) %>%
  mutate(label = row.names(.)) %>%
  filter(label %in% '(Intercept)' | grepl('^Year', label)) %>%
  mutate(int = int + ifelse(label %in% '(Intercept)', 0, int[label %in% '(Intercept)'])) %>%
  mutate(label = as.numeric(ifelse(label %in% '(Intercept)', 2021, gsub('Year', '', label)))) %>%
  mutate(type = 'fixef', mod = 'succ', var = 'Year')


all.efx = rbind(
  grow.year, flow.year, numb.year, recr.year,
  grow.plot, flow.plot, numb.plot, recr.plot,
  surv.plot, seed.plot, succ.plot,
  seed.year, succ.year
)

row.names(all.efx) = NULL

all.efx %>%
  filter(var %in% 'Plot') %>%
  pivot_wider(names_from = mod, values_from = int) %>%
  select(-c(var, type, label)) %>%
  filter(complete.cases(.)) %>%
  cor() # %>%
  #image()

# ah man how do I get rid of the 

all.efx %>%
  filter(var %in% 'Year') %>%
  mutate(label = as.numeric(label)) %>% 
  select(-c(type, var)) %>%
  # I think we can use each of these columns as is... i.e. don't have to
  # add/subtract to align years
  pivot_wider(names_from = mod, values_from = int)
