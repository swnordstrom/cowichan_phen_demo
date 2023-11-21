### Script for trying to fit mark-recapture model
### SN - 17 Nov 2023

### Setup
library(ggplot2)
library(dplyr)
library(tidyr)
library(marked)

# Clear namespace
rm(list = ls())

# Read in data (and merge with treatment dataframe)
demo = merge(
  x = read.csv('01_data_cleaning/out/demo_postcombine.csv'),
  y = read.csv('00_raw_data/plot_treatments.csv') %>% rename(Plot = plot)
)

head(demo)
nrow(demo)

# Select relevant columns only
demo = demo %>%
  select(
    plantid, Plot, trt, Year,
    No.leaves, Leaf.length, No.umbels,
    demo.note, proc.note, edited
  )

head(demo)
str(demo)
length(unique(demo$plantid)) # 1002 plants

### Generate capture history

demo.ch = demo %>%
  arrange(Year) %>%
  mutate(obs.alive = as.numeric(!((!No.leaves | is.na(No.leaves)) & (is.na(No.umbels) | !No.umbels)))) %>%
  mutate(Year = paste0('y', Year)) %>%
  pivot_wider(
    id_cols = c(plantid, Plot, trt),
    names_from = Year, values_from = obs.alive, values_fill = 0
  ) %>%
  unite(col = ch, all_of(starts_with('y')), sep = '') %>%
  mutate(Plot = factor(Plot), trt = factor(trt)) %>%
  # Needed to process data
  as.data.frame()

# Neat!
table(demo.ch$ch)

### Get design matrix?

# Looks like I need to do these separately
# Not sure if I trust how this routine handles more than one group

demo.proc.plt = process.data(
  data = demo.ch,
  begin.time = 2016,
  groups = "Plot"
)

demo.des.plt = make.design.data(
  data = demo.proc.plt,
  parameters = list(Phi = list(static = 'Plot'), p = list(static = 'Plot'))
)

demo.des.plt
demo.des.plt$Phi
demo.des.plt$p

# id is stored in demo.proc.plt

demo.cr_phi.1_p.1 = crm(
  demo.proc.plt,
  demo.des.plt
)

demo.cr_phi.t_p.t = crm(
  demo.proc.plt,
  demo.des.plt,
  model.parameters = list(
    Phi = list(formula = ~ time),
    p   = list(formula = ~ time)
  )
)

demo.cr_phi.t.p_p.t.p = crm(
  demo.proc.plt,
  demo.des.plt,
  model.parameters = list(
    Phi = list(formula = ~ time + Plot),
    p   = list(formula = ~ time + Plot)
  ),
  hessian = TRUE
)

demo.cr_phi.1_p.1
demo.cr_phi.t_p.t
demo.cr_phi.t.p_p.t.p # marginal improvement in AIC

demo.cr_phi.t.p_p.t = crm(
  demo.proc.plt,
  demo.des.plt,
  model.parameters = list(
    Phi = list(formula = ~ time + Plot),
    p   = list(formula = ~ time)
  ),
  hessian = TRUE
)

demo.cr_phi.t_p.t.p = crm(
  demo.proc.plt,
  demo.des.plt,
  model.parameters = list(
    Phi = list(formula = ~ time),
    p   = list(formula = ~ time + Plot)
  ),
  hessian = TRUE
)

demo.cr_phi.tp_p.tp = crm(
  demo.proc.plt,
  demo.des.plt,
  model.parameters = list(
    Phi = list(formula = ~ time * Plot),
    p   = list(formula = ~ time * Plot)
  ),
  hessian = TRUE
) # no convergey :(

demo.cr_phi.t_p.tp = crm(
  demo.proc.plt,
  demo.des.plt,
  model.parameters = list(
    Phi = list(formula = ~ time),
    p   = list(formula = ~ time * Plot)
  ),
  hessian = TRUE
) 


sapply(
  list(demo.cr_phi.1_p.1, demo.cr_phi.t_p.t, demo.cr_phi.t.p_p.t.p, 
       demo.cr_phi.t_p.t.p, demo.cr_phi.t.p_p.t.p,
       demo.cr_phi.tp_p.tp,
       demo.cr_phi.t_p.tp),
  function(x) x$results$AIC
) 

# demo.cr_phi.t.p_p.t.p$results$beta


predict(demo.cr_phi.t_p.tp)$Phi %>%
  mutate(time = as.numeric(as.character(time))) %>%
  filter(time < 2022) %>%
  ggplot(aes(x = time)) +
  geom_line(aes(y = estimate)) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.1) +
  lims(y = c(0.8, 1))

ggsave('02_data_exploration/figs/mrc_surv_est.png', width = 8, height = 3)

merge(
  x = predict(demo.cr_phi.t_p.tp)$p,
  y = distinct(demo.ch, Plot, trt)
) %>% 
  arrange(time, Plot) %>%
  mutate(time = as.numeric(as.character(time))) %>%
  ggplot(aes(x = time, group = Plot)) +
  geom_line(aes(y = estimate, colour = trt)) +
  scale_colour_manual(values = c('black', 'red', 'blue'))

ggsave('02_data_exploration/figs/mrc_detc_est.png', width = 8, height = 3)

pred_phi.t_p.tp = rbind(
  predict(demo.cr_phi.t_p.tp)$Phi %>% mutate(par = 'survival'),
  predict(demo.cr_phi.t_p.tp)$p %>% mutate(par = 'detection')
) %>%
  merge(distinct(demo.ch, Plot, trt))

pred_phi.t.p_p.t.p %>%
  ggplot(aes(x = time, y = estimate, group = Plot, colour = trt)) +
  geom_line() +
  facet_wrap(~ par) +
  labs(x = 'Year', y = 'Probability') +
  scale_colour_manual(values = c('black', 'red', 'blue'))

ggsave('02_data_exploration/figs/mrc_ests.png', width = 8, height = 5)

#######

demo.ch = demo %>%
  arrange(Year) %>%
  mutate(obs.alive = as.numeric(!((!No.leaves) & !is.na(No.leaves) & is.na(No.umbels)))) %>%
  mutate(Year = paste0('y', Year)) %>%
  pivot_wider(
    id_cols = c(plantid, Plot, trt),
    names_from = Year, values_from = obs.alive, values_fill = 0
  ) %>%
  unite(col = ch, all_of(starts_with('y')), sep = '') %>%
  mutate(Plot = factor(Plot), trt = factor(trt), freq = 1) %>%
  # Needed to process data
  as.data.frame()

demo.proc.plt = process.data(
  data = demo.ch,
  begin.time = 2016,
  groups = "Plot"
)

demo.des.plt = make.design.data(
  data = demo.proc.plt,
  parameters = list(Phi = list(static = 'Plot'), p = list(static = 'Plot'))
)

# # This doesn't work (lol)
# abcde = merge_design.covariates(
#   demo.des.plt$Phi,
#   df = demo %>% 
#     mutate(size = ifelse(
#       is.na(Leaf.length) | is.na(No.leaves),
#       NA,
#       log(No.leaves * Leaf.length)
#     )
#   ) %>%
#     select(plantid, size, Year) %>%
#     rename(time = Year) # %>%
#     # group_by(plantid) %>%
#     # filter(!any(is.na(size) & obs.alive) & No.umbels > 0)
# )

demo.des.plt$Phi = merge(
  x = demo.proc.plt$data[,c("plantid", "id")],
  y = demo
) %>%
  mutate(size = ifelse(
    is.na(Leaf.length) | is.na(No.leaves),
    NA,
    log(No.leaves * Leaf.length)
  )
) %>%
  select(id, size, Year) %>%
  rename(time = Year) %>%
  merge(y = demo.des.plt$Phi, all.y = TRUE) %>%
  select(id, occ, time, cohort, age, Plot, size, Time, Cohort, Age, order)

demo.cr_phi.1_p.1 = crm(
  demo.proc.plt,
  demo.des.plt
)

demo.cr_phi.t_p.t = crm(
  demo.proc.plt,
  demo.des.plt,
  model.parameters = list(
    Phi = list(formula = ~ time),
    p   = list(formula = ~ time)
  )
)

demo.cr_phi.t.x_p.t = crm(
  demo.proc.plt,
  demo.des.plt,
  model.parameters = list(
    Phi = list(formula = ~ time + size),
    p   = list(formula = ~ time)
  )
)

# arg....

#####
#####
# Try again, this time taking out plants with NAs for size when alive
#####

demo.rm = demo %>%
  arrange(Year) %>%
  mutate(obs.alive = as.numeric(!((!No.leaves | is.na(No.leaves)) & (is.na(No.umbels) | !No.umbels)))) %>%
  # Remove *all* records for any plant observed alive but without size
  group_by(plantid) %>%
  filter(!any(obs.alive & (is.na(No.leaves) | is.na(Leaf.length)))) %>%
  ungroup()

demo.rm.ch = demo.rm %>%
  mutate(Year = paste0('y', Year)) %>%
  pivot_wider(
    id_cols = c(plantid, Plot, trt),
    names_from = Year, values_from = obs.alive, values_fill = 0
  ) %>%
  unite(col = ch, all_of(starts_with('y')), sep = '') %>%
  mutate(Plot = factor(Plot), trt = factor(trt)) %>%
  # Needed to process data
  as.data.frame()

demo.proc.plt = process.data(
  data = demo.rm.ch,
  begin.time = 2016,
  groups = "Plot"
)

demo.des.plt = make.design.data(
  data = demo.proc.plt,
  parameters = list(Phi = list(static = 'Plot'), p = list(static = 'Plot'))
)

demo.des.plt$Phi = merge(
  x = demo.proc.plt$data[,c("plantid", "id")],
  y = demo.rm
) %>%
  mutate(size = ifelse(
    is.na(Leaf.length) | is.na(No.leaves),
    NA,
    log(No.leaves * Leaf.length)
  )
  ) %>%
  select(id, size, Year) %>%
  rename(time = Year) %>%
  merge(y = demo.des.plt$Phi, all.y = TRUE) %>%
  select(id, occ, time, cohort, age, Plot, size, Time, Cohort, Age, order)

demo.cr_phi.t.x_p.t = crm(
  demo.proc.plt,
  demo.des.plt,
  model.parameters = list(
    Phi = list(formula = ~ time + size),
    p   = list(formula = ~ time)
  )
)
