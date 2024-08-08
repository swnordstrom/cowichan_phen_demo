
### ---------------------------------

library(ggplot2)
library(dplyr)
library(tidyr)

rm(list = ls())

### Observed vitals
# Growth + survival kernel
growsurv = read.csv('03_construct_kernels/out/deterministic_growsurv_kernel.csv')
# Reproductive kernel (all phenology)
reprodct = read.csv('03_construct_kernels/out/deterministic_reprod_kernel_phen.csv')
# Reproductive kernel (for LTRE only)
reprodct.ltre = read.csv('03_construct_kernels/out/determinstic_reprod_kernel_phen_ltre.csv')

head(growsurv)
head(reprodct)

p.germ = .001

kernel.df = merge(
  growsurv %>% select(-c(pred.surv, pred.grow.mean, p.grow.size)),
  reprodct,
  by.x = c('size.prev', 'size.cur', 'trt'), by.y = c('size.prev', 'size.nex', 'trt'),
  suffixes = c('.g', '.r')
) %>%
  mutate(p.size.cur = p.size.cur.g + p.size.cur.r * p.germ) %>%
  select(-c(p.size.cur.g, p.size.cur.r))

head(kernel.df)

all.matr = split(kernel.df, kernel.df[,c("trt", "phen")], sep = '_') %>%
  lapply(
    function(df) {
      df %>%
        arrange(size.prev, size.cur) %>%
        pivot_wider(names_from = size.prev, values_from = p.size.cur) %>%
        select(-c(trt, phen, size.cur)) %>%
        as.matrix()
    }
  )

all.lambda = all.matr %>% 
  sapply(function(m) Re(eigen(m)$values[1])) %>%
  data.frame(lambda = .) %>%
  mutate(trt_phen = row.names(.)) %>%
  separate(trt_phen, into = c('trt', 'phen'), sep = '_')

ltre.kernel.df = merge(
  growsurv %>% select(-c(pred.surv, pred.grow.mean, p.grow.size)),
  reprodct.ltre,
  by.x = c('size.prev', 'size.cur', 'trt'), by.y = c('size.prev', 'size.nex', 'trt'),
  suffixes = c('.g', '.r')
) %>%
  mutate(p.size.cur = p.size.cur.g + p.size.cur.r * p.germ) %>%
  select(-c(p.size.cur.g, p.size.cur.r))

ltre.matr = split(ltre.kernel.df, ltre.kernel.df[,c("trt", "phen")], sep = '_') %>%
  lapply(
    function(df) {
      df %>%
        arrange(size.prev, size.cur) %>%
        pivot_wider(names_from = size.prev, values_from = p.size.cur) %>%
        select(-c(trt, phen, size.cur)) %>%
        as.matrix()
    }
  )

ltre.lambda = ltre.matr %>%
  sapply(function(m) ifelse(nrow(m) > 0, Re(eigen(m)$values[1]), NA)) %>%
  data.frame(lambda = .) %>%
  filter(!is.na(lambda)) %>%
  mutate(trt_phen = row.names(.)) %>%
  separate(trt_phen, into = c('trt', 'phen'), sep = '_') %>%
  arrange(trt) %>%
  mutate(
    phen = as.numeric(phen),
    obs.phen = (trt %in% 'control' & floor(phen) == 125) |
      (trt %in% 'drought' & floor(phen) == 122) | (trt %in% 'irrigated' & floor(phen) == 127),
    phen.date = as.Date(phen, format = '%b-%d')
  )

head(all.lambda)
head(ltre.lambda)

all.lambda %>%
  mutate(phen.date = as.Date(as.numeric((phen)), format = '%b-%d')) %>%
  ggplot(aes(x = phen.date, y = lambda, group = trt, colour = trt)) +
  geom_line() +
  scale_colour_manual(values = c('black', 'red', 'blue')) +
  geom_point(data = ltre.lambda, aes(shape = obs.phen), size = 4) +
  scale_shape_manual(values = c(1, 19)) +
  guides(shape = 'none') +
  theme(panel.background = element_blank())

# --- 

# Perturbations for sensitivity analysis

# Growth sensitivities
growsurv.pert = read.csv('03_construct_kernels/out/deterministic_grow_coef_perturbation_no_phen.csv') %>%
  rename(size.nex = size.cur, param = perturb.param, orig.parval = orig.par.val) %>%
  mutate(param = gsub('\\_', '.', param)) %>%
  # Going to take out the survival terms here because they don't vary by treatment
  # (same with the growth standard deviation)
  filter(!grepl('surv', param), !param %in% 'grow.sigma')
  
reprodct.pert = read.csv('03_construct_kernels/out/deterministic_repr_coef_perturbation_phen.csv')

# For the sake of ease, I'm going to remove the drought-irrigation phen crosses (not necessary)
reprodct.pert = reprodct.pert %>% 
  filter(!(trt %in% 'drought' & floor(mean.phen) == 127)) %>%
  filter(!(trt %in% 'irrigated' & floor(mean.phen) == 122))

reprodct.ltre = reprodct.ltre %>%
  filter(!(trt %in% 'drought' & floor(phen) == 127)) %>%
  filter(!(trt %in% 'irrigated' & floor(phen) == 122))

# Get all of the perturbed kernels together

all.perts = rbind(
  merge(
    growsurv %>% select(-c(pred.surv, p.grow.size, pred.grow.mean)) %>% rename(size.nex = size.cur), 
    reprodct.pert,
    by = c('size.prev', 'size.nex', 'trt'),
    suffixes = c('.g', '.f')
  ) %>%
    select(size.prev, size.nex, trt, p.size.cur.g, p.size.cur.f, param, orig.parval, mean.phen),
  merge(
    growsurv.pert, 
    reprodct.ltre %>% rename(mean.phen = phen),
    by = c('size.prev', 'size.nex', 'trt'),
    suffixes = c('.g', '.f')
  ) %>%
    select(size.prev, size.nex, trt, p.size.cur.g, p.size.cur.f, param, orig.parval, mean.phen)
) %>%
  mutate(p.size.cur = p.size.cur.g + p.size.cur.f * p.germ) %>%
  select(-c(p.size.cur.g, p.size.cur.f))

pert.matr = split(all.perts, all.perts[,c("trt", "param", "mean.phen")], sep = '_', drop = TRUE) %>%
  lapply(
    function(df) {
      df %>%
        arrange(size.prev, size.nex) %>%
        pivot_wider(names_from = size.prev, values_from = p.size.cur) %>%
        select(-c(trt, param, mean.phen, size.nex, orig.parval)) %>%
        as.matrix()
    }
  )

length(pert.matr)

pert.lambda = pert.matr %>%
  sapply(function(m) Re(eigen(m)$values[1])) %>%
  data.frame(lambda = .) %>%
  mutate(trt_phen = row.names(.)) %>%
  separate(trt_phen, into = c('trt', 'param', 'phen'), sep = '_') %>%
  mutate(phen = as.numeric(phen))

head(pert.lambda)

# Merge together to estimate sensitivities

sensitivities = merge(
  ltre.lambda %>% select(-phen.date),
  pert.lambda,
  by = c('trt', 'phen'),
  all.x = FALSE, all.y = TRUE,
  suffixes = c('.orig', '.pert')
) %>%
  mutate(d.lambda = (lambda.pert - lambda.orig) / .0001)


# Get vital rate differences
# I'm going to separate out the phen here...
param.diffs = rbind(growsurv.pert, reprodct.pert %>% select(-mean.phen)) %>%
  filter(!grepl('phen', param)) %>%
  distinct(trt, param, orig.parval) %>%
  pivot_wider(names_from = trt, values_from = orig.parval) %>%
  mutate(d.c = drought - control, i.c = irrigated - control) %>%
  select(-c(drought, irrigated, control)) %>%
  pivot_longer(c(d.c, i.c), names_to = 'contrast', values_to = 'pardiff')

param.diffs

# Now, need to sensitivites from mean matrices

# First step: get mean matrices for perturbed and unperturbed

mean.lambdas = 1

mean.kernels = ltre.kernel.df %>%
  # (remove unnecessary phen-trt combos)
  filter(!(trt %in% 'drought' & floor(phen) == 127)) %>%
  filter(!(trt %in% 'irrigated' & floor(phen) == 122)) %>%
  pivot_wider(names_from = trt, values_from = p.size.cur) %>%
  mutate(d.c = (drought + control) / 2, i.c = (irrigated + control) / 2) %>%
  select(-c(drought, control, irrigated)) %>%
  pivot_longer(c(d.c, i.c), names_to = 'contrast', values_to = 'mentry') %>%
  # (ah - some of these contrasts don't exist because of the now-imbalanced design)
  # (filter out NAs)
  filter(!is.na(mentry))

mean.pert.kernels = all.perts %>%
  select(-orig.parval) %>%
  pivot_wider(names_from = trt, values_from = p.size.cur) %>%
  mutate(d.c = (drought + control) / 2, i.c = (irrigated + control) / 2) %>%
  select(-c(drought, control, irrigated)) %>%
  pivot_longer(c(d.c, i.c), names_to = 'contrast', values_to = 'mentry') %>%
  # filter out NAs
  filter(!is.na(mentry))

mean.lambdas = split(
  mean.kernels,
  mean.kernels[,c("contrast", "phen")],
  sep = '_', drop = TRUE
) %>%
  lapply(
    function(df) {
      df %>%
        arrange(size.prev, size.cur) %>%
        pivot_wider(names_from = size.prev, values_from = mentry) %>%
        select(-c(contrast, phen, size.cur)) %>%
        as.matrix()
    }
  ) %>%
  sapply(function(x) Re(eigen(x)$values[1])) %>%
  data.frame(lambda = .) %>%
  mutate(c_pg = row.names(.)) %>%
  separate(c_pg, into = c('contrast', 'phen'), sep = '_')  

mean.pert.lambdas = split(
  mean.pert.kernels,
  mean.pert.kernels[,c("contrast", 'param', "mean.phen")],
  sep = '_', drop = TRUE
) %>%
  lapply(
    function(df) {
      df %>%
        arrange(size.prev, size.nex) %>%
        pivot_wider(names_from = size.prev, values_from = mentry) %>%
        select(-c(contrast, mean.phen, param, size.nex)) %>%
        as.matrix()
    }
  ) %>%
  sapply(function(x) Re(eigen(x)$values[1])) %>%
  data.frame(lambda = .) %>%
  mutate(c_pg = row.names(.)) %>%
  separate(c_pg, into = c('contrast', 'param', 'phen'), sep = '_')  

mean.sensitivities = merge(
  mean.lambdas, mean.pert.lambdas,
  by = c('contrast', 'phen'), suffixes = c('.orig', '.pert')
) %>%
  mutate(sv = (lambda.pert - lambda.orig) / .0001) %>%
  mutate(phen = as.numeric(phen))

head(mean.sensitivities)
str(mean.sensitivities)

# Now, merge and combine

coef.ltre = merge(
  param.diffs, mean.sensitivities %>% select(-c(lambda.orig, lambda.pert)),
  by = c('contrast', 'param')
) %>%
  mutate(contribution = pardiff * sv)

head(coef.ltre)

coef.ltre.by.rate = coef.ltre %>%
  separate(param, into = c('rate', 'param'), sep = '\\.') %>%
  group_by(contrast, rate, phen) %>%
  summarise(contribution = sum(contribution))

coef.ltre.by.rate %>%
  ggplot(aes(x = rate, y = contribution, fill = contrast)) +
  geom_col() +
  scale_fill_manual(values = c('red', 'blue')) +
  facet_wrap(phen ~ contrast)
# (unsurprisingly, very similar!)

# Okay... let's now compare with the actual lambda differences...

ltre.compare = ltre.lambda %>%
  filter(!(trt %in% 'drought' & floor(phen) == 127)) %>%
  filter(!(trt %in% 'irrigated' & floor(phen) == 122)) %>%
  select(-c(phen.date, obs.phen)) %>%
  pivot_wider(names_from = trt, values_from = lambda) %>%
  mutate(d.c = drought - control, i.c = irrigated - control) %>%
  select(-c(control, drought, irrigated)) %>%
  pivot_longer(c(d.c, i.c), names_to = 'contrast', values_to = 'd.lambda') %>%
  filter(!is.na(d.lambda)) %>%
  merge(coef.ltre.by.rate %>% group_by(contrast, phen) %>% summarise(sum.contrib = sum(contribution)))

ltre.compare
# HAHA HA YES, YES!!

# Next... phen LTRE step

# phen.compare = reprodct.pert %>%
#   distinct(trt, mean.phen) %>%
#   mutate(
#     obs.phen = (trt %in% 'control' & floor(mean.phen) == 125) |
#       (trt %in% 'drought' & floor(mean.phen) == 122) | (trt %in% 'irrigated' & floor(mean.phen) == 127),
#   ) %>%
#   group_by(trt) %>%
#   reframe(phen.diff = diff(mean.phen), sign = diff(obs.phen))

# Differences in lambda and phenology by treatment
# going to do this in an ugly but effective way

phen.lambda.diff = ltre.lambda %>% 
  filter(!(trt %in% 'drought' & floor(phen) %in% 127), !(trt %in% 'irrigated' & floor(phen) %in% 122)) %>%
  select(-phen.date) %>%
  uncount(weights = 1 + as.numeric(trt %in% 'control' & obs.phen)) %>%
  group_by(trt) %>%
  mutate(
    # This column will let me do this easily
    ii = floor((0:(n()-1)) / 2),
    # formatting this int a  readable string for pivoting
    obs.phen = ifelse(obs.phen, 'observed', 'hypothetical')
  ) %>%
  ungroup() %>%
  pivot_wider(id_cols = c(trt, ii), names_from = obs.phen, values_from = c(lambda, phen)) %>%
  mutate(
    lambda.diff = lambda_hypothetical - lambda_observed,
    phen.diff = phen_hypothetical - phen_observed 
  ) %>%
  select(trt, ii, phen.diff, lambda.diff)

# Okay... now, get mean matrices

# # Testing some code below
# testy = all.perts %>%
#   mutate(
#     obs.phen = (trt %in% 'control' & floor(mean.phen) == 125) |
#       (trt %in% 'drought' & floor(mean.phen) == 122) | (trt %in% 'irrigated' & floor(mean.phen) == 127)
#   ) %>%
#   # Duplicating the control @ control phen entries to get 
#   uncount(weight = 1 + as.numeric(trt %in% 'control' & obs.phen)) %>%
#   arrange(mean.phen) %>%
#   select(-c(orig.parval, mean.phen)) %>%
#   group_by(size.prev, size.nex, trt, param, obs.phen) %>%
#   mutate(dupl = duplicated(p.size.cur)) %>%
#   group_by(size.prev, size.nex, trt, param) %>%
#   mutate(ii = cumsum(dupl)) %>%
#   ungroup()
# 
# # these should all be control...
# testy %>% filter(ii > 0) %>% distinct(trt)
# # sweet
# # should be some values here where obs.phen is false
# testy %>% filter(ii > 0) %>% distinct(obs.phen)
# # sweet
# # going to assume this is good

mean.phen.perts = all.perts %>%
  filter(grepl('phen', param)) %>%
  mutate(
    obs.phen = (trt %in% 'control' & floor(mean.phen) == 125) |
      (trt %in% 'drought' & floor(mean.phen) == 122) | (trt %in% 'irrigated' & floor(mean.phen) == 127)
  ) %>%
  # Duplicating the control @ control phen entries to get 
  uncount(weight = 1 + as.numeric(trt %in% 'control' & obs.phen)) %>%
  arrange(mean.phen) %>%
  group_by(size.prev, size.nex, trt, param, mean.phen) %>%
  mutate(dupl = duplicated(p.size.cur)) %>%
  group_by(size.prev, size.nex, trt, param) %>%
  mutate(ii = cumsum(dupl)) %>%
  ungroup() %>%
  select(-c(dupl, orig.parval, mean.phen)) %>%
  # Okay - now for pivoting
  # going to change the obs.phen column into something more legible
  mutate(obs.phen = ifelse(obs.phen, 'observed', 'hypothetical')) %>%
  pivot_wider(names_from = obs.phen, values_from = p.size.cur) %>%
  # Get the mean of each entry
  mutate(mentry = (observed + hypothetical) / 2) %>%
  select(-c(hypothetical, observed))

mean.phen.pert.lambdas = split(
  mean.phen.perts,
  mean.phen.perts[,c("trt", "param", "ii")],
  sep = '_', drop = TRUE
) %>%
  lapply(
    function(df) {
      df %>%
        arrange(size.prev, size.nex) %>%
        pivot_wider(names_from = size.prev, values_from = mentry) %>%
        select(-c(trt, ii, param, size.nex)) %>%
        as.matrix()
    }
  ) %>%
  sapply(function(x) Re(eigen(x)$values[1])) %>%
  data.frame(lambda = .) %>%
  mutate(tpi = row.names(.)) %>%
  separate(tpi, into = c('trt', 'rate', 'ii'), sep = '_')

# Now get means of unperturbed lambdas

mean.phen = ltre.kernel.df %>%
  # Get unnecessary phen-trt combos
  filter(!(trt %in% 'drought' & floor(phen) == 127)) %>%
  filter(!(trt %in% 'irrigated' & floor(phen) == 122)) %>%
  mutate(
    obs.phen = (trt %in% 'control' & floor(phen) == 125) |
      (trt %in% 'drought' & floor(phen) == 122) | (trt %in% 'irrigated' & floor(phen) == 127)
  ) %>%
  # Duplicating the control @ control phen entries to get 
  uncount(weight = 1 + as.numeric(trt %in% 'control' & obs.phen)) %>%
  arrange(phen) %>%
  group_by(size.prev, size.cur, trt, phen) %>%
  mutate(dupl = duplicated(p.size.cur)) %>%
  group_by(size.prev, size.cur, trt) %>%
  mutate(ii = cumsum(dupl)) %>%
  ungroup() %>%
  select(-c(dupl, phen)) %>%
  # Okay - now for pivoting
  # going to change the obs.phen column into something more legible
  mutate(obs.phen = ifelse(obs.phen, 'observed', 'hypothetical')) %>%
  pivot_wider(names_from = obs.phen, values_from = p.size.cur) %>%
  # Get the mean of each entry
  mutate(mentry = (observed + hypothetical) / 2) %>%
  select(-c(hypothetical, observed)) 

mean.phen.lambdas = split(
  mean.phen,
  mean.phen[,c("trt", "ii")],
  sep = '_', drop = TRUE
) %>%
  lapply(
    function(df) {
      df %>%
        arrange(size.prev, size.cur) %>%
        pivot_wider(names_from = size.prev, values_from = mentry) %>%
        select(-c(trt, ii, size.cur)) %>%
        as.matrix()
    }
  ) %>%
  sapply(function(x) Re(eigen(x)$values[1])) %>%
  data.frame(lambda = .) %>%
  mutate(ti = row.names(.)) %>%
  separate(ti, into = c('trt', 'ii'), sep = '_')
 
phen.sensitivities = merge(
  mean.phen.lambdas, mean.phen.pert.lambdas,
  by = c('trt', 'ii'),
  suffixes = c('.orig', '.pert')
) %>%
  mutate(sv = (lambda.pert - lambda.orig) / .0001)

phen.sensitivities

coef.phen.ltre = merge(phen.sensitivities, phen.lambda.diff) %>%
  # somewhere there's a sign change... fix later
  mutate(contribution = -1 * phen.diff * sv)

# Compare
coef.phen.ltre %>% group_by(trt, ii) %>%
  summarise(
    lambda.diff = mean(lambda.diff),
    total.contr = sum(contribution)
  )
# Hell. Fucking. Yeah.

### -----

# Make an LTRE figure

ltre.all.contrs = rbind(
  coef.ltre.by.rate %>%
    filter(floor(phen) == 125) %>%
    select(-phen) %>%
    mutate(ltre.var = 'alpha') %>%
    rename(trt = contrast) %>%
    mutate(trt = ifelse(trt %in% 'd.c', 'drought', 'irrigated')),
  coef.phen.ltre %>%
    select(trt, rate, contribution) %>%
    filter(!trt %in% 'control') %>%
    mutate(rate = gsub('phen\\.', '', rate)) %>%
    mutate(ltre.var = 'beta')
)

ltre.all.contrs %>%
  mutate(
    rate = case_match(
      rate,
      'flow' ~ 'flowering.prob',
      'grow' ~ 'growth',
      'seed' ~ 'seed.production',
      'succ' ~ 'umbel.succ',
      'recr' ~ 'recruit.size'
    )
  ) %>%
  mutate(contr.label = paste0(ltre.var, '[', rate, ']')) %>%
  mutate(trt = paste(trt, 'vs. control')) %>%
  ggplot(aes(x = contr.label, y = contribution)) +
  geom_col_pattern(
    aes(fill = trt, pattern_spacing = ltre.var),
    pattern_fill = 'gray22', pattern_density = 0.05
  ) +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_pattern_spacing_manual(values = c(0.05, 0.1)) +
  scale_fill_manual(values = c('red', 'blue')) +
  labs(x = 'variable', y = 'contribution') +
  facet_wrap(~ trt, nrow = 2) +
  theme(
    legend.position = 'none',
    panel.background = element_blank()
  )

ggsave('04_analysis/figures/draft_figures/ltre_2024-08-07.png', width = 5, height = 5)

all.lambda %>%
  mutate(phen.date = as.Date(as.numeric((phen)), format = '%b-%d')) %>%
  ggplot(aes(x = phen.date, y = lambda, group = trt, colour = trt)) +
  geom_line() +
  scale_colour_manual(values = c('black', 'red', 'blue')) +
  geom_point(
    data = ltre.lambda %>% filter(obs.phen | floor(phen) == 125), 
    aes(shape = obs.phen), 
    size = 4
  ) +
  scale_shape_manual(values = c(1, 19)) +
  guides(shape = 'none') +
  labs(x = 'bud date') +
  theme(
    panel.background = element_blank(),
    legend.position = 'top'
  )

ggsave('04_analysis/figures/draft_figures/lambdas_2024-08-07.png', width = 5, height = 5)
