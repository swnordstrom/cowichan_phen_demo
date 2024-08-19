# Script for doing deterministic LTRE by vital rate regression coefficient
# perturbation with no phenology
# Mostly this is being done as a proof of concept (making sure sums are correct)
# etc.

### ---------------------------------

library(ggplot2)
library(dplyr)
library(tidyr)

rm(list = ls())

# Observed vitals
growsurv = read.csv('03_construct_kernels/out/deterministic_growsurv_kernel.csv')
reprodct = read.csv('03_construct_kernels/out/deterministic_reprod_kernel_no_phen.csv')

# Outputs with perturbations
growsurv.pert = read.csv('03_construct_kernels/out/deterministic_grow_coef_perturbation_no_phen.csv')
reprodct.pert = read.csv('03_construct_kernels/out/deterministic_repr_coef_perturbation_no_phen.csv')

# Get germination rates
trt.germ = expand.grid(
  trt = c('control', 'drought', 'irrigated'),
  p.germ = c(0.001, 0.005, 0.01)
)

# Combine to get kernels
all.kernels = merge(
  x = growsurv %>% select(size.prev, size.cur, trt, p.size.cur),
  y = reprodct %>% select(size.pre , size.cur, trt, p.size.cur),
  by.x = c('size.prev', 'size.cur', 'trt'),
  by.y = c('size.pre' , 'size.cur', 'trt'),
  suffixes = c('.g', '.f')
) %>%
  merge(trt.germ, all.y = TRUE) %>%
  mutate(p.size.cur = p.size.cur.g + p.size.cur.f * p.germ) %>%
  select(-c(p.size.cur.g, p.size.cur.f))

head(all.kernels)

# Get perturbation matrices
all.pert.kernels = rbind(
  merge(
    x = growsurv.pert, 
    y = reprodct %>% select(size.pre, size.cur, trt, p.size.cur),
    by.x = c('size.prev', 'size.cur', 'trt'),
    by.y = c('size.pre' , 'size.cur', 'trt'),
    suffixes = c('.g', '.f')
  ),
  merge(
    x = growsurv %>% select(size.prev , size.cur, trt, p.size.cur),
    y = reprodct.pert,
    by.x = c('size.prev', 'size.cur', 'trt'),
    by.y = c('size.pre' , 'size.cur', 'trt'),
    suffixes = c('.g', '.f')
  )
) %>%
  merge(trt.germ, all.y = TRUE) %>%
  mutate(p.size.cur = p.size.cur.g + p.size.cur.f * p.germ) %>%
  select(-c(p.size.cur.g, p.size.cur.f)) %>%
  mutate(perturb.param = gsub('\\_', '.', perturb.param))

head(all.pert.kernels)

# Split up to make into matrices

all.matr = split(all.kernels, all.kernels[,c("trt", "p.germ")], sep = '_') %>%
  lapply(
    function(df) {
      df %>%
        arrange(size.prev, size.cur) %>%
        pivot_wider(names_from = size.prev, values_from = p.size.cur) %>%
        select(-c(trt, p.germ, size.cur)) %>%
        as.matrix()
    }
  )

all.lambda = all.matr %>% 
  sapply(function(m) Re(eigen(m)$values[1])) %>%
  data.frame(lambda = .) %>%
  mutate(trt_p.germ = row.names(.)) %>%
  separate(trt_p.germ, into = c('trt', 'p.germ'), sep = '_')

all.lambda %>%
  ggplot(aes(x = p.germ, y = lambda, colour = trt, group = trt)) +
  geom_line() +
  geom_point() +
  scale_colour_manual(values = c('black', 'red', 'blue'))

all.pert.matr = split(all.pert.kernels, all.pert.kernels[,c("trt", "p.germ", "perturb.param")], sep = '_') %>%
  lapply(
    function(df) {
      df %>%
        arrange(size.prev, size.cur) %>%
        pivot_wider(names_from = size.prev, values_from = p.size.cur) %>%
        select(-c(trt, p.germ, perturb.param, size.cur, orig.par.val)) %>%
        as.matrix()
    }
  )

all.pert.lambda = all.pert.matr %>% 
  sapply(function(m) Re(eigen(m)$values[1])) %>%
  data.frame(lambda = .) %>%
  mutate(trt_p.germ = row.names(.)) %>%
  separate(trt_p.germ, into = c('trt', 'p.germ', 'parm'), sep = '_')

head(all.pert.lambda)

# Cool, now we mergie

all.sens = merge(
  all.lambda, all.pert.lambda, 
  by = c('trt', 'p.germ'),
  suffixes = c('.orig', '.pert')
) %>%
  mutate(sv = (lambda.pert - lambda.orig) / .0001)

all.sens %>%
  ggplot(aes(x = parm, y = sv, colour = trt, group = trt)) +
  geom_point(position = position_dodge(width = 1/4)) +
  scale_colour_manual(values = c('black', 'red', 'blue')) +
  facet_wrap(~ p.germ, nrow = 3)

# ah... sign error here. the zero-inflation terms need to be done differently
# (ah maybe just keep it that way - signs will cancel out appropriately in LTRE)

### -----------
# Setting up the LTRE now

# Want to get vital rate differences
coef.diffs = rbind(growsurv.pert, reprodct.pert %>% rename(size.prev = size.pre)) %>%
  distinct(trt, perturb.param, orig.par.val) %>%
  pivot_wider(id_cols = perturb.param, names_from = trt, values_from = orig.par.val) %>%
  mutate(d.c = drought - control, i.c = irrigated - control) %>%
  select(-c(control, drought, irrigated)) %>%
  pivot_longer(c(d.c, i.c), names_to = 'contrast', values_to = 'coef.diff') %>%
  mutate(perturb.param = gsub('\\_', '.', perturb.param))

# Vital rate sensitivities

# Estimating lambda at midpoints for perturbed matrices

mean.pert.kernels = all.pert.kernels %>%
  select(-orig.par.val) %>%
  pivot_wider(names_from = trt, values_from = p.size.cur) %>%
  mutate(d.c = (drought + control) / 2, i.c = (irrigated + control) / 2) %>%
  select(-c(control, drought, irrigated)) %>%
  pivot_longer(c(d.c, i.c), names_to = 'contrast', values_to = 'p.size.cur')

mean.pert.lambdas = split(
  mean.pert.kernels, 
  mean.pert.kernels[,c("contrast", "p.germ", "perturb.param")], 
  sep = '_'
) %>%
  lapply(
    function(df) {
      df %>%
        arrange(size.prev, size.cur) %>%
        pivot_wider(names_from = size.prev, values_from = p.size.cur) %>%
        select(-c(contrast, p.germ, perturb.param, size.cur)) %>%
        as.matrix()
    }
  ) %>%
  sapply(function(x) Re(eigen(x)$values[1])) %>%
  data.frame(lambda = .) %>%
  mutate(c_pg_pp = row.names(.)) %>%
  separate(c_pg_pp, into = c('contrast', 'p.germ', 'perturb.param'), sep = '_')

# Estimating lambda at midpoints for unperturbed matrices

mean.kernels = all.kernels %>%
  pivot_wider(names_from = trt, values_from = p.size.cur) %>%
  mutate(d.c = (drought + control) / 2, i.c = (irrigated + control) / 2) %>%
  select(-c(control, drought, irrigated)) %>%
  pivot_longer(c(d.c, i.c), names_to = 'contrast', values_to = 'p.size.cur')

mean.lambdas = split(
  mean.kernels,
  mean.kernels[,c("contrast", "p.germ")],
  sep = '_'
) %>%
  lapply(
    function(df) {
      df %>%
        arrange(size.prev, size.cur) %>%
        pivot_wider(names_from = size.prev, values_from = p.size.cur) %>%
        select(-c(contrast, p.germ, size.cur)) %>%
        as.matrix()
    }
  ) %>%
  sapply(function(x) Re(eigen(x)$values[1])) %>%
  data.frame(lambda = .) %>%
  mutate(c_pg = row.names(.)) %>%
  separate(c_pg, into = c('contrast', 'p.germ'), sep = '_')

mean.sensitivities = merge(
  mean.lambdas, mean.pert.lambdas,
  by = c('contrast', 'p.germ'), suffixes = c('.orig', '.pert')
) %>%
  mutate(sv = (lambda.pert - lambda.orig) / .0001)

head(mean.sensitivities)

# # Above might not have been right...
# # Above was taking mean of perturbed matrices w/in contrast, taking mean of
# # unperturbed matrices w/in contrast, and then comparing those lambdas
# # When instead maybe I should take mean within a contrast and then estimate
# # difference in lambdas between the contrasts...
# # (that's a little less elegant to code I think...)
# 
# # mean.lambdas = merge(
# #   all.kernels, 
# #   all.pert.kernels %>% select(-orig.par.val),
# #   by = c('trt', 'p.germ', 'size.prev', 'size.cur'),
# #   suffixes = c('.orig', '.pert')
# # ) %>%
# #   mutate(p.size.mid = (p.size.cur.orig + p.size.cur.pert) / 2) %>%
# #   select(-c(p.size.cur.orig, p.size.cur.pert)) %>%
# #   split(f = .[,c("trt", "p.germ", "perturb.param")], sep = '_') %>%
# #     lapply(
# #       function(df) {
# #         df %>%
# #           arrange(size.prev, size.cur) %>%
# #           pivot_wider(names_from = size.prev, values_from = p.size.cur) %>%
# #           select(-c(trt, p.germ, perturb.param, size.cur)) %>%
# #           as.matrix()
# #       }
# #     ) %>%
# #   sapply(function(x) Re(eigen(x)$values[1])) %>%
# #   data.frame(lambda = .)
# 
# # wait... this doesn't work either, where would the deltas fit in?

# Okay... now combine the mean sensitivities and the coefficient differences

coef.ltre.df = merge(mean.sensitivities, coef.diffs, by = c('contrast', 'perturb.param')) %>%
  mutate(contrib = coef.diff * sv)

compare.results = coef.ltre.df %>% group_by(contrast, p.germ) %>% summarise(sum.contrib = sum(contrib)) %>%
  merge(
    all.lambda %>%
      pivot_wider(names_from = trt, values_from = lambda) %>%
      mutate(d.c = drought - control, i.c = irrigated - control) %>%
      select(-c(control, drought, irrigated)) %>%
      pivot_longer(c(d.c, i.c), names_to = 'contrast', values_to = 'lambda.diff')
  )

compare.results %>%
  ggplot(aes(x = lambda.diff, y = sum.contrib)) +
  geom_segment(aes(x = 0, xend = 0, y = -.01, yend = .03), linetype = 2) +
  geom_segment(aes(x = -.01, xend = .03, y = 0, yend = 0), linetype = 2) +
  geom_segment(aes(x = -.01, xend = .03, y = -.01, yend = .03), linetype = 3) +
  geom_point(aes(colour = contrast, shape = p.germ), size = 3) +
  scale_colour_manual(values = c('red', 'blue'))
# KING!!!!!!

coef.ltre.df %>% 
  filter(contrib != 0) %>%
  ggplot(aes(x = perturb.param, y = contrib, colour = contrast)) +
  geom_segment(aes(xend = perturb.param, yend = 0)) +
  geom_point(aes(shape = contrib > 0), size = 3) +
  scale_shape_manual(values = c(6, 2), '') +
  scale_colour_manual(values = c('red', 'blue')) +
  facet_wrap(contrast ~ p.germ) +
  theme(
    axis.text.x = element_text(angle = 90),
    legend.position = 'none'
  )

coef.ltre.df %>%
  filter(coef.diff != 0) %>%
  # mutate(coef.sign = coef > 0, sv.sign)
  ggplot(aes(x = perturb.param, y = coef.diff, colour = contrast, fill = contrast)) +
  geom_segment(aes(xend = perturb.param, yend = 0)) +
  geom_point(aes(shape = (coef.diff * sv) > 0), size = 3) +
  scale_colour_manual(values = c('red', 'blue')) +
  scale_fill_manual(values = c('red', 'blue')) +
  scale_shape_manual(values = c(25, 24), '') +
  # facet_wrap(contrast ~ p.germ) +
  facet_wrap(~ contrast) +
  theme(legend.position = 'none')

# Look at our actual lambdas vs matrix mean contrast lambdas
all.lambda %>% 
  ggplot(aes(x = p.germ, y = lambda)) + 
  geom_line(aes(group = trt, colour = trt)) + 
  geom_point(aes(colour = trt), size = 3) + 
  scale_colour_manual(values = c('black', 'red', 'blue')) +
  geom_point(data = mean.lambdas, aes(fill = contrast), shape = 23) +
  scale_fill_manual(values = c('red', 'blue')) + 
  geom_line(data = mean.lambdas, aes(group = contrast), linetype = 2)
# eh I mean the mean lambdas are, unsurprisingly, in between the treatment lambdas!

# still not sure what the issue is
# okay well let's at least get some code to group params into vital rate groupings

coef.ltre.df %>%
  rename(param = perturb.param) %>%
  mutate(
    vital = case_when(
      grepl('seed', param) | grepl('succ', param) | 
        grepl('flow', param) | grepl('numb', param) ~ 'repr',
      grepl('grow', param) | grepl('recr', param) ~ 'grow',
      grepl('surv', param) ~ 'surv'
    )
  ) %>%
  group_by(contrast, p.germ, vital) %>%
  summarise(sum.contrib = sum(contrib)) %>%
  #filter(sum.contrib != 0) %>%
  ggplot(aes(x = p.germ, y = sum.contrib, colour = contrast, shape = vital)) +
  geom_point(position = position_dodge(width = 0.5), size = 4) +
  scale_colour_manual(values = c('red', 'blue'))
