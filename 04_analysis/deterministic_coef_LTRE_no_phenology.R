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


all.pert.matr = split(all.pert.kernels, all.pert.kernels[,c("trt", "p.germ", "perturb.param")], sep = '_') %>%
  lapply(
    function(df) {
      df %>%
        arrange(size.prev, size.cur) %>%
        pivot_wider(names_from = size.prev, values_from = p.size.cur) %>%
        select(-c(trt, p.germ, perturb.param, size.cur)) %>%
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
