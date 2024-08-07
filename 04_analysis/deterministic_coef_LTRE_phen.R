
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
    obs.phen = !duplicated(phen),
    phen = as.Date(as.numeric(phen), format = '%b-%d')
  )

head(all.lambda)
head(ltre.lambda)

all.lambda %>%
  mutate(phen = as.Date(as.numeric((phen)), format = '%b-%d')) %>%
  ggplot(aes(x = phen, y = lambda, group = trt, colour = trt)) +
  geom_line() +
  scale_colour_manual(values = c('black', 'red', 'blue')) +
  geom_point(data = ltre.lambda, aes(shape = obs.phen), size = 4) +
  scale_shape_manual(values = c(1, 19)) +
  guides(shape = 'none') +
  theme(panel.background = element_blank())
