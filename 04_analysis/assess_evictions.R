# Script for assessing eviction (see: Williams, Miller, and Ellner 2012) in our kernels
# (Much of code here is borrowed from other scripts)
# Description of appraoch here:
# https://figshare.com/articles/dataset/Appendix_A_Derivation_of_d_/3554163?file=5623080
# used for deterministic no-phen kernels

library(ggplot2)
library(tidyr)
library(dplyr)

rm(list = ls())

# ------ Read in data -----

# Processed kernel data:
growsurv = read.csv('03_construct_kernels/out/deterministic_growsurv_kernel.csv')
reprod = read.csv('03_construct_kernels/out/deterministic_reprod_kernel_no_phen.csv')

# Demographic data (for growth and survival only - need for fitting models)
all.data = merge(
  x = read.csv('01_data_cleaning/out/demo_phen_seed_2016-2024_final.csv'),
  y = read.csv('00_raw_data/plot_treatments.csv'),
  by.x = 'Plot', by.y = 'plot'
)

# ----- Process data -----

# Kernel first:
# Dataframe of rates
all.rates = merge(
  x = growsurv, y = reprod,
  by.x = c('size.prev', 'size.cur', 'trt'), by.y = c('size.pre', 'size.cur', 'trt'),
  suffixes = c('.g', '.f')
) %>%
  merge(
    y = expand.grid(
      trt = c('control', 'drought', 'irrigated'),
      p.germ = c(0.001, 0.005, 0.01, 0.05)
    ), 
  by = 'trt'
)

# Kernel
kern.df = all.rates %>%
  select(size.prev, size.cur, trt, p.size.cur.g, p.size.cur.f, p.germ) %>%
  mutate(p.size.cur.f = p.size.cur.f * p.germ) %>%
  mutate(p.size.cur = p.size.cur.g + p.size.cur.f) %>%
  select(-c(p.size.cur.g, p.size.cur.f))

# Eigenvectors
all.eigs = split(kern.df, kern.df[,c("trt", "p.germ")], sep = '_') %>%
  lapply(
    function(df) {
      df %>%
        arrange(size.prev, size.cur) %>%
        pivot_wider(names_from = size.prev, values_from = p.size.cur) %>%
        select(-c(trt, p.germ, size.cur))
    }
  ) %>%
  lapply(
    function(x)
      data.frame(
        w = Re(eigen(x)$vectors[,1]),
        v = Re(eigen(t(x))$vectors[,1])
      )
  )


# Now prepare the demographic data
# (needed for vital rate models)

all.demo = all.data %>% 
  filter(in.demo) %>%
  distinct(plantid, Year, .keep_all = TRUE)

# Dataset for survival models
demo.surv = merge(
  # Demo in time step t+1
  x = all.demo %>% 
    mutate(prev.year = Year - 1) %>%
    rename(surv.year = Year) %>%
    select(Plot, plantid, surv.year, prev.year, No.leaves, Leaf.length, surv, trt),
  # Demo in time step t
  y = all.demo %>%
    # we are *only* interested in plants alive in time step t
    filter(surv) %>%
    # Select relevant columns
    select(Plot, plantid, Year, No.leaves, Leaf.length, trt),
  by.x = c('Plot', 'plantid', 'prev.year', 'trt'),
  by.y = c('Plot', 'plantid', 'Year', 'trt'),
  suffixes = c('', '.pre'),
  all.x = FALSE, all.y = FALSE
) %>%
  filter(
    !is.na(Leaf.length.pre) & !is.na(No.leaves.pre) &
      Leaf.length.pre > 0 & No.leaves.pre > 0
  ) %>%
  # Get rid of 2016 records because the sizes are not reliable
  filter(prev.year > 2016) %>%
  # Add size columns
  mutate(size.prev = log(No.leaves.pre * Leaf.length.pre))

# Dataset for growth model
demo.grow = demo.surv %>% 
  filter(surv, !is.na(Leaf.length) & !is.na(No.leaves) & Leaf.length > 0 & No.leaves > 0) %>%
  mutate(size.cur = log(Leaf.length * No.leaves))

# Remove 2024 from survival dataset
demo.surv = demo.surv %>% filter(surv.year < 2024)

# Dataset for recruitment model
demo.recr = all.demo %>%
  arrange(Year) %>%
  # Want demo records (to get size distribution)
  filter(in.demo) %>%
  # Give me only the first sighting for each plant
  distinct(plantid, .keep_all = TRUE) %>%
  # Give me only records after 2017 (these are more likely to be seen for the
  # first time)
  filter(Year > 2017) %>%
  # Give me only records with leaf counts and lengths (for estimating size)
  filter(!is.na(No.leaves), !is.na(Leaf.length)) %>%
  # Give me plants that did not flower (plants unlikely to flower in first year)
  filter(!No.umbels & !is.na(No.umbels)) %>%
  # Okay... going to assume that only plants with one leaf are new recruits
  filter(No.leaves < 2) %>%
  mutate(size = log(Leaf.length))

# ----- Fit vital rate models -----

# Survival model
s_s = glmmTMB(
  formula = surv ~ size.prev + (1 | Plot),
  family = 'binomial',
  data = demo.surv
)

# Growth model
g_st = glmmTMB(
  size.cur ~ size.prev * trt + (1 | prev.year) + (1 | Plot / plantid),
  data = demo.grow
)

# New recruit size model
r_t.y = glmmTMB(size ~ trt + (1 | Year) + (1 | Plot), data = demo.recr)


# ----- Get data frame with rho_L and rho_U -----
# (see supporting material)

eviction.probs = expand.grid(
  size.prev = (5:60)/10,
  size.cur = (c(-10:4, 61:200))/10,
  trt = c('control', 'drought', 'irrigated')
) %>%
  # Predicted survival
  mutate(
    pred.surv = predict(
      newdata = .,
      object = s_s, type = 'response',
      re.form = ~ 0, allow.new.levels = TRUE
    )
  ) %>%
  # Predicted growth
  mutate(
    pred.grow.mean = predict(
      newdata = .,
      object = g_st, type = 'response',
      re.form = ~ 0, allow.new.levels = TRUE
    )
  ) %>%
  # Predicted distribution of sizes in next time step
  mutate(p.grow.size = 0.1 * dnorm(size.cur, mean = pred.grow.mean, sd = summary(g_st)$sigma)) %>%
  # Combine all together to get overall size distribution in next time step
  mutate(p.size.cur = pred.surv * p.grow.size) %>%
  group_by(trt, size.prev, low.hi = ifelse(size.cur < .5, 'l', 'u')) %>%
  summarise(rho.u = sum(p.size.cur)) %>%
  ungroup()

# Formula in SI:
# dlambda_U = v(U) * <rho(U), w> / <v, w>
# dlambda_L = v(L) * <rho(L), w> / <v, w>

# Okay... not sure what the best way to vectorize this is...
eviction.list = eviction.probs %>%
  pivot_wider(names_from = low.hi, values_from = rho.u) %>% 
  split(.$trt) %>% 
  lapply(function(x) x %>% select(l, u)) %>%
  rep(times = 4)

# Not super elegant but we can do it like this
data.frame(a = names(all.eigs), b = names(eviction.list))

for (i in 1:length(all.eigs)) {
  all.eigs[[i]]$l = eviction.list[[i]]$l
  all.eigs[[i]]$u = eviction.list[[i]]$u
}
# (is there a way to merge lists...?)

d.lambda = sapply(
  all.eigs,
  function(x) {
    dl_U = with(x, v[length(v)] * sum(u * w) / sum(v * w))
    dl_L = with(x, v[1] * sum(l * w) / sum(v * w))
    return(c(dl_U = dl_U, dl_L = dl_L))
  }
) %>%
  t() %>%
  data.frame() %>%
  mutate(ctrl.pgerm = row.names(.)) %>%
  separate(ctrl.pgerm, into = c('trt', 'p.germ'), sep = '_')

d.lambda %>%
  pivot_longer(c(dl_U, dl_L), names_to = 'LU', values_to = 'dl') %>%
  ggplot(aes(x = p.germ, y = dl, group = trt)) +
  geom_line(aes(colour = trt)) +
  scale_colour_manual(values = c('black', 'red', 'blue')) +
  facet_wrap(~ LU)

# Looks like the max here (for higher germination probabilities) is 1e-4

# Check of "eviction" caused by lower boundary and recruit size distribution
predict(
  r_t.y, 
  newdata = data.frame(trt = c('control', 'drought', 'irrigated')), 
  allow.new.levels = TRUE, re.form = ~ 0) %>% 
  pnorm(.5, mean = ., sd = summary(r_t.y)$sigma)
# tiny - ~10^-6
