# Script for assessing eviction (see: Williams, Miller, and Ellner 2012) in our kernels
# (Much of code here is borrowed from other scripts)
# Description of approach here:
# https://figshare.com/articles/dataset/Appendix_A_Derivation_of_d_/3554163?file=5623080
# used for full kernels

library(ggplot2)
library(tidyr)
library(dplyr)

rm(list = ls())

# ------ Growth kernels outside boundaries -----

if (!file.exists('04_analysis/out/eviction_growsurv_subkernel.csv')) { 
  # If the growth kernel read-in file does not exist,
  # generate it (and export it)
  
  source('03_construct_kernels/prepare_demo_data_growsurv.R')
  
  # Survival model
  s_s = glmmTMB(
    formula = surv ~ size.prev + (1 | Plot),
    family = 'binomial',
    data = demo.surv.sizes
  )
  
  # Growth model
  g_st.ty = glmmTMB(
    size.cur ~ size.prev + size.prev * trt + (1 | prev.year) + (1 | prev.year:trt) + (1 | Plot / plantid),
    data = demo.grow
  )
  
  # Residual variance in growth models
  gv.sd = summary(g_st.ty)$sigma
  
  # Generate a kernel
  grow.surv.kernel = expand.grid(
    size.prev = (5:60)/10,
    size.cur = (c(-10:4, 61:200))/10,
    trt = 'control'
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
      # Model without phenology
      pred.grow.mean = predict(
        newdata = .,
        object = g_st.ty, type = 'response',
        re.form = ~ 0, allow.new.levels = TRUE
      )
    ) %>%
    # Predicted distribution of sizes in next time step
    mutate(pv.grow.size = 0.1 * dnorm(size.cur, mean = pred.grow.mean, sd = gv.sd))
  
  write.csv(
    grow.surv.kernel, row.names = FALSE, 
    file = '04_analysis/out/eviction_growsurv_subkernel.csv'
  )
  
}

# ------ Read in all kernel data -----

# Our subkernels:
growsurv = read.csv('03_construct_kernels/out/deterministic_growsurv_kernel_phen_ltre.csv') %>%
  # Give us only the control kernels
  filter(trt.phen.idx %in% 7) %>%
  select(-trt.phen.idx)
reprod = read.csv('03_construct_kernels/out/deterministic_reprod_kernel_phen_ltre.csv') %>%
  # Give us only the control kernels
  filter(trt.phen.idx %in% 7) %>%
  select(-trt.phen.idx)

# Eviction subkernels
ev_growsurv = read.csv('04_analysis/out/eviction_growsurv_subkernel.csv')


# ----- Process data (build kernels) -----

# Set a germination rate
p.germ = 0.004632985

# Kernel first:
# Dataframe of rates
all.kernels = merge(
  x = growsurv, y = reprod,
  by.x = c('size.prev', 'size.cur'), 
  by.y = c('size.prev', 'size.nex'),
) %>%
  mutate(
    p.size.cur = pred.surv * (pv.grow.size * (1 - prob.flower) + pf.grow.size * prob.flower) + (p.size.cur * p.germ)
  ) %>%
  select(-c(pred.surv, pv.grow.size, pf.grow.size, prob.flower)) %>%
  # Sort df rows/kernel entries (just in case)
  arrange(size.prev, size.cur)

# Eigenvectors for our main kernel
all.eigs = all.kernels %>%
  (
    function(df) {
      df %>%
        arrange(size.prev, size.cur) %>%
        pivot_wider(names_from = size.prev, values_from = p.size.cur) %>%
        select(-c(size.cur))
    }
  ) %>%
  (
    function(x)
      data.frame(
        l = Re(eigen(x)$values[1]),
        w = Re(eigen(x)$vectors[,1]),
        v = Re(eigen(t(x))$vectors[,1])
      )
  )

(lambda = all.eigs$l[1])

# ----- Get data frame with rho_L and rho_U -----
# (see supporting material)

eviction.probs = ev_growsurv %>%
  # Get the probability of 
  mutate(p.size.cur = pred.surv * pv.grow.size) %>%
  select(size.prev, size.cur, p.size.cur) %>%
  # Sort df rows/kernel entries (just in case)
  arrange(size.prev, size.cur) %>%
  # Break down transition probabilities into below/above kernel dimensions
  group_by(size.prev, low.hi = ifelse(size.cur < 0.5, 'l', 'u')) %>%
  # take the sum and multiply by 0.1 (the binwidth; dx in the Riemann sum)
  summarise(rho = sum(p.size.cur) * 0.1) %>%
  ungroup()

# Formula in SI:
# dlambda_U = v(U) * <rho(U), w> / <v, w>
# dlambda_L = v(L) * <rho(L), w> / <v, w>

# Pivot to wide form to get l and u as vectors (cols in data frame)
eviction.probs = eviction.probs %>%
  pivot_wider(names_from = low.hi, values_from = rho, names_prefix = 'rho_') %>% 
  select(size.prev, rho_l, rho_u) %>%
  arrange(size.prev)

# all.eigs (the eigenvectors) should already be sorted by size.prev
all.eigs = cbind(eviction.probs, w = all.eigs$w, v = all.eigs$v)

dl.vals = all.eigs %>%
  arrange(size.prev) %>%
  mutate(
    vstar = v + (1 / lambda) * (last(v)*rho_u + first(v)*rho_l)
  ) %>%
  summarise(
    dl_U = last(vstar) * sum(rho_u * w) / sum(vstar * w),
    dl_L = first(vstar) * sum(rho_l * w) / sum(vstar * w)
  )

dl.vals
#           dl_U         dl_L
# 1 8.956528e-07 1.462964e-07
sum(dl.vals)

# Check of "eviction" caused by lower boundary and recruit size distribution

source('03_construct_kernels/prepare_demo_data_repr.R')

r_t.y = glmmTMB(size ~ trt + (1 | Year) + (1 | Plot), data = demo.recr)

predict(
  r_t.y, 
  newdata = data.frame(trt = c('control', 'drought', 'irrigated')), 
  allow.new.levels = TRUE, re.form = ~ 0
) %>% 
  pnorm(.5, mean = ., sd = summary(r_t.y)$sigma)
# tiny less than ~10^-6 in all cases
