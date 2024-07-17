# Deterministic LTRE analysis 
# So far we have entrywise contributions and growth rate-group (G, S, R)
# contributions

# Matrices: should have rows corresponding to size t (size.cur), 
# columns corresponding to size t-1 (size.prev)
# This means pivoting_wider (DF to matrix) should take names_from size.prev
# And on ggplots of the matrices the x should be size.prev, y should be size.cur

# ------------------------------------------------------------------------------

# ----- Setup -----

# Packages
library(ggplot2)
library(tidyr)
library(dplyr)

# Clear namespace
rm(list = ls())

# Observed vital rates
growsurv = read.csv('03_construct_kernels/out/deterministic_growsurv_kernel.csv')
reprod = read.csv('03_construct_kernels/out/deterministic_reprod_kernel_no_phen.csv')

# # Perturbed vital rate estimates for sensitivity analysis
# perturb.grow = read.csv('03_construct_kernels/out/deterministic_grow_perturbation.csv')


# ----- Processing/combining -----

# Matrix of treatment-germination probability combinations
germin = expand.grid(
  trt = c('control', 'drought', 'irrigated'),
  p.germ = c(0.001, 0.005, 0.01, 0.05)
)

# Data frame of all estimated vital rates
all.rates = merge(
  x = growsurv, y = reprod,
  by.x = c('size.prev', 'size.cur', 'trt'), by.y = c('size.pre', 'size.cur', 'trt'),
  suffixes = c('.g', '.f')
) %>%
  merge(germin, by = 'trt')

# Kernels in data frame form - one row per matrix entry per matrix
kern.df = all.rates %>%
  select(size.prev, size.cur, trt, p.size.cur.g, p.size.cur.f, p.germ) %>%
  # Number of establishing individuals per size bin (number of individuals times
  # germination/establishment rate)
  mutate(p.size.cur.f = p.size.cur.f * p.germ) %>%
  # Number of individuals in size class is number that grow into size glass plus
  # number born/establishing into size class
  mutate(p.size.cur = p.size.cur.g + p.size.cur.f) %>%
  select(-c(p.size.cur.g, p.size.cur.f))

head(kern.df)

# Visualize kernels
# (note log scale - otherwise the large individuals blow up the scale)
kern.df %>%
  ggplot(aes(x = size.prev, y = size.cur, fill = log(p.size.cur, base = 10))) +
  geom_raster() +
  scale_y_reverse() +
  scale_fill_viridis_c() +
  facet_wrap(trt ~ p.germ)

# # Data frame of all estimated vital rates
# perturbed.kern.df = growsurv %>%
#   select(-p.size.cur) %>%
#   merge(y = perturb.grow, by = c('size.prev', 'size.cur', 'trt')) %>%
#   merge(
#     y = reprod,
#     by.x = c('size.prev', 'size.cur', 'trt'), by.y = c('size.pre', 'size.cur', 'trt'),
#     suffixes = c('.g', '.f')
# ) %>%
#   merge(germin, by = 'trt') %>%
#   mutate(p.size.cur.f = p.size.cur * p.germ) %>%
#   # Number of individuals in size class is number that grow into size glass plus
#   # number born/establishing into size class
#   mutate(p.size.cur.perturb.growth.mu = p.size.cur.perturb.growth.mu + p.size.cur.f) %>%
#   select(size.prev, size.cur, trt, p.germ, p.size.cur.perturb.growth.mu)

# ----- Generating matrices and matrix products (e.g., eigenvectors)

all.matr = split(kern.df, kern.df[,c("trt", "p.germ")], sep = '_') %>%
  lapply(
    function(df) {
      df %>%
        arrange(size.prev, size.cur) %>%
        pivot_wider(names_from = size.prev, values_from = p.size.cur) %>%
        select(-c(trt, p.germ, size.cur))
    }
  )

length(all.matr)

# Estimate lambda from matrices
all.lambda = sapply(all.matr, function(x) Re(eigen(x)$values[1])) %>%
  data.frame(lambda = .) %>%
  mutate(t.germ = row.names(.)) %>%
  separate(t.germ, into = c('trt', 'p.germ'), sep = '_')

all.lambda %>%
  ggplot(aes(x = p.germ, y = lambda, group = trt, colour = trt)) +
  geom_point(size = 3) +
  geom_line() +
  labs(x = 'probability of germination') +
  scale_colour_manual(values = c('black', 'red', 'blue'))
# LMAOOOOOOOOO

# Seed production
all.rates %>% 
  select(size.prev, trt, seeds.total) %>% 
  ggplot(aes(x = size.prev, y = seeds.total, group = trt, colour = trt)) + 
  geom_line() + 
  scale_y_log10() +
  labs(x = 'plant size', y = 'seeds per plant') +
  scale_colour_manual(values = c('black', 'red', 'blue'))
# Looks like control plants produce more seeds *when they are large*
# Irrigated plants produce more seeds when small...

# Growth
all.rates %>%
  filter(p.germ < .005) %>%
  select(size.prev, pred.grow.mean, trt) %>%
  ggplot(aes(x = size.prev, y = pred.grow.mean, group = trt, colour = trt)) + 
  geom_segment(aes(x = .5, xend = 6, y = .5, yend = 6), inherit.aes = FALSE, linetype = 2, colour = 'gray77') +
  geom_line() + 
  labs(x = 'plant size (pre-growth)', y = 'plant size (post-growth)') +
  scale_colour_manual(values = c('black', 'red', 'blue'))
# Here's what's up - control plots have the worst growth
# (also why does this plot take so long to generate?)

# (There's no survival difference)

# For LTRE

all.vector = all.matr %>%
  lapply(
    function(x)
      data.frame(
        w = Re(eigen(x)$vectors[,1]),
        v = Re(eigen(t(x))$vectors[,1])
      )
  )

all.sensit = lapply(all.vector, function(x) x$v %*% t(x$w) / (sum(x$v * x$w)))

all.sensit

# Plot of sensitivities?
all.sensit.df = all.sensit %>%
  lapply(
    function(df) {
      data.frame(df) %>%
        mutate(size.cur = (5:60)/10) %>%
        pivot_longer(-size.cur, names_to = 'size.pre', values_to = 'sens') %>%
        mutate(size.pre = as.numeric(gsub('X', '', size.pre)))
    }
  )

for (i in 1:length(all.sensit.df)) all.sensit.df[[i]]$pct = names(all.sensit.df)[i]

all.sensit.df = all.sensit.df %>%
  do.call(what = rbind) %>%
  separate(pct, into = c('trt', 'p.germ'), sep = '_')

all.sensit.df %>%
  mutate(log.sens = log(sens, base = 10)) %>%
  ggplot(aes(x = size.pre, y = size.cur, fill = log.sens)) +
  geom_tile() +
  scale_fill_viridis_c() +
  scale_y_reverse() +
  facet_wrap(trt ~ p.germ)

all.sensit.df %>%
  group_by(size.pre, trt, p.germ) %>%
  summarise(sum.sens = sum(sens)) %>%
  merge(data.frame(size.pre = 1:56, y = (5:60)/10)) %>%
  select(-size.pre) %>%
  rename(size.pre = y) %>%
  ggplot(aes(x = size.pre, y = sum.sens, group = interaction(trt, p.germ))) +
  geom_line(aes(colour = trt, linetype = factor(p.germ))) +
  scale_colour_manual(values = c('black', 'red', 'blue'))

# --- Okay... do LTRE just on matrix entries

diff.df = kern.df %>%
  pivot_wider(names_from = trt, values_from = p.size.cur) %>%
  mutate(
    d.c = drought - control,
    i.c = irrigated - control
  ) %>%
  select(-c(control, drought, irrigated)) %>%
  pivot_longer(c(d.c, i.c), names_to = 'contrast', values_to = 'diff')

diff.matr = split(diff.df, diff.df[,c("contrast", "p.germ")], sep = '_') %>%
  lapply(
    function(df) {
      df %>%
        arrange(size.prev, size.cur) %>%
        pivot_wider(names_from = size.prev, values_from = diff) %>%
        select(-c(contrast, p.germ, size.cur))
    }
  )

diff.matr[[1]]
names(diff.matr)

# I need the MATRIX MEANS for sensitivities for the LTRE

mean.df = kern.df %>%
  pivot_wider(names_from = trt, values_from = p.size.cur) %>%
  mutate(
    d.c = (drought + control) / 2,
    i.c = (irrigated + control) / 2
  ) %>%
  select(-c(control, drought, irrigated)) %>%
  pivot_longer(c(d.c, i.c), names_to = 'contrast', values_to = 'mean')

mean.vectors = split(mean.df, mean.df[,c("contrast", "p.germ")], sep = '_') %>%
  lapply(
    function(df) {
      df %>%
        arrange(size.prev, size.cur) %>%
        pivot_wider(names_from = size.prev, values_from = mean) %>%
        select(-c(contrast, p.germ, size.cur))
    }
  ) %>%
  lapply(
    function(x)
      data.frame(
        w = Re(eigen(x)$vectors[,1]),
        v = Re(eigen(t(x))$vectors[,1])
      )
  )

mean.sens = lapply(mean.vectors, function(x) x$v %*% t(x$w) / (sum(x$v * x$w)))

(names(mean.sens) == names(diff.matr))

# Okay... let's just try to get one of these up and running

ltre1 = diff.matr[[1]] * mean.sens[[1]]

data.frame(ltre1) %>%
  mutate(size.cur = (5:60)/10) %>%
  pivot_longer(-size.cur, names_to = 'size.pre', values_to = 'entry') %>%
  mutate(size.pre = as.numeric(gsub('X', '', size.pre))) %>%
  ggplot(aes(x = size.pre, y = size.cur, fill = entry)) +
  geom_raster() +
  geom_segment(aes(x = 0.5, xend = 6, y = 0.5, yend = 6), linetype = 2, linewidth = 0.25) +
  scale_y_reverse() +
  scale_fill_gradient2(low = 'red', high = 'blue', midpoint = 0, mid = 'white')
# Huh
# Ooh... much more interesting!

# Row sums will give summed contributions from each size
colSums(ltre1) %>% plot(type = 'l')
abline(h = 0)

sum(ltre1) # I think this is it!
  
# --- Do all treatments

ltre.list = vector(mode = 'list', length = length(mean.sens))
names(ltre.list) = names(mean.sens)

for (i in 1:length(ltre.list)) ltre.list[[i]] = diff.matr[[i]] * mean.sens[[i]]

ltre.as.matr = ltre.list %>%
  lapply(
    function(lentry) {
      data.frame(lentry) %>%
        mutate(size.cur = (5:60)/10) %>%
        pivot_longer(-size.cur, names_to = 'size.pre', values_to = 'entry') %>%
        mutate(size.pre = as.numeric(gsub('X', '', size.pre)))
    }
  )
  
# There's probably a slicker way to do this, but, add treatment column here
for (i in 1:length(ltre.list)) ltre.as.matr[[i]]$c.trt = names(ltre.list)[i]

ltre.as.df = do.call(rbind, ltre.as.matr) %>%
  separate(c.trt, into = c('contr', 'p.germ'), sep = '_')

head(ltre.as.df)

ltre.as.df %>%
  ggplot(aes(x = size.pre, y = size.cur, fill = entry)) +
  geom_tile() +
  geom_segment(
    aes(x = .5, xend = 6, y = 0, yend = 6), linetype = 2, colour = 'gray88'
  ) +
  scale_y_reverse() +
  scale_fill_gradient2(high = 'royalblue', low = 'red', mid = 'white', midpoint = 0) +
  facet_wrap(contr ~ p.germ, nrow = 2)
# Good!

ltre.as.df %>%
  group_by(contr, p.germ, size.pre) %>%
  summarise(total.contr = sum(entry)) %>%
  ggplot(aes(x = size.pre, y = total.contr, colour = p.germ, linetype = contr)) +
  geom_line()
# Positive contributions largely from small-intermediate plants
# contributions grow larger (in both directions but moreso for positive) with higher germ rates

ltre.as.df %>% group_by(contr, p.germ) %>% summarise(dlambda = sum(entry))
# Good!

all.lambda %>% 
  pivot_wider(names_from = trt, values_from = lambda) %>% 
  mutate(d.c = drought - control, i.c = irrigated - control)
# Should be getting sums that look like this (and we do!)

# ---

# Okay. Now, attempt to do this but for vital rates.

# Want rate differences for each treatment
# Rates: reproduction (including germination), survival, growth

rate.compare.wide = all.rates %>%
  mutate(germn.total = p.size.cur.f * p.germ) %>%
  select(
    trt, size.prev, size.cur, p.germ,
    surv = pred.surv, grow = p.grow.size, repr = germn.total
  ) %>%
  pivot_longer(cols = c(surv, grow, repr), names_to = 'vital', values_to = 'val') %>%
  pivot_wider(names_from = trt, values_from = val) %>%
  mutate(
    d.c = drought - control, 
    i.c = irrigated - control,
    d.c.mean = (drought + control) / 2,
    i.c.mean = (irrigated + control) / 2
  ) %>%
  select(-c(control, drought, irrigated))

rate.diff = rate.compare.wide %>%
  select(-c(d.c.mean, i.c.mean)) %>%
  pivot_longer(cols = c(d.c, i.c), names_to = 'contrast', values_to = 'diff') %>%
  pivot_wider(names_from = vital, values_from = diff) %>%
  arrange(size.cur, size.prev)

rate.mean = rate.compare.wide %>%
  select(-c(d.c, i.c)) %>%
  pivot_longer(cols = c(d.c.mean, i.c.mean), names_to = 'contrast', values_to = 'mean') %>%
  mutate(contrast = gsub('\\.mean', '', contrast)) %>%
  pivot_wider(names_from = vital, values_from = mean) %>%
  arrange(size.cur, size.prev)

# The sensitivities are the matrix sensitivities (found above) with a multiplier
# (chain rule)
# But, the growth sensitivity is not straightforward to estimate or interpret analytically
# So do it numerically with perturbations
# (for now - will only do growth numerically)
# NOTE that this does *not* give an entry-wise sensitivity but an overall
# sensitivity to a change in the *mean* of the growth kernel
# (but we are ultimately interested in relative effects of each vital rate on
# the whole, not necessarily entrywise, so while mildly constraining this is
# probably okay on the whole)

# # Code below gives sensitivity to shift in mean of growth (mu) but not
# entry-wise sensitivities which are useful for LTRE...
# g_mu.sens = split(perturbed.kern.df, f = perturbed.kern.df[,c("trt", "p.germ")], sep = '_') %>%
#   sapply(
#     function(df) {
#       m = df %>%
#         arrange(size.prev, size.cur) %>%
#         pivot_wider(names_from = size.prev, values_from = p.size.cur.perturb.growth.mu) %>%
#         select(-c(trt, p.germ, size.cur)) %>%
#         as.matrix()
#       return(Re(eigen(m)$values[1]))
#     }
#   ) %>%
#   data.frame(lambda_h = .) %>%
#   mutate(trt_p.germ = row.names(.)) %>%
#   separate(trt_p.germ, into = c('trt', 'p.germ'), sep = '_') %>%
#   merge(all.lambda) %>%
#   mutate(g_mu.sens = (lambda_h - lambda)/.001)
# 
# # A plot of growth sensitivities as a function of p-germ
# g_mu.sens %>% 
#   ggplot(aes(x = p.germ, y = g_mu.sens, group = trt)) + 
#   geom_line(aes(colour = trt)) + 
#   scale_colour_manual(values = c('black', 'red', 'blue'))

# okay let's just try getting one of these to work

rate.mean.matr = rate.mean %>%
  pivot_longer(c(surv, grow, repr), names_to = 'vital', values_to = 'mean') %>%
  split(.[,c("contrast", "p.germ")], sep = '_') %>%
  lapply(
    function(df) {
      df %>%
        arrange(size.prev, size.cur) %>%
        pivot_wider(names_from = size.prev, values_from = mean) %>%
        select(-c(contrast, p.germ, size.cur)) %>%
        split(.$vital) %>%
        lapply(function(x) x %>% select(-vital))
    }
  )

rate.diff.matr = rate.diff %>%
  pivot_longer(c(surv, grow, repr), names_to = 'vital', values_to = 'diff') %>%
  split(.[,c("contrast", "p.germ")], sep = '_') %>%
  lapply(
    function(df) {
      df %>%
        arrange(size.prev, size.cur) %>%
        pivot_wider(names_from = size.prev, values_from = diff) %>%
        select(-c(contrast, p.germ, size.cur)) %>%
        split(.$vital) %>%
        lapply(function(x) x %>% select(-vital))
    }
  )

names(rate.mean.matr) == names(mean.sens)
names(rate.diff.matr) == names(mean.sens)
# great

rate.ltre.list = vector('list', length = length(mean.sens))
names(rate.ltre.list) = names(mean.sens)

for (i in 1:length(rate.ltre.list)) {
  
  # growth sensitivities with the compensation term
  gs.num = (mean.sens[[i]] * rate.mean.matr[[i]]$grow * rate.mean.matr[[i]]$surv)
  gs.num = matrix(colSums(gs.num), nrow = nrow(gs.num), ncol = ncol(gs.num), byrow = TRUE) - gs.num
  gs.den = matrix(colSums(rate.mean.matr[[i]]$grow), nrow = nrow(gs.num), ncol = ncol(gs.num), byrow = TRUE) -
    rate.diff.matr[[i]]$grow
  growth.sens = (rate.mean.matr[[i]]$surv * mean.sens[[i]]) - gs.num / gs.den
  
  rate.ltre.list[[i]]$surv = rate.diff.matr[[i]]$surv * mean.sens[[i]] * rate.mean.matr[[i]]$grow
  rate.ltre.list[[i]]$grow = rate.diff.matr[[i]]$grow * growth.sens
  # rate.ltre.list[[i]]$grow = rate.diff.matr[[i]]$grow * mean.sens[[i]] * rate.mean.matr[[i]]$surv
  rate.ltre.list[[i]]$repr = rate.diff.matr[[i]]$repr * mean.sens[[i]]
}

sapply(rate.ltre.list[[1]], sum) #%>% sum()

summed.rate.ltre.contrs = sapply(rate.ltre.list, function(x) sapply(x, sum)) %>% 
  data.frame() %>%
  mutate(rate = row.names(.)) %>%
  pivot_longer(-rate, names_to = 'contr_p.germ', values_to = 'contrib') %>%
  separate(contr_p.germ, into = c('contrast', 'p.germ'), sep = '_')

summed.rate.ltre.contrs %>%
  mutate(contrast = paste0(ifelse(contrast %in% 'd.c', 'drought', 'irrigation'), ' vs. control')) %>%
  ggplot(aes(x = p.germ, y = contrib, shape = rate, group = rate)) +
  geom_point(position = position_dodge(width = .5), size = 4) +
  scale_colour_manual(values = c('red', 'blue')) +
  facet_wrap(~ contrast)

