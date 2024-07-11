# Deterministic LTRE analysis 
# BUT - I think this script isn't accounting for adult survival when doing the
# seedling recruitment

# CHECK LIFE CYCLE DIAGRAM
# and then GO THROUGH SCRIPT and ADD SURVIVAL TO RECRUITMENT as needed


library(ggplot2)
library(tidyr)
library(dplyr)

growsurv = read.csv('03_construct_kernels/out/deterministic_growsurv_kernel.csv')
reprod = read.csv('03_construct_kernels/out/deterministic_reprod_kernel_no_phen.csv')

germin = expand.grid(
  trt = c('control', 'drought', 'irrigated'),
  p.germ = c(0.001, 0.005, 0.01, 0.05)
)

all.rates = merge(
  x = growsurv, y = reprod,
  by.x = c('size.prev', 'size.cur', 'trt'), by.y = c('size.pre', 'size.cur', 'trt'),
  suffixes = c('.g', '.f')
) %>%
  merge(germin, by = 'trt')

kern.df = all.rates %>%
  select(size.prev, size.cur, trt, p.size.cur.g, p.size.cur.f, p.germ) %>%
  mutate(p.size.cur.f = p.size.cur.f * p.germ) %>%
  mutate(p.size.cur = p.size.cur.g + p.size.cur.f) %>%
  select(-c(p.size.cur.g, p.size.cur.f))

head(kern.df)

kern.df %>%
  ggplot(aes(x = size.prev, y = size.cur, fill = p.size.cur)) +
  geom_raster() +
  scale_y_reverse() +
  scale_fill_viridis_c() +
  facet_wrap(p.germ ~ trt)

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

all.lambda = sapply(all.matr, function(x) Re(eigen(x)$values[1])) %>%
  data.frame(lambda = .) %>%
  mutate(t.germ = row.names(.)) %>%
  separate(t.germ, into = c('trt', 'p.germ'), sep = '_')

all.lambda %>%
  ggplot(aes(x = p.germ, y = lambda, group = trt, colour = trt)) +
  geom_point(size = 3) +
  geom_line() +
  scale_colour_manual(values = c('black', 'red', 'blue'))
# LMAOOOOOOOOO

# Seed production
all.rates %>% 
  select(size.prev, trt, seeds.total) %>% 
  ggplot(aes(x = size.prev, y = seeds.total, group = trt, colour = trt)) + 
  geom_line() + 
  scale_y_log10() +
  scale_colour_manual(values = c('black', 'red', 'blue'))
# Looks like contorl plants produce more seeds *when they are large*
# Irrigated plants produce more seeds when small...

# Growth
all.rates %>% 
  ggplot(aes(x = size.prev, y = pred.grow.mean, group = trt, colour = trt)) + 
  geom_line() + 
  scale_y_log10() +
  scale_colour_manual(values = c('black', 'red', 'blue'))
# Here's what's up - control plots have the worst growth

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
  mutate(size.pre = (5:60)/10) %>%
  pivot_longer(-size.pre, names_to = 'size.cur', values_to = 'entry') %>%
  mutate(size.cur = as.numeric(gsub('X', '', size.cur))) %>%
  ggplot(aes(x = size.cur, y = size.pre, fill = entry)) +
  geom_raster() +
  scale_y_reverse() +
  scale_fill_gradient2(low = 'red', high = 'blue', midpoint = 0, mid = 'white')
# Huh
# Ooh... much more interesting!

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

ltre.as.df %>% group_by(contr, p.germ) %>% summarise(dlambda = sum(entry))
# Good!

all.lambda %>% 
  pivot_wider(names_from = trt, values_from = lambda) %>% 
  mutate(d.c = drought - control, i.c = irrigated - control)
# Should be getting sums that look like this...

# ---

# Okay. Now, attempt to do this but for vital rates.

# Want rate differences for each treatment
# Rates: reproduction (including germination), survival, growth

rate.compare.wide = all.rates %>%
  mutate(germn.total = seeds.total * p.germ) %>%
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
  pivot_wider(names_from = vital, values_from = diff)

rate.mean = rate.compare.wide %>%
  select(-c(d.c, i.c)) %>%
  pivot_longer(cols = c(d.c.mean, i.c.mean), names_to = 'contrast', values_to = 'mean') %>%
  mutate(contrast = gsub('\\.mean', '', contrast)) %>%
  pivot_wider(names_from = vital, values_from = mean)

# The sensitivities are the matrix sensitivities (found above) with a multiplier
# (chain rule)

# Might be easiest to do this with a hadamard product?

rate.mean.matr = rate.mean %>%
  pivot_longer(c(surv, grow, repr), names_to = 'vital', values_to = 'mean') %>%
  split(.[,c("vital", "contrast", "p.germ")], sep = '_') %>%
  lapply(
    function(df) {
      df %>%
        arrange(size.prev, size.cur) %>%
        pivot_wider(names_from = size.prev, values_from = mean) %>%
        select(-c(contrast, p.germ, size.cur, vital))
    }
  )

rate.diff.matr = rate.diff %>%
  pivot_longer(c(surv, grow, repr), names_to = 'vital', values_to = 'diff') %>%
  split(.[,c("vital", "contrast", "p.germ")], sep = '_') %>%
  lapply(
    function(df) {
      df %>%
        arrange(size.prev, size.cur) %>%
        pivot_wider(names_from = size.prev, values_from = diff) %>%
        select(-c(contrast, p.germ, size.cur, vital))
    }
  )

# Okay... this is not hte most hygenic way to do this, but we can just replicate
# the list

rate.mean.sens = rep(mean.sens, each = 3)
# the `each` argument here means elements are duplicated in the same order as the lists above

rate.ltre.list = vector('list', length = length(rate.mean.sens))

names(rate.diff.matr) == names(rate.mean.matr)

# OH... lmao this won't work out like I hoped

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
names(ltre.list) = names(mean.sens)

for (i in 1:length(rate.ltre.list)) {
  rate.ltre.list[[i]]$surv = rate.diff.matr[[i]]$surv * mean.sens[[i]] *
    with(rate.mean.matr[[i]], grow + repr)
  rate.ltre.list[[i]]$grow = rate.diff.matr[[i]]$grow * mean.sens[[i]] * rate.mean.matr[[i]]$surv
  rate.ltre.list[[i]]$repr = rate.diff.matr[[i]]$repr * mean.sens[[i]] * rate.mean.matr[[i]]$surv
}

# Okay... well at the very least I can do this:
sapply(rate.ltre.list[[1]], sum) #%>% sum()
# reproduction value looks too big
# oh and also not sure how to get the VR column in here...

sapply(rate.ltre.list, function(x) sapply(x, sum))
# reproduction values