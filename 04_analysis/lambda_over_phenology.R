# ---------------------------
# Estimating lambda across a range of mean population bud dates
# Reads in some *very large files*
# Primarily producing a figure
# ---------------------------

# --- Setup ---------------------------------------------------
library(ggplot2)
library(tidyr)
library(dplyr)
library(cowplot)

rm(list = ls())

# # Get point estimates for kernels estimated across a phenology range

# Growth + survival kernel
growsurv = read.csv('03_construct_kernels/out/deterministic_growsurv_kernel.csv') %>%
  # (note: at some point I should go back and re-export this file without these columns,
  # after which this line of code should be deleted)
  select(-c(pred.surv, pred.grow.mean, p.grow.size))
# Reproductive kernel (all phenology)
reprodct.all = read.csv('03_construct_kernels/out/deterministic_reprod_kernel_phen.csv')
# Reproductive kernel (for LTRE only)
reprodct.ltre = read.csv('03_construct_kernels/out/determinstic_reprod_kernel_phen_ltre.csv')

# # Get bootstrapped intervals
# Growth + survival 
gs.boot = read.csv('03_construct_kernels/out/deterministic_growsurv_bootstrap.csv')
# Reproductive (all phenology)
fr.boot.all = read.csv('03_construct_kernels/out/deterministic_reprod_bootstrap_allphen.csv')
# Reproductive (for LTRE only)
fr.boot.ltre = read.csv('03_construct_kernels/out/deterministic_reprod_bootstrap_ltre.csv')

head(growsurv)
head(reprodct.all)

# Germination probability
p.germ = .001

# --- All-phenology kernels

# Data frame to produce kernels (with point estimates) in data frame form

kernel.all.df = merge(
  # Survival + growth subkernel
  growsurv,
  # Reproductive subkernels (for *each phenology*)
  reprodct.all,
  by.x = c('size.prev', 'size.cur', 'trt'), by.y = c('size.prev', 'size.nex', 'trt'),
  suffixes = c('.g', '.r')
) %>%
  # Combine growth and survival entries into single kernel entry
  mutate(p.size.cur = p.size.cur.g + p.size.cur.r * p.germ) %>%
  select(-c(p.size.cur.g, p.size.cur.r))

head(kernel.all.df)

# Get lambda estimates for each treatment on each mean buddate
# note: these will be only point estimates of lambda; uncertainty will come from
# bootstrapped estimates
all.lambda = split(kernel.all.df, kernel.all.df[,c("trt", "phen")], sep = '_', drop = TRUE) %>%
  # Split the kernel data up by phenology/treatment and convert each into matrix form
  lapply(
    function(df) {
      df %>%
        arrange(size.prev, size.cur) %>%
        pivot_wider(names_from = size.prev, values_from = p.size.cur) %>%
        select(-c(trt, phen, size.cur)) %>%
        as.matrix()
    }
  ) %>% 
  # Get lambda (maximum eigenvalue of each matrix)
  sapply(function(m) Re(eigen(m)$values[1])) %>%
  # Convert to data frame with treatrment/phen info
  data.frame(lambda = .) %>%
  mutate(trt_phen = row.names(.)) %>%
  separate(trt_phen, into = c('trt', 'phen'), sep = '_') %>%
  # Convert phen column into a date type
  mutate(phen.date = as.Date(as.numeric(phen), format = '%b-%d'))

head(all.lambda)
# good

# Merge together bootstrapped subkernels
# slow - takes about a minute
kernel.all.boot.df = merge(
  gs.boot     %>% pivot_longer(starts_with('b'), names_to = 'boot', values_to = 'p.size.cur'),
  fr.boot.all %>% pivot_longer(starts_with('b'), names_to = 'boot', values_to = 'p.size.cur'),
  by.x = c('size.prev', 'size.cur', 'trt', 'boot'), by.y = c('size.prev', 'size.nex', 'trt', 'boot'),
  suffixes = c('.g', '.f')
) %>%
  mutate(p.size.cur = p.size.cur.g + p.germ * p.size.cur.f) %>%
  select(-c(p.size.cur.g, p.size.cur.f))

# Split this and convert to matrix form, then estimate lambda from matrices
# (also slowish)
all.boot.lambda = split(
  kernel.all.boot.df, kernel.all.boot.df[,c("trt", "boot", "mean.phen")], 
  sep = '_', drop = TRUE
) %>%
  lapply(
    function(df) {
      df %>%
        arrange(size.prev, size.cur) %>%
        pivot_wider(names_from = size.prev, values_from = p.size.cur) %>%
        select(-c(trt, mean.phen, boot, size.cur)) %>%
        as.matrix()
    }
  ) %>%
  # Get eigenvalues
  sapply(function(x) Re(eigen(x)$values[1])) %>%
  data.frame(lambda = .) %>%
  mutate(ttbp = row.names(.)) %>%
  separate(ttbp, into = c('trt', 'boot', 'phen'), sep = '_') %>%
  # Convert phen column into a date type
  mutate(phen.date = as.Date(as.numeric(phen), format = '%b-%d'))

# --- LTRE (observed trt-phen combos)

# Data frame for point estimate kernels
kernel.ltre.df = merge(
  growsurv, reprodct.ltre,
  by.x = c('size.prev', 'size.cur', 'trt'), by.y = c('size.prev', 'size.nex', 'trt'),
  suffixes = c('.g', '.r')
) %>%
  # Combining growth/surv and reproduction subkernels 
  mutate(p.size.cur = p.size.cur.g + p.size.cur.r * p.germ) %>%
  select(-c(p.size.cur.g, p.size.cur.r))

# Convert to get matrices and eigenvalues
ltre.lambda = split(
  kernel.ltre.df, kernel.ltre.df[,c("trt", "trt.phen", "phen")], sep = '_', drop = TRUE
) %>%
  # Split the kernel data up by phenology/treatment and convert each into matrix form
  lapply(
    function(df) {
      df %>%
        arrange(size.prev, size.cur) %>%
        pivot_wider(names_from = size.prev, values_from = p.size.cur) %>%
        select(-c(trt, trt.phen, phen, size.cur)) %>%
        as.matrix()
    }
  ) %>% 
  # Get lambda (maximum eigenvalue of each matrix)
  sapply(function(m) Re(eigen(m)$values[1])) %>%
  # Convert to data frame with treatrment/phen info
  data.frame(lambda = .) %>%
  mutate(trt_phen = row.names(.)) %>%
  separate(trt_phen, into = c('trt', 'trt.phen', 'phen'), sep = '_') %>%
  # Convert phen column into a date type
  mutate(phen.date = as.Date(as.numeric(phen), format = '%b-%d'))

# Get bootstrapped LTRE kernels
# (will take a sec to run)
kernel.boot.ltre.df = merge(
  gs.boot      %>% pivot_longer(starts_with('b'), names_to = 'boot', values_to = 'p.size.cur'),
  fr.boot.ltre %>% pivot_longer(starts_with('b'), names_to = 'boot', values_to = 'p.size.cur'),
  by.x = c('size.prev', 'size.cur', 'trt', 'boot'), by.y = c('size.prev', 'size.nex', 'trt', 'boot'),
  suffixes = c('.g', '.f')
) %>%
  mutate(p.size.cur = p.size.cur.g + p.germ * p.size.cur.f) %>%
  select(-c(p.size.cur.g, p.size.cur.f))

ltre.boot.lambda = split(
  kernel.boot.ltre.df, kernel.boot.ltre.df[,c("trt", 'trt.phen', "boot")], sep = '_', drop = TRUE
) %>%
  lapply(
    function(df) {
      df %>%
        arrange(size.prev, size.cur) %>%
        pivot_wider(names_from = size.prev, values_from = p.size.cur) %>%
        select(-c(trt, trt.phen, boot, size.cur)) %>%
        as.matrix()
    }
  ) %>%
  # Get eigenvalues
  sapply(function(x) Re(eigen(x)$values[1])) %>%
  data.frame(lambda = .) %>%
  mutate(ttb = row.names(.)) %>%
  separate(ttb, into = c('trt', 'trt.phen', 'boot'), sep = '_')

# # Summaries for plots

# Get 95% bootstrapped intervals for all bootstrapped datasets

# Bootstrapped intervals at each date for each treatment
all.boot.intervals = all.boot.lambda %>%
  group_by(trt, phen.date) %>%
  reframe(
    cibound = quantile(lambda, c(0.025, 0.975)),
    lohi = c('lo', 'hi')
  ) %>%
  pivot_wider(names_from = lohi, values_from = cibound)


# Bootstrapped intervals just for just the observed LTRE days
ltre.boot.intervals = ltre.boot.lambda %>%
  group_by(trt, trt.phen) %>%
  reframe(
    cibound = quantile(lambda, c(0.025, 0.975)),
    lohi = c('lo', 'hi')
  ) %>%
  pivot_wider(names_from = lohi, values_from = cibound)

# Read in phen dates for mean phen among treatments:
source('03_construct_kernels/prepare_demo_data_repr.R')

# Fit phen model
d_t = glmmTMB(
  phen.julian ~ trt + Year + (1 | Plot / plantid),
  data = phen
)

annual.phen.dates = expand.grid(
  trt = c('control', 'drought', 'irrigated'),
  Year = factor(2021:2024)
) %>%
  mutate(
    mean.phen = predict(d_t, newdata = ., allow.new.levels = TRUE, re.form = ~ 0)
  ) %>%
  # Convert to date format (for plotting)
  mutate(phen.date = as.Date(mean.phen, format = '%b-%d')) %>%
  mutate(
    ydodge = case_when(
      trt %in% 'control' ~ 0.915,
      trt %in% 'drought' ~ 0.9175,
      trt %in% 'irrigated' ~ 0.9125
    )
  )

# Make plot

lambda.trt.pan = all.lambda %>%
  ggplot(aes(x = phen.date, group = trt)) +
  geom_ribbon(
    data = all.boot.intervals,
    aes(x = phen.date, ymin = lo, ymax = hi, fill = trt),
    alpha = 0.125
  ) +
  geom_line(aes(y = lambda, colour = trt), linewidth = 1.2) +
  # geom_segment(
  #   data = ltre.boot.intervals,
  #   aes(xend = phen.date, y = lo, yend = hi, colour = trt)
  # ) +
  geom_point(
    data = annual.phen.dates,
    aes(x = phen.date, y = ydodge, colour = trt),
    size = 10, shape = '*'
  ) +
  geom_point(
    data = ltre.lambda %>% mutate(phen.date = as.Date(as.numeric(phen), format = '%b-%d')),
    aes(y = lambda, colour = trt, shape = (trt == trt.phen)), size = 4
  ) +
  scale_shape_manual(values = c(1, 19)) +
  scale_colour_manual(values = c('black', 'red', 'blue')) +
  scale_fill_manual(values = c('black', 'red', 'blue')) +
  guides(shape = 'none') +
  labs(x = 'Mean bud date') +
  theme(
    panel.background = element_blank(),
    legend.position = 'none'
  )

lambda.trt.pan

# But... plot bootstrapped treatment differences over time
# (first need to assemble these)

boot.lambda.diff = all.boot.lambda %>%
  pivot_wider(names_from = trt, values_from = lambda) %>%
  mutate(d.c = drought - control, i.c = irrigated - control) %>%
  select(-c(drought, control, irrigated)) %>%
  pivot_longer(c(d.c, i.c), names_to = 'contrast', values_to = 'd.lambda')

boot.lambda.diff.interval = boot.lambda.diff %>%
  group_by(phen.date, contrast) %>%
  reframe(
    d.lambda = quantile(d.lambda, probs = c(0.025, 0.975)),
    lohi = c('lo', 'hi')
  ) %>%
  pivot_wider(names_from = lohi, values_from = d.lambda) %>%
  # Do a merge to get the observed mean d.lambda
  merge(
    y = all.lambda %>%
      pivot_wider(names_from = trt, values_from = lambda) %>%
      mutate(d.c = drought - control, i.c = irrigated - control) %>%
      select(-c(drought, control, irrigated)) %>%
      pivot_longer(c(d.c, i.c), names_to = 'contrast', values_to = 'mean.d.lambda')
  )

ltre.lambda.diff = all.boot.lambda %>%
  pivot_wider(names_from = trt, values_from = lambda) %>%
  mutate(d.c = drought - control, i.c = irrigated - control) %>%
  select(-c(drought, control, irrigated)) %>%
  pivot_longer(c(d.c, i.c), names_to = 'contrast', values_to = 'd.lambda') %>%
  filter(!is.na(d.lambda))

lambda.contr.pan = boot.lambda.diff %>%
  mutate(
    contrast = paste(ifelse(contrast %in% 'd.c', 'drought', 'irrigated'), 'vs. control')
  ) %>%
  ggplot(aes(x = phen.date, group = contrast)) +
  annotate(
    'segment',
    x = as.Date('1970-04-08'), xend = as.Date('1970-06-03'),
    y = 0, yend = 0,
    linetype = 2, colour = 'gray'
  ) +
  geom_point(
    aes(y = d.lambda),
    position = position_jitter(width = 1), alpha = 0.25
  ) +
  geom_point(
    data = boot.lambda.diff.interval %>%
      mutate(
        contrast = paste(ifelse(contrast %in% 'd.c', 'drought', 'irrigated'), 'vs. control')
      ),
    aes(y = mean.d.lambda),
    size = 4, shape = 21, stroke = 2
  ) +
  geom_point(
    # y limits here found by:
    # ggplot_build(lambda.contr.pan)$layout$panel_scales_y
    data = annual.phen.dates %>%
      filter(!(trt %in% 'control')) %>%
      mutate(contrast = paste(trt, 'vs. control')),
    aes(x = phen.date, y = -0.007, colour = contrast),
    shape = '*', size = 10
  ) +
  geom_segment(
    data = boot.lambda.diff.interval %>%
      mutate(
        contrast = paste(ifelse(contrast %in% 'd.c', 'drought', 'irrigated'), 'vs. control')
      ),
    aes(xend = phen.date, y = lo, yend = hi)
  ) +
  # scale_shape_manual(values = c(1, 19)) +
  labs(x = 'Mean bud date', y = expression(Delta~lambda)) +
  guides(colour = 'none', fill = 'none') +
  scale_colour_manual(values = c('red', 'blue')) +
  # scale_fill_manual(values = c('red', 'blue')) +
  facet_wrap(~ contrast, nrow = 2) +
  theme(panel.background = element_blank())

plot_grid(
  lambda.trt.pan, lambda.contr.pan, nrow = 1
)

# good start
