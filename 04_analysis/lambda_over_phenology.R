# ---------------------------
# Estimating lambda across a range of mean population bud dates
# Reads in some *very large files*
# Primarily producing a figure
# ---------------------------

# --- Setup ---------------------------------------------------
library(ggplot2)
library(ggh4x)
library(tidyr)
library(dplyr)
library(glmmTMB)
library(cowplot)

rm(list = ls())

# # Get point estimates for kernels estimated across a phenology range

# Growth + survival kernel
growsurv.all = read.csv('03_construct_kernels/out/deterministic_growsurv_kernel_phen.csv') %>%
  # filter to only two weeks before/after mean
  filter(phen.c %in% -14:14)
# Reproductive kernel (all phenology)
reprodct.all = read.csv('03_construct_kernels/out/deterministic_reprod_kernel_phen.csv') %>%
  # filter to only two weeks before/after mean
  filter(phen.c %in% -14:14)

# Growth/survival kernel for LTRE only
growsurv.ltre = read.csv('03_construct_kernels/out/deterministic_growsurv_kernel_phen_ltre.csv')
# Reproductive kernel for LTRE only
reprodct.ltre = read.csv('03_construct_kernels/out/determinstic_reprod_kernel_phen_ltre.csv')

# # Get bootstrapped intervals
# Growth + survival (all phenology)
gs.boot.all = read.csv('03_construct_kernels/out/deterministic_growsurv_bootstrap_allphen.csv') %>%
  filter(phen.c %in% -14:14)
# Growth + survival (for LTRE only)
gs.boot.ltre = read.csv('03_construct_kernels/out/deterministic_growsurv_bootstrap_ltre.csv')
# Reproductive (all phenology)
fr.boot.all = read.csv('03_construct_kernels/out/deterministic_reprod_bootstrap_allphen.csv') %>%
  # filter to only two weeks before/after mean
  filter(phen.c %in% -14:14)
# Reproductive (for LTRE only)
fr.boot.ltre = read.csv('03_construct_kernels/out/deterministic_reprod_bootstrap_ltre.csv')

head(growsurv.all)
head(reprodct.all)

# Read in LTRE treatment-phenology info
trt.phen.ltre.key = merge(
  x = read.csv('03_construct_kernels/ltre_treatment_key.csv'),
  y = read.csv('03_construct_kernels/phen_treatment_means.csv'),
  by.x = 'trt.phen', by.y = 'trt'
) %>%
  arrange(trt.phen.idx) %>%
  select(trt.phen.idx, everything())

# Control mean for re-centering phenology
phen.ctrl.mean = read.csv('03_construct_kernels/phen_treatment_means.csv') %>%
  filter(trt %in% 'control') %>%
  pull(mean.phen)

# Germination probability
p.germ = .001
# p.germ = 0.0058007812


# --- All-phenology kernels

# Data frame to produce kernels (with point estimates) in data frame form

kernel.all.df = merge(
  # Survival + growth subkernel
  growsurv.all,
  # Reproductive subkernels (for *each phenology*)
  reprodct.all,
  by.x = c('size.prev', 'size.cur', 'trt', 'phen.c'), by.y = c('size.prev', 'size.nex', 'trt', 'phen.c'),
  suffixes = c('.g', '.r')
) %>%
  # Combine growth and survival entries into single kernel entry
  mutate(
    p.size.cur = pred.surv * (pv.grow.size * (1 - prob.flower) + pf.grow.size * prob.flower) + (p.size.cur * p.germ)
  ) %>%
  # Re-center phenology around control mean
  mutate(phen = phen.c + phen.ctrl.mean) %>%
  select(-c(pred.surv, pv.grow.size, pf.grow.size, prob.flower, phen.c))

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
  gs.boot.all,    # %>% pivot_longer(starts_with('b'), names_to = 'boot', values_to = 'p.size.cur'),
  fr.boot.all,    # %>% pivot_longer(starts_with('b'), names_to = 'boot', values_to = 'p.size.cur'),
  by.x = c('size.prev', 'size.cur', 'trt', 'phen.c', 'boot'), by.y = c('size.prev', 'size.nex', 'trt', 'phen.c', 'boot'),
  suffixes = c('.g', '.f')
) %>%
  # mutate(p.size.cur = p.size.cur.g + p.germ * p.size.cur.f) %>%
  # Get size for matrix entries
  mutate(p.size.cur = pred.surv * (pv.grow.size * (1 - prob.flower) + pf.grow.size * prob.flower) + (p.size.cur * p.germ)) %>%
  # Re-center phenology to observed control mean
  mutate(phen = phen.c + phen.ctrl.mean) %>%
  select(-c(pred.surv, pv.grow.size, pf.grow.size, prob.flower, phen.c))

# Split this and convert to matrix form, then estimate lambda from matrices
# (also slowish)
all.boot.lambda = split(
  kernel.all.boot.df, kernel.all.boot.df[,c("trt", "boot", "phen")], 
  sep = '_', drop = TRUE
) %>%
  lapply(
    function(df) {
      df %>%
        arrange(size.prev, size.cur) %>%
        pivot_wider(names_from = size.prev, values_from = p.size.cur) %>%
        select(-c(trt, phen, boot, size.cur)) %>%
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
  growsurv.ltre, reprodct.ltre,
  by.x = c('size.prev', 'size.cur', 'trt.phen.idx'), 
  by.y = c('size.prev', 'size.nex', 'trt.phen.idx'),
  suffixes = c('.g', '.r')
) %>%
  # Combining growth/surv and reproduction subkernels 
  mutate(
    p.size.cur = pred.surv * (pv.grow.size * (1 - prob.flower) + pf.grow.size * prob.flower) + (p.size.cur * p.germ)
  ) %>%
  select(-c(pred.surv, pv.grow.size, pf.grow.size, prob.flower))

# Convert to get matrices and eigenvalues
ltre.lambda = split(
  kernel.ltre.df, kernel.ltre.df$trt.phen.idx, sep = '_', drop = TRUE
) %>%
  # Split the kernel data up by phenology/treatment and convert each into matrix form
  lapply(
    function(df) {
      df %>%
        arrange(size.prev, size.cur) %>%
        pivot_wider(names_from = size.prev, values_from = p.size.cur) %>%
        select(-c(trt.phen.idx, size.cur)) %>%
        as.matrix()
    }
  ) %>% 
  # Get lambda (maximum eigenvalue of each matrix)
  sapply(function(m) Re(eigen(m)$values[1])) %>%
  # Convert to data frame with treatrment/phen info
  data.frame(lambda = .) %>%
  mutate(trt.phen.idx = row.names(.)) %>%
  merge(trt.phen.ltre.key) %>%
  # Convert phen column into a date type
  mutate(phen.date = as.Date(as.numeric(mean.phen)))

# Get bootstrapped LTRE kernels
# (will take a sec to run)
kernel.boot.ltre.df = merge(
  gs.boot.ltre, # %>% pivot_longer(starts_with('b'), names_to = 'boot', values_to = 'p.size.cur'),
  fr.boot.ltre %>% mutate(boot = gsub('b', '', boot)), # %>% pivot_longer(starts_with('b'), names_to = 'boot', values_to = 'p.size.cur'),
  by.x = c('size.prev', 'size.cur', 'trt.phen.idx', 'boot'), 
  by.y = c('size.prev', 'size.nex', 'trt.phen.idx', 'boot'),
  suffixes = c('.g', '.f')
) %>%
  mutate(
    p.size.cur = pred.surv * (pv.grow.size * (1 - prob.flower) + pf.grow.size * prob.flower) + (p.size.cur * p.germ)
  ) %>%
  select(-c(pred.surv, pv.grow.size, pf.grow.size, prob.flower))

ltre.boot.lambda = split(
  kernel.boot.ltre.df, kernel.boot.ltre.df[,c("trt.phen.idx", "boot")], sep = '_', drop = TRUE
) %>%
  lapply(
    function(df) {
      df %>%
        arrange(size.prev, size.cur) %>%
        pivot_wider(names_from = size.prev, values_from = p.size.cur) %>%
        select(-c(trt.phen.idx, boot, size.cur)) %>%
        as.matrix()
    }
  ) %>%
  # Get eigenvalues
  sapply(function(x) Re(eigen(x)$values[1])) %>%
  data.frame(lambda = .) %>%
  mutate(idx.boot = row.names(.)) %>%
  separate(idx.boot, into = c('trt.phen.idx', 'boot'), sep = '_') %>%
  merge(trt.phen.ltre.key) %>%
  # Convert phen column into a date type
  mutate(phen.date = as.Date(as.numeric(mean.phen)))

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
  group_by(trt.phen.idx) %>%
  reframe(
    cibound = quantile(lambda, c(0.025, 0.975)),
    lohi = c('lo', 'hi')
  ) %>%
  merge(trt.phen.ltre.key) %>%
  rename(trt = trt.rate) %>%
  mutate(phen.date = as.Date(as.numeric(mean.phen), format = '%b-%d')) %>%
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
  ggplot(aes(x = phen.date)) +
  annotate(
    'segment', linetype = 2, colour = 'gray',
    x = min(all.lambda$phen.date), xend = max(all.lambda$phen.date),
    y = 1, yend = 1
  ) +
  geom_ribbon(
    data = all.boot.intervals,
    aes(x = phen.date, ymin = lo, ymax = hi, fill = trt, group = trt),
    alpha = 0.125
  ) +
  geom_line(aes(y = lambda, colour = trt, group = trt), linewidth = 1.2) +
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
    data = ltre.lambda %>% filter(trt.rate == trt.phen),
    aes(y = lambda, colour = trt.rate), size = 4, shape = 19
  ) +
  # scale_shape_manual(values = c(NA, 19)) +
  scale_colour_manual(values = c('black', 'goldenrod', 'dodgerblue'), '') +
  scale_fill_manual(values = c('black', 'goldenrod', 'dodgerblue'), '') +
  guides(shape = 'none') +
  labs(x = 'Mean emergence date', y = expression(lambda)) +
  theme(
    panel.background = element_blank(),
    legend.position = 'top'
    # legend.position = 'inside',
    # legend.position.inside = c(0.8, 0.8)
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
    # contrast = paste(ifelse(contrast %in% 'd.c', 'drought', 'irrigated'), 'vs. control')
    contrast = ifelse(
      contrast %in% 'd.c',
      'i) drought vs. control',
      'ii) irrigated vs. control'
    )
  ) %>%
  ggplot(aes(x = phen.date, group = contrast)) +
  annotate(
    'segment',
    x = as.Date('1970-04-22'), xend = as.Date('1970-05-20'),
    y = 0, yend = 0,
    linetype = 2, colour = 'gray'
  ) +
  geom_point(
    aes(y = d.lambda),
    position = position_jitter(width = 1), alpha = 0.125
  ) +
  geom_point(
    data = boot.lambda.diff.interval %>%
      mutate(
        # contrast = paste(ifelse(contrast %in% 'd.c', 'drought', 'irrigated'), 'vs. control')
        contrast = ifelse(
          contrast %in% 'd.c',
          'i) drought vs. control',
          'ii) irrigated vs. control'
        )
      ),
    aes(y = mean.d.lambda),
    size = 4, shape = 21, stroke = 2
  ) +
  geom_point(
    # y limits here found by:
    # ggplot_build(lambda.contr.pan)$layout$panel_scales_y
    data = annual.phen.dates %>%
      filter(!(trt %in% 'control')) %>%
      mutate(
        # contrast - paste(trt, 'vs. control')),
        contrast = ifelse(
          trt %in% 'drought',
          'i) drought vs. control',
          'ii) irrigated vs. control'
        )
      ),
    aes(x = phen.date, y = -0.007, colour = contrast),
    shape = '*', size = 10
  ) +
  geom_segment(
    data = boot.lambda.diff.interval %>%
      mutate(
        # contrast = paste(ifelse(contrast %in% 'd.c', 'drought', 'irrigated'), 'vs. control')
        contrast = ifelse(
          contrast %in% 'd.c',
          'i) drought vs. control',
          'ii) irrigated vs. control'
        )
      ),
    aes(xend = phen.date, y = lo, yend = hi)
  ) +
  # scale_shape_manual(values = c(1, 19)) +
  labs(x = 'Mean emergence date', y = expression(Delta~lambda)) +
  guides(colour = 'none', fill = 'none') +
  # scale_colour_manual(values = c('red', 'blue')) +
  # scale_fill_manual(values = c('red', 'blue')) +
  scale_colour_manual(values = c('goldenrod', 'dodgerblue')) +
  scale_x_continuous(
    breaks = as.Date(c('1970-04-27', '1970-05-04', '1970-05-11', '1970-05-18')),
    labels = format(as.Date(c('1970-04-27', '1970-05-04', '1970-05-11', '1970-05-18')), '%b %d')
  ) +
  facet_wrap(~ contrast, nrow = 2) +
  theme(
    panel.background = element_blank(),
    strip.background = element_part_rect(fill = 'white', side = 'b', colour = 'gray22')
  )

# lambda.legend = get_plot_component(
#   lambda.trt.pan + theme(legend.position = 'top'),
#   pattern = 'guide-box', return_all = TRUE
# )[[4]]

plot_grid(
  lambda.trt.pan, lambda.contr.pan, 
  labels = c('a)', 'b)'), rel_widths = c(1, 0.5),
  nrow = 1
) %>%
  save_plot(
    filename = '04_analysis/figures/draft_figures/lambdas_phen.png',
    base_height = 5, base_width = 8
  )

# legend needs to be smaller... now sure how to do this and keep size consistent...

# Trying something...

pp = all.lambda %>%
  filter(phen.date %in% as.Date(121:129)) %>%
  ggplot(aes(x = phen.date)) +
  geom_line(aes(y = lambda, colour = trt, group = trt), linewidth = 1.2) +
  geom_point(
    data = ltre.lambda %>% filter(trt.rate == trt.phen),
    aes(y = lambda, colour = trt.rate), size = 4, shape = 19
  ) +
  geom_point(
    data = ltre.lambda %>% filter(trt.rate %in% 'control'),
    aes(y = lambda, colour = trt.phen), size = 4, shape = 21
  ) +
  scale_shape_manual(values = c(NA, 19)) +
  scale_colour_manual(values = c('black', 'goldenrod', 'dodgerblue'), '') +
  scale_fill_manual(values = c('black', 'goldenrod', 'dodgerblue'), '') +
  guides(shape = 'none') +
  labs(x = 'Mean bud date', y = expression(lambda)) +
  theme(
    panel.background = element_blank(),
    legend.position = 'none'
    # legend.position = 'inside',
    # legend.position.inside = c(0.8, 0.8)
  )

p2 = pp +
  # Vertical lines for beta sums
  geom_line(
    data = ltre.lambda %>% filter(!trt.phen %in% 'control'),
    aes(y = lambda, colour = trt.phen),
    linetype = 2
  ) +
  # Vertical lines for alpha sums
  geom_line(
    data = ltre.lambda %>% filter(trt.rate %in% 'control'),
    aes(x = as.Date(phen.ctrl.mean), y = lambda),
    linetype = 2
  ) +
  # Horizontal lines for alpha (phen shift) - drought
  geom_line(
    data = ltre.lambda %>% filter(trt.rate %in% 'control', !trt.phen %in% 'irrigated'),
    aes(
      x = phen.date, 
      y = ltre.lambda$lambda[ltre.lambda$trt.rate %in% 'control' & ltre.lambda$trt.phen %in% 'drought']
    ),
    linetype = 2
  ) +
  # Horizontal lines for alpha (phen shift) - irrigated
  geom_line(
    data = ltre.lambda %>% filter(trt.rate %in% 'control', !trt.phen %in% 'drought'),
    aes(
      x = phen.date, 
      y = ltre.lambda$lambda[ltre.lambda$trt.rate %in% 'control' & ltre.lambda$trt.phen %in% 'irrigated']
    ),
    linetype = 2
  ) +
  annotate(
    'text', x = as.Date(124.75), y = 0.9405,
    label = expression(Sigma ~ alpha), colour = 'dodgerblue',
    vjust = 'center', hjust = 'center'
  ) +
  annotate(
    'text', x = as.Date(125.75), y = 0.942,
    label = expression(Sigma ~ alpha), colour = 'goldenrod', 
    vjust = 'center', hjust = 'center'
  ) +
  annotate(
    'text', x = as.Date(121.5), y = 0.946, 
    label = expression(Sigma ~ beta), colour = 'goldenrod',
    vjust = 'center', hjust = 'center'
  ) +
  annotate(
    'text', x = as.Date(127.5), y = 0.9425,
    label = expression(Sigma ~ beta), colour = 'dodgerblue',
    vjust = 'center', hjust = 'center'
  )

plot_grid(
  get_plot_component(lambda.trt.pan, 'guide-box', return_all = TRUE)[[4]],
  plot_grid(
    lambda.trt.pan + labs(x = '') + theme(legend.position = 'none'), 
    p2 + labs(x = '', y = ''), 
    labels = c('a)', 'b)'),
    align = 'v', nrow = 1
  ),
  rel_heights = c(0.1, 1), nrow = 2
) %>%
  save_plot(filename = '~/Desktop/eg_figfig.png', base_height = 5, base_width = 8)
