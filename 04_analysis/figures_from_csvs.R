library(ggplot2)
library(cowplot)
library(ggh4x)
library(dplyr)
library(tidyr)

rm(list = ls())

# ======================================================= #
# Lambda-phenology figure
# ======================================================= #

all.lambda = read.csv('04_analysis/out/all_dates_lambda.csv') %>%
  mutate(phen.date = as.Date(phen))

all.boot.lambda = read.csv('04_analysis/out/all_dates_bootstrapped_lambda.csv') %>%
  mutate(phen.date = as.Date(phen))

all.boot.intervals = all.boot.lambda %>%
  group_by(trt, phen.date) %>%
  reframe(
    cibound = quantile(lambda, c(0.025, 0.975)),
    lohi = c('lo', 'hi')
  ) %>%
  pivot_wider(names_from = lohi, values_from = cibound)

ltre.lambda = read.csv('04_analysis/out/ltre_design_lambda.csv') %>%
  mutate(phen.date = as.Date(phen)) %>%
  rename(trt.rate = trt)

annual.phen.dates = read.csv('03_construct_kernels/out/phen_annual_treatment_means.csv') %>%
  # Convert to date format (for plotting)
  mutate(phen.date = as.Date(mean.bud, format = '%b-%d')) %>%
  mutate(
    ydodge = case_when(
      trt %in% 'control' ~ 0.9725,
      trt %in% 'drought' ~ 0.975,
      trt %in% 'irrigated' ~ 0.97
    )
  )

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

lambda.trt.pan = all.lambda %>%
  ggplot(aes(x = phen.date)) +
  annotate(
    'segment', linetype = 2, colour = 'gray',
    x = min(all.lambda$phen.date), xend = max(all.lambda$phen.date),
    y = 1, yend = 1, linewidth = 0.25
  ) +
  geom_ribbon(
    data = all.boot.intervals,
    aes(x = phen.date, ymin = lo, ymax = hi, fill = trt, group = trt),
    alpha = 0.125, show.legend = FALSE
  ) +
  geom_line(aes(y = lambda, colour = trt, group = trt)) +
  geom_point(
    data = annual.phen.dates,
    aes(x = phen.date, y = ydodge, colour = trt),
    size = 4, shape = '*', show.legend = FALSE
  ) +
  geom_point(
    data = ltre.lambda %>% filter(trt.rate == trt.phen),
    aes(y = lambda, fill = trt.rate, shape = trt.rate), 
    size = 1.2
  ) +
  scale_shape_manual(values = c(21, 24, 25), '') +
  scale_colour_manual(values = c('black', 'goldenrod1', 'dodgerblue'), '') +
  scale_fill_manual(values = c('black', 'goldenrod1', 'dodgerblue'), '') +
  labs(x = '', y = expression(lambda)) +
  theme(
    panel.background = element_blank(),
    legend.position = 'inside',
    legend.position.inside = c(0.5, 0.9),
    legend.direction = 'horizontal',
    legend.background = element_blank(),
    legend.key.size = unit(0.2, 'cm'),
    legend.key.spacing = unit(0.1, 'cm'),
    axis.ticks = element_line(linewidth = 0.125),
    text = element_text(size = 6)
  )

# lambda.trt.pan

lambda.contr.pan = merge(boot.lambda.diff, boot.lambda.diff.interval) %>%
  filter(d.lambda < lo | d.lambda > hi) %>%
  mutate(
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
    linetype = 2, colour = 'gray', linewidth = 0.25
  ) +
  geom_point(
    data = boot.lambda.diff.interval %>%
      mutate(
        contrast = ifelse(
          contrast %in% 'd.c',
          'i) drought vs. control',
          'ii) irrigated vs. control'
        )
      ),
    aes(y = mean.d.lambda),
    size = 1.2, shape = 21
  ) +
  geom_line(
    data = boot.lambda.diff.interval %>%
      mutate(
        contrast = ifelse(
          contrast %in% 'd.c',
          'i) drought vs. control',
          'ii) irrigated vs. control'
        )
      ),
    aes(y = mean.d.lambda),
    linewidth = 0.25
  ) +
  geom_point(
    # y limits here found by:
    # ggplot_build(lambda.contr.pan)$layout$panel_scales_y
    data = annual.phen.dates %>%
      filter(!(trt %in% 'control')) %>%
      mutate(
        contrast = ifelse(
          trt %in% 'drought',
          'i) drought vs. control',
          'ii) irrigated vs. control'
        )
      ),
    aes(x = phen.date, y = -0.020, colour = contrast),
    shape = '*', size = 4
  ) +
  geom_ribbon(
    data = boot.lambda.diff.interval %>%
          mutate(
            contrast = ifelse(
              contrast %in% 'd.c',
              'i) drought vs. control',
              'ii) irrigated vs. control'
            )
          ),
    aes(x = phen.date, ymin = lo, ymax = hi, fill = contrast),
    alpha = 0.125
  ) +
  labs(x = '', y = expression(Delta~lambda)) +
  guides(colour = 'none', fill = 'none') +
  scale_colour_manual(values = c('goldenrod1', 'dodgerblue')) +
  scale_fill_manual(values = c('goldenrod1', 'dodgerblue')) +
  scale_x_continuous(
    breaks = as.Date(c('1970-04-27', '1970-05-04', '1970-05-11', '1970-05-18')),
    labels = format(as.Date(c('1970-04-27', '1970-05-04', '1970-05-11', '1970-05-18')), '%b %d')
  ) +
  facet_wrap(~ contrast, nrow = 2) +
  theme(
    panel.background = element_blank(),
    strip.background = element_part_rect(fill = 'white', side = 'b', colour = 'gray22'),
    legend.position = 'none',
    axis.ticks = element_line(linewidth = 0.125),
    text = element_text(size = 6)
  )

x.ax.lab = ggdraw() +
  draw_label('Flowering date', size = 7, vjust = 0) +
  theme(plot.margin = margin(0, 0, 10, 0))

plot_grid(
  plot_grid(
    lambda.trt.pan, lambda.contr.pan,
    labels = c('a', 'b'), label_size = 6,
    rel_widths = c(1, 0.75),
    nrow = 1
  ),
  x.ax.lab,
  ncol = 1, rel_heights = c(1, 0.01)
) %>%
  save_plot(
    filename = '04_analysis/figures/Fig3.tiff',
    base_width = 10, base_height = 6, units = 'cm'
  )

rm(list = ls())

# ======================================================= #
# LTRE figure
# ======================================================= #

# control.ltre.summ = read.csv('04_analysis/out/overall_ltre_summary.csv')
# 
# obsv.ltre = read.csv('04_analysis/out/rate-type-combo_ltre_summary.csv')
# 
# obsv.by.demo.type = read.csv('04_analysis/out/rate_combo_alone_ltre_summary.csv')
# obsv.by.demo = obsv.by.demo.type %>% filter(group %in% c('grow', 'repr')) # %>% rename(demo = group)
# obsv.by.type = obsv.by.demo.type %>% filter(group %in% c('phen', 'trt')) # %>% rename(type = group)

### Reading in plot components:

ltre.trt.contrib  = read.csv('04_analysis/out/trt_ltre_contribs.csv')
ltre.phen.contrib = read.csv('04_analysis/out/phen_ltre_contribs.csv')

ltre.boot.trt.contrib  = read.csv('04_analysis/out/boot_trt_ltre_contribs.csv')
ltre.boot.phen.contrib = read.csv('04_analysis/out/boot_phen_ltre_contribs.csv')

### Combine into a single data frame:

# Mean contributions
ltre.contribs = rbind(
  # Treatment contributions:
  ltre.trt.contrib %>%
    # Main results: date used is flowering date of control
    # Mirrored results: date used is the non-control treatments
    mutate(result.set = ifelse(trt.phen %in% 'control', 'main', 'mirrored')) %>%
    # Label these as the treatment effects
    mutate(effect = 'trt') %>%
    select(result.set, effect, contrast, rate, contrib),
  # Phenology contributions:
  ltre.phen.contrib %>%
    # Main results: treatment is non-control treatments
    # Mirrored results: uses control
    mutate(result.set = ifelse(trt.rate %in% 'control', 'mirrored', 'main')) %>%
    # Label these as the phenology effects
    mutate(effect = 'phen') %>%
    select(result.set, effect, contrast = contrast.phen, rate, contrib)
)

# Bootstrapped contributions
boot.contribs = rbind(
  # Treatment contributions:
  ltre.boot.trt.contrib %>%
    # Main results: date used is flowering date of control
    # Mirrored results: date used is the non-control treatments
    mutate(result.set = ifelse(trt.phen %in% 'control', 'main', 'mirrored')) %>%
    # Label these as the treatment effects
    mutate(effect = 'trt') %>%
    select(result.set, effect, contrast, rate, samp, contrib),
  # Phenology contributions:
  ltre.boot.phen.contrib %>%
    # Main results: treatment is non-control treatments
    # Mirrored results: uses control
    mutate(result.set = ifelse(trt.rate %in% 'control', 'mirrored', 'main')) %>%
    # Label these as the phenology effects
    mutate(effect = 'phen') %>%
    select(result.set, effect, contrast = contrast.phen, rate, samp, contrib)
)


### Panel a: all vital rates 

# Merge together the bootstrapped confidence intervals with the estimates:
ltre.all.rates.summary = merge(
  ltre.contribs,
  boot.contribs %>% 
    group_by(result.set, effect, contrast, rate) %>% 
    reframe(qq = quantile(contrib, probs = c(0.025, 0.975)), lohi = c('lo', 'hi')) %>% 
    pivot_wider(names_from = lohi, values_from = qq)
)  %>%
  # Add some features to make the plot prettier
  mutate(
    ltre.varb = paste0(ifelse(effect %in% 'phen', 'phi', 'psi'), '[', rate, ']'),
    contr.pretty = paste(ifelse(contrast %in% 'd.c', 'drought', 'irrigated'), 'vs. control')
  )

pa = ltre.all.rates.summary %>%
  filter(result.set %in% 'main') %>%
  ggplot(aes(x = ltre.varb)) +
  geom_col(
    aes(y = contrib, fill = contrast), colour = 'gray22'
  ) +
  geom_segment(aes(xend = ltre.varb, y = lo, yend = hi), linewidth = 1.2) +
  scale_x_discrete(
    labels = scales::label_parse(),
    limits = c(
      'phi[grow]', 'phi[succ]', 'phi[seed]',
      'psi[grow]', 'psi[flow]', 'psi[seed]', 'psi[recr]'
    ),
    guide = guide_axis(n.dodge = 2)
  ) +
  # scale_pattern_manual(values = c('stripe', 'crosshatch')) +
  scale_fill_manual(values = c('goldenrod1', 'dodgerblue')) +
  facet_wrap(~ contr.pretty) +
  labs(x = '', y = expression(paste('Contribution to ', Delta, lambda))) +
  guides(fill = 'none', pattern = 'none') +
  theme(
    panel.background = element_blank(),
    strip.background = element_part_rect(fill = 'white', side = 'b', colour = 'gray22'),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 10),
    strip.text = element_text(size = 10)
  )

pa
ggsave('04_analysis/figures/ltre_panel_a.png', width = 8, height = 5)


### Panel b: 

# Aggregate vital rates into groups:
ltre.contribs.growth.repr = ltre.contribs %>%
  group_by(result.set, effect, contrast, demo = ifelse(rate %in% c('grow', 'recr'), 'grow', 'repr')) %>%
  summarise(contrib = sum(contrib)) %>%
  ungroup()

boot.contribs.growth.repr = boot.contribs %>%
  group_by(result.set, effect, contrast, demo = ifelse(rate %in% c('grow', 'recr'), 'grow', 'repr'), samp) %>%
  summarise(contrib = sum(contrib)) %>%
  ungroup()

# Combine into single data frame with confidence interavls:
ltre.growth.repr.summary = merge(
  ltre.contribs.growth.repr,
  boot.contribs.growth.repr %>% 
    group_by(result.set, effect, contrast, demo) %>% 
    reframe(qq = quantile(contrib, probs = c(0.025, 0.975)), lohi = c('lo', 'hi')) %>% 
    pivot_wider(names_from = lohi, values_from = qq)
)  %>%
  # Add some features to make the plot prettier
  mutate(
    effect.pretty = ifelse(effect %in% 'phen', 'phenology', 'treatment'),
    # demo.pretty = ifelse(demo %in% 'grow', 'growth', 'reproduction'),
    contr.pretty = paste(ifelse(contrast %in% 'd.c', 'drought', 'irrigated'), 'vs. control')
  )

pb = ltre.growth.repr.summary %>%
  filter(result.set %in% 'main') %>%
  ggplot(aes(x = demo, y = contrib)) +
  geom_col(aes(fill = contr.pretty), linewidth = 0.25, colour = 'gray22') +
  geom_segment(aes(xend = demo, y = lo, yend = hi), linewidth = 0.25) + # ,linewidth = 1.2) +
  scale_fill_manual(values = c('goldenrod1', 'dodgerblue')) +
  scale_x_discrete(
    limits = c('grow', 'repr'), labels = c('growth', 'reproduction'), guide = guide_axis(n.dodge = 2)
  ) +
  guides(fill = 'none') +
  labs(x = '', y = expression(paste('Contribution to ', Delta, lambda))) +
  facet_nested( ~ contr.pretty + effect.pretty) +
  theme(
    strip.background = element_part_rect(fill = 'white', side = 'b', colour = 'gray22'),
    panel.background = element_blank(),
    text = element_text(size = 6),
    axis.ticks = element_line(linewidth = 0.125),
    # axis.text = element_text(size = 4),
    # strip.text = element_text(size = 4)
  )


pb
ggsave('04_analysis/figures/Fig4.tiff', width = 8, height = 5, units = 'cm')


### Panel c: 

# Summed into cumulative effects of effect type (i.e., aggregating over demographic pathways)
# (used for left panel)

ltre.contribs.by.effect = ltre.contribs.growth.repr %>%
  group_by(result.set, effect, contrast) %>%
  summarise(contrib = sum(contrib))

boot.contribs.by.effect = boot.contribs.growth.repr %>%
  group_by(result.set, effect, contrast, samp) %>%
  summarise(contrib = sum(contrib))

# Single data frame with confidence intervals:
ltre.by.effect.summary = merge(
  ltre.contribs.by.effect,
  boot.contribs.by.effect %>% 
    group_by(result.set, effect, contrast) %>% 
    reframe(qq = quantile(contrib, probs = c(0.025, 0.975)), lohi = c('lo', 'hi')) %>% 
    pivot_wider(names_from = lohi, values_from = qq)
)  %>%
  # Add some features to make the plot prettier
  mutate(
    effect.pretty = ifelse(effect %in% 'phen', 'phenology', 'treatment'),
    # demo.pretty = ifelse(demo %in% 'grow', 'growth', 'reproduction'),
    contr.pretty = paste(ifelse(contrast %in% 'd.c', 'drought', 'irrigated'), 'vs. control')
  )

# Summed into cumulative effects of pathway (i.e., aggregating over drivers)

ltre.contribs.by.demo = ltre.contribs.growth.repr %>%
  group_by(result.set, demo, contrast) %>%
  summarise(contrib = sum(contrib))

boot.contribs.by.demo = boot.contribs.growth.repr %>%
  group_by(result.set, demo, contrast, samp) %>%
  summarise(contrib = sum(contrib))

# Single data frame with confidence intervals:
ltre.by.demo.summary = merge(
  ltre.contribs.by.demo,
  boot.contribs.by.demo %>% 
    group_by(result.set, demo, contrast) %>% 
    reframe(qq = quantile(contrib, probs = c(0.025, 0.975)), lohi = c('lo', 'hi')) %>% 
    pivot_wider(names_from = lohi, values_from = qq)
)  %>%
  # Add some features to make the plot prettier
  mutate(
    # effect.pretty = ifelse(effect %in% 'phen', 'phenology', 'treatment'),
    # demo.pretty = ifelse(demo %in% 'grow', 'growth', 'reproduction'),
    contr.pretty = paste(ifelse(contrast %in% 'd.c', 'drought', 'irrigated'), 'vs. control')
  )


pc = ltre.by.effect.summary %>%
  filter(result.set %in% 'main') %>%
  ggplot(aes(x = effect, y = contrib)) +
  geom_col(aes(fill = contrast), colour = 'gray22') +
  geom_segment(
    aes(xend = effect, y = lo, yend = hi),
    linewidth = 1.2
  ) +
  scale_fill_manual(values = c('goldenrod1', 'dodgerblue')) +
  scale_colour_manual(values = c('goldenrod1', 'dodgerblue')) +
  scale_x_discrete(
    limits = c('phen', 'trt'), labels = c('phenology', 'treatment'),
    guide = guide_axis(n.dodge = 2)
  ) +
  scale_y_continuous(limits = c(-0.03, 0.05)) +
  # scale_y_continuous(limits = c(-0.025, 0.0375)) +
  # labs(x = '', y = '') +
  labs(x = '', y = expression(paste('Contribution to ', Delta, lambda))) +
  guides(fill = 'none') +
  facet_wrap( ~ contr.pretty, nrow = 1) +
  theme(
    strip.background = element_part_rect(fill = 'white', side = 'b', colour = 'gray22'),
    panel.background = element_blank(),
    # axis.text.x = element_text(angle = 45),
    # axis.text.y = element_blank(),
    # axis.ticks.y = element_blank(),
    axis.text = element_text(size = 7),
    strip.text = element_text(size = 7)# ,
    # plot.margin = margin(l = 5, r = 0)
  )

pd = ltre.by.demo.summary %>%
  filter(result.set %in% 'main') %>%
  ggplot(aes(x = demo, y = contrib)) +
  geom_col(aes(fill = contr.pretty), colour = 'gray22') +
  geom_segment(
    aes(xend = demo, y = lo, yend = hi),
    linewidth = 1.2
  ) +
  scale_fill_manual(values = c('goldenrod1', 'dodgerblue')) +
  scale_colour_manual(values = c('goldenrod1', 'dodgerblue')) +
  scale_x_discrete(labels = c('growth', 'reproduction'), guide = guide_axis(n.dodge = 2)) +
  scale_y_continuous(limits = c(-0.03, 0.05)) +
  # scale_y_continuous(limits = c(-0.025, 0.0375)) +
  # labs(x = '', y = expression(paste('Contribution to ', Delta, lambda))) +
  labs(x = '', y = '') +
  guides(fill = 'none') +
  facet_wrap( ~ contr.pretty, nrow = 1) +
  theme(
    strip.background = element_part_rect(fill = 'white', side = 'b', colour = 'gray22'),
    panel.background = element_blank(),
    axis.text = element_text(size = 7),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.text = element_text(size = 7)# ,
    # plot.margin = margin(l = 0, r = 5)
  )

# pc
# pd

plot_grid(pc, pd, labels = c('a', 'b'), rel_widths = c(1, 1), align = 'vh') %>%
  save_plot(filename = '04_analysis/figures/ltre_panel_c.png', base_width = 5, base_height = 3)


# ======================================================= #
# Mirrored LTRE figures
# ======================================================= #

pma = ltre.all.rates.summary %>%
  filter(result.set %in% 'mirrored') %>%
  ggplot(aes(x = ltre.varb)) +
  geom_col(
    aes(y = contrib, fill = contrast), colour = 'gray22'
  ) +
  geom_segment(aes(xend = ltre.varb, y = lo, yend = hi), linewidth = 1.2) +
  scale_x_discrete(
    labels = scales::label_parse(),
    limits = c(
      'phi[grow]', 'phi[succ]', 'phi[seed]',
      'psi[grow]', 'psi[flow]', 'psi[seed]', 'psi[recr]'
    ),
    guide = guide_axis(n.dodge = 2)
  ) +
  # scale_pattern_manual(values = c('stripe', 'crosshatch')) +
  scale_fill_manual(values = c('goldenrod1', 'dodgerblue')) +
  facet_wrap(~ contr.pretty) +
  labs(x = '', y = expression(paste('Contribution to ', Delta, lambda))) +
  guides(fill = 'none', pattern = 'none') +
  theme(
    panel.background = element_blank(),
    strip.background = element_part_rect(fill = 'white', side = 'b', colour = 'gray22'),
    axis.text.x = element_text(size = 9),
    axis.text.y = element_text(size = 7),
    strip.text = element_text(size = 7)
  )

pma
ggsave('04_analysis/figures/ltre_mirror_a.png', width = 8, height = 5)


pmb = ltre.growth.repr.summary %>%
  filter(result.set %in% 'mirrored') %>%
  ggplot(aes(x = demo, y = contrib)) +
  geom_col(aes(fill = contr.pretty), colour = 'gray22') +
  geom_segment(
    aes(xend = demo, y = lo, yend = hi),
    linewidth = 1.2
  ) +
  scale_fill_manual(values = c('goldenrod1', 'dodgerblue')) +
  scale_x_discrete(
    limits = c('grow', 'repr'), labels = c('growth', 'reproduction'), guide = guide_axis(n.dodge = 2)
  ) +
  guides(fill = 'none') +
  labs(x = '', y = expression(paste('Contribution to ', Delta, lambda))) +
  facet_nested( ~ contr.pretty + effect.pretty) +
  theme(
    strip.background = element_part_rect(fill = 'white', side = 'b', colour = 'gray22'),
    panel.background = element_blank(),
    axis.text = element_text(size = 7),
    strip.text = element_text(size = 7)
  )


pmb
ggsave('04_analysis/figures/ltre_mirror_b.png', width = 5, height = 3)


# ======================================================= #
# Comparing LTRE contributions with observed Delta lambda
# ======================================================= #

rm(list = ls())

# Observed lambdas
ltre.lambda = read.csv('04_analysis/out/ltre_design_lambda.csv') %>%
  rename(trt.rate = trt)

# LTRE contributions
ltre.trt.contrib  = read.csv('04_analysis/out/trt_ltre_contribs.csv')
ltre.phen.contrib = read.csv('04_analysis/out/phen_ltre_contribs.csv')

### Individual trt/phen estimates
trt.delta.lambda = ltre.lambda %>%
  select(trt.phen, trt.rate, lambda) %>%
  pivot_wider(names_from = trt.rate, values_from = lambda) %>%
  mutate(d.c = drought - control, i.c = irrigated - control) %>%
  mutate(result.set = ifelse(trt.phen %in% 'control', 'main', 'mirrored')) %>%
  select(result.set, d.c, i.c) %>%
  pivot_longer(
    c(d.c, i.c), 
    names_to = 'contrast', values_to = 'Dlambda', 
    values_drop_na = TRUE
  ) %>%
  mutate(contrast.type = 'trt')
phen.delta.lambda = ltre.lambda %>%
  select(trt.phen, trt.rate, lambda) %>%
  pivot_wider(names_from = trt.phen, values_from = lambda) %>%
  mutate(d.c = drought - control, i.c = irrigated - control) %>%
  mutate(result.set = ifelse(trt.rate %in% 'control', 'mirrored', 'main')) %>%
  select(result.set, d.c, i.c) %>%
  pivot_longer(
    c(d.c, i.c),
    names_to = 'contrast', values_to = 'Dlambda',
    values_drop_na = TRUE
  ) %>%
  mutate(contrast.type = 'phen')

# Get the sums of LTRE contributions across all vital rates
ltre.trt.contribs = ltre.trt.contrib %>%
  mutate(result.set = ifelse(trt.phen %in% 'control', 'main', 'mirrored')) %>%
  group_by(contrast, result.set) %>%
  summarise(trt.contrib = sum(contrib)) %>%
  ungroup()
ltre.phen.contribs = ltre.phen.contrib %>%
  mutate(result.set = ifelse(trt.rate %in% 'control', 'mirrored', 'main')) %>%
  group_by(contrast = contrast.phen, result.set) %>%
  summarise(phen.contrib = sum(contrib)) %>%
  ungroup()

compare.trt.dlambda = merge(trt.delta.lambda, ltre.trt.contribs) %>%
  # Relative error
  mutate(reltve.error = (Dlambda - trt.contrib) / Dlambda)
compare.phen.dlambda = merge(phen.delta.lambda, ltre.phen.contribs) %>%
  # Relative error
  mutate(reltve.error = (Dlambda - phen.contrib) / Dlambda)

### Overall contrast estimates

# Get Delta lambdas
ltre.delta.lambda = ltre.lambda %>%
  filter(trt.phen == trt.rate) %>%
  select(trt = trt.phen, lambda) %>%
  pivot_wider(names_from = trt, values_from = lambda) %>%
  mutate(d.c = drought - control, i.c = irrigated - control) %>%
  select(d.c, i.c) %>%
  pivot_longer(c(d.c, i.c), names_to = 'contrast', values_to = 'Dlambda')

compare.contribs.dlambda = merge(ltre.trt.contribs, ltre.phen.contribs) %>%
  merge(ltre.delta.lambda) %>%
  mutate(
    ltre.Dlambda = trt.contrib + phen.contrib,
    reltve.error = (ltre.Dlambda - Dlambda) / Dlambda
  )

compare.contribs.dlambda %>%
  mutate(across(where(is.double), ~ round(., 6)))
