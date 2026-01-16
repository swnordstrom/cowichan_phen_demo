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
  mutate(phen.date = as.Date(mean.phen))

annual.phen.dates = read.csv('03_construct_kernels/out/phen_annual_treatment_means.csv') %>%
  # Convert to date format (for plotting)
  mutate(phen.date = as.Date(mean.bud, format = '%b-%d')) %>%
  mutate(
    ydodge = case_when(
      trt %in% 'control' ~ 0.9225,
      trt %in% 'drought' ~ 0.925,
      trt %in% 'irrigated' ~ 0.92
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
    aes(y = lambda, fill = trt.rate), size = 4, shape = 21
  ) +
  # scale_shape_manual(values = c(NA, 19)) +
  scale_colour_manual(values = c('black', 'goldenrod1', 'dodgerblue'), '') +
  scale_fill_manual(values = c('black', 'goldenrod1', 'dodgerblue'), '') +
  guides(shape = 'none') +
  labs(x = '', y = expression(lambda)) +
  theme(
    panel.background = element_blank(),
    # legend.position = 'top'
    legend.position = 'inside',
    legend.position.inside = c(0.5, 0.85),
    legend.direction = 'horizontal'
  )

# lambda.trt.pan

lambda.contr.pan = boot.lambda.diff %>%
  filter(as.numeric(gsub('b', '', boot)) < 101) %>%
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
  labs(x = '', y = expression(Delta~lambda)) +
  guides(colour = 'none', fill = 'none') +
  # scale_colour_manual(values = c('red', 'blue')) +
  # scale_fill_manual(values = c('red', 'blue')) +
  scale_colour_manual(values = c('goldenrod1', 'dodgerblue')) +
  scale_x_continuous(
    breaks = as.Date(c('1970-04-27', '1970-05-04', '1970-05-11', '1970-05-18')),
    labels = format(as.Date(c('1970-04-27', '1970-05-04', '1970-05-11', '1970-05-18')), '%b %d')
  ) +
  facet_wrap(~ contrast, nrow = 2) +
  theme(
    panel.background = element_blank(),
    strip.background = element_part_rect(fill = 'white', side = 'b', colour = 'gray22'),
    legend.position = 'none'
  )

lambda.legend = get_plot_component(
  lambda.trt.pan + theme(legend.position = 'top'),
  pattern = 'guide-box', return_all = TRUE
)[[4]]

x.ax.lab = ggdraw() +
  draw_label('Flowering date', vjust = 0) +
  theme(plot.margin = margin(0, 0, 10, 0))

plot_grid(
  plot_grid(
    lambda.trt.pan, lambda.contr.pan,
    labels = c('a', 'b'), rel_widths = c(1, 0.5),
    nrow = 1
  ),
  x.ax.lab,
  ncol = 1, rel_heights = c(1, 0.01)
) %>%
  save_plot(
    filename = '04_analysis/figures/draft_figures/lambdas_phen.png',
    base_height = 5, base_width = 8
  )

rm(list = ls())


# ======================================================= #
# LTRE figure
# ======================================================= #

control.ltre.summ = read.csv('04_analysis/out/overall_ltre_summary.csv')

obsv.ltre = read.csv('04_analysis/out/rate-type-combo_ltre_summary.csv')

obsv.by.demo.type = read.csv('04_analysis/out/rate_combo_alone_ltre_summary.csv')
obsv.by.demo = obsv.by.demo.type %>% filter(group %in% c('grow', 'repr')) # %>% rename(demo = group)
obsv.by.type = obsv.by.demo.type %>% filter(group %in% c('phen', 'trt')) # %>% rename(type = group)

pa = control.ltre.summ %>%
  mutate(
    contr.pretty = paste(ifelse(contrast %in% 'd.c', 'drought', 'irrigated'), 'vs. control')
  ) %>%
  ggplot(aes(x = ltre.varb)) +
  # geom_col_pattern(
  #   aes(y = contrib, fill = contrast, pattern = varb),
  #   colour = 'gray22',
  #   pattern_colour = 'gray22', pattern_fill = 'gray22',
  #   pattern_density = 0.025
  # ) +
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

pa
ggsave('04_analysis/figures/ltre_panel_a.png', width = 8, height = 5)


pb = obsv.ltre %>%
  mutate(
    contr.pretty = paste(ifelse(contrast %in% 'd.c', 'drought', 'irrigated'), 'vs. control'),
    demo = ifelse(demo %in% 'grow', 'growth', 'reproduction')
  ) %>%
  ggplot(aes(x = type, y = contrib)) +
  geom_col(aes(fill = contr.pretty), colour = 'gray22') +
  geom_segment(
    aes(xend = type, y = lo, yend = hi),
    linewidth = 1.2
  ) +
  scale_fill_manual(values = c('goldenrod1', 'dodgerblue')) +
  scale_x_discrete(
    limits = c('phen', 'trt'), labels = c('phenology', 'treatment'), guide = guide_axis(n.dodge = 2)
  ) +
  guides(fill = 'none') +
  labs(x = '', y = expression(paste('Contribution to ', Delta, lambda))) +
  facet_nested( ~ contr.pretty + demo) +
  theme(
    strip.background = element_part_rect(fill = 'white', side = 'b', colour = 'gray22'),
    panel.background = element_blank(),
    axis.text = element_text(size = 7),
    strip.text = element_text(size = 7)
  )


pb
ggsave('04_analysis/figures/ltre_panel_b.png', width = 5, height = 3)

pc = obsv.by.type %>%
  rename(type = group) %>%
  ggplot(aes(x = type, y = contrib)) +
  geom_col(aes(fill = contrast), colour = 'gray22') +
  geom_segment(
    aes(xend = type, y = lo, yend = hi),
    linewidth = 1.2
  ) +
  scale_fill_manual(values = c('goldenrod1', 'dodgerblue')) +
  scale_colour_manual(values = c('goldenrod1', 'dodgerblue')) +
  scale_x_discrete(
    limits = c('phen', 'trt'), labels = c('phenology', 'treatment'),
    guide = guide_axis(n.dodge = 2)
  ) +
  scale_y_continuous(limits = c(-0.008, 0.0215)) +
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

pd = obsv.by.demo %>%
  rename(demo = group) %>%
  ggplot(aes(x = demo, y = contrib)) +
  geom_col(aes(fill = contr.pretty), colour = 'gray22') +
  geom_segment(
    aes(xend = demo, y = lo, yend = hi),
    linewidth = 1.2
  ) +
  scale_fill_manual(values = c('goldenrod1', 'dodgerblue')) +
  scale_colour_manual(values = c('goldenrod1', 'dodgerblue')) +
  scale_x_discrete(labels = c('growth', 'reproduction'), guide = guide_axis(n.dodge = 2)) +
  scale_y_continuous(limits = c(-0.008, 0.0215)) +
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


rm(list = ls())


# ======================================================= #
# Mirrored LTRE figures
# ======================================================= #

mirrored.ltre = read.csv('04_analysis/out/out/mirror_ltre_results.csv')

mirrored.ltre.summ = merge(
  mirrored.ltre %>% filter(samp %in% 'obsv') %>% select(-c(rate, samp)),
  mirrored.ltre %>%
    filter(samp %in% 'boot') %>%
    group_by(contrast, ltre.varb, varb) %>%
    reframe(
      cilim = quantile(contrib, probs = c(0.025, 0.975)),
      lohi = c('lo', 'hi')
    ) %>%
    pivot_wider(names_from = lohi, values_from = cilim)
) %>%
  ungroup()

# LTRE Panel a) (all vital rate contributions)

pa = mirrored.ltre.summ %>%
  mutate(
    contr.pretty = paste(ifelse(contrast %in% 'd.c', 'drought', 'irrigated'), 'vs. control')
  ) %>%
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

pa
ggsave('04_analysis/figures/ltre_mirror_panel_a.png', width = 8, height = 5)


obsv.by.demo.type = mirrored.ltre %>%
  filter(samp %in% 'obsv') %>%
  group_by(contrast, demo = ifelse(rate %in% c('grow', 'recr'), 'grow', 'repr'), type) %>%
  summarise(contrib = sum(contrib)) %>%
  ungroup()

boot.by.demo.type = mirrored.ltre %>%
  filter(samp %in% 'boot') %>%
  group_by(contrast, samp, demo = ifelse(rate %in% c('grow', 'recr'), 'grow', 'repr'), type) %>%
  summarise(contrib = sum(contrib)) %>%
  ungroup()

boot.by.demo.type.summ = boot.by.demo.type %>%
  group_by(contrast, demo, type) %>%
  reframe(
    cilim = quantile(contrib, probs = c(0.025, 0.975)),
    lohi = c('lo', 'hi')
  ) %>%
  pivot_wider(names_from = lohi, values_from = cilim) %>%
  mutate(contr.pretty = paste(ifelse(contrast %in% 'd.c', 'drought', 'irrigated'), 'vs. control'))

pb = obsv.by.demo.type %>%
  mutate(
    contr.pretty = paste(ifelse(contrast %in% 'd.c', 'drought', 'irrigated'), 'vs. control'),
    demo = ifelse(demo %in% 'grow', 'growth', 'reproduction')
  ) %>%
  ggplot(aes(x = type, y = contrib)) +
  geom_col(aes(fill = contr.pretty), colour = 'gray22') +
  geom_segment(
    data = boot.by.demo.type.summ %>% mutate(demo = ifelse(demo %in% 'grow', 'growth', 'reproduction')),
    aes(xend = type, y = lo, yend = hi),
    linewidth = 1.2
  ) +
  scale_fill_manual(values = c('goldenrod1', 'dodgerblue')) +
  scale_x_discrete(
    limits = c('phen', 'trt'), labels = c('phenology', 'treatment'),
    guide = guide_axis(n.dodge = 2)
  ) +
  guides(fill = 'none') +
  labs(x = '', y = expression(paste('Contribution to ', Delta, lambda))) +
  facet_nested( ~ contr.pretty + demo) +
  theme(
    strip.background = element_part_rect(fill = 'white', side = 'b', colour = 'gray22'),
    panel.background = element_blank(),
    axis.text = element_text(size = 7),
    strip.text = element_text(size = 7)
  )


# pb
# ggsave('04_analysis/figures/ltre_mirror_panel_b.png', width = 5, height = 3)


obsv.by.demo = obsv.by.demo.type %>%
  group_by(demo, contrast) %>%
  summarise(contrib = sum(contrib)) %>%
  ungroup()

obsv.by.type = obsv.by.demo.type %>%
  group_by(type, contrast) %>%
  summarise(contrib = sum(contrib)) %>%
  ungroup()

boot.by.demo = boot.by.demo.type %>%
  group_by(demo, contrast, samp) %>%
  summarise(contrib = sum(contrib)) %>%
  ungroup()

boot.by.type = boot.by.demo.type %>%
  group_by(type, contrast, samp) %>%
  summarise(contrib = sum(contrib)) %>%
  ungroup()

boot.by.demo.summ = boot.by.demo %>%
  group_by(contrast, demo) %>%
  reframe(
    cilim = quantile(contrib, probs = c(0.025, 0.975)),
    lohi = c('lo', 'hi')
  ) %>%
  pivot_wider(names_from = lohi, values_from = cilim) %>%
  mutate(
    contr.pretty = paste(ifelse(contrast %in% 'd.c', 'drought', 'irrigated'), '\nvs. control')
  )

boot.by.type.summ = boot.by.type %>%
  group_by(contrast, type) %>%
  reframe(
    cilim = quantile(contrib, probs = c(0.025, 0.975)),
    lohi = c('lo', 'hi')
  ) %>%
  pivot_wider(names_from = lohi, values_from = cilim) %>%
  mutate(
    contr.pretty = paste(ifelse(contrast %in% 'd.c', 'drought', 'irrigated'), '\nvs. control')
  )

pc = obsv.by.type %>%
  mutate(
    contr.pretty = paste(ifelse(contrast %in% 'd.c', 'drought', 'irrigated'), '\nvs. control')
  ) %>%
  ggplot(aes(x = type, y = contrib)) +
  geom_col(aes(fill = contrast), colour = 'gray22') +
  geom_segment(
    data = boot.by.type.summ,
    aes(xend = type, y = lo, yend = hi),
    linewidth = 1.2
  ) +
  scale_fill_manual(values = c('goldenrod1', 'dodgerblue')) +
  scale_colour_manual(values = c('goldenrod1', 'dodgerblue')) +
  scale_x_discrete(
    limits = c('phen', 'trt'), labels = c('phenology', 'treatment'),
    guide = guide_axis(n.dodge = 2)
  ) +
  scale_y_continuous(limits = c(-0.008, 0.0215)) +
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

pd = obsv.by.demo %>%
  mutate(
    contr.pretty = paste(ifelse(contrast %in% 'd.c', 'drought', 'irrigated'), '\nvs. control')
  ) %>%
  ggplot(aes(x = demo, y = contrib)) +
  geom_col(aes(fill = contr.pretty), colour = 'gray22') +
  geom_segment(
    data = boot.by.demo.summ,
    aes(xend = demo, y = lo, yend = hi),
    linewidth = 1.2
  ) +
  scale_fill_manual(values = c('goldenrod1', 'dodgerblue')) +
  scale_colour_manual(values = c('goldenrod1', 'dodgerblue')) +
  scale_x_discrete(labels = c('growth', 'reproduction'), guide = guide_axis(n.dodge = 2)) +
  scale_y_continuous(limits = c(-0.008, 0.0215)) +
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

plot_grid(pc, pd, labels = c('i', 'ii'), rel_widths = c(1, 1), align = 'vh')
ggsave('04_analysis/figures/ltre_mirror_panel_c.png', width = 5, height = 3)
