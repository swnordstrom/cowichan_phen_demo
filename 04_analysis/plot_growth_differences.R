# Script for visualizing differences among treatments in mean inter-annual
# growth rates
# This script exports Fig. S2 

library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyr)

rm(list = ls())

grow.boots = read.csv('~/Desktop/deterministic_growsurv_bootstrapped_perturbed_params.csv') %>%
  pivot_longer(-boot, names_to = 'param', values_to = 'estimate') %>%
  separate_wider_delim(param, delim = '_', names = c('param', 'trt'))

grow.preds = grow.boots %>%
  pivot_wider(names_from = param, values_from = estimate) %>%
  merge(data.frame(size = (5:60)/10)) %>%
  mutate(p.g  = grow.int + grow.slope * size)

# Plot expectations of growth for each bootstrapped model
grow.preds %>%
  ggplot(aes(x = size, y = p.g, colour = trt, group = interaction(boot, trt))) +
  geom_line(linewidth = 0.01) +
  scale_colour_manual(values = c('black', 'goldenrod1', 'dodgerblue')) +
  theme(panel.background = element_blank())

# Get the predicted growth for each treatment for each bootstrap,
# and then get differences between treatments and control
pred.diffs = grow.preds %>%
  select(boot, trt, size, p.g) %>%
  # Getting differences between treatment and control
  pivot_wider(names_from = trt, values_from = p.g) %>%
  mutate(d.c = drought - control, i.c = irrigated - control) %>%
  select(boot, size, d.c, i.c) %>%
  pivot_longer(c(d.c, i.c), names_to = 'contrast', values_to = 'grow.diff') %>%
  mutate(contr.pretty = paste(ifelse(contrast %in% 'd.c', 'drought', 'irrigated'), 'vs. control'))

pan.a = pred.diffs %>%
  ggplot(aes(x = size, y = grow.diff, colour = contr.pretty, group = interaction(boot, contrast))) +
  annotate('segment', x = 0.5, xend = 6, y = 0, yend = 0, linetype = 2, colour = 'gray44') +
  geom_line(linewidth = 0.025) +
  scale_colour_manual(values = c('goldenrod1', 'dodgerblue')) +
  labs(
    x = '', y = 'difference in size, year t\n'
  ) +
  theme(
    panel.background = element_blank(),
    legend.position = 'none',
    axis.title = element_text(size = 6),
    axis.text = element_text(size = 6)
  )

pan.b = pred.diffs %>%
  group_by(contr.pretty, size) %>%
  summarise(p.larger = mean(grow.diff > 0)) %>%
  ungroup() %>%
  mutate(signif = p.larger < 0.025 | p.larger > 0.975) %>%
  ggplot(aes(x = size, y = p.larger, group = contr.pretty, colour = contr.pretty)) +
  geom_line(linewidth = 0.5) +
  geom_point(aes(shape = signif), size = 2) +
  scale_shape_manual(values = c(1, 19)) +
  scale_colour_manual(values = c('goldenrod1', 'dodgerblue')) +
  labs(
    x = 'size, year t-1', 
    y = 'fraction of bootstrapps where\n treatment exceeds control'
  ) +
  theme(
    panel.background = element_blank(),
    legend.position = 'none',
    axis.title = element_text(size = 6),
    axis.text = element_text(size = 6)
  )


plt.legend = get_legend(
  pan.a + 
    guides(col = guide_legend(title = '', override.aes = list(linewidth = 0.5))) + 
    theme(legend.position = 'top', legend.text = element_text(size = 6))
)


plot_grid(
  plt.legend, pan.a, pan.b, nrow = 3, rel_heights = c(0.1, 1, 1),
  align = 'h', axis = 'l'
) %>%
  save_plot(
    filename = '04_analysis/figures/fig_s2_growth_differences.png', 
    base_width = 3, base_height = 5
  )
