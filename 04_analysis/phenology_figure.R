# Script for creating Fig. 2 of manuscript,, summarizing treatment effects on
# flowering phenology.

# Load in packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(glmmTMB)

# Clear namespace
rm(list = ls())

# This script will load the phenology data as `phen`
source('03_construct_kernels/prepare_demo_data_repr.R')

# Read in the phenology effects and bootstrapped estimates:
phen_boot = read.csv('03_construct_kernels/out/phenology_bootstrapped_means.csv')

# Aggregate by plant (because that is the unit of replication in analysis)
phen_by_plant = phen |>
  group_by(plantid, Year = as.factor(Year), Plot = as.factor(Plot), trt) |>
  summarise(mean.phen = mean(phen.julian, na.rm = TRUE)) |>
  ungroup()

# Manually set breaks for date axis
plot_breaks = as.Date(
  # paste0(c('03-15', '04-01', '04-15', '05-01', '05-15', '06-01', '06-15'), '-1970'),
  paste0(c('04-01', '04-15', '05-01', '05-15', '06-01', '06-15'), '-1970'),
  format = '%m-%d-%Y'
)

# Get estimated mean treatment effects on phenology
# (object phen.treatment.means is loaded in from `prepare_demo_data_repr.R`
# script above)
phen_trt_effects = phen.treatment.means |>
  pivot_wider(names_from = trt, values_from = mean.phen) |>
  mutate(across(everything(), ~ . - control)) |>
  select(-control) |>
  pivot_longer(everything(), names_to = 'contrast', values_to = 'true_d_phen')

# Get bootstrapped treatment-phenology effect sizes:
phen_boot_contrasts = phen_boot |>
  mutate(across(c(control, drought, irrigated), ~ . - control)) |>
  select(-control) |>
  pivot_longer(c(drought, irrigated), names_to = 'contrast', values_to = 'boot_d_phen') |>
  # Merge in the true effect size:
  merge(phen_trt_effects)

# Get intervals on the treatment effects:
phen_boot_contrast_intervals = phen_boot_contrasts |>
  group_by(contrast, true_d_phen) |>
  reframe(
    ci = quantile(boot_d_phen, probs = c(0.025, 0.975)),
    hilo = c('lo', 'hi')
  ) |>
  pivot_wider(names_from = hilo, values_from = ci)

# Just in case we also want annual means:
d_t = glmmTMB(
  mean.phen ~ trt + Year + (1 | Plot / plantid),
  data = phen_by_plant
)

# Backbone for generating model predictions
annual_means = expand.grid(
  Year = factor(2021:2024), 
  trt = c('control', 'drought', 'irrigated')
)

# Add model predictions
annual_means = annual_means |>
  mutate(
    mean.phen = predict(
      d_t, newdata = annual_means, allow.new.levels = TRUE, re.form = ~ 0
    ),
    # modifying the factor order for plot aesthetics
    trt = factor(trt, levels = c('irrigated', 'control', 'drought')),
    Year = factor(Year, levels = 2024:2021)
  )


### Make the plot

# Panel with raw phenology data (plus annual means)
pan_raw = phen_by_plant |>
  mutate(
    # modifying the factor order for plot aesthetics
    trt = factor(trt, levels = c('irrigated', 'control', 'drought')),
    Year = factor(Year, levels = 2024:2021)
  ) |>
  ggplot(aes(x = Year, y = mean.phen, group = trt)) +
  geom_point(
    # aes(shape = trt),
    aes(colour = trt),
    size = 0.5, stroke = 0.25, alpha = 0.25,
    position = position_jitterdodge(
      jitter.height = 0.125, jitter.width = 0.5, dodge.width = 0.75
    )
  ) +
  geom_point(
    data = annual_means, inherit.aes = TRUE,
    aes(shape = trt, fill = trt),
    size = 2, position = position_dodge(width = 0.75)
  ) +
  scale_shape_manual(
    # values = c(1, 2, 6),
    values = c(21, 24, 25),
    breaks = c('control', 'drought', 'irrigated'),
    'treatment'
  ) +
  # scale_colour_manual(values = c('goldenrod1', 'black', 'dodgerblue'), 'treatment') +
  scale_colour_manual(
    # values = c('dodgerblue', 'black', 'goldenrod1'), 
    values = c('black', 'goldenrod1', 'dodgerblue'), 
    breaks = c('control', 'drought', 'irrigated'),
    'treatment'
  ) +
  scale_fill_manual(
    # values = c('dodgerblue', 'black', 'goldenrod1'), 
    values = c('black', 'goldenrod1', 'dodgerblue'), 
    breaks = c('control', 'drought', 'irrigated'),
    'treatment'
  ) +
  scale_y_continuous(
    breaks = plot_breaks,
    labels = format.Date(plot_breaks, '%b %d')# ,
    # limits = as.numeric(c(as.Date('03-01-1970', format = '%m-%d-%Y'), max(plot_breaks)))
    
  ) +
  labs(x = 'Year', y = 'Day of year') +
  theme(
    panel.background = element_blank(),
    legend.position = 'none',
    axis.ticks = element_line(linewidth = 0.125),
    text = element_text(size = 6)
  ) +
  coord_flip()

# Consider flipping axes (would then want to totally reverse factor order)

# Panel with summary statistcs for treatment effects
pan_stat = phen_boot_contrast_intervals |>
  mutate(contrast = factor(contrast, levels = c('irrigated', 'drought'))) |>
  ggplot(aes(x = contrast)) +
  geom_segment(aes(xend = contrast, y = lo, yend = hi), linewidth = 0.25) +
  # geom_point(
  #   aes(y = boot_d_phen, colour = contrast, shape = contrast), 
  #   position = position_jitter(width = 0.2), alpha = 0.5
  # ) +
  geom_point(aes(y = true_d_phen, fill = contrast, shape = contrast), size = 3) +
  scale_shape_manual(
    values = 24:25,
    breaks = c('drought', 'irrigated'),
  ) +
  scale_fill_manual(
    values = c('goldenrod1', 'dodgerblue'),
    breaks = c('drought', 'irrigated')
  ) +
  scale_colour_manual(
    values = c('goldenrod1', 'dodgerblue'), 
    breaks = c('drought', 'irrigated'),
  ) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.25, colour = 'gray88') +
  labs(
    x = 'Treatment',
    y = 'Effect on mean flowering date (days)',
  ) +
  theme(
    panel.background = element_blank(),
    legend.position = 'none',
    axis.ticks = element_line(linewidth = 0.125),
    text = element_text(size =  6)
  ) +
  coord_flip()

# Legend (for treatment)
trt_legend = get_legend(pan_raw + theme(legend.position = 'top'))

# Flag for whether to include the picture of Lomatium bud as inset
# (default will be FALSE because I don't want to include the photo in the repo)

plot_img = FALSE

if (plot_img) {
  
  plot_grid(
    trt_legend,
    plot_grid(
      ggdraw(pan_raw) + draw_image('9.png', halign = 0.2, scale = 0.48, valign = 0.8),
      pan_stat,
      labels = c('a)', 'b)'),
      rel_widths = c(1, .75),
      label_size = 6, nrow = 1
    ),
    nrow = 2, rel_heights = c(0.1, 1)
  ) |>
    save_plot(
      # filename = '04_analysis/figures/Fig2_trt_phenology.tiff',
      filename = '04_analysis/figures/Fig2_trt_phenology.pdf',
      base_width = 11, base_height = 6.75, units = 'cm'
    )
  
} else {
  plot_grid(
    trt_legend,
    plot_grid(pan_raw, pan_stat, nrow = 1, align = 'v'),
    nrow = 2, rel_heights = c(0.1, 1)
  ) # |>
    # save_plot(
    #   filename = '04_analysis/figures/trt_phen_effects.png',
    #   base_width = 11, base_height = 6, units = 'cm'
    # )
}
