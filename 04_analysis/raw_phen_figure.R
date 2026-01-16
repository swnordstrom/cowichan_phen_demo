### Supplemental figure for displaying observed phenology data

# Load in packages
library(ggplot2)
library(dplyr)
library(tidyr)

# Clear namespace
rm(list = ls())

# This script will load the phenology data as `phen`
source('03_construct_kernels/prepare_demo_data_repr.R')

# Aggregate by plant (because that is the unit of replication in analysis)
phen_by_plant = phen |>
  group_by(plantid, Year, trt) |>
  summarise(mean.phen = mean(phen.julian, na.rm = TRUE)) |>
  ungroup() |>
  mutate(trt = factor(trt, levels = c('drought', 'control', 'irrigated')))

# Manually set breaks for date axis
plot_breaks = as.Date(
  paste0(c('03-15', '04-01', '04-15', '05-01', '05-15', '06-01', '06-15'), '-1970'),
  format = '%m-%d-%Y'
)

# Make the plot
phen_by_plant |>
  ggplot(aes(x = Year, y = mean.phen, group = trt, colour = trt)) +
  geom_point(
    aes(shape = trt),
    size = 2.75, stroke = 0.5,
    position = position_jitterdodge(
      jitter.height = 0.5, jitter.width = 0.5, dodge.width = 0.75
    )
  ) +
  scale_shape_manual(values = c(6, 1, 2), 'treatment') +
  scale_colour_manual(values = c('goldenrod1', 'black', 'dodgerblue'), 'treatment') +
  scale_y_continuous(
    breaks = plot_breaks,
    labels = format.Date(plot_breaks, '%b %d')
  ) +
  labs(x = 'Year', y = 'Day of year') +
  theme(
    panel.background = element_blank(),
    legend.position = 'top'
  ) 

# Consider flipping axes (would then want to totally reverse factor order)

ggsave('04_analysis/figures/fig_supp_phen_raw.png', width = 8, height = 5)  



