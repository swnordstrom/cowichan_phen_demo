library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpattern)

# using ggpattern package for fill
# https://cran.r-project.org/web/packages/ggpatstern/vignettes/patterns-stripes.html

backbone = expand.grid(
  contrast = c('drought', 'irrigated'),
  contrate = c(
    "alpha[surv]",
    "alpha[grow]",
    "alpha[umbl]",
    "alpha[succ]",
    "alpha[seed]",
    "beta[succ]",
    "beta[seed]"
  )
)

set.seed(920)

plt.backbone = backbone %>%
  mutate(contvalu = ifelse(grepl('surv', contrate), 0, rnorm(nrow(.), sd = 0.1)))
         
plt.backbone %>%
  mutate(
    ltre = ifelse(grepl('alpha', contrate), 'trt', 'phen'),
    contrast = paste(contrast, 'vs. control')
  ) %>%
  ggplot(aes(x = contrate, y = contvalu)) +
  geom_col_pattern(
    aes(fill = contrast, pattern_spacing = ltre), # pattern_angle = ltre), 
    pattern_fill = 'gray33', pattern_density = 0.05
  ) +
  # scale_pattern_angle_manual(values = c(30, -30)) +
  # scale_x_reverse() +
  scale_fill_manual(values = c('red', 'blue')) +
  scale_pattern_spacing_manual(values = c(0.05, 0.1)) +
  scale_x_discrete(labels = scales::label_parse()) +
  facet_wrap(~ contrast, nrow = 2) +
  labs(x = 'vital rate contribution', y = 'contribution value') +
  theme(panel.background = element_blank(), legend.position = 'none')

ggsave(
  '04_analysis/figures/draft_figures/eg_ltre_fig.png',
  width = 5, height = 5
)
