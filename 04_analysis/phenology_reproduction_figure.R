library(ggplot2)
library(dplyr)
library(tidyr)
library(glmmTMB)
library(cowplot)

rm(list = ls())

# Read in data with wrapper script
source('03_construct_kernels/prepare_demo_data_repr.R')


# --------------------------------------
# ----- Fit model ----------------------
# --------------------------------------

# === Seed model ===
# Response: number of seeds (negative binomial distribution) of an umbel
# Predictors: treatment (categorical), year (factor) , mean budding date of
# plant (centered, continuous), plant size (continuous)
s_st.p_s.u.p2 = glmmTMB(
  no.seeds ~ trt * size + Year + phen.c + (1 | Plot / plantid),
  family = 'nbinom2',
  ziformula = ~ size + phen.umbels + Year + poly(phen.c, 2) + (1 | Plot / plantid),
  data = seed
)

# --------------------------------------
# ----- Get model predictions ----------
# --------------------------------------

# First - what is the IQR for sizes of flowering plants?
seed %>%
  distinct(plantid, Year, .keep_all = TRUE) %>%
  reframe(iqrs = quantile(size, probs = c(0.25, 0.5, 0.75)))

seed.lin.preds = expand.grid(
  Year = factor(2021:2024),
  # size = c(3.5, 3.9, 4.3),
  # phen.umbels = c(1, 2, 5),
  phen.umbels = 1,
  size = 3.9,
  trt = c('control', 'drought', 'irrigated'),
  phen.c = (-4:4)*7
) %>%
  mutate(
    zlnk = predict(s_st.p_s.u.p2, newdata = ., allow.new.levels = TRUE, re.form = ~ 0, type = 'zlink'),
    link = predict(s_st.p_s.u.p2, newdata = ., allow.new.levels = TRUE, re.form = ~ 0, type = 'link'),
  )

seed.preds = seed.lin.preds %>%
  group_by(size = factor(size), phen.umbels, trt, phen.c) %>%
  summarise(across(c(zlnk, link), mean)) %>%
  ungroup() %>%
  mutate(
    # Probability of an umbel surviving to make seed
    p.succ = 1 - (1 / (1 + exp(-zlnk))),
    # Expected number of seeds (given that umbel survives)
    s.cond = exp(link),
    # Expected number of seeds()
    n.seed = p.succ * s.cond
  ) %>%
  mutate(
    size = relevel(size, ref = '3.9'),
    phen = phen.c + round(mean(seed$mean.phen)),
    phen.date = as.Date(phen, format = '%b-%d')
  )

# Add a formatted phen date to seed data
seed = seed %>% mutate(phen.date = as.Date(mean.phen, format = '%b-%d'))

# Figure panel a)
# - probability of umbel success/survival
# - data points vertically jittered, coloured by treatment
# - curves giving estimates for each treatment

pan.a = seed.preds %>%
  # filter(phen.umbels < 2) %>%
  ggplot(aes(x = phen.date)) +
  geom_point(
    data = seed %>% 
      filter(!is.na(no.seeds), phen.c > -40) %>%
      mutate(umbel.succ = as.numeric(no.seeds > 0)),
    aes(y = umbel.succ, colour = trt, shape = no.seeds > 0),
    position = position_jitter(height = 0.0625),
    alpha = 0.125, size = 2
  ) +
  geom_line(
    aes(y = p.succ, group = trt)
  ) +
  scale_y_continuous(breaks = (0:4)/4) +
  scale_shape_manual(values = c(4, 19)) +
  # scale_linewidth_manual(values = c(1, 0.25, 0.25)) +
  scale_colour_manual(values = c('black', 'goldenrod', 'dodgerblue')) +
  labs(x = '', y = 'Probability of umbel surviving') +
  theme(
    # axis.text.x = element_blank(),
    legend.position = 'none',
    panel.background = element_blank()
  )

# Figure panel b)
# - Seeds per surviving umbel
# - Data points, lines

pan.b = seed.preds %>%
  # filter(phen.umbels < 2) %>%
  ggplot(aes(x = phen.date)) +
  geom_point(
    data = seed %>% filter(!is.na(no.seeds), phen.c > -40, no.seeds > 0),
    aes(y = no.seeds, colour = trt),
    alpha = 0.125, size = 2
  ) +
  geom_line(
    aes(
      y = s.cond, group = trt, colour = trt
    )
  ) +
  scale_y_log10() +
  # scale_linewidth_manual(values = c(1, 0.5, 0.5)) +
  scale_colour_manual(values = c('black', 'goldenrod', 'dodgerblue')) +
  labs(x = '', y = 'Seeds per surviving umbel') +
  theme(
    # axis.text.x = element_blank(),
    legend.position = 'none',
    panel.background = element_blank()
  )

# wow this fucking blows, fucking awful, fuck

pan.c = seed.preds %>%
  # filter(phen.umbels < 2) %>%
  ggplot(aes(x = phen.date)) +
  geom_point(
    data = seed %>% 
      filter(!is.na(no.seeds), phen.c > -40) %>%
      mutate(
        n.seeds = ifelse(no.seeds > 0, no.seeds, 1/2),
        n.seeds.lab = ifelse(n.seeds < 1, 'failed', 'surviving')
      ),
    aes(y = n.seeds, colour = trt, shape = n.seeds.lab),
    alpha = 0.125, size = 2
  ) +
  geom_line(aes(y = n.seed, group = trt, colour = trt)) +
  scale_y_log10() +
  # scale_linewidth_manual(values = c(1, 0.5, 0.5)) +
  scale_shape_manual(values = c(4, 19)) +
  scale_colour_manual(values = c('black', 'goldenrod', 'dodgerblue')) +
  labs(x = '', y = 'Seeds per umbel') +
  theme(
    legend.position = 'none',
    panel.background = element_blank()
  )

leg.z = get_plot_component(
  pan.c +
    guides(shape = guide_legend('umbel fate')) +
    theme(legend.position = 'top'),
  'guide-box',
  return_all = TRUE
)[[4]]

# okay - something in cowplot must have changed...
# look for the non-empty element of get_plot_component()

plot_grid(
  NULL, leg.z, NULL, pan.a, pan.b, pan.c, byrow = TRUE,
  nrow = 2, rel_heights = c(0.1, 1)
) %>%
  save_plot(
    filename = '04_analysis/figures/draft_figures/phen_reproduction.png',
    base_width = 8, base_height = 5
  )

