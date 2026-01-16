library(ggplot2)
library(dplyr)
library(tidyr)
library(glmmTMB)
library(cowplot)

rm(list = ls())

# Read in data with wrapper script
source('03_construct_kernels/prepare_demo_data_growsurv.R')
source('03_construct_kernels/prepare_demo_data_repr.R')

# --------------------------------------
# ----- Fit model ----------------------
# --------------------------------------

# === Seed model ===
# Response: number of seeds (negative binomial distribution) of an umbel
# Predictors: treatment (categorical), year (factor) , mean budding date of
# plant (centered, continuous), plant size (continuous)
s_st.p_s.u.p = glmmTMB(
  no.seeds ~ trt * size + Year + phen.c + (1 | Plot / plantid),
  family = 'nbinom2',
  ziformula = ~ size + phen.umbels + Year + phen.c + (1 | Plot / plantid),
  data = seed
)

# === Growth model (vegetative) ===
# Response: size of plant after transition (gaussian link)
# Predictors: treatment (categorical), size prior to transition
g_st.ty = glmmTMB(
  size.cur ~ size.prev + size.prev * trt + (1 | prev.year) + (1 | prev.year:trt) + (1 | Plot / plantid),
  data = demo.grow
)

# === Growth model (flowering) ===
# Response: size of plant after transition (gaussian link)
# Predictors: treatment (categorical), size prior to transition, phenology (centered)
g_phen = glmmTMB(
  size.cur ~ size.prev * prev.year + trt * prev.year + phen.c + (1 | Plot / plantid),
  data = demo.grow %>% filter(!is.na(phen.c))
)

# --------------------------------------
# ----- Reproduction  panels -----------
# --------------------------------------

# # First - what is the IQR for sizes of flowering plants?
# seed %>%
#   distinct(plantid, Year, .keep_all = TRUE) %>%
#   reframe(iqrs = quantile(size, probs = c(0.25, 0.5, 0.75)))

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
    zlnk = predict(s_st.p_s.u.p, newdata = ., allow.new.levels = TRUE, re.form = ~ 0, type = 'zlink'),
    link = predict(s_st.p_s.u.p, newdata = ., allow.new.levels = TRUE, re.form = ~ 0, type = 'link'),
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
    aes(y = p.succ, group = trt),
    colour = 'gray77', linewidth = 1.2
  ) +
  scale_y_continuous(breaks = (0:4)/4) +
  scale_shape_manual(values = c(4, 19)) +
  # scale_linewidth_manual(values = c(1, 0.25, 0.25)) +
  scale_colour_manual(values = c('black', 'goldenrod1', 'dodgerblue')) +
  labs(x = '', y = 'Probability of umbel success') +
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
  geom_line(aes(y = s.cond, group = trt, colour = trt), linewidth = 1.2) +
  scale_y_log10() +
  # scale_linewidth_manual(values = c(1, 0.5, 0.5)) +
  scale_colour_manual(values = c('black', 'goldenrod1', 'dodgerblue')) +
  labs(x = '', y = 'Seeds per successful umbel') +
  theme(
    # axis.text.x = element_blank(),
    legend.position = 'none',
    panel.background = element_blank()
  )

# Not great.

pan.c = seed.preds %>%
  # filter(phen.umbels < 2) %>%
  ggplot(aes(x = phen.date)) +
  geom_point(
    data = seed %>% 
      filter(!is.na(no.seeds), phen.c > -40) %>%
      mutate(
        n.seeds = ifelse(no.seeds > 0, no.seeds, 1/2),
        n.seeds.lab = ifelse(n.seeds < 1, 'failed', 'successful')
      ),
    aes(y = n.seeds, colour = trt, shape = n.seeds.lab),
    alpha = 0.125, size = 2
  ) +
  geom_line(aes(y = n.seed, group = trt, colour = trt), linewidth = 1.2) +
  scale_y_log10() +
  # scale_linewidth_manual(values = c(1, 0.5, 0.5)) +
  scale_shape_manual(values = c(4, 19)) +
  scale_colour_manual(values = c('black', 'goldenrod1', 'dodgerblue'), 'treatment') +
  labs(x = '', y = 'Seeds per umbel') +
  theme(
    legend.position = 'none',
    panel.background = element_blank()
  )

leg.z = get_legend(
  pan.c +
    guides(shape = guide_legend('umbel fate'), colour = guide_legend('')) +
    theme(legend.position = 'top')
)

# # okay - something in cowplot must have changed...
# # look for the non-empty element of get_plot_component()
# 
# plot_grid(
#   NULL, leg.z, NULL, pan.a, pan.b, pan.c, byrow = TRUE,
#   labels = c('', '', '', 'a', 'b', 'c'),
#   nrow = 2, rel_heights = c(0.1, 1)
# ) %>%
#   save_plot(
#     filename = '04_analysis/figures/draft_figures/phen_reproduction.png',
#     base_width = 8, base_height = 5
#   )


# --------------------------------------
# ---------- Growth panels -------------
# --------------------------------------


# Figure panel d)
# - vegetative growth kernel
# - data points colored by treatment
# - lines giving growth estimates for each treatment

growth.veg.pred = expand.grid(
  size.prev = (5:60)/10, 
  trt = c('control', 'drought', 'irrigated')
) %>%
  mutate(
    size.cur = predict(g_st.ty, newdata = ., allow.new.levels = TRUE, re.form = ~ 0)
  )

# Plot

# set.seed(3007)

pan.d = demo.grow %>%
  ggplot(aes(x = size.prev, y = size.cur, colour = trt, fill = trt)) +
  geom_point(size = 2, alpha = 0.1) +
  geom_line(
    data = growth.veg.pred,
    aes(group = trt),
    linewidth = 1.2
  ) +
  scale_colour_manual(values = c('black', 'goldenrod1', 'dodgerblue'), 'treatment') +
  labs(x = '', y = 'Size, year t+1') +
  theme(
    legend.position = 'none',
    panel.background = element_blank()
  )

# Figure panel e)
# - flowering growth kernel with phenology
# - data points colored by phenology
# - lines giving growth estimates for each treatment

growth.flow.pred = expand.grid(
  size.prev = (5:60)/10, 
  phen.c = (-2:2) * 7,
  prev.year = factor(2021:2023),
  trt = 'control'
) %>%
  mutate(
    size.cur = predict(g_phen, newdata = ., allow.new.levels = TRUE, re.form = ~ 0)
  ) %>%
  group_by(size.prev, trt, phen.c) %>%
  summarise(size.cur = mean(size.cur))

growth.flow.pred

pan.e = demo.grow %>%
  filter(!is.na(phen.c)) %>%
  ggplot(aes(x = size.prev, y = size.cur)) +
  geom_point(size = 2, alpha = 0.1) +
  geom_line(
    data = growth.flow.pred,
    aes(group = phen.c, colour = phen.c),
    linewidth = 1.2
  ) +
  scale_colour_gradient2(
    low = 'yellow', high = 'magenta', mid = 'black', midpoint = 0,
    breaks = c(-1:1) * 14,
    labels = format(as.Date((-1:1) * 14 + phen.ctrl.mean), '%b %d')
  ) +
  lims(x = c(0.5, 6), y = c(0.5, 6)) +
  labs(x = '', y = '') +
  theme(
    legend.position = 'none',
    panel.background = element_blank()
  )

leg.trt = get_legend(
  pan.c +
    guides(shape = guide_legend('umbel fate'), colour = guide_legend('treatment')) +
    theme(legend.position = 'top')
)

leg.phen = get_legend(
  pan.e + 
    guides(colour = guide_legend('Flowering date, year t')) +
    theme(legend.position = 'top')
)

# --------------------------------------
# ---------- Combining panels ----------
# --------------------------------------

plot_grid(
  plot_grid(leg.trt, leg.phen, nrow = 2),
  plot_grid(
    pan.a, pan.b, pan.c, nrow = 1,
    labels = c('a', 'b', 'c'), label_x = -0.005
  ), 
  grid::textGrob("Flowering date",  vjust = 0),
  plot_grid(
    pan.d, pan.e, align = 'v',
    labels = c('d', 'e'), label_x = -0.005
  ), 
  grid::textGrob("Size, year t", vjust = 0),
  nrow = 5, rel_heights = c(0.2, 1, 0.025, 1, 0.025)
) %>%
  save_plot(
    filename = '04_analysis/figures/draft_figures/reproduction_growth_phen_trt.png',
    base_height = 8, base_width = 8,
  )
