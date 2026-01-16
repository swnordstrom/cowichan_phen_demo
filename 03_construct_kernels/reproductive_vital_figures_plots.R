# Script for generating figures and estimates for Lomatium manuscript
# Here working only with reproductive vital rates
library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyr)
library(glmmTMB)


# =====================================
#  ------------------------------------
# Data read in and processing/reshaping
#  ------------------------------------
# =====================================

rm(list = ls())

# Read in data with wrapper script
source('03_construct_kernels/prepare_demo_data_repr.R')

# =======================================================
# -------------------------------------------------------
# Fit models (model selection performed in other scripts)
# -------------------------------------------------------
# =======================================================

# === Flowering/umbel count model ===
u_s_s.ty = glmmTMB(
  No.umbels ~ size + (1 | Year) + (1 | Plot / plantid),
  family = 'truncated_poisson',
  ziformula = ~ size + trt + (1 | Year) + (1 | Year:trt) + (1 | Plot / plantid),
  data = demo.flow
)

# Mean effect of irrigation on flowering:
# coefficient 0.55808, irrigation reduces odds of flowering by 1 - exp(-.55808) = ~ 43%
# drought effect on flowering is quite small

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

# === Recruit size model ===
# r_t.y = glmmTMB(size ~ trt + (1 | Year) + (1 | Plot), data = demo.recr)

# =======================================================
# -------------------------------------------------------
# Quantities reported in the manuscript
# -------------------------------------------------------
# =======================================================

# === Estimate of phenology effects on seed production === 
expand.grid(phen.c = 0:-7, Year = 2021:2024) %>%
  mutate(size = 3.9, phen.umbels = 1, trt = 'control') %>%
  mutate(
    lin.zinf = predict(s_st.p_s.u.p, newdata = ., allow.new.levels = TRUE, re.form = ~ 0, type = 'zlink'),
    lin.cond = predict(s_st.p_s.u.p, newdata = ., allow.new.levels = TRUE, re.form = ~ 0, type = 'link')
  ) %>%
  group_by(phen.c, trt) %>%
  summarise(across(c(lin.zinf, lin.cond), mean)) %>%
  mutate(
    # Transform from linear scale to response scale
    seeds.per.umbel = (1 / (1 + exp(lin.zinf))) * exp(lin.cond),
  )

expand.grid(phen.c = c(-3.29, 0, 1.31), Year = 2021:2024) %>%
  mutate(size = 3.9, phen.umbels = 1, trt = 'control') %>%
  mutate(
    lin.zinf = predict(s_st.p_s.u.p, newdata = ., allow.new.levels = TRUE, re.form = ~ 0, type = 'zlink'),
    lin.cond = predict(s_st.p_s.u.p, newdata = ., allow.new.levels = TRUE, re.form = ~ 0, type = 'link')
  ) %>%
  group_by(phen.c, trt) %>%
  summarise(across(c(lin.zinf, lin.cond), mean)) %>%
  mutate(
    # Transform from linear scale to response scale
    seeds.per.umbel = (1 / (1 + exp(lin.zinf))) * exp(lin.cond),
  )



# =======================================================
# -------------------------------------------------------
# Plots for supplement
# -------------------------------------------------------
# =======================================================

# Remove NAs for plot
seed = seed %>% filter(!is.na(no.seeds))

pred.umbel.counts = expand.grid(
  size = (5:60)/10,
  trt = factor(c('control', 'drought', 'irrigated'), levels = c('irrigated', 'drought', 'control'))
  # trt  = factor(c('irrigated', 'drought', 'control'))
) %>%
  mutate(
    pred.zinf = predict(u_s_s.ty, newdata = ., re.form = ~ 0, allow.new.levels = TRUE, type = 'zprob'),
    pred.cond = predict(u_s_s.ty, newdata = ., re.form = ~ 0, allow.new.levels = TRUE, type = 'conditional'),
    p.flower = 1 - pred.zinf,
    total.umbels = p.flower * pred.cond
  )

pred.seed.counts = merge(
  pred.umbel.counts,
  data.frame(Year = 2022:2024)
) %>%
  rename(phen.umbels = total.umbels) %>%
  mutate(phen.c = 0) %>%
  mutate(
    pred.zinf = predict(s_st.p_s.u.p, newdata = ., re.form = ~ 0, allow.new.levels = TRUE, type = 'zlink'),
    pred.cond = predict(s_st.p_s.u.p, newdata = ., re.form = ~ 0, allow.new.levels = TRUE, type = 'link'),
  ) %>%
  group_by(size, trt, phen.umbels) %>%
  summarise(across(c(pred.zinf, pred.cond), mean)) %>%
  mutate(
    p.surv = 1 / (1 + exp(pred.zinf)),
    n.seed = exp(pred.cond),
    seed.per.umbel = p.surv * n.seed,
    total.seed = seed.per.umbel * phen.umbels
  )  %>%
  ungroup()

pred.seed.counts %>%
  select(size, trt, p.surv, n.seed, seed.per.umbel, total.seed) %>%
  pivot_longer(-c(size, trt), names_to = 'varb', values_to = 'value') %>%
  ggplot(aes(x = size, y = value, group = trt, colour = trt)) +
  geom_line() +
  scale_colour_manual(
    values = c('black', 'goldenrod', 'dodgerblue'), 
    breaks = c('control', 'drought', 'irrigated')
  ) +
  facet_wrap(~ varb, scale = 'free_y')

# Quadratic-looking effect for pr(survival) and seeds per umbel is due to the
# umbel count influencing probability of umbel survival.

# Get seeds per plant for figure
seeds.per.plant = seed %>%
  group_by(plantid, Year, trt) %>%
  summarise(
    no.seeds = sum(no.seeds),
    size = mean(size)
  ) %>%
  ungroup()

### Plot elements

p.a = pred.umbel.counts %>%
  ggplot(aes(x = size, y = p.flower, group = trt, colour = trt)) +
  geom_point(
    data = demo.flow %>% mutate(p.flower = as.numeric(No.umbels > 0)),
    aes(shape = as.logical(p.flower)),
    position = position_jitter(height = 0.1), size = 2, alpha = 0.05
  ) +
  geom_line(linewidth = 1.2) +
  scale_colour_manual(
    values = c('black', 'goldenrod', 'dodgerblue'), 
    breaks = c('control', 'drought', 'irrigated')
  ) +
  scale_shape_manual(values = c(4, 19), labels = c('no', 'yes')) +
  labs(x = '', y = 'probability of flowering') +
  theme(
    legend.position = 'none',
    panel.background = element_blank()
  )

p.b = pred.umbel.counts %>%
  ggplot(aes(x = size, y = pred.cond, group = trt)) +
  geom_point(
    data = demo.flow %>% filter(No.umbels > 0) %>% rename(pred.cond = No.umbels),
    aes(colour = trt),
    position = position_jitter(height = 0.2), size = 2, alpha = 0.05
  ) +
  geom_line(linewidth = 1.2, colour = 'gray77') +
  labs(x = '', y = 'umbels / flowering plant') +
  scale_colour_manual(
    values = c('black', 'goldenrod', 'dodgerblue'), 
    breaks = c('control', 'drought', 'irrigated')
  ) +
  # scale_y_log10() +
  theme(
    legend.position = 'none',
    panel.background = element_blank()
  )

p.c = pred.umbel.counts %>%
  ggplot(aes(x = size, y = total.umbels, group = trt, colour = trt)) +
  geom_point(
    data = demo.flow %>% rename(total.umbels = No.umbels),
    aes(shape = as.logical(total.umbels > 0)),
    position = position_jitter(height = 0.2), size = 2, alpha = 0.05
  ) +
  geom_line(linewidth = 1.2) +
  scale_colour_manual(
    values = c('black', 'goldenrod', 'dodgerblue'), 
    breaks = c('control', 'drought', 'irrigated')
  ) +
  scale_shape_manual(values = c(4, 19)) +
  labs(x = '', y = 'umbels / plant') +
  # scale_y_log10() +
  theme(
    legend.position = 'none',
    panel.background = element_blank()
  )

umbel.legend = get_legend(
  p.a +
    guides(
      colour = guide_legend('treatment', order = 1), 
      shape = guide_legend('flowering', order = 2, override.aes = list(alpha = 1))
    ) +
    theme(legend.position = 'top', legend.box = 'vertical', legend.spacing.y = unit(0.1, 'points'))
)

# plot_grid(
#   NULL, umbel.legend, NULL, p.a, p.b, p.c, 
#   rel_heights = c(0.1, 1), nrow = 2
# )

p.d = pred.seed.counts %>%
  ggplot(aes(x = size, y = p.surv, group = trt)) +
  geom_point(
    data = seed,
    aes(
      x = size, y = as.numeric(no.seeds > 0), colour = trt, shape = as.logical(no.seeds > 0)
    ),
    position = position_jitter(height = 0.1), size = 2, alpha = 0.05
  ) +
  geom_line(linewidth = 1.2, colour = 'gray77') +
  scale_shape_manual(values = c(4, 19), labels = c('unsuccessful', 'successful')) +
  scale_colour_manual(
    values = c('black', 'goldenrod', 'dodgerblue'), 
    breaks = c('control', 'drought', 'irrigated')
  ) +
  labs(x = '', y = 'prob. of umbel success') +
  theme(
    legend.position = 'none',
    panel.background = element_blank()
  )

p.e = pred.seed.counts %>%
  ggplot(aes(x = size, y = n.seed, group = trt, colour = trt)) +
  geom_point(
    data = seed %>% filter(no.seeds > 0),
    aes(x = size, y = no.seeds, colour = trt),
    size = 2, alpha = 0.05
  ) +
  geom_line(linewidth = 1.2) +
  scale_colour_manual(
    values = c('black', 'goldenrod', 'dodgerblue'), 
    breaks = c('control', 'drought', 'irrigated')
  ) +
  scale_y_log10() +
  labs(x = '', y = 'seeds / successful umbel') +
  theme(
    legend.position = 'none',
    panel.background = element_blank()
  )

p.f = pred.seed.counts %>%
  ggplot(aes(x = size, y = seed.per.umbel, group = trt, colour = trt)) +
  geom_point(
    data = seed %>% mutate(no.seeds = ifelse(no.seeds < 1, 0.5, no.seeds)),
    aes(x = size, y = no.seeds, shape = as.logical(no.seeds > 0.5), colour = trt),
    size = 2, alpha = 0.05
  ) +
  geom_line(linewidth = 1.2) +
  scale_shape_manual(values = c(4, 19)) +
  scale_colour_manual(
    values = c('black', 'goldenrod', 'dodgerblue'), 
    breaks = c('control', 'drought', 'irrigated')
  ) +
  scale_y_log10() +
  labs(x = '', y = 'seeds / umbel') +
  theme(
    legend.position = 'none',
    panel.background = element_blank()
  )

seed.legend = get_legend(
  p.d +
    guides(shape = guide_legend('umbel fate', override.aes = list(alpha = 1)), colour = 'none') +
    theme(legend.position = 'top', legend.box = 'vertical', legend.spacing.y = unit(0.1, 'points'))
)

# plot_grid(
#   NULL, umbel.legend, NULL, p.a, p.b, p.c, 
#   NULL, seed.legend,  NULL, p.d, p.e, p.f,
#   rel_heights = c(0.2, 1, 0.1, 1), nrow = 4
# )

p.g = pred.seed.counts %>%
  ggplot(aes(x = size, y = total.seed, group = trt, colour = trt)) +
  geom_point(
    data = seeds.per.plant,
    aes(x = size, y = no.seeds, colour = trt),
    size = 2, alpha = 0.1
  ) +
  geom_line(linewidth = 1.2) +
  scale_colour_manual(
    values = c('black', 'goldenrod', 'dodgerblue'), 
    breaks = c('control', 'drought', 'irrigated')
  ) +
  labs(x = '', y = 'seeds / plant') +
  theme(
    legend.position = 'none',
    panel.background = element_blank()
  )

p.h = pred.seed.counts %>%
  filter(size > 3.1, size < 4.68) %>%
  ggplot(aes(x = size, y = total.seed, group = trt, colour = trt)) +
  geom_point(
    data = seeds.per.plant %>% 
      filter(size > 3.1, size < 4.68), # %>%
      # mutate(no.seeds = ifelse(no.seeds < 1, 0.5, no.seeds)),
    aes(x = size, y = no.seeds, colour = trt),
    size = 2, alpha = 0.1
  ) +
  geom_line(linewidth = 1.2) +
  scale_colour_manual(
    values = c('black', 'goldenrod', 'dodgerblue'), 
    breaks = c('control', 'drought', 'irrigated')
  ) +
  # scale_y_log10() +
  labs(x = '', y = '') +
  theme(
    legend.position = 'none',
    panel.background = element_blank()
  )


plot_grid(
  plot_grid(
    NULL, umbel.legend, NULL, p.a, p.b, p.c, 
    NULL, seed.legend,  NULL, p.d, p.e, p.f,
    labels = c('', '', '', 'a)', 'b)', 'c)', '', '', '', 'd)', 'e)', 'f)'),
    rel_heights = c(0.2, 1, 0.1, 1), nrow = 4
  ),
  ggdraw(plot_grid(p.g, p.h, nrow = 1, align = 'v', labels = c('g)', 'h)'))) + 
    draw_label('size', x = 0.5, y = 0, vjust = -1),
  rel_heights = c(3, 1), nrow = 2
) %>%
  save_plot(
    '04_analysis/figures/draft_figures/fig_supp_reproduction.png', 
    plot = ., base_height = 8, base_asp = 1.1
  )
