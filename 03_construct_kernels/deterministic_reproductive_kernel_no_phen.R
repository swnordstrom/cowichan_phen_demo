# Script for re-running reproductive vital rate models and performing model
# selection and exporting a deterministic, no-phenology reproductive kernel for
# IPM analysis.
# Reads in processed demo/seed (including phenology) data, 2016-2024
# (sn init july 2024)

# --- Setup ---------------------------------------------------------

library(ggplot2)
library(dplyr)
library(tidyr)
library(glmmTMB)

rm(list = ls())

# --- Read in and prepare data  -------------------------------------

# All data (demo, phenology, seed set)
all.data = merge(
  x = read.csv('01_data_cleaning/out/demo_phen_seed_2016-2024_final.csv'),
  y = read.csv('00_raw_data/plot_treatments.csv'),
  by.x = 'Plot', by.y = 'plot'
)

# Dataset for flowering (additional processing below)
demo.flow = all.data %>%
  # Get only plants that are:
  # - known/estimated to be alive
  # - Has size measurements (from demography) and has non-zero leaf counts and lengths (for size measurements)
  # - Is not missing umbel counts (because we're subsetting this for umbel count analysis)
  filter(
    in.demo, Year > 2016, surv,
    !is.na(No.leaves) & !is.na(Leaf.length) & No.leaves > 0 & Leaf.length > 0,
    !is.na(No.umbels)
  ) %>%
  # Get a size column
  mutate(size = log(No.leaves * Leaf.length))

head(demo.flow)
nrow(demo.flow)
table(demo.flow$Year)

# Dataset for seed analysis
demo.seed = demo.flow %>%
  # Give only the records for which we have a record of seed set
  filter(Year > 2020, phen.umbels > 0) %>%
  filter(in.seed)

# Need to add in rows for umbels that died before the seed counting
# This is necessary for estimating the probability of an umbel producing zero seeds
demo.seed = rbind(
  demo.seed,
  demo.seed %>%
    group_by(Year, plantid) %>%
    mutate(miss.umbel = ifelse(phen.umbels < n(), 0, phen.umbels - n())) %>%
    ungroup() %>%
    uncount(miss.umbel) %>%
    mutate(no.seeds = 0)
) %>%
  mutate(Year = factor(Year))

# Dataset for estimating recruit size distribution
demo.recr = all.data %>%
  arrange(Year) %>%
  # Want demo records (to get size distribution)
  filter(in.demo) %>%
  # Give me only the first sighting for each plant
  distinct(plantid, .keep_all = TRUE) %>%
  # Give me only records after 2017 (these are more likely to be seen for the
  # first time)
  filter(Year > 2017) %>%
  # Give me only records with leaf counts and lengths (for estimating size)
  filter(!is.na(No.leaves), !is.na(Leaf.length)) %>%
  # Give me plants that did not flower (plants unlikely to flower in first year)
  filter(!No.umbels & !is.na(No.umbels)) %>%
  # Okay... going to assume that only plants with one leaf are new recruits
  filter(No.leaves < 2) %>%
  mutate(size = log(Leaf.length))

# Finalize the flowering demo dataset 
demo.flow = demo.flow %>%
  # (currently there is one row per umbel for plants in the seed set dataset -
  # give me just one record per plant)
  distinct(Year, plantid, .keep_all = TRUE) %>%
  # Remove plot 2 in 2024 (missing a lot of vegetative plant sizes)
  filter(!(Plot %in% 2 & Year %in% 2024))


# # ----- MODEL SELECTION ------------------------------------------
# 
# # --- Umbel count model ------------------------------------------
# 
# # Strategy here: use a hurdle model
# # - Response is number of umbels per plant
# # - Hurdle portion (through zero-inflation) is the probability of (not) flowering
# # - Umbel count is modeled with a zero-inflated Poisson distribution
# 
# # Test first for size-effects on both parts of the model
# # Then test for treatment effects
# # Then test for potential interactions
# 
# u_0 = glmmTMB(
#   No.umbels ~ (1 | Year) + (1 | Plot / plantid),
#   family = 'truncated_poisson',
#   ziformula = ~ (1 | Year) + (1 | Plot / plantid),
#   data = demo.flow
# )
# 
# u_s_0 = glmmTMB(
#   No.umbels ~ size + (1 | Year) + (1 | Plot / plantid),
#   family = 'truncated_poisson',
#   ziformula = ~ (1 | Year) + (1 | Plot / plantid),
#   data = demo.flow
# )
# 
# u_0_s = glmmTMB(
#   No.umbels ~ (1 | Year) + (1 | Plot / plantid),
#   family = 'truncated_poisson',
#   ziformula = ~ size + (1 | Year) + (1 | Plot / plantid),
#   data = demo.flow
# )
# 
# u_s_s = glmmTMB(
#   No.umbels ~ size + (1 | Year) + (1 | Plot / plantid),
#   family = 'truncated_poisson',
#   ziformula = ~ size + (1 | Year) + (1 | Plot / plantid),
#   data = demo.flow
# )
# 
# AIC(u_0, u_0_s, u_s_0, u_s_s) %>% mutate(daic = round(AIC - min(AIC), 2))
# # Keep size in the models
# 
# # Treatment effects
# 
# u_s.t_s = glmmTMB(
#   No.umbels ~ trt + size + (1 | Year) + (1 | Plot / plantid),
#   family = 'truncated_poisson',
#   ziformula = ~ size + (1 | Year) + (1 | Plot / plantid),
#   data = demo.flow
# )
# 
# u_s_s.t = glmmTMB(
#   No.umbels ~ size + (1 | Year) + (1 | Plot / plantid),
#   family = 'truncated_poisson',
#   ziformula = ~ trt * size + (1 | Year) + (1 | Plot / plantid),
#   data = demo.flow
# )
# 
# u_st_s = glmmTMB(
#   No.umbels ~ size * trt + (1 | Year) + (1 | Plot / plantid),
#   family = 'truncated_poisson',
#   ziformula = ~ size + (1 | Year) + (1 | Plot / plantid),
#   data = demo.flow
# )
# 
# u_s_st = glmmTMB(
#   No.umbels ~ size + (1 | Year) + (1 | Plot / plantid),
#   family = 'truncated_poisson',
#   ziformula = ~ trt * size + (1 | Year) + (1 | Plot / plantid),
#   data = demo.flow
# )
# 
# u_s.ty_s = glmmTMB(
#   No.umbels ~ size + trt + (1 | Year) + (1 | Year:trt) + (1 | Plot / plantid),
#   family = 'truncated_poisson',
#   ziformula = ~ size + (1 | Year) + (1 | Plot / plantid),
#   data = demo.flow
# )
# 
# u_s_s.ty = glmmTMB(
#   No.umbels ~ size + (1 | Year) + (1 | Plot / plantid),
#   family = 'truncated_poisson',
#   ziformula = ~ size + trt + (1 | Year) + (1 | Year:trt) + (1 | Plot / plantid),
#   data = demo.flow
# )
# 
# u_s.ty_s.ty = glmmTMB(
#   No.umbels ~ size + trt + (1 | Year) + (1 | Year:trt) + (1 | Plot / plantid),
#   family = 'truncated_poisson',
#   ziformula = ~ size + trt + (1 | Year) + (1 | Year:trt) + (1 | Plot / plantid),
#   data = demo.flow
# )
# 
# AIC(u_s.t_s, u_s_s.t, u_st_s, u_s_st, u_s_s, u_s.ty_s, u_s_s.ty, u_s.ty_s.ty) %>% 
#   mutate(daic = round(AIC - min(AIC), 2))
# # Best performing model has treatment-year effect for probability of flowering
# # But how large is the effect size?
# 
# summary(u_s_s.ty)
# 
# u_st_s.ty = glmmTMB(
#   No.umbels ~ size * trt + (1 | Year) + (1 | Plot / plantid),
#   family = 'truncated_poisson',
#   ziformula = ~ size + trt + (1 | Year) + (1 | Year:trt) + (1 | Plot / plantid),
#   data = demo.flow
# )
# 
# u_s_st.ty = glmmTMB(
#   No.umbels ~ size + (1 | Year) + (1 | Plot / plantid),
#   family = 'truncated_poisson',
#   ziformula = ~ size * trt + (1 | Year) + (1 | Year:trt) + (1 | Plot / plantid),
#   data = demo.flow
# )
# 
# AIC(u_st_s.ty, u_s_st.ty, u_s_s.ty) %>% mutate(daic = round(AIC - min(AIC), 2))
# # No size-treatment effects
# 
# # Okay... compare model with year-varying effects with non-varying effects
# 
# expand.grid(size = (5:60)/10, trt = c('control', 'drought', 'irrigated')) %>%
#   mutate(
#     u.nrf = predict(u_s_s, newdata = ., allow.new.levels = TRUE, re.form = ~ 0, type = 'response'),
#     u.yrf = predict(u_s_s.ty, newdata = ., allow.new.levels = TRUE, re.form = ~ 0, type = 'response')
#   ) %>%
#   pivot_longer(c(u.nrf, u.yrf), names_to = 'model', values_to = 'pred') %>%
#   ggplot(aes(x = size, y = pred, linetype = model, colour = trt)) +
#   geom_line() +
#   scale_colour_manual(values = c('black', 'red', 'blue')) +
#   scale_y_log10()
# # interesting - the year-varying effects model looks like it always favors the controls
# # i.e. control plants always produce more umbels
# # presumably this is probability of flowering-related - smaller plants less likely to flower
# 
# expand.grid(size = (5:60)/10, trt = c('control', 'drought', 'irrigated')) %>%
#   mutate(
#     pred.umbl = predict(
#       u_s_s.ty, newdata = ., allow.new.levels = TRUE, re.form = ~ 0, type = 'response'
#     )
#   ) %>%
#   ggplot(aes(x = size, y = pred.umbl, colour = trt)) +
#   geom_point(
#     data = demo.flow,
#     aes(x = size, y = No.umbels),
#     size = 3, position = position_jitter(height = 0.25), alpha = 0.2
#   ) +
#   geom_line() + # scale_y_log10()
#   scale_colour_manual(values = c('black', 'red', 'blue'))
# 
# 
# # --- Seed set model  
# 
# # Here, we'll also use a zero-inflated model
# # - Probability of producing zero seed estimated with zero inflation term
# # - Distribution of seeds after the zero inflation is modeled by a Negative Binomial
# #   - doing this because it accurately captures overdispersion and because the 
# #     data-generating process can also account for umbels that were simply not 
# #     sufficiently pollinated
# 
# s_0 = glmmTMB(
#   no.seeds ~ (1 | Plot / plantid),
#   ziformula = ~ (1 | Plot / plantid),
#   family = 'nbinom2',
#   data = demo.seed
# )
# 
# s_0_y = glmmTMB(
#   no.seeds ~ (1 | Plot / plantid),
#   ziformula = ~ Year + (1 | Plot / plantid),
#   family = 'nbinom2',
#   data = demo.seed
# )
# 
# s_y_0 = glmmTMB(
#   no.seeds ~ Year + (1 | Plot / plantid),
#   ziformula = ~ (1 | Plot / plantid),
#   family = 'nbinom2',
#   data = demo.seed
# )
# 
# s_y_y = glmmTMB(
#   no.seeds ~ Year + (1 | Plot / plantid),
#   ziformula = ~ Year + (1 | Plot / plantid),
#   family = 'nbinom2',
#   data = demo.seed
# )
# 
# AIC(s_0, s_y_0, s_0_y, s_y_y) %>% arrange(AIC)
# # Unsurprisingly, year effects supported in both models
# 
# s_s.y_y = glmmTMB(
#   no.seeds ~ size + Year + (1 | Plot / plantid),
#   ziformula = ~ Year + (1 | Plot / plantid),
#   family = 'nbinom2',
#   data = demo.seed
# )
# 
# s_y_s.y = glmmTMB(
#   no.seeds ~ Year + (1 | Plot / plantid),
#   ziformula = ~ size + Year + (1 | Plot / plantid),
#   family = 'nbinom2',
#   data = demo.seed
# )
# 
# s_s.y_s.y = glmmTMB(
#   no.seeds ~ size + Year + (1 | Plot / plantid),
#   ziformula = ~ size + Year + (1 | Plot / plantid),
#   family = 'nbinom2',
#   data = demo.seed
# )
# 
# AIC(s_y_y, s_s.y_y, s_y_s.y, s_s.y_s.y) %>% mutate(daic = AIC - min(AIC))
# 
# # Include size in both
# 
# s_s.t.y_s.y = glmmTMB(
#   no.seeds ~ trt + size + Year + (1 | Plot / plantid),
#   ziformula = ~ size + Year + (1 | Plot / plantid),
#   family = 'nbinom2',
#   data = demo.seed
# )
# 
# s_s.y_s.t.y = glmmTMB(
#   no.seeds ~ size + Year + (1 | Plot / plantid),
#   ziformula = ~ size + Year + trt + (1 | Plot / plantid),
#   family = 'nbinom2',
#   data = demo.seed
# )
# 
# s_s.t.y_s.t.y = glmmTMB(
#   no.seeds ~ trt + size + Year + (1 | Plot / plantid),
#   ziformula = ~ trt + size + Year + (1 | Plot / plantid),
#   family = 'nbinom2',
#   data = demo.seed
# )
# 
# AIC(s_s.y_s.y, s_s.t.y_s.y, s_s.y_s.t.y, s_s.t.y_s.t.y) %>% mutate(daic = round(AIC - min(AIC), 2))
# # Eh actually not seeing evidence of a treatment effect here...
# # Weird! it did show up for the poisson model...
# 
# s_st.y_s.y = glmmTMB(
#   no.seeds ~ trt * size + Year + (1 | Plot / plantid),
#   ziformula = ~ size + Year + (1 | Plot / plantid),
#   family = 'nbinom2',
#   data = demo.seed
# )
# 
# s_s.ty_s.y = glmmTMB(
#   no.seeds ~ trt * Year + size + (1 | Plot / plantid),
#   ziformula = ~ size + Year + (1 | Plot / plantid),
#   family = 'nbinom2',
#   data = demo.seed
# )
# 
# s_sy.ty_s.y = glmmTMB(
#   no.seeds ~ trt * Year + trt * size + (1 | Plot / plantid),
#   ziformula = ~ size + Year + (1 | Plot / plantid),
#   family = 'nbinom2',
#   data = demo.seed
# )
# 
# s_s.y_st.y = glmmTMB(
#   no.seeds ~ size + Year + (1 | Plot / plantid),
#   ziformula = ~ trt * size + Year + (1 | Plot / plantid),
#   family = 'nbinom2',
#   data = demo.seed
# )
# 
# s_s.y_s.ty = glmmTMB(
#   no.seeds ~ size + Year + (1 | Plot / plantid),
#   ziformula = ~ trt * Year + size + (1 | Plot / plantid),
#   family = 'nbinom2',
#   data = demo.seed
# )
# 
# s_s.y_st.ty = glmmTMB(
#   no.seeds ~ size + Year + (1 | Plot / plantid),
#   ziformula = ~ trt * Year + trt * size + (1 | Plot / plantid),
#   family = 'nbinom2',
#   data = demo.seed
# )
# 
# AIC(s_st.y_s.y, s_s.ty_s.y, s_sy.ty_s.y, s_s.y_st.y, s_s.y_s.ty, s_s.y_st.ty, s_s.y_s.y) %>%
#   arrange(AIC)
# # Won't even bother with any of the models with the treatment interactions in
# # the zero-inflation formula (their models look bad in this table)
# 
# 
# summary(s_st.y_s.y)
# 
# # Number of umbels
# 
# s_st.y.u_s.y = glmmTMB(
#   no.seeds ~ trt * size + phen.umbels + Year + (1 | Plot / plantid),
#   ziformula = ~ size + Year + (1 | Plot / plantid),
#   family = 'nbinom2',
#   data = demo.seed
# )
# 
# s_st.y_s.y.u = glmmTMB(
#   no.seeds ~ trt * size + Year + (1 | Plot / plantid),
#   ziformula = ~ size + Year + phen.umbels + (1 | Plot / plantid),
#   family = 'nbinom2',
#   data = demo.seed
# )
# 
# s_st.y.u_s.y.u = glmmTMB(
#   no.seeds ~ trt * size + phen.umbels + Year + (1 | Plot / plantid),
#   ziformula = ~ size + Year + phen.umbels + (1 | Plot / plantid),
#   family = 'nbinom2',
#   data = demo.seed
# )
# 
# AIC(s_st.y_s.y, s_st.y.u_s.y, s_st.y_s.y.u, s_st.y.u_s.y.u) %>% arrange(AIC)
# 
# # Look for size-year interactions
# 
# s_sty_s.y.u = glmmTMB(
#   no.seeds ~ trt * size * Year + (1 | Plot / plantid),
#   ziformula = ~ size + Year + phen.umbels + (1 | Plot / plantid),
#   family = 'nbinom2',
#   data = demo.seed
# )
# 
# s_st.y_sy.u = glmmTMB(
#   no.seeds ~ trt * size + Year + (1 | Plot / plantid),
#   ziformula = ~ size * Year + phen.umbels + (1 | Plot / plantid),
#   family = 'nbinom2',
#   data = demo.seed
# )
# 
# s_sty_sy.u = glmmTMB(
#   no.seeds ~ trt * size * Year + (1 | Plot / plantid),
#   ziformula = ~ size * Year + phen.umbels + (1 | Plot / plantid),
#   family = 'nbinom2',
#   data = demo.seed
# )
# 
# AIC(s_st.y_s.y.u, s_sty_s.y.u, s_st.y_sy.u, s_sty_sy.u) %>% arrange(AIC)
# 
# 
# # Best mod is s_st.y_sy.u
# 
# summary(s_st.y_sy.u)
# 
# pred.seeds = expand.grid(
#   Year = factor(2021:2024),
#   size = (5:60)/10,
#   trt = c('control', 'drought', 'irrigated'),
#   phen.umbels = c(1, 5)
# ) %>%
#   mutate(
#     pred.seed = predict(
#       s_st.y_sy.u, newdata = ., 
#       re.form = ~ 0, type = 'response',
#       allow.new.levels = TRUE
#     )
#   )
# 
# pred.seeds %>% 
#   mutate(phen.umbels = factor(phen.umbels)) %>%
#   ggplot(aes(x = size, y = pred.seed, group = interaction(trt, phen.umbels))) +
#   geom_line(aes(colour = trt, linetype = phen.umbels)) +
#   scale_colour_manual(values = c('black', 'red', 'blue')) +
#   scale_y_log10() +
#   facet_wrap(~ Year)
# 
# # --- Recruit size model  ---------------------------------------
# 
# # Null model
# r_0 = glmmTMB(size ~ (1 | Plot), data = demo.recr)
# 
# # A model with a year random effect
# r_y = glmmTMB(size ~ (1 | Year) + (1 | Plot), data = demo.recr)
# 
# # A model with a treatment effect
# r_t = glmmTMB(size ~ trt + (1 | Plot), data = demo.recr)
# 
# # A model with both year efect and treatment effect
# r_t.y = glmmTMB(size ~ trt + (1 | Year) + (1 | Plot), data = demo.recr)
# 
# # A model with both year efect and treatment effect
# r_ty = glmmTMB(size ~ trt + (1 | Year) + (1 | Year:trt) + (1 | Plot), data = demo.recr)
# 
# AIC(r_0, r_y, r_t, r_t.y, r_ty) %>% mutate(daic = round(AIC - min(AIC), 2))
# # Interesting. Definite evidence for a year effect. Weak evidence of a treatment effect.
# 
# summary(r_t.y)
# # The effect size is kinda marginal...
# # Wait... drought plants are slightly larger than control plants
# # Huh.
# 
# ggplot(demo.recr, aes(x = size, group = trt, fill = trt)) +
#   geom_histogram(position = 'identity', alpha = 0.5, binwidth = 0.2) + 
#   scale_fill_manual(values = c('black', 'red', 'blue'))
# # hmm... 

# --- Estimate predictions from model forming backbone of kernel -------

u_s_s.ty = glmmTMB(
  No.umbels ~ size + (1 | Year) + (1 | Plot / plantid),
  family = 'truncated_poisson',
  ziformula = ~ size + trt + (1 | Year) + (1 | Year:trt) + (1 | Plot / plantid),
  data = demo.flow
)

s_st.y_sy.u = glmmTMB(
  no.seeds ~ trt * size + Year + (1 | Plot / plantid),
  ziformula = ~ size * Year + phen.umbels + (1 | Plot / plantid),
  family = 'nbinom2',
  data = demo.seed
)

r_t.y = glmmTMB(size ~ trt + (1 | Year) + (1 | Plot), data = demo.recr)

# Standard deviation (residual error) term from recruitment model
sigma.recr = summary(r_t.y)$sigma

backbone = expand.grid(
  size = (5:60)/10,
  size.cur = (5:60)/10,
  trt = c('control', 'drought', 'irrigated'),
  Year = factor(2021:2024)
) %>%
  mutate(
    # Umbel count
    phen.umbels = predict(
      u_s_s.ty, newdata = ., 
      allow.new.levels = TRUE, re.form = ~ 0, type = 'response'
    )
  ) %>%
  mutate(
    # Use this umbel count to get seeds/umbel
    seeds.per.umbel = predict(
      s_st.y_sy.u, newdata = .,
      allow.new.levels = TRUE, re.form = ~ 0
    )
  ) %>%
  group_by(size, size.cur, trt) %>%
  summarise(across(c(phen.umbels, seeds.per.umbel), mean)) %>%
  ungroup() %>%
  mutate(
    # Transform from linear scale to response scale
    seeds.per.umbel = exp(seeds.per.umbel),
    # Total umbels per plant
    seeds.total = seeds.per.umbel * phen.umbels
  ) %>%
  mutate(
    # Mean recruit size
    recr.mean = predict(
      r_t.y, allow.new.levels = TRUE, re.form = ~ 0, newdata = .
    ),
    # Get the number of seeds produced for each size grouping
    p.size.cur = 0.1 * seeds.total * dnorm(x = size.cur, mean = recr.mean, sd = sigma.recr)
  ) %>%
  # Rename column
  rename(size.pre = size)

head(backbone)
tail(backbone)

backbone %>%
  ggplot(aes(x = size.pre, y = seeds.total, group = trt)) +
  geom_line(aes(colour = trt)) +
  scale_colour_manual(values = c('black', 'red', 'blue'))
# Plausible numbers

backbone %>%
  ggplot(aes(x = size.pre, y = size.cur, fill = p.size.cur)) + # fill = log(p.size.cur, base = 10))) +
  geom_raster() +
  scale_fill_viridis_c() +
  scale_y_reverse() +
  facet_wrap(~ trt)
# Lol, cool


# Export kernel in CSV form
write.csv(
  backbone,
  na = '', row.names = FALSE,
  '03_construct_kernels/out/deterministic_reprod_kernel_no_phen.csv'
)

