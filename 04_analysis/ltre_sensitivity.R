# CHECK DELTA VALUE

### ---------------------------------

library(ggplot2)
library(ggh4x)
library(dplyr)
library(tidyr)
library(cowplot)

rm(list = ls())

cat('Building kernels for LTRE and Figure 3... ')

# ------------------------------------------------------                  
# ------ Read in all data ------------------------------
# ------------------------------------------------------                  

# --- Read in observed data kernels (subkernels) 
# I estimated reproductive and growth/survival subkernels separately
# Naming convention in files/object names:
# - gs is growth + survival
# - fr is flowering + reproduction

# Growth + survival subkernel
gs.obsv = read.csv('03_construct_kernels/out/deterministic_growsurv_kernel_phen_ltre.csv')
# Flowering + reproduction subkernel
fr.obsv = read.csv('03_construct_kernels/out/deterministic_reprod_kernel_phen_ltre.csv')


# --- Read in perturbed subkernels (on observed data)
# (these get used for sensitivity analysis)

# Growth perturbed kernels
# (why is the filename 'no_phen'... growth coefficients don't have any phen at all...)
gs.obsv.pert = read.csv('03_construct_kernels/out/deterministic_grow_coef_perturbation_phen.csv') %>%
  rename(size.nex = size.cur, param = perturb.param, orig.parval = orig.par.val) %>%
  mutate(param = gsub('\\_', '.', param))

# Reproductive perturbed kernels
fr.obsv.pert = read.csv('03_construct_kernels/out/deterministic_repr_coef_perturbation_phen.csv')


# --- Read in parameters used in bootstrapping
# (these give the differences in beta in the LTRE)

gs.pert.pars = read.csv('03_construct_kernels/out/deterministic_growsurv_bootstrapped_perturbed_params.csv') # %>%
# pivot_longer(-boot, names_to = 'rate_trt', values_to = 'parval') %>%
# separate(rate_trt, into = c('rate', 'trt'), sep = '_')

fr.pert.pars = read.csv('03_construct_kernels/out/deterministic_reprod_bootstrapped_perturbed_params.csv') 
# Do this in two steps because this file also contains the phen dates used in the bootstarp

# Read in LTRE treatment-phenology info
# NOTE: in some cases we will want to crop out the treatment means 
trt.phen.ltre.key = merge(
  x = read.csv('03_construct_kernels/ltre_treatment_key.csv'),
  y = read.csv('03_construct_kernels/out/phen_treatment_means.csv'),
  by.x = 'trt.phen', by.y = 'trt'
) %>%
  arrange(trt.phen.idx) %>%
  select(trt.phen.idx, everything())

# ------------------------------------------------------                  
# ------ Build kernels ---------------------------------
# ------------------------------------------------------

# Idea here:
# Each one of these imported CSVs has columns for size before and after the
# transition, as well as for treatment.
# We can create kernels by merging the data frames by these sizes (plus
# treatment, phenology, etc.) to get the kernel entries in DF form.
# Then we can use wrapper scripts to convert these into matrices and estimate
# lambdas.

# Germination probabilities
p.germ = as.vector((c(.5, 1) %o% 10^(-(4:2))))
# p.germ = 0.0058007812

# --- Observed kernel
obsv.kernel.df = merge(
  gs.obsv, fr.obsv,
  by.x = c('size.prev', 'size.cur', 'trt.phen.idx'), 
  by.y = c('size.prev', 'size.nex', 'trt.phen.idx'),
) %>%
  merge(y = data.frame(p.germ = p.germ)) %>%
  mutate(
    p.size.cur = pred.surv * (pv.grow.size * (1 - prob.flower) + pf.grow.size * prob.flower) + (p.size.cur * p.germ)
  ) %>%
  select(-c(pred.surv, pv.grow.size, pf.grow.size, prob.flower))


# --- Perturbed kernels (just point estimates, not bootstrapped)
# Idea here is to merge the perturbed kernels (where only *one* parameter is
# perturbed at a time) with the un-perturbed kernels to make sure we're getting
# estimates where only one parameter is being perturbed at a time.
# We have the perturbed kernels in two separate data frames, so we'll do two
# different merges and then rbind them together
obsv.pert.kernel.df = rbind(
  merge(
    # need to remove some columns
    gs.obsv %>% rename(size.nex = size.cur), 
    fr.obsv.pert,
    by = c('size.prev', 'size.nex', 'trt.phen.idx')
  ) %>%
    select(
      size.prev, size.nex, trt.phen.idx, pv.grow.size, pf.grow.size, 
      pred.surv, prob.flower, p.size.cur, param, orig.parval
    ),
  merge(
    gs.obsv.pert, 
    fr.obsv,
    by = c('size.prev', 'size.nex', 'trt.phen.idx'),
  ) %>%
    select(
      size.prev, size.nex, trt.phen.idx, pv.grow.size, pf.grow.size, 
      pred.surv, prob.flower, p.size.cur, param, orig.parval
    )
) %>%
  merge(y = data.frame(p.germ =  p.germ)) %>%
  mutate(
    p.size.cur = pred.surv * (pv.grow.size * (1 - prob.flower) + pf.grow.size * prob.flower) + (p.size.cur * p.germ)
  ) %>%
  select(-c(pred.surv, pv.grow.size, pf.grow.size, prob.flower))

# Merge in treatment info and dates

# Observed kernels
obsv.kernel.df = merge(obsv.kernel.df, trt.phen.ltre.key) %>%
  select(-trt.phen.idx) %>%
  rename(trt = trt.rate)

# Perturbed kernels
obsv.pert.kernel.df = merge(obsv.pert.kernel.df, trt.phen.ltre.key) %>%
  select(-trt.phen.idx) %>%
  rename(trt = trt.rate)

# ------------------------------------------------------                  
# ------ Build midpoint kernels ------------------------
# ------------------------------------------------------

# These are used for the sensitivity analysis
# The sensitivity of lambda to a given vital rate will be estimated by getting
# the midpoint kernel of the two treatments for each level (combo of phen,
# perturbed vital rate, and bootstrap sample).

# Here the `contrast` column has levels `d.c` (comparing or averaging drought
# and control) and `i.c` (comparing or averaging irrigation and control)

# Midpoints of observed data
midp.obsv.kernel.df = obsv.kernel.df %>%
  pivot_wider(names_from = trt, values_from = p.size.cur) %>%
  mutate(d.c = (drought + control) / 2, i.c = (irrigated + control) / 2) %>%
  select(-c(drought, control, irrigated)) %>%
  # 'mentry' is just mean entry
  pivot_longer(c(d.c, i.c), names_to = 'contrast', values_to = 'mentry') %>%
  # remove empty rows - this is for differences that we don't estimate
  # differences for
  filter(!is.na(mentry))

# Midpoints from perturbed kernels
midp.pert.kernel.df = obsv.pert.kernel.df %>%
  select(-orig.parval) %>%
  # remove phen perturbations
  filter(!grepl('phen', param)) %>%
  pivot_wider(names_from = trt, values_from = p.size.cur) %>%
  mutate(d.c = (drought + control) / 2, i.c = (irrigated + control) / 2) %>%
  select(-c(drought, control, irrigated)) %>%
  pivot_longer(c(d.c, i.c), names_to = 'contrast', values_to = 'mentry') %>%
  # filter out NAs
  filter(!is.na(mentry))

# Get midpoints for PHENOLOGY
# Here - to estimate the sensitivity of lambda to phenology, getting a midpoint
# kernel between phen levels

midp.phen.kernel.df = obsv.kernel.df %>%
  select(-mean.phen) %>%
  # NOTE: pivoting out by trt.phen instead of trt (because we're averaging
  # across different phens)
  pivot_wider(names_from = trt.phen, values_from = p.size.cur) %>%
  mutate(d.c = (drought + control) / 2, i.c = (irrigated + control) / 2) %>%
  select(-c(control, drought, irrigated)) %>%
  pivot_longer(c(d.c, i.c), names_to = 'contrast.phen', values_to = 'mentry', values_drop_na = TRUE)

midp.phen.pert.kernel.df = obsv.pert.kernel.df %>%
  select(-c(orig.parval,  mean.phen)) %>%
  # subsetting out ONLY the phenology-related vital rates
  filter(grepl('phen', param)) %>%
  # NOTE: pivoting out by trt.phen instead of trt (because we're averaging
  # across different phens)
  pivot_wider(names_from = trt.phen, values_from = p.size.cur) %>%
  mutate(d.c = (drought + control) / 2, i.c = (irrigated + control) / 2) %>%
  select(-c(control, drought, irrigated)) %>%
  pivot_longer(c(d.c, i.c), names_to = 'contrast.phen', values_to = 'mentry', values_drop_na = TRUE)

# ------------------------------------------------------                  
# ------ Estimate lambdas ------------------------------
# ------------------------------------------------------

# Here - convert the data frames into matrices and estimate lambda for those
# matrices
# Easy, quick way to do this is using split() and lapply()
# Procedure here is:
# - split the data frame into a bunch of lists (split())
# - use lapply to convert the data frame in each list into a matrix
# - use sapply to estimate lambda
# - do slight data frame manipulation for the rest

# Build data frame with observed lambdas
obsv.lambda = split(
  obsv.kernel.df, obsv.kernel.df[,c("trt", "trt.phen", "mean.phen", "p.germ")], 
  sep = '_', drop = TRUE
) %>%
  # split() splits the kernel df into a list where each entry is a data frame
  # for one phen-treatment kernel
  lapply(
    function(df) {
      df %>%
        arrange(size.prev, size.cur) %>%
        pivot_wider(names_from = size.prev, values_from = p.size.cur) %>%
        select(-c(trt, trt.phen, mean.phen, size.cur, p.germ)) %>%
        as.matrix()
    }
  ) %>%
  # estimate the dominant eigenvalue 
  sapply(function(m) Re(eigen(m)$values[1])) %>%
  data.frame(lambda = .) %>%
  # get treatment and phen from rownames
  mutate(ttp = row.names(.)) %>%
  separate(ttp, into = c('trt', 'trt.phen', 'phen', 'p.germ'), sep = '_') %>%
  arrange(trt) %>%
  # convert phen to date (for plotting)
  mutate(
    mean.phen = as.numeric(phen),
    phen.date = as.Date(mean.phen, format = '%b-%d')
  )

# Get lambda for the midpoint of observed matrices
midp.obsv.lambda = split(
  midp.obsv.kernel.df, midp.obsv.kernel.df[,c("contrast", 'trt.phen', "mean.phen", "p.germ")],
  sep = '_', drop = TRUE
) %>%
  lapply(
    function(df) {
      df %>%
        arrange(size.prev, size.cur) %>%
        pivot_wider(names_from = size.prev, values_from = mentry) %>%
        select(-c(contrast, trt.phen, mean.phen, size.cur, p.germ)) %>%
        as.matrix()
    }
  ) %>%
  sapply(function(x) Re(eigen(x)$values[1])) %>%
  data.frame(lambda = .) %>%
  mutate(c_tp = row.names(.)) %>%
  separate(c_tp, into = c('contrast', 'trt.phen', 'phen', 'p.germ'), sep = '_')  

# Get lambdas for the midpoint perturbed matrices
midp.pert.lambda = split(
  midp.pert.kernel.df,
  midp.pert.kernel.df[,c("contrast", 'trt.phen', 'mean.phen', 'param', 'p.germ')],
  sep = '_', drop = TRUE
) %>%
  lapply(
    function(df) {
      df %>%
        arrange(size.prev, size.nex) %>%
        pivot_wider(names_from = size.prev, values_from = mentry) %>%
        select(-c(contrast, trt.phen, param, mean.phen, size.nex, p.germ)) %>%
        as.matrix()
    }
  ) %>%
  sapply(function(x) Re(eigen(x)$values[1])) %>%
  data.frame(lambda = .) %>%
  mutate(c_pg = row.names(.)) %>%
  separate(c_pg, into = c('contrast', 'trt.phen', 'phen', 'param', 'p.germ'), sep = '_')  

# Phenology midpoint lambdas
midp.phen.lambda = split(
  midp.phen.kernel.df,
  midp.phen.kernel.df[,c("trt", "contrast.phen", "p.germ")],
  sep = '_', drop = TRUE
) %>%
  lapply(
    function(df) {
      df %>%
        arrange(size.prev, size.cur) %>%
        pivot_wider(names_from = size.prev, values_from = mentry) %>%
        select(-c(trt, contrast.phen, size.cur, p.germ)) %>%
        as.matrix()
    }
  ) %>%
  sapply(function(x) Re(eigen(x)$values[1])) %>%
  data.frame(lambda = .) %>%
  mutate(tc = row.names(.)) %>%
  separate(tc, into = c('trt', 'contrast.phen', 'p.germ'), sep = '_')

# Phenology midpoint lambdas after perturbation
midp.phen.pert.lambda = split(
  midp.phen.pert.kernel.df,
  midp.phen.pert.kernel.df[,c("trt", "param", "contrast.phen", "p.germ")],
  sep = '_', drop = TRUE
) %>%
  lapply(
    function(df) {
      df %>%
        arrange(size.prev, size.nex) %>%
        pivot_wider(names_from = size.prev, values_from = mentry) %>%
        select(-c(trt, contrast.phen, param, size.nex, p.germ)) %>%
        as.matrix()
    }
  ) %>%
  sapply(function(x) Re(eigen(x)$values[1])) %>%
  data.frame(lambda = .) %>%
  mutate(tpi = row.names(.)) %>%
  separate(tpi, into = c('trt', 'rate', 'contrast.phen', 'p.germ'), sep = '_')

# # Clear some more memory
# rm(
#   boot.kernel.df, boot.pert.kernel.df, midp.boot.kernel.df,
#   midp.boot.pert.kernel.df, midp.obsv.kernel.df, midp.pert.kernel.df,
#   midp.phen.boot.kernel.df, midp.phen.boot.pert.kernel.df, 
#   midp.phen.kernel.df, midp.phen.pert.kernel.df, obsv.pert.kernel.df
# )

# ------------------------------------------------------                  
# ------ Estimate sensitivities ------------------------
# ------------------------------------------------------

midp.obsv.sens = merge(
  midp.obsv.lambda, midp.pert.lambda,
  by = c('contrast', 'trt.phen', 'phen', 'p.germ'), suffixes = c('.orig', '.pert', 'p.germ')
) %>%
  # NOTE the delta value is hard-coded in here
  mutate(sv = (lambda.pert - lambda.orig) / 0.001) %>%
  mutate(phen = as.numeric(phen))

midp.phen.sens = merge(
  midp.phen.lambda, midp.phen.pert.lambda,
  by = c('trt', 'contrast.phen', 'p.germ'), suffixes = c('.orig', '.pert', 'p.germ')
) %>%
  rename(trt.rate = trt) %>%
  mutate(sv = ((lambda.pert - lambda.orig) / 0.001))

# # Need a -1 in here for when the contrast in phenology is positive or negative
# mutate(
#   sv = ((lambda.pert - lambda.orig) / 0.0001) * ifelse(grepl('^d', trt) | grepl('^d', contrast.phen), -1, 1)
# )

# ------------------------------------------------------                  
# ------ Get parameter differences ---------------------
# ------------------------------------------------------

# Here, just getting the paramter differences (to multiply by the sensitivities)

# For the observed datasets, these were stored in the observed perturbed data frames
# For the bootstrap dataset, these were stored in separate CSVs

obsv.param.diffs = rbind(
  gs.obsv.pert %>% distinct(trt.phen.idx, param, orig.parval), 
  fr.obsv.pert %>% distinct(trt.phen.idx, param, orig.parval)
) %>%
  merge(trt.phen.ltre.key %>% select(-mean.phen)) %>%
  select(-c(trt.phen.idx)) %>%
  rename(trt = trt.rate) %>%
  # get rid of the phen 
  filter(!grepl('phen', param)) %>%
  # Get differences between treatments
  pivot_wider(names_from = trt, values_from = orig.parval) %>%
  mutate(d.c = drought - control, i.c = irrigated - control) %>%
  select(-c(drought, irrigated, control)) %>%
  pivot_longer(c(d.c, i.c), names_to = 'contrast', values_to = 'pardiff')


# Phenology parameter differences

obsv.phen.diffs = rbind(
  gs.obsv.pert %>% distinct(trt.phen.idx, param, orig.parval), 
  fr.obsv.pert %>% distinct(trt.phen.idx, param, orig.parval)
) %>%
  merge(trt.phen.ltre.key %>% select(-mean.phen)) %>%
  select(-c(trt.phen.idx)) %>%
  rename(trt = trt.phen) %>%
  # Now, select only the phen
  filter(grepl('phen', param)) %>%
  rename(rate = param) %>%
  # Get differences between treatments
  pivot_wider(names_from = trt, values_from = orig.parval) %>%
  mutate(d.c = drought - control, i.c = irrigated - control) %>%
  select(-c(drought, irrigated, control)) %>%
  pivot_longer(c(d.c, i.c), names_to = 'contrast.phen', values_to = 'pardiff')


# ------------------------------------------------------                  
# ------ Estimate LTRE contributions -------------------
# ------------------------------------------------------

# Multipling the sensitivities by the rate differences
# So, merging and then making a column for the distinct contribution

# From observed data
obsv.ltre = merge(midp.obsv.sens, obsv.param.diffs) %>%
  mutate(contrib = pardiff * sv)

# Phenology effects
phen.ltre = merge(midp.phen.sens, obsv.phen.diffs) %>%
  mutate(contrib = pardiff * sv)

# ------------------------------------------------------                  
# ------ Check that lambda differences match LTRE sums -
# ------------------------------------------------------

obsv.dlambda.compare = merge(
  obsv.ltre %>% 
    group_by(contrast, trt.phen, p.germ) %>% 
    summarise(csum = sum(contrib)),
  obsv.lambda %>% 
    pivot_wider(names_from = trt, values_from = lambda) %>% 
    mutate(d.c = drought - control, i.c = irrigated - control) %>% 
    select(-c(control, drought, irrigated)) %>% 
    pivot_longer(c(d.c, i.c), names_to = 'contrast', values_to = 'd.lambda', values_drop_na = TRUE) %>%
    select(-c(phen.date, mean.phen))
)

obsv.dlambda.compare %>% mutate(relerr = (csum - d.lambda) / d.lambda)
# Okay better than before! 3/4 are <1% and the final one is at 2.2%...


phen.dlambda.compare = merge(
  phen.ltre %>% 
    group_by(contrast.phen, trt.rate, p.germ) %>% 
    summarise(csum = sum(contrib)),
  obsv.lambda %>% 
    select(-c(phen, mean.phen, phen.date)) %>%
    rename(trt.rate = trt) %>%
    pivot_wider(names_from = trt.phen, values_from = lambda) %>% 
    mutate(d.c = drought - control, i.c = irrigated - control) %>% 
    select(-c(control, drought, irrigated)) %>% 
    pivot_longer(c(d.c, i.c), names_to = 'contrast.phen', values_to = 'd.lambda', values_drop_na = TRUE)
)

phen.dlambda.compare %>% mutate(relerr = (csum - d.lambda) / d.lambda)
# uh...

merge(
  obsv.ltre %>% 
    group_by(contrast, trt.phen, p.germ) %>% 
    summarise(csum = sum(contrib)),
  phen.ltre %>% 
    group_by(contrast = contrast.phen, trt.rate, p.germ) %>% 
    summarise(csum = sum(contrib)),
  by = 'contrast'
) %>%
  # (want to have estimates that don't have the same treatment for phen and rate)
  filter(trt.phen != trt.rate) %>%
  mutate(csum = csum.x + csum.y) %>%
  merge(
    obsv.lambda %>% filter(trt == trt.phen) %>% select(-contains('phen')) %>%
      pivot_wider(names_from = trt, values_from = lambda) %>% 
      mutate(d.c = drought - control, i.c = irrigated - control) %>% 
      select(-c(control, drought, irrigated)) %>% 
      pivot_longer(c(d.c, i.c), names_to = 'contrast', values_to = 'd.lambda', values_drop_na = TRUE)
  ) %>%
  mutate(relerr = (csum - d.lambda) / d.lambda)


####################

obsv.lambda %>%
  mutate(p.germ = as.numeric(p.germ)) %>%
  filter(trt == trt.phen) %>%
  ggplot(aes(x = p.germ, y = lambda, group = trt)) +
  geom_line() +
  geom_point() +
  scale_x_log10()

obsv.d.lambda = obsv.lambda %>% 
  pivot_wider(names_from = trt, values_from = lambda) %>% 
  mutate(d.c = drought - control, i.c = irrigated - control) %>% 
  select(-c(control, drought, irrigated)) %>% 
  pivot_longer(c(d.c, i.c), names_to = 'contrast', values_to = 'd.lambda', values_drop_na = TRUE) %>%
  select(-c(phen.date, mean.phen))

obsv.d.lambda %>%
  filter(trt.phen %in% 'control') %>%
  mutate(p.germ = as.numeric(p.germ)) %>%
  ggplot(aes(x = p.germ, y = d.lambda, group = contrast)) +
  annotate('segment', x = min(p.germ), xend = max(p.germ), y = 0, yend = 0, colour = 'gray77', linetype = 2) +
  geom_line(aes(colour = contrast)) +
  geom_point(aes(fill = contrast), size = 3, shape = 21) +
  scale_x_log10() +
  # scale_y_continuous(limits = c(0, NA)) +
  scale_colour_manual(
    values = c('goldenrod1', 'dodgerblue'),
    labels = c('drought vs. control', 'irrigated vs. control'),
    ''
  ) +
  scale_fill_manual(
    values = c('goldenrod1', 'dodgerblue'),
    labels = c('drought vs. control', 'irrigated vs. control'),
    ''
  ) +
  labs(x = 'Probability of germination', y = expression(paste(Delta, lambda))) +
  theme(
    panel.background = element_blank(),
    legend.position = 'top'
  )

ggsave('04_analysis/figures/ltre_germ_delta_lambda.png', width = 5, height = 5)

obsv.ltre %>%
  mutate(p.germ = as.numeric(p.germ)) %>%
  ggplot(
    aes(
      x = p.germ, y = contrib, colour = contrast, 
      group = interaction(contrast, trt.phen, param)
    )
  ) +
  geom_line() +
  geom_point() +
  scale_x_log10() +
  facet_wrap(~ param)

# hmm...

# obsv.ltre %>%
#   mutate(p.germ = as.numeric(p.germ)) %>%
#   ggplot(
#     aes(
#       x = param, y = contrib, fill = contrast, 
#       group = interaction(contrast, p.germ)
#     )
#   ) +
#   geom_col(position = position_dodge()) +
#   facet_wrap(~ trt.phen + phen)
# # oh...
# 
# phen.ltre %>%
#   mutate(p.germ = as.numeric(p.germ)) %>%
#   ggplot(
#     aes(
#       x = rate, y = contrib, fill = contrast.phen, 
#       group = interaction(contrast.phen, p.germ)
#     )
#   ) +
#   geom_col(position = position_dodge()) +
#   facet_wrap(~ trt.rate + contrast.phen)

obsv.trt.ltre = obsv.ltre %>%
  separate(param, into = c('rate', 'param'), sep = '\\.') %>%
  select(-param) %>%
  group_by(contrast, trt.phen, rate, p.germ) %>%
  summarise(contrib = sum(contrib)) %>%
  ungroup()

obsv.phen.ltre = phen.ltre %>%
  mutate(rate = gsub('phen\\.', '', rate)) %>%
  group_by(trt.rate, contrast.phen, rate, p.germ) %>%
  summarise(contrib = sum(contrib)) %>%
  ungroup()

obsv.trt.ltre %>%
  mutate(p.germ = as.numeric(p.germ)) %>%
  ggplot(aes(x = rate, y = contrib, fill = contrast)) +
  geom_col() +
  facet_wrap(~ p.germ + contrast)

obsv.phen.ltre %>%
  mutate(p.germ = as.numeric(p.germ)) %>%
  ggplot(aes(x = rate, y = contrib, fill = contrast.phen)) +
  geom_col(position = position_dodge()) +
  facet_wrap(~ p.germ + contrast.phen)


control.ltre.all = rbind(
  # --- Observed treatment effects
  obsv.trt.ltre %>%
    # give me LTRE values for the control dates and remove column
    filter(trt.phen %in% 'control') %>%
    select(-trt.phen) %>%
    # marker for type of observation
    mutate(varb = 'psi', type = 'trt'),
  # --- Observed phenology effects (within treatment)
  obsv.phen.ltre %>%
    # give me LTRE values where the reference date is the control
    # and remove unnecessary column
    filter(trt.rate %in% 'control') %>%
    select(-trt.rate) %>%
    # Rename column for column agreement
    rename(contrast = contrast.phen) %>%
    mutate(varb = 'phi', type = 'phen')
) %>%
  mutate(ltre.varb = paste0(varb, '[', rate, ']'))

control.ltre.all %>%
  mutate(p.germ = as.numeric(p.germ)) %>%
  ggplot(aes(x = ltre.varb, y = contrib)) +
  geom_col(aes(fill = contrast)) +
  facet_wrap(~ p.germ + contrast)

control.ltre.all %>%
  mutate(p.germ = as.factor(as.numeric(p.germ))) %>%
  ggplot(aes(x = ltre.varb, y = contrib, group = p.germ)) +
  geom_col(aes(fill = p.germ), colour = 'gray22', position = position_dodge()) +
  scale_x_discrete(
    labels = scales::label_parse(),
    limits = c(
      'phi[grow]', 'phi[succ]', 'phi[seed]',
      'psi[grow]', 'psi[flow]', 'psi[seed]', 'psi[recr]'
    ),
    guide = guide_axis(n.dodge = 2)
  ) +
  scale_fill_viridis_d(option = 'A') +
  facet_wrap(~ contrast) +
  theme(
    panel.background = element_blank(),
    strip.background = element_part_rect(fill = 'white', side = 'b', colour = 'gray22'),
    axis.text.x = element_text(size = 9),
    axis.text.y = element_text(size = 7),
    strip.text = element_text(size = 7)
  )


obsv.contribs = rbind(
  # --- Observed treatment effects
  obsv.trt.ltre %>%
    # give me LTRE values for the control dates and remove column
    filter(trt.phen %in% 'control') %>%
    select(-trt.phen) %>%
    # marker for type of observation
    mutate(varb = 'psi', type = 'trt'),
  # --- Observed phenology effects (within treatment)
  obsv.phen.ltre %>%
    # give me LTRE values where the reference date is the control
    # and remove unnecessary column
    filter(trt.rate %in% 'control') %>%
    select(-trt.rate) %>%
    # Rename column for column agreement
    rename(contrast = contrast.phen) %>%
    mutate(varb = 'phi', type = 'phen')
)

obsv.by.demo.type = obsv.contribs %>%
  group_by(contrast, demo = ifelse(rate %in% c('grow', 'recr'), 'grow', 'repr'), type, p.germ) %>%
  summarise(contrib = sum(contrib)) %>%
  ungroup()

obsv.by.demo.type %>%
  mutate(
    contr.pretty = paste(ifelse(contrast %in% 'd.c', 'drought', 'irrigated'), 'vs. control'),
    demo = ifelse(demo %in% 'grow', 'growth', 'reproduction'),
    p.germ = as.factor(as.numeric(p.germ))
  ) %>%
  ggplot(aes(x = type, y = contrib, group = p.germ)) +
  geom_col(aes(fill = p.germ), position = position_dodge(), colour = 'gray22') +
  scale_x_discrete(
    limits = c('phen', 'trt'), labels = c('phenology', 'treatment'),
    guide = guide_axis(n.dodge = 2)
  ) +
  scale_fill_viridis_d(option = 'A', 'germintion rate') +
  labs(x = '', y = expression(paste('Contribution to ', Delta, lambda))) +
  guides(fill = guide_legend(nrow = 1)) +
  facet_nested( ~ contr.pretty + demo) +
  theme(
    strip.background = element_part_rect(fill = 'white', side = 'b', colour = 'gray22'),
    panel.background = element_blank(),
    axis.text = element_text(size = 7),
    strip.text = element_text(size = 7),
    legend.position = 'top',
    legend.direction = 'horizontal'
  )

ggsave('04_analysis/figures/ltre_germ_panb.png', width = 8, height = 5)
