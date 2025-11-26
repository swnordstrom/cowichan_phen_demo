# Script for doing two-way LTRE design with the following features
# - Phenology treatment
# - Sensitivites determined by regression coefficient perturbation
# - Bootstrapped samples for uncertainty assessment

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


# --- Read in kernels from bootstrapped resampling

gs.boot = read.csv('03_construct_kernels/out/deterministic_growsurv_bootstrap_ltre.csv')
fr.boot = read.csv('03_construct_kernels/out/deterministic_reprod_bootstrap_ltre.csv') %>%
  mutate(boot = gsub('b', '', boot))


# --- Read in perturbed kernels from bootstrapped sampling

gs.boot.pert = read.csv('03_construct_kernels/out/deterministic_growsurv_perturb_bootstraps.csv')
fr.boot.pert = read.csv('03_construct_kernels/out/deterministic_reprod_perturb_bootstraps.csv')


# --- Read in parameters used in bootstrapping
# (these give the differences in beta in the LTRE)

gs.pert.pars = read.csv('03_construct_kernels/out/deterministic_growsurv_bootstrapped_perturbed_params.csv') # %>%
  # pivot_longer(-boot, names_to = 'rate_trt', values_to = 'parval') %>%
  # separate(rate_trt, into = c('rate', 'trt'), sep = '_')

fr.pert.pars = read.csv('03_construct_kernels/out/deterministic_reprod_bootstrapped_perturbed_params.csv') 
# Do this in two steps because this file also contains the phen dates used in the bootstarp

# Dates used in phen bootstrapping
phen.boots = read.csv('03_construct_kernels/out/phenology_bootstrapped_means.csv') # %>%
  # select(c(boot, contains('phen'))) %>%
  # pivot_longer(-boot, names_to = 'trt.phen', values_to = 'phen') %>%
  # mutate(
  #   trt.phen = gsub('phen\\_', '', trt.phen),
  #   phen = as.Date(as.numeric(phen), format = '%b-%d')
  # )

# # Reproductive bootstrapped parameters
# fr.pert.pars = fr.pert.pars %>%
#   select(-contains('phen')) %>%
#   pivot_longer(-boot, names_to = 'rate_trt', values_to = 'parval') %>%
#   separate(rate_trt, into = c('rate', 'trt'), sep = '_')

# Read in LTRE treatment-phenology info
# NOTE: in some cases we will want to crop out the treatment means 
# (e.g. when using bootstrapped means)
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

# Germination probability
p.germ = .001
# p.germ = 0.0058007812

# --- Observed kernel
obsv.kernel.df = merge(
  gs.obsv, fr.obsv,
  by.x = c('size.prev', 'size.cur', 'trt.phen.idx'), 
  by.y = c('size.prev', 'size.nex', 'trt.phen.idx'),
) %>%
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
  mutate(
    p.size.cur = pred.surv * (pv.grow.size * (1 - prob.flower) + pf.grow.size * prob.flower) + (p.size.cur * p.germ)
  ) %>%
  select(-c(pred.surv, pv.grow.size, pf.grow.size, prob.flower))

  
# --- Bootstrapped kernel, unperturbed
# (this takes some time - also merging by bootstrap sample number)
boot.kernel.df = merge(
  gs.boot, # %>% pivot_longer(starts_with('b'), names_to = 'boot', values_to = 'p.size.cur'),
  fr.boot, # %>% pivot_longer(starts_with('b'), names_to = 'boot', values_to = 'p.size.cur'),
  by.x = c('size.prev', 'size.cur', 'trt.phen.idx', 'boot'), by.y = c('size.prev', 'size.nex', 'trt.phen.idx', 'boot'),
  suffixes = c('.g', '.f')
) %>%
  mutate(
    p.size.cur = pred.surv * (pv.grow.size * (1 - prob.flower) + pf.grow.size * prob.flower) + (p.size.cur * p.germ)
  ) %>%
  select(-c(pred.surv, pv.grow.size, pf.grow.size, prob.flower)) %>%
  rename(samp = boot)
  # mutate(samp = as.numeric(gsub('b', '', boot))) %>%
  # select(-c(p.size.cur.g, p.size.cur.f, boot))

# --- Perturbed bootstrap kernel
# (this takes even more time, because it's merging by bootstrap and perturbed
# parameter)
boot.pert.kernel.df = rbind(
  merge(
    gs.boot.pert, # %>% 
      # pivot_longer(starts_with('b'), names_to = 'samp', values_to = 'p.size.cur') %>%
      # mutate(samp = gsub('b', '', samp)),
    fr.boot, # %>% 
      # pivot_longer(starts_with('b'), names_to = 'samp', values_to = 'p.size.cur') %>%
      # mutate(samp = gsub('b', '', samp)),
    by.x = c('size.prev', 'size.cur', 'trt.phen.idx', 'boot'), by.y = c('size.prev', 'size.nex', 'trt.phen.idx', 'boot'),
    suffixes = c('.g', '.f')
  ) %>%
    select(
      size.prev, size.cur, trt.phen.idx, boot, param = perturb.param, 
      pred.surv, pv.grow.size, pf.grow.size, prob.flower, p.size.cur
    ),
  merge(
    gs.boot, #%>%
      # pivot_longer(starts_with('b'), names_to = 'samp', values_to = 'p.size.cur') %>%
      # mutate(samp = gsub('b', '', samp)),
    fr.boot.pert, # %>%
      # pivot_longer(starts_with('b'), names_to = 'samp', values_to = 'p.size.cur') %>%
      # mutate(samp = gsub('b', '', samp)),
    by.x = c('size.prev', 'size.cur', 'trt.phen.idx', 'boot'), by.y = c('size.prev', 'size.nex', 'trt.phen.idx', 'boot'),
    suffixes = c('.g', '.f')
  ) %>%
    select(
      size.prev, size.cur, trt.phen.idx, boot, param, 
      pred.surv, pv.grow.size, pf.grow.size, prob.flower, p.size.cur
    )
) %>%
  mutate(
    p.size.cur = pred.surv * (pv.grow.size * (1 - prob.flower) + pf.grow.size * prob.flower) + (p.size.cur * p.germ)
  ) %>%
  select(-c(pred.surv, pv.grow.size, pf.grow.size, prob.flower))

# Remove input files
rm(fr.obsv, gs.obsv, fr.boot, gs.boot, fr.boot.pert, gs.boot.pert)

# Merge in treatment info and dates

# Observed kernels
obsv.kernel.df = merge(obsv.kernel.df, trt.phen.ltre.key) %>%
  select(-trt.phen.idx) %>%
  rename(trt = trt.rate)

# Bootstrapped kernels
boot.kernel.df = merge(boot.kernel.df, trt.phen.ltre.key %>% select(-mean.phen)) %>%
  select(-trt.phen.idx) %>%
  rename(trt = trt.rate)

# Perturbed kernels
obsv.pert.kernel.df = merge(obsv.pert.kernel.df, trt.phen.ltre.key) %>%
  select(-trt.phen.idx) %>%
  rename(trt = trt.rate)

# Bootstrapped perturbed kernels
boot.pert.kernel.df = merge(boot.pert.kernel.df, trt.phen.ltre.key %>% select(-mean.phen)) %>%
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

# Bootstrapped midpoints (not perturbed)
midp.boot.kernel.df = boot.kernel.df %>%
  pivot_wider(names_from = trt, values_from = p.size.cur) %>%
  mutate(d.c = (drought + control) / 2, i.c = (irrigated + control) / 2) %>%
  select(-c(drought, control, irrigated)) %>%
  pivot_longer(c(d.c, i.c), names_to = 'contrast', values_to = 'p.size.cur') %>%
  filter(!is.na(p.size.cur))

# Bootstrapped midpoints (perturbed)
midp.boot.pert.kernel.df = boot.pert.kernel.df %>%
  # remove phen perturbations
  filter(!grepl('phen', param)) %>%
  pivot_wider(names_from = trt, values_from = p.size.cur) %>%
  mutate(d.c = (drought + control) / 2, i.c = (irrigated + control) / 2) %>%
  select(-c(drought, control, irrigated)) %>%
  pivot_longer(c(d.c, i.c), names_to = 'contrast', values_to = 'p.size.cur') %>%
  filter(!is.na(p.size.cur))

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

midp.phen.boot.kernel.df = boot.kernel.df %>%
  # NOTE: pivoting out by trt.phen instead of trt (because we're averaging
  # across different phens)
  pivot_wider(names_from = trt.phen, values_from = p.size.cur) %>%
  mutate(d.c = (drought + control) / 2, i.c = (irrigated + control) / 2) %>%
  select(-c(control, drought, irrigated)) %>%
  pivot_longer(c(d.c, i.c), names_to = 'contrast.phen', values_to = 'mentry', values_drop_na = TRUE)

midp.phen.boot.pert.kernel.df = boot.pert.kernel.df %>%
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
  obsv.kernel.df, obsv.kernel.df[,c("trt", "trt.phen", "mean.phen")], 
  sep = '_', drop = TRUE
) %>%
  # split() splits the kernel df into a list where each entry is a data frame
  # for one phen-treatment kernel
  lapply(
    function(df) {
      df %>%
        arrange(size.prev, size.cur) %>%
        pivot_wider(names_from = size.prev, values_from = p.size.cur) %>%
        select(-c(trt, trt.phen, mean.phen, size.cur)) %>%
        as.matrix()
    }
  ) %>%
  # estimate the dominant eigenvalue 
  sapply(function(m) Re(eigen(m)$values[1])) %>%
  data.frame(lambda = .) %>%
  # get treatment and phen from rownames
  mutate(ttp = row.names(.)) %>%
  separate(ttp, into = c('trt', 'trt.phen', 'phen'), sep = '_') %>%
  arrange(trt) %>%
  # convert phen to date (for plotting)
  mutate(
    mean.phen = as.numeric(phen),
    phen.date = as.Date(mean.phen, format = '%b-%d')
  )

# Build data frame with bootstrapped (non-perturbed) lambdas
boot.lambda = split(
  boot.kernel.df, boot.kernel.df[,c("trt", "trt.phen", "samp")],
  sep = '_', drop = TRUE
) %>%
  lapply(
    function(df) {
      df %>%
        arrange(size.prev, size.cur) %>%
        pivot_wider(names_from = size.prev, values_from = p.size.cur) %>%
        select(-c(trt, trt.phen, samp, size.cur)) %>%
        as.matrix()
    }
  ) %>%
  sapply(function (m) Re(eigen(m)$values[1])) %>%
  data.frame(lambda = .) %>%
  mutate(tps = row.names(.)) %>%
  separate(tps, into = c('trt', 'trt.phen', 'samp'), sep = '_')

# Get lambda for the midpoint of observed matrices
midp.obsv.lambda = split(
  midp.obsv.kernel.df, midp.obsv.kernel.df[,c("contrast", 'trt.phen', "mean.phen")],
  sep = '_', drop = TRUE
) %>%
  lapply(
    function(df) {
      df %>%
        arrange(size.prev, size.cur) %>%
        pivot_wider(names_from = size.prev, values_from = mentry) %>%
        select(-c(contrast, trt.phen, mean.phen, size.cur)) %>%
        as.matrix()
    }
  ) %>%
  sapply(function(x) Re(eigen(x)$values[1])) %>%
  data.frame(lambda = .) %>%
  mutate(c_tp = row.names(.)) %>%
  separate(c_tp, into = c('contrast', 'trt.phen', 'phen'), sep = '_')  

# Get lambdas for the midpoint perturbed matrices
midp.pert.lambda = split(
  midp.pert.kernel.df,
  midp.pert.kernel.df[,c("contrast", 'trt.phen', 'mean.phen', 'param')],
  sep = '_', drop = TRUE
) %>%
  lapply(
    function(df) {
      df %>%
        arrange(size.prev, size.nex) %>%
        pivot_wider(names_from = size.prev, values_from = mentry) %>%
        select(-c(contrast, trt.phen, param, mean.phen, size.nex)) %>%
        as.matrix()
    }
  ) %>%
  sapply(function(x) Re(eigen(x)$values[1])) %>%
  data.frame(lambda = .) %>%
  mutate(c_pg = row.names(.)) %>%
  separate(c_pg, into = c('contrast', 'trt.phen', 'phen', 'param'), sep = '_')  

# Get lambdas for the midpoint of the bootstrapped observed matrices
# (takes a sec to run - mclapply may be useful...)
midp.boot.lambda = split(
  midp.boot.kernel.df,
  midp.boot.kernel.df[,c("contrast", "trt.phen", "samp")],
  sep = '_', drop = TRUE
) %>%
  lapply(
    function(df) {
      df %>%
        arrange(size.prev, size.cur) %>%
        pivot_wider(names_from = size.prev, values_from = p.size.cur) %>%
        select(-c(contrast, samp, trt.phen, size.cur)) %>%
        as.matrix()
    }
  ) %>%
  sapply(function (m) Re(eigen(m)$values[1])) %>%
  data.frame(lambda = .) %>%
  mutate(cps = row.names(.)) %>%
  separate(cps, into = c('contrast', 'trt.phen', 'samp'), sep = '_')

# Get lambdas for midpoint of perturbed bootstrap matrices
# (also slow - slower than the above, takes about a minute)
midp.boot.pert.lambda = split(
  midp.boot.pert.kernel.df, 
  midp.boot.pert.kernel.df[,c("contrast", "trt.phen", "param", "boot")],
  sep = '_', drop = TRUE
) %>%
  lapply(
    function(df) {
      df %>%
        arrange(size.prev, size.cur) %>%
        pivot_wider(names_from = size.prev, values_from = p.size.cur) %>%
        select(-c(contrast, boot, param, trt.phen, size.cur)) %>%
        as.matrix()
    }
  ) %>%
  sapply(function(x) Re(eigen(x)$values[1])) %>%
  data.frame(lambda = .) %>%
  mutate(tpps = row.names(.)) %>%
  separate(tpps, into = c('contrast', 'trt.phen', "param", 'samp'), sep = '_')

# Phenology midpoint lambdas
midp.phen.lambda = split(
  midp.phen.kernel.df,
  midp.phen.kernel.df[,c("trt", "contrast.phen")],
  sep = '_', drop = TRUE
) %>%
  lapply(
    function(df) {
      df %>%
        arrange(size.prev, size.cur) %>%
        pivot_wider(names_from = size.prev, values_from = mentry) %>%
        select(-c(trt, contrast.phen, size.cur)) %>%
        as.matrix()
    }
  ) %>%
  sapply(function(x) Re(eigen(x)$values[1])) %>%
  data.frame(lambda = .) %>%
  mutate(tc = row.names(.)) %>%
  separate(tc, into = c('trt', 'contrast.phen'), sep = '_')

# Phenology midpoint lambdas after perturbation
midp.phen.pert.lambda = split(
  midp.phen.pert.kernel.df,
  midp.phen.pert.kernel.df[,c("trt", "param", "contrast.phen")],
  sep = '_', drop = TRUE
) %>%
  lapply(
    function(df) {
      df %>%
        arrange(size.prev, size.nex) %>%
        pivot_wider(names_from = size.prev, values_from = mentry) %>%
        select(-c(trt, contrast.phen, param, size.nex)) %>%
        as.matrix()
    }
  ) %>%
  sapply(function(x) Re(eigen(x)$values[1])) %>%
  data.frame(lambda = .) %>%
  mutate(tpi = row.names(.)) %>%
  separate(tpi, into = c('trt', 'rate', 'contrast.phen'), sep = '_')

midp.phen.boot.lambda = split(
  midp.phen.boot.kernel.df,
  midp.phen.boot.kernel.df[,c("trt", "samp", "contrast.phen")],
  sep = '_', drop = TRUE
) %>%
  lapply(
    function(df) {
      df %>%
        arrange(size.prev, size.cur) %>%
        pivot_wider(names_from = size.prev, values_from = mentry) %>%
        select(-c(trt, contrast.phen, samp, size.cur)) %>%
        as.matrix()
    }
  ) %>%
  sapply(function(x) Re(eigen(x)$values[1])) %>%
  data.frame(lambda = .) %>%
  mutate(tpi = row.names(.)) %>%
  separate(tpi, into = c('trt', 'samp', 'contrast.phen'), sep = '_')

# Lambda for phenology midpoint bootstraps
midp.phen.boot.pert.lambda = split(
  midp.phen.boot.pert.kernel.df,
  midp.phen.boot.pert.kernel.df[,c("trt", "param", "contrast.phen", "boot")],
  sep = '_', drop = TRUE
) %>%
  lapply(
    function(df) {
      df %>%
        arrange(size.prev, size.cur) %>%
        pivot_wider(names_from = size.prev, values_from = mentry) %>%
        select(-c(trt, contrast.phen, param, boot, size.cur)) %>%
        as.matrix()
    }
  ) %>%
  sapply(function(x) Re(eigen(x)$values[1])) %>%
  data.frame(lambda = .) %>%
  mutate(tpi = row.names(.)) %>%
  separate(tpi, into = c('trt', 'rate', 'contrast.phen', 'samp'), sep = '_')

# Clear some more memory
rm(
  boot.kernel.df, boot.pert.kernel.df, midp.boot.kernel.df,
  midp.boot.pert.kernel.df, midp.obsv.kernel.df, midp.pert.kernel.df,
  midp.phen.boot.kernel.df, midp.phen.boot.pert.kernel.df, 
  midp.phen.kernel.df, midp.phen.pert.kernel.df, obsv.pert.kernel.df
)

# ------------------------------------------------------                  
# ------ Estimate sensitivities ------------------------
# ------------------------------------------------------

midp.obsv.sens = merge(
  midp.obsv.lambda, midp.pert.lambda,
  by = c('contrast', 'trt.phen', 'phen'), suffixes = c('.orig', '.pert')
) %>%
  # NOTE the delta value is hard-coded in here
  mutate(sv = (lambda.pert - lambda.orig) / .0001) %>%
  mutate(phen = as.numeric(phen))

midp.boot.sens = merge(
  midp.boot.lambda, midp.boot.pert.lambda,
  by = c('contrast', 'trt.phen', 'samp'), suffixes = c('.orig', '.pert')
) %>%
  # NOTE the delta value is hard-coded here too
  mutate(sv = (lambda.pert - lambda.orig) / .0001)

midp.phen.sens = merge(
   midp.phen.lambda, midp.phen.pert.lambda,
   by = c('trt', 'contrast.phen'), suffixes = c('.orig', '.pert')
) %>%
  rename(trt.rate = trt) %>%
  mutate(sv = ((lambda.pert - lambda.orig) / 0.0001))

  # # Need a -1 in here for when the contrast in phenology is positive or negative
  # mutate(
  #   sv = ((lambda.pert - lambda.orig) / 0.0001) * ifelse(grepl('^d', trt) | grepl('^d', contrast.phen), -1, 1)
  # )

midp.phen.boot.sens = merge(
  midp.phen.boot.lambda, midp.phen.boot.pert.lambda,
  by = c('trt', 'contrast.phen', 'samp'), suffixes = c('.orig', '.pert')
) %>%
  rename(trt.rate = trt) %>%
  mutate(sv = ((lambda.pert - lambda.orig) / 0.0001))
  # Do NOT need a -1 in here because 
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

boot.param.diffs = merge(gs.pert.pars, fr.pert.pars, by = 'boot') %>%
  pivot_longer(-boot, names_to = 'partrt', values_to = 'parval') %>%
  separate_wider_delim(partrt, names = c('rate', 'trt'), delim = '_') %>%
  pivot_wider(names_from = trt, values_from = parval) %>%
  mutate(d.c = drought - control, i.c = irrigated - control) %>%
  select(-c(control, drought, irrigated)) %>%
  pivot_longer(c(d.c, i.c), names_to = 'contrast', values_to = 'pardiff', values_drop_na = TRUE) %>%
  rename(
    samp = boot,
    param = rate
  )

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

# boot.phen.diffs = phen.boots %>%
#   rename(samp = boot) %>%
#   # looks like I already exported this as wide-format
#   # mutate(phen = as.numeric(phen)) %>%
#   # pivot_wider(names_from = trt.phen, values_from = phen) %>%
#   mutate(d.c = drought - control, i.c = irrigated - control) %>%
#   select(-c(control, drought, irrigated)) %>%
#   pivot_longer(c(d.c, i.c), names_to = 'contrast.phen', values_to = 'phen.diff')
# 
# # The code above gives only the differences *from the control perspective*
# # But we will also want to get differences from the drought/irrigated perspective
# # which will require negating the phen diff
# boot.phen.diffs = rbind(
#     boot.phen.diffs %>% mutate(trt = 'control'),
#     boot.phen.diffs %>% 
#       mutate(trt = ifelse(contrast.phen %in% 'd.c', 'drought', 'irrigated')) %>% 
#       mutate(phen.diff = -1 * phen.diff)
#   )

# Checks:
# head(boot.phen.diffs)
# sum(boot.phen.diffs$phen.diff) # good - should be zero

# # Phenology: just take from the observed lambda data frame
# phen.param.diffs = obsv.lambda %>% 
#   # Get the observed phenology column
#   mutate(
#     trt.obs.phen = case_when(
#       trt %in% 'control' & floor(mean.phen) == 125 ~ 'observed',
#       trt %in% 'drought' & floor(mean.phen) == 122 ~ 'observed',
#       trt %in% 'irrigated' & floor(mean.phen) == 127 ~ 'observed',
#       .default = 'hypothetical'
#     )
#   ) %>%
#   select(-c(phen.date, phen)) %>%
#   uncount(weights = 1 + as.numeric(trt %in% 'control' & trt.obs.phen %in% 'observed')) %>%
#   group_by(trt) %>%
#   # ii column will help us distinguish between the two different control-on-control dates
#   mutate(ii = cumsum(duplicated(mean.phen))) %>%
#   ungroup() %>%
#   pivot_wider(id_cols = c(trt, ii), names_from = trt.obs.phen, values_from = c(lambda, mean.phen)) %>%
#   mutate(
#     lambda.diff = lambda_hypothetical - lambda_observed,
#     phen.diff = mean.phen_hypothetical - mean.phen_observed 
#   ) %>%
#   select(trt, ii, phen.diff, lambda.diff)

# ------------------------------------------------------                  
# ------ Estimate LTRE contributions -------------------
# ------------------------------------------------------

# Multipling the sensitivities by the rate differences
# So, merging and then making a column for the distinct contribution

# From observed data
obsv.ltre = merge(midp.obsv.sens, obsv.param.diffs) %>%
  mutate(contrib = pardiff * sv)

# From bootstraps
boot.ltre = merge(midp.boot.sens, boot.param.diffs %>% filter(!grepl('phen', param))) %>%
  mutate(contrib = pardiff * sv)

# Phenology effects
phen.ltre = merge(midp.phen.sens, obsv.phen.diffs) %>%
  mutate(contrib = pardiff * sv)

phen.boot.ltre = merge(
  midp.phen.boot.sens, 
  boot.param.diffs %>% filter(grepl('phen', param)) %>% rename(contrast.phen = contrast, rate = param)
) %>%
  mutate(contrib = pardiff * sv)

# ------------------------------------------------------                  
# ------ Check that lambda differences match LTRE sums -
# ------------------------------------------------------

obsv.dlambda.compare = merge(
    obsv.ltre %>% 
      group_by(contrast, trt.phen) %>% 
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
    group_by(contrast.phen, trt.rate) %>% 
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
# even more accurate!

merge(
  obsv.ltre %>% 
    group_by(contrast, trt.phen) %>% 
    summarise(csum = sum(contrib)),
  phen.ltre %>% 
    group_by(contrast = contrast.phen, trt.rate) %>% 
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

# Also has max error of 2.5%, everything else <1%

# ------------------------------------------------------                  
# ------ Combine contributions by rate (not param) -----
# ------------------------------------------------------

obsv.trt.ltre = obsv.ltre %>%
  separate(param, into = c('rate', 'param'), sep = '\\.') %>%
  select(-param) %>%
  group_by(contrast, trt.phen, rate) %>%
  summarise(contrib = sum(contrib)) %>%
  ungroup()
  
boot.trt.ltre = boot.ltre %>%
  separate(param, into = c('rate', 'param'), sep = '\\.') %>%
  select(-param) %>%
  group_by(contrast, trt.phen, samp, rate) %>%
  summarise(contrib = sum(contrib)) %>%
  ungroup()

obsv.phen.ltre = phen.ltre %>%
  mutate(rate = gsub('phen\\.', '', rate)) %>%
  group_by(trt.rate, contrast.phen, rate) %>%
  summarise(contrib = sum(contrib)) %>%
  ungroup()

boot.phen.ltre = phen.boot.ltre %>%
  mutate(rate = gsub('phen\\.', '', rate)) %>%
  group_by(trt.rate, contrast.phen, samp, rate) %>%
  summarise(contrib = sum(contrib)) %>%
  ungroup()

# ------------------------------------------------------                  
# ------ Crude plots -----------------------------------
# ------------------------------------------------------

# obsv.trt.ltre %>%
#   # I want to do this on the control buddate
#   filter(trt.phen %in% 'control') %>%
#   ggplot(aes(x = rate, y = contrib, fill = contrast)) +
#   geom_col(position = 'dodge')
# 
# # Picture here: more growth in treatments, less flowering
# # differing treatment effects on seed production's influence
# # differing treatment effects on recruit size (other dir.)
# 
# boot.trt.ltre %>%
#   # I want to do this on the control buddate
#   filter(trt.phen %in% 'control') %>%
#   ggplot(aes(x = rate, y = contrib, colour = contrast)) +
#   geom_point(position = position_dodge(width = 0.25), alpha = 0.5)
# 
# obsv.phen.ltre %>%
#   # think about which trt we want...
#   filter(trt.rate %in% 'control') %>%
#   ggplot(aes(x = rate, y = contrib, fill = contrast.phen)) +
#   geom_col(position = 'dodge')
# 
# boot.phen.ltre %>%
#   # We'll do the drought/irrigated differences for these
#   filter(trt.rate %in% 'control') %>%
#   ggplot(aes(x = rate, y = contrib, colour = contrast.phen)) +
#   geom_point(position = position_dodge(width = 0.25), alpha = 0.5)

# Combining...

control.ltre.all = rbind(
  # --- Observed treatment effects
  obsv.trt.ltre %>%
    # give me LTRE values for the control dates and remove column
    filter(trt.phen %in% 'control') %>%
    select(-trt.phen) %>%
    # marker for type of observation
    mutate(varb = 'psi', samp = 'obsv', type = 'trt'),
  # --- Bootstrapped treatment effects
  boot.trt.ltre %>%
    # give me LTRE values for the control dates and remove unneeded columns
    filter(trt.phen %in% 'control') %>%
    select(-c(trt.phen, samp)) %>%
    # marker for type of observation
    mutate(varb = 'psi', samp = 'boot', type = 'trt'),
  # --- Observed phenology effects (within treatment)
  obsv.phen.ltre %>%
    # give me LTRE values where the reference date is the control
    # and remove unnecessary column
    filter(trt.rate %in% 'control') %>%
    select(-trt.rate) %>%
    # Rename column for column agreement
    rename(contrast = contrast.phen) %>%
    mutate(varb = 'phi', samp = 'obsv', type = 'phen'),
  # --- Bootstrapped phenology effects
  boot.phen.ltre %>%
    # give me LTRE values where the reference date is the control
    # and remove unnecessary columns
    filter(trt.rate %in% 'control') %>%
    select(-c(trt.rate, samp)) %>%
    # Rename column for column agreement
    rename(contrast = contrast.phen) %>%
    mutate(varb = 'phi', samp = 'boot', type = 'phen')
) %>%
  mutate(ltre.varb = paste0(varb, '[', rate, ']'))

# control.ltre.all %>%
#   mutate(ltre.varb = paste0(varb, '[', rate, ']')) %>%
#   ggplot(aes(x = ltre.varb, y = contrib)) +
#   geom_point(aes(shape = samp, alpha = samp, size = samp, colour = type)) +
#   scale_alpha_manual(values = c(0.5, 1)) +
#   scale_shape_manual(values = c(1, 19)) +
#   scale_size_manual(values = c(1, 4)) +
#   scale_colour_manual(values = c('gray11', 'gray66')) +
#   scale_x_discrete(labels = scales::label_parse()) +
#   facet_wrap(~ contrast, nrow = 2)

# ugly.

control.ltre.summ = merge(
  control.ltre.all %>% filter(samp %in% 'obsv') %>% select(-c(rate, samp)),
  control.ltre.all %>%
    filter(samp %in% 'boot') %>%
    group_by(contrast, ltre.varb, varb) %>%
    reframe(
      cilim = quantile(contrib, probs = c(0.025, 0.975)),
      lohi = c('lo', 'hi')
    ) %>%
    pivot_wider(names_from = lohi, values_from = cilim)
) %>%
  ungroup()

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

# distribution of bootstrap estimates - normal?
control.ltre.all %>% 
  filter(samp %in% 'boot') %>% 
  group_by(ltre.varb, contrast) %>%
  mutate(std.contrib = (contrib - mean(contrib)) / sd(contrib)) %>%
  ggplot(aes(x = std.contrib, group = interaction(contrast, ltre.varb), colour = contrast)) + 
  geom_density(aes(colour = contrast))
# looks normal to me

# Combinations across treatments

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

boot.contribs = rbind(
  # --- Bootstrapped treatment effects
  boot.trt.ltre %>%
    # give me LTRE values for the control dates and remove unneeded columns
    filter(trt.phen %in% 'control') %>%
    select(-trt.phen) %>%
    # marker for type of observation
    mutate(varb = 'psi', type = 'trt'),
  # --- Bootstrapped phenology effects
  boot.phen.ltre %>%
    # give me LTRE values where the reference date is the control
    # and remove unnecessary columns
    filter(trt.rate %in% 'control') %>%
    select(-trt.rate) %>%
    # Rename column for column agreement
    rename(contrast = contrast.phen) %>%
    mutate(varb = 'phi', type = 'phen')
) 

head(boot.contribs)
# Aggregations to do:
# - Growth vs. reproduction
# - Treatment vs. phenology (alpha vs. beta)

obsv.by.demo.type = obsv.contribs %>%
  group_by(contrast, demo = ifelse(rate %in% c('grow', 'recr'), 'grow', 'repr'), type) %>%
  summarise(contrib = sum(contrib)) %>%
  ungroup()

boot.by.demo.type = boot.contribs %>%
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

head(boot.by.demo.type)

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


pb
ggsave('04_analysis/figures/ltre_panel_b.png', width = 5, height = 3)

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

ggsave('04_analysis/figures/ltre_panel_c.png', width = 5, height = 3)



# plot limits
# ggplot_build(pa)$layout$panel_scales_y[[1]]$range$range

# ------------------------------------------------------                  
# ------ Export plot and csvs --------------------------
# ------------------------------------------------------

### Plot

# R panel of plot (panel a will be left panel)
right.panel = plot_grid(
  pb, plot_grid(pc, pd, labels = c('ci', 'cii'), rel_widths = c(1, 1), align = 'vh'), 
  nrow = 2, labels = c('b', '')
)

# left.panel

# Export
plot_grid(pa, right.panel, ncol = 2, labels = c('a', '')) # %>%
  save_plot(filename = '04_analysis/figures/ltre_fig_allpanels.png', base_width = 8, base_height = 5)
          
### Export CSVs

# Export observed ltre contributions
write.csv(
  obsv.contribs, row.names = FALSE,
  file = '04_analysis/out/all_ltre_observed_contributions.csv'
)

# Export all bootstrapped contributions
write.csv(
  boot.contribs, row.names = FALSE,
  file = '04_analysis/out/all_ltre_bootstrapped_contributions.csv'
)

# LTRE summary with finest terms (individual alpha-beta terms)
write.csv(
  control.ltre.summ, row.names = FALSE,
  file = '04_analysis/out/overall_ltre_summary.csv'
)

# LTRE summary binned by rate (growth vs. reproduction) and effect type (treatment vs. phen)
write.csv(
  merge(obsv.by.demo.type, boot.by.demo.type.summ), row.names = FALSE,
  file = '04_analysis/out/rate-type-combo_ltre_summary.csv'
)

# LTRE summary binned by rate OR effect type (each individual contribution is counted twice here)
write.csv(
  rbind(
    merge(obsv.by.demo, boot.by.demo.summ) %>% rename(group = demo),
    merge(obsv.by.type, boot.by.type.summ) %>% rename(group = type)
  ),
  row.names = FALSE,
  file = '04_analysis/out/rate_combo_alone_ltre_summary.csv'
)

# Bootstrapped lambda values for LTRE design
write.csv(
  boot.lambda, row.names = FALSE,
  file = '04_analysis/out/ltre_design_bootstrapped_lambdas.csv'
)

cat('Done.\n')

# ------------------------------------------------------                  
# ------ Summary statistics for MS ---------------------
# ------------------------------------------------------

obsv.lambda %>% 
  filter(trt == trt.phen) %>% 
  select(trt, lambda) %>% 
  pivot_wider(names_from = trt, values_from = lambda) %>% 
  mutate(across(everything(), ~ . - control))

boot.lambda %>% 
  filter(trt == trt.phen) %>% 
  select(-trt.phen) %>% 
  pivot_wider(names_from = trt, values_from = lambda) %>% 
  mutate(across(c(control, drought, irrigated), ~ . - control)) %>%
  reframe(
    ddrought = quantile(drought, probs = c(0.025, 0.975)),
    dirrigat = quantile(irrigated, probs = c(0.025, 0.975))
  ) %>%
  mutate(across(everything(), ~ round(., 3)))

