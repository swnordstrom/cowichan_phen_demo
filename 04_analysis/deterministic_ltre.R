# Script for doing two-way LTRE design with the following features
# - Phenology treatment
# - Sensitivites determined by regression coefficient perturbation
# - Bootstrapped samples for uncertainty assessment

### ---------------------------------

library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)

rm(list = ls())

# ------------------------------------------------------                  
# ------ Read in all data ------------------------------
# ------------------------------------------------------                  

# --- Read in observed data kernels (subkernels) 
# I estimated reproductive and growth/survival subkernels separately
# Naming convention in files/object names:
# - gs is growth + survival
# - fr is flowering + reproduction

# Growth + survival subkernel
gs.obsv = read.csv('03_construct_kernels/out/deterministic_growsurv_kernel.csv') %>% 
  # (will likely remove the following call when the data frame is next updated/exported)
  select(-c(pred.surv, pred.grow.mean, p.grow.size))
# Flowering + reproduction subkernel
fr.obsv = read.csv('03_construct_kernels/out/determinstic_reprod_kernel_phen_ltre.csv') %>%  
  # Remove phen-treatment levels we haven't observed
  filter(!(trt %in% 'drought' & floor(phen) == 127)) %>%
  filter(!(trt %in% 'irrigated' & floor(phen) == 122))


# --- Read in perturbed subkernels (on observed data)
# (these get used for sensitivity analysis)

# Growth perturbed kernels
# (why is the filename 'no_phen'... growth coefficients don't have any phen at all...)
gs.obsv.pert = read.csv('03_construct_kernels/out/deterministic_grow_coef_perturbation_no_phen.csv') %>%
  rename(size.nex = size.cur, param = perturb.param, orig.parval = orig.par.val) %>%
  mutate(param = gsub('\\_', '.', param)) %>%
  # Going to take out the survival terms here because they don't vary by treatment
  # (same with the growth standard deviation)
  filter(!grepl('surv', param), !param %in% 'grow.sigma')

# Reproductive perturbed kernels
fr.obsv.pert = read.csv('03_construct_kernels/out/deterministic_repr_coef_perturbation_phen.csv') %>%
  # For the sake of ease, I'm going to remove the drought-irrigation phen
  # treatments we don't observe (these will not be part of the LTRE)
  filter(!(trt %in% 'drought' & floor(mean.phen) == 127)) %>%
  filter(!(trt %in% 'irrigated' & floor(mean.phen) == 122))


# --- Read in kernels from bootstrapped resampling

gs.boot = read.csv('03_construct_kernels/out/deterministic_growsurv_bootstrap.csv')
fr.boot = read.csv('03_construct_kernels/out/deterministic_reprod_bootstrap_ltre.csv')


# --- Read in perturbed kernels from bootstrapped sampling

gs.boot.pert = read.csv('03_construct_kernels/out/deterministic_growsurv_perturb_bootstraps.csv')
fr.boot.pert = read.csv('03_construct_kernels/out/deterministic_reprod_perturb_bootstraps.csv')


# --- Read in parameters used in bootstrapping
# (these give the differences in beta in the LTRE)

gs.pert.pars = read.csv('03_construct_kernels/out/growsurv_bootstrapped_perturbed_params.csv') %>%
  pivot_longer(-boot, names_to = 'rate_trt', values_to = 'parval') %>%
  separate(rate_trt, into = c('rate', 'trt'), sep = '_')

fr.pert.pars = read.csv('03_construct_kernels/out/reprod_bootstrapped_perturbed_params.csv') %>%
  pivot_longer(-boot, names_to = 'rate_trt', values_to = 'parval') %>%
  separate(rate_trt, into = c('rate', 'trt'), sep = '_')


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

# --- Observed kernel
obsv.kernel.df = merge(
  gs.obsv, fr.obsv,
  by.x = c('size.prev', 'size.cur', 'trt'), by.y = c('size.prev', 'size.nex', 'trt'),
  suffixes = c('.g', '.r')
) %>%
  mutate(p.size.cur = p.size.cur.g + p.size.cur.r * p.germ) %>%
  select(-c(p.size.cur.g, p.size.cur.r))


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
    by = c('size.prev', 'size.nex', 'trt'),
    suffixes = c('.g', '.f')
  ) %>%
    select(size.prev, size.nex, trt, p.size.cur.g, p.size.cur.f, param, orig.parval, mean.phen),
  merge(
    gs.obsv.pert, 
    fr.obsv %>% rename(mean.phen = phen),
    by = c('size.prev', 'size.nex', 'trt'),
    suffixes = c('.g', '.f')
  ) %>%
    select(size.prev, size.nex, trt, p.size.cur.g, p.size.cur.f, param, orig.parval, mean.phen)
) %>%
  mutate(p.size.cur = p.size.cur.g + p.size.cur.f * p.germ) %>%
  select(-c(p.size.cur.g, p.size.cur.f))


# --- Bootstrapped kernel, unperturbed
# (this takes some time - also merging by bootstrap sample number)
boot.kernel.df = merge(
  gs.boot %>% pivot_longer(starts_with('b'), names_to = 'boot', values_to = 'p.size.cur'),
  fr.boot %>% pivot_longer(starts_with('b'), names_to = 'boot', values_to = 'p.size.cur'),
  by.x = c('size.prev', 'size.cur', 'trt', 'boot'), by.y = c('size.prev', 'size.nex', 'trt', 'boot'),
  suffixes = c('.g', '.f')
) %>%
  mutate(p.size.cur = p.size.cur.g + p.germ * p.size.cur.f) %>%
  mutate(samp = gsub('b', '', boot)) %>%
  select(-c(p.size.cur.g, p.size.cur.f, boot))

# --- Perturbed bootstrap kernel
# (this takes even more time, because it's merging by bootstrap and perturbed
# parameter)
boot.pert.kernel.df = rbind(
  merge(
    gs.boot.pert %>% 
      pivot_longer(starts_with('b'), names_to = 'samp', values_to = 'p.size.cur') %>%
      mutate(samp = gsub('b', '', samp)),
    fr.boot %>% 
      pivot_longer(starts_with('b'), names_to = 'samp', values_to = 'p.size.cur') %>%
      mutate(samp = gsub('b', '', samp)),
    by.x = c('size.prev', 'size.cur', 'trt', 'samp'), by.y = c('size.prev', 'size.nex', 'trt', 'samp'),
    suffixes = c('.g', '.f')
  ) %>%
    select(size.prev, size.cur, trt, samp, param = perturb.param, mean.phen, p.size.cur.g, p.size.cur.f),
  merge(
    gs.boot %>%
      pivot_longer(starts_with('b'), names_to = 'samp', values_to = 'p.size.cur') %>%
      mutate(samp = gsub('b', '', samp)),
    fr.boot.pert %>%
      pivot_longer(starts_with('b'), names_to = 'samp', values_to = 'p.size.cur') %>%
      mutate(samp = gsub('b', '', samp)),
    by.x = c('size.prev', 'size.cur', 'trt', 'samp'), by.y = c('size.prev', 'size.nex', 'trt', 'samp'),
    suffixes = c('.g', '.f')
  ) %>%
    select(size.prev, size.cur, trt, samp, param, mean.phen, p.size.cur.g, p.size.cur.f)
) %>%
  mutate(p.size.cur = p.size.cur.g + p.size.cur.f * p.germ) %>%
  select(-c(p.size.cur.g, p.size.cur.f))

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
  pivot_wider(names_from = trt, values_from = p.size.cur) %>%
  mutate(d.c = (drought + control) / 2, i.c = (irrigated + control) / 2) %>%
  select(-c(drought, control, irrigated)) %>%
  pivot_longer(c(d.c, i.c), names_to = 'contrast', values_to = 'p.size.cur') %>%
  filter(!is.na(p.size.cur))

# Get midpoints for PHENOLOGY
# Here - to estimate the sensitivity of lambda to phenology, getting a midpoint
# kernel between phen levels

midp.phen.kernel.df = obsv.kernel.df %>%
  # Okay... going to need to pivot wide to get an average across two days.
  # With the averaging across treatments, we pivoted out the treatment columns
  # and took the averages/differences of those columns.
  # We'll want to do that here with the different phenology columns.
  # A complication is that we want the control treatment on the control date
  # twice (to contrast with each of the treatment dates).
  # It would be helpful to the observed treatment bud dates stored in a column
  mutate(
    trt.obs.phen = case_when(
      trt %in% 'control' & floor(phen) == 125 ~ 'observed',
      trt %in% 'drought' & floor(phen) == 122 ~ 'observed',
      trt %in% 'irrigated' & floor(phen) == 127 ~ 'observed',
      .default = 'hypothetical'
    )
  ) %>%
  # Start by duplicating the control @ the mean control date
  uncount(weight = 1 + as.numeric(trt %in% 'control' & trt.obs.phen %in% 'observed')) %>%
  # sort rows to make sure the pivoting works correctly
  arrange(phen) %>%
  group_by(size.prev, size.cur, trt, phen) %>%
  # this 'dupl' column will just denote the 
  mutate(dupl = duplicated(p.size.cur)) %>%
  group_by(size.prev, size.cur, trt) %>%
  # this 'ii' column will just be used for distinguishing the different
  # control-observed parts of the data frame separable from each other
  mutate(ii = cumsum(dupl)) %>%
  ungroup() %>%
  select(-c(dupl, phen)) %>%
  # Okay - now for pivoting
  # going to change the obs.phen column into something more legible
  pivot_wider(names_from = trt.obs.phen, values_from = p.size.cur) %>%
  # Get the mean of each entry
  mutate(mentry = (observed + hypothetical) / 2) %>%
  select(-c(hypothetical, observed))

# Estimate with phenology perturbed
midp.phen.pert.kernel.df = obsv.pert.kernel.df %>%
  # Give us only the phenology-related kernels
  filter(grepl('phen', param)) %>%
  # Repeat the approach used above
  mutate(
    trt.obs.phen = case_when(
      trt %in% 'control' & floor(mean.phen) == 125 ~ 'observed',
      trt %in% 'drought' & floor(mean.phen) == 122 ~ 'observed',
      trt %in% 'irrigated' & floor(mean.phen) == 127 ~ 'observed',
      .default = 'hypothetical'
    )
  ) %>%
  # Start by duplicating the control @ the mean control date
  uncount(weight = 1 + as.numeric(trt %in% 'control' & trt.obs.phen %in% 'observed')) %>%
  # sort rows to make sure the pivoting works correctly
  arrange(mean.phen) %>%
  group_by(size.prev, size.nex, trt, param, mean.phen) %>%
  # this 'dupl' column will just denote the 
  mutate(dupl = duplicated(p.size.cur)) %>%
  group_by(size.prev, size.nex, trt, param) %>%
  # this 'ii' column will just be used for distinguishing the different
  # control-observed parts of the data frame separable from each other
  mutate(ii = cumsum(dupl)) %>%
  ungroup() %>%
  select(-c(dupl, orig.parval, mean.phen)) %>%
  # Okay - now for pivoting
  # going to change the obs.phen column into something more legible
  pivot_wider(names_from = trt.obs.phen, values_from = p.size.cur) %>%
  # Get the mean of each entry
  mutate(mentry = (observed + hypothetical) / 2) %>%
  select(-c(hypothetical, observed))

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
obsv.lambda = split(obsv.kernel.df, obsv.kernel.df[,c("trt", "phen")], sep = '_', drop = TRUE) %>%
  # split() splits the kernel df into a list where each entry is a data frame
  # for one phen-treatment kernel
  lapply(
    function(df) {
      df %>%
        arrange(size.prev, size.cur) %>%
        pivot_wider(names_from = size.prev, values_from = p.size.cur) %>%
        select(-c(trt, phen, size.cur)) %>%
        as.matrix()
    }
  ) %>%
  # estimate the dominant eigenvalue 
  sapply(function(m) Re(eigen(m)$values[1])) %>%
  data.frame(lambda = .) %>%
  # get treatment and phen from rownames
  mutate(trt_phen = row.names(.)) %>%
  separate(trt_phen, into = c('trt', 'phen'), sep = '_') %>%
  arrange(trt) %>%
  # convert phen to date (for plotting)
  mutate(
    mean.phen = as.numeric(phen),
    phen.date = as.Date(mean.phen, format = '%b-%d')
  )

# Build data frame with bootstrapped (non-perturbed) lambdas
boot.lambda = split(
  boot.kernel.df,
  boot.kernel.df[,c("trt", "mean.phen", "samp")],
  sep = '_', drop = TRUE
) %>%
  lapply(
    function(df) {
      df %>%
        arrange(size.prev, size.cur) %>%
        pivot_wider(names_from = size.prev, values_from = p.size.cur) %>%
        select(-c(trt, samp, mean.phen, size.cur)) %>%
        as.matrix()
    }
  ) %>%
  sapply(function (m) Re(eigen(m)$values[1])) %>%
  data.frame(lambda = .) %>%
  mutate(tps = row.names(.)) %>%
  separate(tps, into = c('trt', 'mean.phen', 'samp'), sep = '_')

# Get lambda for the midpoint of observed matrices
midp.obsv.lambda = split(
  midp.obsv.kernel.df,
  midp.obsv.kernel.df[,c("contrast", "phen")],
  sep = '_', drop = TRUE
) %>%
  lapply(
    function(df) {
      df %>%
        arrange(size.prev, size.cur) %>%
        pivot_wider(names_from = size.prev, values_from = mentry) %>%
        select(-c(contrast, phen, size.cur)) %>%
        as.matrix()
    }
  ) %>%
  sapply(function(x) Re(eigen(x)$values[1])) %>%
  data.frame(lambda = .) %>%
  mutate(c_pg = row.names(.)) %>%
  separate(c_pg, into = c('contrast', 'phen'), sep = '_')  

# Get lambdas for the midpoint perturbed matrices
midp.pert.lambda = split(
  midp.pert.kernel.df,
  midp.pert.kernel.df[,c("contrast", 'param', "mean.phen")],
  sep = '_', drop = TRUE
) %>%
  lapply(
    function(df) {
      df %>%
        arrange(size.prev, size.nex) %>%
        pivot_wider(names_from = size.prev, values_from = mentry) %>%
        select(-c(contrast, mean.phen, param, size.nex)) %>%
        as.matrix()
    }
  ) %>%
  sapply(function(x) Re(eigen(x)$values[1])) %>%
  data.frame(lambda = .) %>%
  mutate(c_pg = row.names(.)) %>%
  separate(c_pg, into = c('contrast', 'param', 'phen'), sep = '_')  

# Get lambdas for the midpoint of the bootstrapped observed matrices
# (takes a sec to run - mclapply may be useful...)
midp.boot.lambda = split(
  midp.boot.kernel.df,
  midp.boot.kernel.df[,c("contrast", "mean.phen", "samp")],
  sep = '_', drop = TRUE
) %>%
  lapply(
    function(df) {
      df %>%
        arrange(size.prev, size.cur) %>%
        pivot_wider(names_from = size.prev, values_from = p.size.cur) %>%
        select(-c(contrast, samp, mean.phen, size.cur)) %>%
        as.matrix()
    }
  ) %>%
  sapply(function (m) Re(eigen(m)$values[1])) %>%
  data.frame(lambda = .) %>%
  mutate(cps = row.names(.)) %>%
  separate(cps, into = c('contrast', 'mean.phen', 'samp'), sep = '_')

# Get lambdas for midpoint of perturbed bootstrap matrices
# (also slow - slower than the above, takes about a minute)
midp.boot.pert.lambda = split(
  midp.boot.pert.kernel.df, 
  midp.boot.pert.kernel.df[,c("contrast", "mean.phen", "param", "samp")],
  sep = '_', drop = TRUE
) %>%
  lapply(
    function(df) {
      df %>%
        arrange(size.prev, size.cur) %>%
        pivot_wider(names_from = size.prev, values_from = p.size.cur) %>%
        select(-c(contrast, samp, param, mean.phen, size.cur)) %>%
        as.matrix()
    }
  ) %>%
  sapply(function(x) Re(eigen(x)$values[1])) %>%
  data.frame(lambda = .) %>%
  mutate(tpps = row.names(.)) %>%
  separate(tpps, into = c('contrast', 'mean.phen', "param", 'samp'), sep = '_') %>%
  mutate(mean.phen = as.numeric(mean.phen))

# Phenology midpoint lambdas
midp.phen.lambda = split(
  midp.phen.kernel.df,
  midp.phen.kernel.df[,c("trt", "ii")],
  sep = '_', drop = TRUE
) %>%
  lapply(
    function(df) {
      df %>%
        arrange(size.prev, size.cur) %>%
        pivot_wider(names_from = size.prev, values_from = mentry) %>%
        select(-c(trt, ii, size.cur)) %>%
        as.matrix()
    }
  ) %>%
  sapply(function(x) Re(eigen(x)$values[1])) %>%
  data.frame(lambda = .) %>%
  mutate(ti = row.names(.)) %>%
  separate(ti, into = c('trt', 'ii'), sep = '_')

# Phenology midpoint lambdas after perturbation
midp.phen.pert.lambda = split(
  midp.phen.pert.kernel.df,
  midp.phen.pert.kernel.df[,c("trt", "param", "ii")],
  sep = '_', drop = TRUE
) %>%
  lapply(
    function(df) {
      df %>%
        arrange(size.prev, size.nex) %>%
        pivot_wider(names_from = size.prev, values_from = mentry) %>%
        select(-c(trt, ii, param, size.nex)) %>%
        as.matrix()
    }
  ) %>%
  sapply(function(x) Re(eigen(x)$values[1])) %>%
  data.frame(lambda = .) %>%
  mutate(tpi = row.names(.)) %>%
  separate(tpi, into = c('trt', 'rate', 'ii'), sep = '_')


# ------------------------------------------------------                  
# ------ Estimate sensitivities ------------------------
# ------------------------------------------------------

midp.obsv.sens = merge(
  midp.obsv.lambda, midp.pert.lambda,
  by = c('contrast', 'phen'), suffixes = c('.orig', '.pert')
) %>%
  # NOTE the delta value is hard-coded in here
  mutate(sv = (lambda.pert - lambda.orig) / .0001) %>%
  mutate(phen = as.numeric(phen))

midp.boot.sens = merge(
  midp.boot.lambda, midp.boot.pert.lambda,
  by = c('contrast', 'mean.phen', 'samp'), suffixes = c('.orig', '.pert')
) %>%
  # NOTE the delta value is hard-coded here too
  mutate(sv = (lambda.pert - lambda.orig) / .0001) %>%
  mutate(mean.phen = as.numeric(mean.phen))

midp.phen.sens = merge(
   midp.phen.lambda, midp.phen.pert.lambda,
   by = c('trt', 'ii'), suffixes = c('.orig', '.pert')
) %>%
  mutate(sv = (lambda.pert - lambda.orig) / 0.0001)


# ------------------------------------------------------                  
# ------ Get parameter differences ---------------------
# ------------------------------------------------------

# Here, just getting the paramter differences (to multiply by the sensitivities)

# For the observed datasets, these were stored in the observed perturbed data frames
# For the bootstrap dataset, these were stored in separate CSVs

obsv.param.diffs = rbind(gs.obsv.pert, fr.obsv.pert %>% select(-mean.phen)) %>%
  # get rid of the phen 
  filter(!grepl('phen', param)) %>%
  # distinct (because we only need original parameter values once)
  distinct(trt, param, orig.parval) %>%
  # Get differences between treatments
  pivot_wider(names_from = trt, values_from = orig.parval) %>%
  mutate(d.c = drought - control, i.c = irrigated - control) %>%
  select(-c(drought, irrigated, control)) %>%
  pivot_longer(c(d.c, i.c), names_to = 'contrast', values_to = 'pardiff')

boot.param.diffs = rbind(gs.pert.pars, fr.pert.pars) %>%
  pivot_wider(names_from = trt, values_from = parval) %>%
  mutate(d.c = drought - control, i.c = irrigated - control) %>%
  select(-c(control, drought, irrigated)) %>%
  pivot_longer(c(d.c, i.c), names_to = 'contrast', values_to = 'pardiff') %>%
  filter(!is.na(pardiff)) %>%
  rename(
    samp = boot,
    param = rate
  )

# Phenology: just take from the observed lambda data frame
phen.param.diffs = obsv.lambda %>% 
  # Get the observed phenology column
  mutate(
    trt.obs.phen = case_when(
      trt %in% 'control' & floor(mean.phen) == 125 ~ 'observed',
      trt %in% 'drought' & floor(mean.phen) == 122 ~ 'observed',
      trt %in% 'irrigated' & floor(mean.phen) == 127 ~ 'observed',
      .default = 'hypothetical'
    )
  ) %>%
  select(-c(phen.date, phen)) %>%
  uncount(weights = 1 + as.numeric(trt %in% 'control' & trt.obs.phen %in% 'observed')) %>%
  group_by(trt) %>%
  # ii column will help us distinguish between the two different control-on-control dates
  mutate(ii = cumsum(duplicated(mean.phen))) %>%
  ungroup() %>%
  pivot_wider(id_cols = c(trt, ii), names_from = trt.obs.phen, values_from = c(lambda, mean.phen)) %>%
  mutate(
    lambda.diff = lambda_hypothetical - lambda_observed,
    phen.diff = mean.phen_hypothetical - mean.phen_observed 
  ) %>%
  select(trt, ii, phen.diff, lambda.diff)

# ------------------------------------------------------                  
# ------ Estimate LTRE contributions -------------------
# ------------------------------------------------------

# Multipling the sensitivities by the rate differences
# So, merging and then making a column for the distinct contribution

# From observed data
obsv.ltre = merge(midp.obsv.sens, obsv.param.diffs) %>%
  mutate(contrib = pardiff * sv)

# From bootstraps
boot.ltre = merge(midp.boot.sens, boot.param.diffs) %>%
  mutate(contrib = pardiff * sv)

# Phenology effects
phen.ltre = merge(midp.phen.sens, phen.param.diffs) %>%
  mutate(contrib = phen.diff * sv)


# ------------------------------------------------------                  
# ------ OLD CODE BELOW --------------------------------
# ------------------------------------------------------

# will be going through here nad cleaning some up over the next day or two

coef.ltre = merge(
  param.diffs, mean.sensitivities %>% select(-c(lambda.orig, lambda.pert)),
  by = c('contrast', 'param')
) %>%
  mutate(contribution = pardiff * sv)

head(coef.ltre)

coef.ltre.by.rate = coef.ltre %>%
  separate(param, into = c('rate', 'param'), sep = '\\.') %>%
  group_by(contrast, rate, phen) %>%
  summarise(contribution = sum(contribution))

coef.ltre.by.rate %>%
  ggplot(aes(x = rate, y = contribution, fill = contrast)) +
  geom_col() +
  scale_fill_manual(values = c('red', 'blue')) +
  facet_wrap(phen ~ contrast)
# (unsurprisingly, very similar!)

# Okay... let's now compare with the actual lambda differences...

ltre.compare = ltre.lambda %>%
  filter(!(trt %in% 'drought' & floor(phen) == 127)) %>%
  filter(!(trt %in% 'irrigated' & floor(phen) == 122)) %>%
  select(-c(phen.date, obs.phen)) %>%
  pivot_wider(names_from = trt, values_from = lambda) %>%
  mutate(d.c = drought - control, i.c = irrigated - control) %>%
  select(-c(control, drought, irrigated)) %>%
  pivot_longer(c(d.c, i.c), names_to = 'contrast', values_to = 'd.lambda') %>%
  filter(!is.na(d.lambda)) %>%
  merge(coef.ltre.by.rate %>% group_by(contrast, phen) %>% summarise(sum.contrib = sum(contribution)))

ltre.compare
# HAHA HA YES, YES!!

# Next... phen LTRE step

# phen.compare = reprodct.pert %>%
#   distinct(trt, mean.phen) %>%
#   mutate(
#     obs.phen = (trt %in% 'control' & floor(mean.phen) == 125) |
#       (trt %in% 'drought' & floor(mean.phen) == 122) | (trt %in% 'irrigated' & floor(mean.phen) == 127),
#   ) %>%
#   group_by(trt) %>%
#   reframe(phen.diff = diff(mean.phen), sign = diff(obs.phen))

# Differences in lambda and phenology by treatment
# going to do this in an ugly but effective way

phen.lambda.diff = ltre.lambda %>% 
  filter(!(trt %in% 'drought' & floor(phen) %in% 127), !(trt %in% 'irrigated' & floor(phen) %in% 122)) %>%
  select(-phen.date) %>%
  uncount(weights = 1 + as.numeric(trt %in% 'control' & obs.phen)) %>%
  group_by(trt) %>%
  mutate(
    # This column will let me do this easily
    ii = floor((0:(n()-1)) / 2),
    # formatting this int a  readable string for pivoting
    obs.phen = ifelse(obs.phen, 'observed', 'hypothetical')
  ) %>%
  ungroup() %>%
  pivot_wider(id_cols = c(trt, ii), names_from = obs.phen, values_from = c(lambda, phen)) %>%
  mutate(
    lambda.diff = lambda_hypothetical - lambda_observed,
    phen.diff = phen_hypothetical - phen_observed 
  ) %>%
  select(trt, ii, phen.diff, lambda.diff)

# Okay... now, get mean matrices

# # Testing some code below
# testy = all.perts %>%
#   mutate(
#     obs.phen = (trt %in% 'control' & floor(mean.phen) == 125) |
#       (trt %in% 'drought' & floor(mean.phen) == 122) | (trt %in% 'irrigated' & floor(mean.phen) == 127)
#   ) %>%
#   # Duplicating the control @ control phen entries to get 
#   uncount(weight = 1 + as.numeric(trt %in% 'control' & obs.phen)) %>%
#   arrange(mean.phen) %>%
#   select(-c(orig.parval, mean.phen)) %>%
#   group_by(size.prev, size.nex, trt, param, obs.phen) %>%
#   mutate(dupl = duplicated(p.size.cur)) %>%
#   group_by(size.prev, size.nex, trt, param) %>%
#   mutate(ii = cumsum(dupl)) %>%
#   ungroup()
# 
# # these should all be control...
# testy %>% filter(ii > 0) %>% distinct(trt)
# # sweet
# # should be some values here where obs.phen is false
# testy %>% filter(ii > 0) %>% distinct(obs.phen)
# # sweet
# # going to assume this is good

mean.phen.perts = all.perts %>%
  filter(grepl('phen', param)) %>%
  mutate(
    obs.phen = (trt %in% 'control' & floor(mean.phen) == 125) |
      (trt %in% 'drought' & floor(mean.phen) == 122) | (trt %in% 'irrigated' & floor(mean.phen) == 127)
  ) %>%
  # Duplicating the control @ control phen entries to get 
  uncount(weight = 1 + as.numeric(trt %in% 'control' & obs.phen)) %>%
  arrange(mean.phen) %>%
  group_by(size.prev, size.nex, trt, param, mean.phen) %>%
  mutate(dupl = duplicated(p.size.cur)) %>%
  group_by(size.prev, size.nex, trt, param) %>%
  mutate(ii = cumsum(dupl)) %>%
  ungroup() %>%
  select(-c(dupl, orig.parval, mean.phen)) %>%
  # Okay - now for pivoting
  # going to change the obs.phen column into something more legible
  mutate(obs.phen = ifelse(obs.phen, 'observed', 'hypothetical')) %>%
  pivot_wider(names_from = obs.phen, values_from = p.size.cur) %>%
  # Get the mean of each entry
  mutate(mentry = (observed + hypothetical) / 2) %>%
  select(-c(hypothetical, observed))

mean.phen.pert.lambdas = split(
  mean.phen.perts,
  mean.phen.perts[,c("trt", "param", "ii")],
  sep = '_', drop = TRUE
) %>%
  lapply(
    function(df) {
      df %>%
        arrange(size.prev, size.nex) %>%
        pivot_wider(names_from = size.prev, values_from = mentry) %>%
        select(-c(trt, ii, param, size.nex)) %>%
        as.matrix()
    }
  ) %>%
  sapply(function(x) Re(eigen(x)$values[1])) %>%
  data.frame(lambda = .) %>%
  mutate(tpi = row.names(.)) %>%
  separate(tpi, into = c('trt', 'rate', 'ii'), sep = '_')

# Now get means of unperturbed lambdas

mean.phen = ltre.kernel.df %>%
  # Get unnecessary phen-trt combos
  filter(!(trt %in% 'drought' & floor(phen) == 127)) %>%
  filter(!(trt %in% 'irrigated' & floor(phen) == 122)) %>%
  mutate(
    obs.phen = (trt %in% 'control' & floor(phen) == 125) |
      (trt %in% 'drought' & floor(phen) == 122) | (trt %in% 'irrigated' & floor(phen) == 127)
  ) %>%
  # Duplicating the control @ control phen entries to get 
  uncount(weight = 1 + as.numeric(trt %in% 'control' & obs.phen)) %>%
  arrange(phen) %>%
  group_by(size.prev, size.cur, trt, phen) %>%
  mutate(dupl = duplicated(p.size.cur)) %>%
  group_by(size.prev, size.cur, trt) %>%
  mutate(ii = cumsum(dupl)) %>%
  ungroup() %>%
  select(-c(dupl, phen)) %>%
  # Okay - now for pivoting
  # going to change the obs.phen column into something more legible
  mutate(obs.phen = ifelse(obs.phen, 'observed', 'hypothetical')) %>%
  pivot_wider(names_from = obs.phen, values_from = p.size.cur) %>%
  # Get the mean of each entry
  mutate(mentry = (observed + hypothetical) / 2) %>%
  select(-c(hypothetical, observed)) 

mean.phen.lambdas = split(
  mean.phen,
  mean.phen[,c("trt", "ii")],
  sep = '_', drop = TRUE
) %>%
  lapply(
    function(df) {
      df %>%
        arrange(size.prev, size.cur) %>%
        pivot_wider(names_from = size.prev, values_from = mentry) %>%
        select(-c(trt, ii, size.cur)) %>%
        as.matrix()
    }
  ) %>%
  sapply(function(x) Re(eigen(x)$values[1])) %>%
  data.frame(lambda = .) %>%
  mutate(ti = row.names(.)) %>%
  separate(ti, into = c('trt', 'ii'), sep = '_')

phen.sensitivities = merge(
  mean.phen.lambdas, mean.phen.pert.lambdas,
  by = c('trt', 'ii'),
  suffixes = c('.orig', '.pert')
) %>%
  mutate(sv = (lambda.pert - lambda.orig) / .0001)

phen.sensitivities

coef.phen.ltre = merge(phen.sensitivities, phen.lambda.diff) %>%
  # somewhere there's a sign change... fix later
  mutate(contribution = -1 * phen.diff * sv)

# Compare
coef.phen.ltre %>% group_by(trt, ii) %>%
  summarise(
    lambda.diff = mean(lambda.diff),
    total.contr = sum(contribution)
  )
# Hell. Fucking. Yeah.

### -----

# Make an LTRE figure

ltre.all.contrs = rbind(
  coef.ltre.by.rate %>%
    filter(floor(phen) == 125) %>%
    select(-phen) %>%
    mutate(ltre.var = 'alpha') %>%
    rename(trt = contrast) %>%
    mutate(trt = ifelse(trt %in% 'd.c', 'drought', 'irrigated')),
  coef.phen.ltre %>%
    select(trt, rate, contribution) %>%
    filter(!trt %in% 'control') %>%
    mutate(rate = gsub('phen\\.', '', rate)) %>%
    mutate(ltre.var = 'beta')
)

ltre.all.contrs %>%
  mutate(
    rate = case_match(
      rate,
      'flow' ~ 'flowering.prob',
      'grow' ~ 'growth',
      'seed' ~ 'seed.production',
      'succ' ~ 'umbel.succ',
      'recr' ~ 'recruit.size'
    )
  ) %>%
  mutate(contr.label = paste0(ltre.var, '[', rate, ']')) %>%
  mutate(trt = paste(trt, 'vs. control')) %>%
  ggplot(aes(x = contr.label, y = contribution)) +
  geom_col_pattern(
    aes(fill = trt, pattern_spacing = ltre.var),
    pattern_fill = 'gray22', pattern_density = 0.05
  ) +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_pattern_spacing_manual(values = c(0.05, 0.1)) +
  scale_fill_manual(values = c('red', 'blue')) +
  labs(x = 'variable', y = 'contribution') +
  facet_wrap(~ trt, nrow = 2) +
  theme(
    legend.position = 'none',
    panel.background = element_blank()
  )

ggsave('04_analysis/figures/draft_figures/ltre_2024-08-07.png', width = 5, height = 5)

all.lambda %>%
  mutate(phen.date = as.Date(as.numeric((phen)), format = '%b-%d')) %>%
  ggplot(aes(x = phen.date, y = lambda, group = trt, colour = trt)) +
  geom_line() +
  scale_colour_manual(values = c('black', 'red', 'blue')) +
  geom_point(
    data = ltre.lambda %>% filter(obs.phen | floor(phen) == 125), 
    aes(shape = obs.phen), 
    size = 4
  ) +
  scale_shape_manual(values = c(1, 19)) +
  guides(shape = 'none') +
  labs(x = 'bud date') +
  theme(
    panel.background = element_blank(),
    legend.position = 'top'
  )

ggsave('04_analysis/figures/draft_figures/lambdas_2024-08-07.png', width = 5, height = 5)


# -----------
# ----------- Look at bootstrapped LTRE contributions
# -----------


pert.lambda = split(
  all.pert.kernel, all.pert.kernel[,c("trt", "param", "mean.phen", "samp")], 
  sep = "_", drop = TRUE
) %>%
  lapply(
    function(df) {
      df %>%
        arrange(size.prev, size.cur) %>%
        pivot_wider(names_from = size.prev, values_from = p.size.cur) %>%
        select(-c(trt, param, samp, mean.phen, size.cur)) %>%
        as.matrix()
    }
  ) %>%
  sapply(function(x) Re(eigen(x)$values[1])) %>%
  data.frame(lambda = .) %>%
  mutate(tpps = row.names(.)) %>%
  separate(tpps, into = c('trt', 'param', 'mean.phen', 'samp'), sep = '_') %>%
  mutate(mean.phen = as.numeric(mean.phen))

# Wait... no we don't want the lambda differences here...
# pert.lambda.diff = all.pert.lambda %>%
#   pivot_wider(names_from = trt, values_from = lambda) %>%
#   mutate(d.c = drought - control, i.c = irrigated - control) %>%
#   filter(!is.na(d.c), !is.na(i.c)) %>%
#   pivot_longer(c(d.c, i.c), names_to = 'contrast', values_to = 'd.lambda')

# Estimating sensitivities of lambda to vital rates for each treatment
boot.sensitivities = merge(
  ltre.lambda, pert.lambda, 
  by = c('trt', 'mean.phen', 'samp'), suffixes = c('.orig', '.pert')
) %>%
  mutate(sv = (lambda.pert - lambda.orig) / 0.0001)
# HOWEVER, this is not what we use for the LTRE... need the midpoints for our
# treatment contrasts

# (requires re-doing several things lmaooo)

midp.boot = ltre.boot %>%
  pivot_wider(names_from = trt, values_from = p.size.cur) %>%
  mutate(d.c = (drought + control) / 2, i.c = (irrigated + control) / 2) %>%
  select(-c(drought, control, irrigated)) %>%
  pivot_longer(c(d.c, i.c), names_to = 'contrast', values_to = 'p.size.cur') %>%
  filter(!is.na(p.size.cur))

midp.boot.lambda = split(
  midp.boot, midp.boot[,c("contrast", "mean.phen", "samp")],
  sep = '_', drop = TRUE
) %>%
  lapply(
    function(df) {
      df %>%
        arrange(size.prev, size.cur) %>%
        pivot_wider(names_from = size.prev, values_from = p.size.cur) %>%
        select(-c(contrast, samp, mean.phen, size.cur)) %>%
        as.matrix()
    }
  ) %>%
  sapply(function(x) Re(eigen(x)$values[1])) %>%
  data.frame(lambda = .) %>%
  mutate(tpps = row.names(.)) %>%
  separate(tpps, into = c('contrast', 'mean.phen', 'samp'), sep = '_') %>%
  mutate(mean.phen = as.numeric(mean.phen))

midp.pert.boot = pert.kernel %>%
  pivot_wider(names_from = trt, values_from = p.size.cur) %>%
  mutate(d.c = (drought + control) / 2, i.c = (irrigated + control) / 2) %>%
  select(-c(drought, control, irrigated)) %>%
  pivot_longer(c(d.c, i.c), names_to = 'contrast', values_to = 'p.size.cur') %>%
  filter(!is.na(p.size.cur))

midp.pert.boot.lambda = split(
  midp.pert.boot, midp.pert.boot[,c("contrast", "mean.phen", "param", "samp")],
  sep = '_', drop = TRUE
) %>%
  lapply(
    function(df) {
      df %>%
        arrange(size.prev, size.cur) %>%
        pivot_wider(names_from = size.prev, values_from = p.size.cur) %>%
        select(-c(contrast, samp, param, mean.phen, size.cur)) %>%
        as.matrix()
    }
  ) %>%
  sapply(function(x) Re(eigen(x)$values[1])) %>%
  data.frame(lambda = .) %>%
  mutate(tpps = row.names(.)) %>%
  separate(tpps, into = c('contrast', 'mean.phen', "param", 'samp'), sep = '_') %>%
  mutate(mean.phen = as.numeric(mean.phen))

midp.boot.sensitivities = merge(
  midp.boot.lambda, midp.pert.boot.lambda,
  by = c('contrast', 'mean.phen', 'samp'), suffixes = c('.orig', '.pert')
) %>%
  mutate(sv = (lambda.pert - lambda.orig) / 0.0001)

head(midp.boot.sensitivities)
# very cool!

# Now... to get those rate differences...

gs.pert.pars = read.csv('03_construct_kernels/out/growth_bootstrapped_perturbed_params.csv') %>%
  pivot_longer(-boot, names_to = 'rate_trt', values_to = 'parval') %>%
  separate(rate_trt, into = c('rate', 'trt'), sep = '_')

fr.pert.pars = read.csv('03_construct_kernels/out/reprod_bootstrapped_perturbed_params.scv') %>%
  pivot_longer(-boot, names_to = 'rate_trt', values_to = 'parval') %>%
  separate(rate_trt, into = c('rate', 'trt'), sep = '_')

boot.pert.pars = rbind(gs.pert.pars, fr.pert.pars) %>%
  pivot_wider(names_from = trt, values_from = parval) %>%
  mutate(d.c = drought - control, i.c = irrigated - control) %>%
  select(-c(control, drought, irrigated)) %>%
  pivot_longer(c(d.c, i.c), names_to = 'contrast', values_to = 'pardiff') %>%
  filter(!is.na(pardiff)) %>%
  rename(samp = boot)

boot.ltre = merge(
  boot.pert.pars, midp.boot.sensitivities,
  by.x = c('contrast', 'rate', 'samp'), 
  by.y = c('contrast', 'param', 'samp')
) %>%
  mutate(ltre.contrib = pardiff * sv)

boot.ltre

boot.ltre %>%
  mutate(phen.date = as.Date(mean.phen, format = '%b-%d')) %>%
  ggplot(aes(x = rate, y = ltre.contrib, colour = ltre.contrib > 0)) + 
  geom_point(position = position_jitter(width = 0.125), alpha = 0.5) +
  geom_point(
    data = coef.ltre %>% mutate(phen.date = as.Date(phen)),
    inherit.aes = FALSE,
    aes(x = param, y = contribution),
    shape = 4, size = 4
  ) +
  facet_wrap(contrast ~ phen.date)
# Nice
# wait... lmao this looks

# Combine by vital rate/process

boot.ltre.by.rate = boot.ltre %>%
  separate(rate, into = c('rate', 'rate.param'), sep = '\\.') %>%
  group_by(contrast, mean.phen, rate, samp) %>%
  summarise(ltre.contrib = sum(ltre.contrib)) %>%
  group_by(contrast, rate, mean.phen) %>%
  # not the sleekest way to do this but whatever
  mutate(
    ci.hi = quantile(ltre.contrib, 0.975),
    ci.lo = quantile(ltre.contrib, 0.025),
    medin = median(ltre.contrib)
  ) %>%
  ungroup()

boot.ltre.by.rate %>%
  mutate(phen.date = as.Date(mean.phen, format = '%b-%d')) %>%
  ggplot(aes(x = rate)) +
  geom_point(
    aes(y = ltre.contrib, colour = ltre.contrib > 0), 
    position = position_jitter(width = 0.125), alpha = 0.5
  ) +
  # The mean contributions look a lot larger than the bootstraps... not good.
  geom_point(
    data = coef.ltre.by.rate %>%
      # rename(phen.)
      filter(!(contrast %in% 'd.c' & floor(phen) == 127)) %>%
      filter(!(contrast %in% 'i.c' & floor(phen) == 122)) %>%
      mutate(phen.date = as.Date(phen)),
    aes(x = rate, y = contribution, fill = contribution > 0),
    shape = 23, size = 3
  ) +
  geom_point(aes(y = medin, fill = medin > 0), shape = 21, size = 3) +
  geom_segment(aes(xend = rate, y = ci.lo, yend = ci.hi), colour = 'gray33') +
  facet_wrap(contrast ~ phen.date)

# AHH it's because I was doing the summing incorrectly (doing means instead of sums lmaoo)
# BUT IT SEEMS TO WORK!

boot.ltre.by.rate %>%
  group_by(contrast, mean.phen, rate) %>%
  summarise(p.sign = max(mean(ltre.contrib > 0), mean(ltre.contrib < 0)))

# next up - check to make sure these sum up to the correct differences in lambda?

# Compare our lambdas...
compare.boot.lambda = merge(
  boot.lambda, ltre.lambda,
  by.x = c('trt', 'mean.phen'), by.y = c('trt', 'phen'),
  suffixes = c('.boot', '.ptes')
)

compare.boot.lambda %>%
  ggplot(aes(x = phen.date, colour = trt, group = trt)) +
  geom_point(aes(y = lambda.boot), position = position_dodge(width = 0.5)) +
  geom_point(
    aes(y = lambda.ptes), position = position_dodge(width = 0.5), 
    shape = 4, size = 8
  ) +
  scale_colour_manual(values = c('black', 'red', 'blue'))

# okay this looks more symmetric now

compare.boot.contrib = merge(
  coef.ltre.by.rate, boot.ltre.by.rate,
  by.x = c('contrast', 'phen', 'rate'), by.y = c('contrast', 'mean.phen', 'rate')
) %>%
  mutate(boot.contribution = ltre.contrib)

compare.boot.contrib %>%
  ggplot(aes(x = rate)) +
  geom_point(aes(y = boot.contribution), position = position_jitter(width = 0.25), alpha = 0.25) +
  geom_point(aes(y = contribution), size = 8, pch = 4) +
  facet_wrap(phen ~ contrast)
# Growth contributions are still off but otherwise good
# (actually growth and seed are now looking biased)

d.lambda.compare = merge(
  ltre.lambda %>%
    select(-obs.phen) %>%
    pivot_wider(names_from = trt, values_from = lambda) %>%
    mutate(d.c = drought - control, i.c = irrigated - control) %>%
    select(-c(drought, control, irrigated)) %>%
    pivot_longer(c(d.c, i.c), names_to = 'contrast', values_to = 'd.lambda') %>%
    filter(!is.na(d.lambda)),
  boot.lambda %>%
    pivot_wider(names_from = trt, values_from = lambda) %>%
    mutate(d.c = drought - control, i.c = irrigated - control) %>%
    select(-c(drought, control, irrigated)) %>%
    pivot_longer(c(d.c, i.c), names_to = 'contrast', values_to = 'd.lambda.boot') %>%
    filter(!is.na(d.lambda.boot)),
  by.x = c('contrast', 'phen'), by.y = c('contrast', 'mean.phen')
)

d.lambda.compare %>%
  ggplot(aes(x = phen, group = contrast)) +
  geom_point(
    aes(y = d.lambda.boot, colour = contrast), 
    position = position_dodge(width = 0.25), alpha = 0.25
  ) +
  geom_point(
    aes(y = d.lambda, fill = contrast),
    shape = 21, size = 4, position = position_dodge(width = 0.25)
  )
# ehhh these still look maybe biased?

d.lambda.compare %>%
  mutate(boot.error = d.lambda.boot - d.lambda) %>%
  ggplot(aes(x = boot.error, fill = boot.error > 0)) +
  geom_histogram(binwidth = 0.001) +
  facet_wrap(contrast ~ phen)
# okay getting some bias here as well.

# eh I'm not convinced this isn't biased...

# parameter variation
merge(
  boot.pert.pars, param.diffs %>% rename(truediff = pardiff), 
  by.x = c('contrast', 'rate'), by.y = c('contrast', 'param')
) %>% 
  group_by(contrast, rate) %>%
  summarise(
    se.diff.mean = sd(pardiff) / sqrt(n() - 1),
    param.diff.t = (mean(pardiff) - mean(truediff)) / se.diff.mean
  ) %>%
  ungroup() %>%
  ggplot(aes(x = rate, y = param.diff.t, colour = contrast, group = contrast)) +
  geom_point(position = position_dodge(width = 0.5)) +
  scale_colour_manual(values = c('red', 'blue'))
# not seeing evidence of bias in the bootstrapped rate differences

boot.pert.pars %>% 
  ggplot(aes(x = rate, y = pardiff, colour = contrast)) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_point(
    data = param.diffs,
    aes(x = param, y = pardiff, colour = contrast),
    position = position_dodge(width = 0.5),
    shape = 4, size = 8
  ) +
  scale_colour_manual(values = c('red', 'blue'))

# what about the sensitivities

compare.sensitivities = merge(
  mean.sensitivities, midp.boot.sensitivities,
  by.x = c('contrast', 'phen', 'param'),
  by.y = c('contrast', 'mean.phen', 'param'),
  suffixes = c('.ptes', '.boot')
)

compare.sensitivities %>%
  group_by(contrast, phen, param) %>%
  mutate(
    se.boot.mean = sd(sv.boot) / sqrt(n() - 1),
    boot.scaled = (sv.boot - mean(sv.boot)) / se.boot.mean,
    ptes.scaled = (sv.ptes - mean(sv.boot)) / se.boot.mean,
  ) %>%
  ungroup() %>%
  ggplot(aes(x = param, group = phen)) +
  # geom_point(
  #   aes(y = boot.scaled), 
  #   alpha = 0.5, position = position_dodge(width = 0.5)
  # ) +
  geom_point(
    aes(y = ptes.scaled),
    size = 4, pch = 4, position = position_dodge(width = 0.5)
  ) +
  scale_colour_manual(values = c('red', 'blue')) +
  facet_wrap(~ contrast) +
  theme(axis.text.x = element_text(angle = 90))
# not seeing any bias here...
# our boostrapped sensitivities are unbiased (at least for growth they are!)

compare.sensitivities %>%
  group_by(contrast, phen, param) %>%
  mutate(
    se.boot.mean = sd(sv.boot) / sqrt(n() - 1),
    boot.t = (sv.boot - sv.ptes) / se.boot.mean
  ) %>%
  ggplot(aes(x = param, group = phen)) +
  geom_point(
    aes(y = boot.t),
    alpha = 0.5, position = position_dodge(width = 0.5)
  ) +
  scale_colour_manual(values = c('red', 'blue')) +
  facet_wrap(~ contrast) +
  theme(axis.text.x = element_text(angle = 90))
# yeah I think the bootstrapped sensitivities are unbiased

all.ltre = merge(
  coef.ltre, boot.ltre, 
  by.x = c('contrast', 'param', 'phen'), by.y = c('contrast', 'rate', 'mean.phen'),
  suffixes = c('.ptes', '.boot')
)

head(all.ltre)

all.ltre %>%
  ggplot() +
  geom_point(aes(x = sv.ptes, y = pardiff.ptes), shape = 4, size = 4) +
  geom_point(aes(x = sv.boot, y = pardiff.boot), alpha = 0.25) +
  facet_wrap(param ~ paste(contrast, phen), scales = 'free')
