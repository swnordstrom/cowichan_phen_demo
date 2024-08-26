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
fr.obsv = read.csv('03_construct_kernels/out/determinstic_reprod_kernel_phen_ltre.csv')


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
fr.obsv.pert = read.csv('03_construct_kernels/out/deterministic_repr_coef_perturbation_phen.csv')


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

fr.pert.pars = read.csv('03_construct_kernels/out/reprod_bootstrapped_perturbed_params.csv') 
# Do this in two steps because this file also contains the phen dates used in the bootstarp

# Dates used in phen bootstrapping
phen.boots = fr.pert.pars %>%
  select(c(boot, contains('phen'))) %>%
  pivot_longer(-boot, names_to = 'trt.phen', values_to = 'phen') %>%
  mutate(
    trt.phen = gsub('phen\\_', '', trt.phen),
    phen = as.Date(as.numeric(phen), format = '%b-%d')
  )

# Reproductive bootstrapped parameters
fr.pert.pars = fr.pert.pars %>%
  select(-contains('phen')) %>%
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
    select(size.prev, size.nex, trt, trt.phen, p.size.cur.g, p.size.cur.f, param, orig.parval),
  merge(
    gs.obsv.pert, 
    fr.obsv %>% rename(mean.phen = phen),
    by = c('size.prev', 'size.nex', 'trt'),
    suffixes = c('.g', '.f')
  ) %>%
    select(size.prev, size.nex, trt, trt.phen, p.size.cur.g, p.size.cur.f, param, orig.parval)
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
    select(size.prev, size.cur, trt, trt.phen, samp, param = perturb.param, p.size.cur.g, p.size.cur.f),
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
    select(size.prev, size.cur, trt, trt.phen, samp, param, p.size.cur.g, p.size.cur.f)
) %>%
  mutate(p.size.cur = p.size.cur.g + p.size.cur.f * p.germ) %>%
  select(-c(p.size.cur.g, p.size.cur.f))

# Remove input files
rm(fr.obsv, gs.obsv, fr.boot, gs.boot, fr.boot.pert, gs.boot.pert)

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
  pivot_wider(names_from = trt, values_from = p.size.cur) %>%
  mutate(d.c = (drought + control) / 2, i.c = (irrigated + control) / 2) %>%
  select(-c(drought, control, irrigated)) %>%
  pivot_longer(c(d.c, i.c), names_to = 'contrast', values_to = 'p.size.cur') %>%
  filter(!is.na(p.size.cur))

# Get midpoints for PHENOLOGY
# Here - to estimate the sensitivity of lambda to phenology, getting a midpoint
# kernel between phen levels

midp.phen.kernel.df = obsv.kernel.df %>%
  select(-phen) %>%
  # NOTE: pivoting out by trt.phen instead of trt (because we're averaging
  # across different phens)
  pivot_wider(names_from = trt.phen, values_from = p.size.cur) %>%
  mutate(d.c = (drought + control) / 2, i.c = (irrigated + control) / 2) %>%
  select(-c(control, drought, irrigated)) %>%
  pivot_longer(c(d.c, i.c), names_to = 'contrast.phen', values_to = 'mentry', values_drop_na = TRUE)

midp.phen.pert.kernel.df = obsv.pert.kernel.df %>%
  select(-orig.parval) %>%
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
  obsv.kernel.df, obsv.kernel.df[,c("trt", "trt.phen", "phen")], 
  sep = '_', drop = TRUE
) %>%
  # split() splits the kernel df into a list where each entry is a data frame
  # for one phen-treatment kernel
  lapply(
    function(df) {
      df %>%
        arrange(size.prev, size.cur) %>%
        pivot_wider(names_from = size.prev, values_from = p.size.cur) %>%
        select(-c(trt, trt.phen, phen, size.cur)) %>%
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
  boot.kernel.df,
  boot.kernel.df[,c("trt", "trt.phen", "samp")],
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
  midp.obsv.kernel.df,
  midp.obsv.kernel.df[,c("contrast", 'trt.phen', "phen")],
  sep = '_', drop = TRUE
) %>%
  lapply(
    function(df) {
      df %>%
        arrange(size.prev, size.cur) %>%
        pivot_wider(names_from = size.prev, values_from = mentry) %>%
        select(-c(contrast, trt.phen, phen, size.cur)) %>%
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
  midp.pert.kernel.df[,c("contrast", 'trt.phen', 'param')],
  sep = '_', drop = TRUE
) %>%
  lapply(
    function(df) {
      df %>%
        arrange(size.prev, size.nex) %>%
        pivot_wider(names_from = size.prev, values_from = mentry) %>%
        select(-c(contrast, trt.phen, param, size.nex)) %>%
        as.matrix()
    }
  ) %>%
  sapply(function(x) Re(eigen(x)$values[1])) %>%
  data.frame(lambda = .) %>%
  mutate(c_pg = row.names(.)) %>%
  separate(c_pg, into = c('contrast', 'trt.phen', 'param'), sep = '_')  

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
  midp.boot.pert.kernel.df[,c("contrast", "trt.phen", "param", "samp")],
  sep = '_', drop = TRUE
) %>%
  lapply(
    function(df) {
      df %>%
        arrange(size.prev, size.cur) %>%
        pivot_wider(names_from = size.prev, values_from = p.size.cur) %>%
        select(-c(contrast, samp, param, trt.phen, size.cur)) %>%
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
  midp.phen.boot.pert.kernel.df[,c("trt", "param", "contrast.phen", "samp")],
  sep = '_', drop = TRUE
) %>%
  lapply(
    function(df) {
      df %>%
        arrange(size.prev, size.cur) %>%
        pivot_wider(names_from = size.prev, values_from = mentry) %>%
        select(-c(trt, contrast.phen, param, samp, size.cur)) %>%
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
  by = c('contrast', 'trt.phen'), suffixes = c('.orig', '.pert')
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

# Initution check: sensitivites to an increase in phenology should all be negative
# (because lambda is always declining with bud date, in our models)

midp.phen.sens = merge(
   midp.phen.lambda, midp.phen.pert.lambda,
   by = c('trt', 'contrast.phen'), suffixes = c('.orig', '.pert')
) %>%
  mutate(sv = ((lambda.pert - lambda.orig) / 0.0001))

  # # Need a -1 in here for when the contrast in phenology is positive or negative
  # mutate(
  #   sv = ((lambda.pert - lambda.orig) / 0.0001) * ifelse(grepl('^d', trt) | grepl('^d', contrast.phen), -1, 1)
  # )

midp.phen.boot.sens = merge(
  midp.phen.boot.lambda, midp.phen.boot.pert.lambda,
  by = c('trt', 'contrast.phen', 'samp'), suffixes = c('.orig', '.pert')
) %>%
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

obsv.param.diffs = rbind(gs.obsv.pert, fr.obsv.pert %>% select(-c(trt.phen, mean.phen))) %>%
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
  pivot_longer(c(d.c, i.c), names_to = 'contrast', values_to = 'pardiff', values_drop_na = TRUE) %>%
  rename(
    samp = boot,
    param = rate
  )

# Okay... here we will need to introduce some negative ones

obsv.phen.diffs = obsv.kernel.df %>% 
  distinct(trt, trt.phen, phen) %>%
  pivot_wider(names_from = trt.phen, values_from = phen) %>%
  mutate(
    # Need sign corrections here:
    # (d/i - control) gives the difference *from the perspective of the control*
    # but in the drought/irrig treatments, it should be control - d/i
    d.c = (drought - control) * ifelse(trt %in% 'control', 1, -1), 
    i.c = (irrigated - control) * ifelse(trt %in% 'control', 1, -1)
  ) %>%
  select(-c(control, drought, irrigated)) %>%
  pivot_longer(c(d.c, i.c), names_to = 'contrast.phen', values_to = 'phen.diff', values_drop_na = TRUE)

boot.phen.diffs = phen.boots %>%
  rename(samp = boot) %>%
  mutate(phen = as.numeric(phen)) %>%
  pivot_wider(names_from = trt.phen, values_from = phen) %>%
  mutate(d.c = drought - control, i.c = irrigated - control) %>%
  select(-c(control, drought, irrigated)) %>%
  pivot_longer(c(d.c, i.c), names_to = 'contrast.phen', values_to = 'phen.diff')

# The code above gives only the differences *from the control perspective*
# But we will also want to get differences from the drought/irrigated perspective
# which will require negating the phen diff
boot.phen.diffs = rbind(
    boot.phen.diffs %>% mutate(trt = 'control'),
    boot.phen.diffs %>% 
      mutate(trt = ifelse(contrast.phen %in% 'd.c', 'drought', 'irrigated')) %>% 
      mutate(phen.diff = -1 * phen.diff)
  )

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
boot.ltre = merge(midp.boot.sens, boot.param.diffs) %>%
  mutate(contrib = pardiff * sv)

# Phenology effects
phen.ltre = merge(midp.phen.sens, obsv.phen.diffs) %>%
  mutate(contrib = phen.diff * sv)

phen.boot.ltre = merge(midp.phen.boot.sens, boot.phen.diffs) %>%
  mutate(contrib = phen.diff * sv)

# ------------------------------------------------------                  
# ------ Check that lambda differences match LTRE sums -
# ------------------------------------------------------

ltre.dlambda.compare = merge(
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

ltre.dlambda.compare %>%
  mutate(relerr = (csum - d.lambda) / d.lambda)
# Slightly negatively biased, but all by <1% of the true lambda difference

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
  group_by(trt, contrast.phen, rate) %>%
  summarise(contrib = sum(contrib)) %>%
  ungroup()

boot.phen.ltre = phen.boot.ltre %>%
  mutate(rate = gsub('phen\\.', '', rate)) %>%
  group_by(trt, contrast.phen, samp, rate) %>%
  summarise(contrib = sum(contrib)) %>%
  ungroup()

# ------------------------------------------------------                  
# ------ Crude plots -----------------------------------
# ------------------------------------------------------

obsv.trt.ltre %>%
  # I want to do this on the control buddate
  filter(trt.phen %in% 'control') %>%
  ggplot(aes(x = rate, y = contrib, fill = contrast)) +
  geom_col(position = 'dodge')

# Picture here: more growth in treatments, less flowering
# differing treatment effects on seed production's influence
# differing treatment effects on recruit size (other dir.)

boot.trt.ltre %>%
  # I want to do this on the control buddate
  filter(trt.phen %in% 'control') %>%
  ggplot(aes(x = rate, y = contrib, colour = contrast)) +
  geom_point(position = position_dodge(width = 0.25), alpha = 0.5)

# This is backwards...

obsv.phen.ltre %>%
  # think about which trt we want...
  filter(trt %in% 'control') %>%
  ggplot(aes(x = rate, y = contrib, fill = contrast.phen)) +
  geom_col(position = 'dodge')

boot.phen.ltre %>%
  # We'll do the drought/irrigated differences for these
  filter(trt %in% 'control') %>%
  ggplot(aes(x = rate, y = contrib, colour = contrast.phen)) +
  geom_point(position = position_dodge(width = 0.25), alpha = 0.5)

# Combining...

control.ltre.all = rbind(
  # --- Observed treatment effects
  obsv.trt.ltre %>%
    # give me LTRE values for the control dates and remove column
    filter(trt.phen %in% 'control') %>%
    select(-trt.phen) %>%
    # marker for type of observation
    mutate(varb = 'alpha', samp = 'obsv', type = 'trt'),
  # --- Bootstrapped treatment effects
  boot.trt.ltre %>%
    # give me LTRE values for the control dates and remove unneeded columns
    filter(trt.phen %in% 'control') %>%
    select(-c(trt.phen, samp)) %>%
    # marker for type of observation
    mutate(varb = 'alpha', samp = 'boot', type = 'trt'),
  # --- Observed phenology effects (within treatment)
  obsv.phen.ltre %>%
    # give me LTRE values where the reference date is the control
    # and remove unnecessary column
    filter(trt %in% 'control') %>%
    select(-trt) %>%
    # Rename column for column agreement
    rename(contrast = contrast.phen) %>%
    mutate(varb = 'beta', samp = 'obsv', type = 'phen'),
  # --- Bootstrapped phenology effects
  boot.phen.ltre %>%
    # give me LTRE values where the reference date is the control
    # and remove unnecessary columns
    filter(trt %in% 'control') %>%
    select(-c(trt, samp)) %>%
    # Rename column for column agreement
    rename(contrast = contrast.phen) %>%
    mutate(varb = 'beta', samp = 'boot', type = 'phen')
) %>%
  mutate(ltre.varb = paste0(varb, '[', rate, ']'))

control.ltre.all %>%
  mutate(ltre.varb = paste0(varb, '[', rate, ']')) %>%
  ggplot(aes(x = ltre.varb, y = contrib)) +
  geom_point(aes(shape = samp, size = samp, colour = type)) +
  scale_shape_manual(values = c(1, 19)) +
  scale_size_manual(values = c(1, 4)) +
  scale_colour_manual(values = c('gray11', 'gray66')) +
  scale_x_discrete(labels = scales::label_parse()) +
  facet_wrap(~ contrast, nrow = 2)

# ugly.

control.ltre.summ = merge(
  control.ltre.all %>% filter(samp %in% 'obsv') %>% select(-c(varb, rate, samp)),
  control.ltre.all %>%
    filter(samp %in% 'boot') %>%
    group_by(contrast, ltre.varb) %>%
    reframe(
      cilim = quantile(contrib, probs = c(0.025, 0.975)),
      lohi = c('lo', 'hi')
    ) %>%
    pivot_wider(names_from = lohi, values_from = cilim)
) %>%
  ungroup()

control.ltre.summ %>%
  mutate(
    contr.pretty = paste(ifelse(contrast %in% 'd.c', 'drought', 'irrigated'), 'vs. control')
  ) %>%
  ggplot(aes(x = ltre.varb)) +
  geom_col(aes(y = contrib, fill = contrast)) +
  geom_segment(aes(xend = ltre.varb, y = lo, yend = hi)) +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_fill_manual(values = c('red', 'blue')) +
  facet_wrap(~ contr.pretty) +
  guides(fill = 'none')

# Okay... something funny is going on with the growth
# Asymmetric, contributions perhaps are non-normal...
# (plot above also kind of suggested these bootstraps were not normally distributed...)

control.ltre.all %>% 
  filter(samp %in% 'boot') %>% 
  group_by(ltre.varb, contrast) %>%
  mutate(std.contrib = (contrib - mean(contrib)) / sd(contrib)) %>%
  ggplot(aes(x = std.contrib, group = interaction(contrast, ltre.varb), colour = contrast)) + 
  geom_density(aes(colour = contrast))
# ah... lack of normality looks to be common
# may be a result of small sample size

