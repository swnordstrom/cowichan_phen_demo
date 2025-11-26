# Make directory for storing outputs if it's not already created
if (!dir.exists('03_construct_kernels/out')) {
  dir.create('03_construct_kernels/out')
}
if (!dir.exists('03_construct_kernels/bootstrapped_model_coefs/')) {
  dir.create('03_construct_kernels/bootstrapped_model_coefs/')
}
if (!dir.exists('04_analysis/figures')) {
  dir.create('04_analysis/figures')
}
if (!dir.exists('04_analysis/out')) {
  dir.create('04_analysis/out')
}

# Get rid of this super annoying feature
options(
  dplyr.summarise.inform = FALSE
)

### Scripts to get subkernels and related values

# Run + export phenology treatment estimates 
source('03_construct_kernels/phenology_mean_estimates.R')
rm(list = ls())

# Run + export growth+survival subkernel (including perturbations)
source('03_construct_kernels/deterministic_growthsurv_kernel_phen.R')
rm(list = ls())

# Run + export reproduction kernel components (including perturbations)
source('03_construct_kernels/deterministic_reproductive_kernel_phen.R')
rm(list = ls())

# Run phenology bootstrapping
source('03_construct_kernels/phenology_bootstrapped_estimates.R')
rm(list = ls())

# Run growth+survival bootstrapping (including perturbations)
source('03_construct_kernels/deterministic_growthsurv_bootstrapping_phen.R')
rm(list = ls())

# Run reproduction bootstrapping (including perturbations)
source('03_construct_kernels/deterministic_reproductive_boostrap_phen.R')
rm(list = ls())

### Scripts to construct kernels and run analysis

# Figure 2 (Lambda across growing season)
source('04_analysis/lambda_over_phenology.R')
rm(list = ls())

# Figure 3 (LTRE)
source('04_analysis/deterministic_ltre.R')
