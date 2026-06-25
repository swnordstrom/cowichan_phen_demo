# cowichan_phen_demo

Data processing and analysis for manuscript "Earlier flowering explains only a small part of experimental drought’s effects on wildflower’s population growth" by Nordstrom et al.

For questions, please contact Scott Nordstrom (scottwatsonnordstrom&gmail:com).

For more information about the Williams research group and related projects, visit https://williamslabubc.weebly.com/.

### General notes

Analysis was performed entirely in `R`, in both versions 4.5.1 and 4.5.2. The following packages are used:

* `dplyr` v.1.1.4 - used extensively
* `tidyr` v.1.3.1 - used extensively
* `purrr` v.1.1.0 - used primarily in bootstrapping scripts
* `glmmTMB` v1.1.13 - used to fit vital rate models and generate predictions
* `ggplot2` v4.0.0 - used for all plots (save for Fig. 1, made in draw.io)
* `cowplot` v1.2.0 - used for various multi-paneled plots
* `ggh4x` v0.3.1 - used for double-faceting feature in Fig. 4

Scripts may include calls to load `parallel` but this package is not used in final analysis.

The IPM kernel construction and analysis can be run by running the wrapper script, `04_analysis/run_all_analysis.R`. 
Scripts are run in serial. 
The default option is to run all bootstrapping scripts with 500 intervals. 
This can be modified by going into any script with `parallel` in the title and changing the call to `n.straps` to the desired number.

Because of the time- and computational-intensity involved, the script was run on Portland State University's `coeus` server.
As written, it took approximately 1.5 days to run.
It was run with 60Gb of memory requested.

### Repository structure

The directory has the following structure:

* `00_raw_data` - raw demographic data and experimental design data
* `01_data_cleaning` - scripts for cleaning *Lomatium utriculatum* demographic and phenological data
* `02_data_exploration` - scripts for exploring data
* `03_construct_kernels` - scripts for vital rate modeling and constructing subkernels
* `04_analysis` - scripts for combining subkernels into kernels, performing analysis on kernels, and generating main figures

### Files of interest

##### `03_` scripts

* `model_selection_survival_growth.R` - model selection for growth and survival models
* `model_selection_reproduction.R` - model selection for reproductive vital rates
* `prepare_demo_data_growsurv.R` - setup growth and survival vital rate data sources and models
* `prepare_demo_data_repr.R` - setup reproductive vital rates and recruit size data sources and models
* `write_ltre_design_key.R` - writes a CSV (`ltre_treatment_key.csv` in the same subdirectory) assigning a numeric value to phenology and treatment combinations (used for reducing space in export files)
* `phenology_mean_estimates.R` - includes model selection for phenology, exports mean phenology estimates per treatment, and includes code to estimate pseudo-R^2
* `phenology_bootstrapped_estimates.R` - performs bootstrapping on phenology models and exports estimates; also includes code to generate and export Fig. S1 (annual phenology estimates in model with treatment effects varying by year)
* `deterministic_growthsurv_kernel_phen.R` - fits and exports growth + survival subkernels for each phenology+treatment combination and runs parameter perturbations
* `deterministic_reproductive_kernel_phen.R` - fits and exports reproduction subkernels for each phenology+treatment combination and runs parameter perturbations
* `deterministic_growthsurv_bootstrapping_phen.R` - performs bootstrapping of relevant vital rate parameters for growth+survival subkernels and exports bootstrapped subkernels
* `deterministic_reproductive_kernel_bootstrap_phen.R` - performs bootstrapping of relevant vital rate parameters for reproductive subkernels and exports bootstrapped subkernels
* `find_optimal_p_germ.R` - numerically solving for the seedling establishment rate used in analysis (*note: requires inputs from the `deterministic_<>_kernel_phen.R` scripts)
* subdirectory `analyze_ffd` - model fitting with first flowering day rather than mean flowering day
* `reproductive_vital_rates_figures_plots.R` - generates Fig. S3 (reproductive components as a function of size and treatment)

###### `04_` scripts

Main analysis scripts:

* `deterministic_ltre.R` - generates kernels to fit two-way LTRE; exports LTRE contributions in CSV form to create Fig. 5 (LTRE figure) and Table S19 (accuracy of LTRE and mirrored LTRE)
* `lambda_over_phenology.R` - generates kernels to estimate lambda over weekly intervals; exports CSVs with lambda estimates to create Fig. 4 (estimates of lambda estimated at different mean flowering dates)

`lambda_over_phenology.R` and `deterministic_ltre.R` export CSVs that are imported and processed in the script `figures_from_csvs.R` to export Figs. 4 and 5 as well as several supporting figures.

Additional scripts:

* `assess_evictions.R` - quantifying the eviction rate for IPMs
* `ltre_sensitivity.R` - estimates Delta lambda and LTRE effects at different establishment rates and exports Figs. S4 (Delta lambda between treatments at different establishment rates) and S5 (LTRE contributions at different establishment rates)
* `reproduction_growth_phen_trt_figure.R` - combines growth and reproduction data to export Fig. 3 (phenology effects on reproduction, phenology effects on growth, treatment effects on growth)
* `phenology_figure.R` - creates and exports Fig. 2 (treatment effects on phenology)
* `plot_growth_differences.R` - creates and exports Fig. S2 (treatment effects on growth rates)

###### `04_analysis/run_all_analysis.R`

This script is a wrapper to run the entire analysis from phenology estimates and subkernel construction through to the analysis of the IPMs. 
Be warned that it takes over a day to run. I suggest running it remotely.

Here is the body of the script:

```
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
```

The individual scripts sourced in this file can also be run or examined on their own (provided that the input files exist).

*Note*: Fig. 1 (conceptual diagram) was not created with any materials in this repository. It was created with [draw.io](https://www.drawio.com/).

### Changelog

* 2026-03-26: updated the `.gitignore` and updated the README
* 2026-06-25: updated repository and README for resubmission
  + Add phenology ~ treatment figure to main text (`03_construct_kernels/raw_phen_fig.R` is now `04_analysis/phenology_figure.R` and exports an expanded plot, Fig. 2)
  + Correspondingly, Figs. 2-4 in main text are now Figs. 3-5
  + Fig. 5 (in resubmission; this is the LTRE figure) is now two combined panels and the corresponding script does not plot the mirrored or double-aggregated LTRE components.
  + Removal of some SI figures means that SI figure labels have changed
  + pushed some missing files (e.g., growth differences script)

