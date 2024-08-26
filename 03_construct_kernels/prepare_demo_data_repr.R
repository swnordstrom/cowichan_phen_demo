##### Wrapper script to prepare demographic data for reproduction models


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

# # Need to add in rows for umbels that died before the seed counting
# # This is necessary for estimating the probability of an umbel producing zero seeds
# demo.seed = rbind(
#   demo.seed,
#   demo.seed %>%
#     group_by(Year, plantid) %>%
#     mutate(miss.umbel = ifelse(phen.umbels < n(), 0, phen.umbels - n())) %>%
#     ungroup() %>%
#     distinct(Year, plantid, .keep_all = TRUE) %>%
#     uncount(miss.umbel) %>%
#     mutate(no.seeds = 0)
# ) %>%
#   mutate(Year = factor(Year))

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

# Umbel-level budding phenology data
phen = all.data %>% 
  filter(in.phen) %>%
  # Split out the bud dates for bud date models; the most umbels seen in a
  # plant is 12, so use separate() to kick these out and then pivot_long to get
  # one row per umbel
  # (first - need to get one row per plant - do a distinct())
  distinct(Year, plantid, .keep_all = TRUE) %>%
  separate_wider_delim(phen.julis, names = paste0('uu', 1:12), delim = ';', too_few = 'align_start') %>%
  pivot_longer(starts_with('uu'), names_to = 'umbel.number', values_to = 'phen.julian') %>%
  filter(!is.na(phen.julian)) %>%
  mutate(phen.julian = as.numeric(gsub('\\s', '', phen.julian))) %>%
  mutate(Year = factor(Year))

# All plants in phen, seed, and demo
seed = all.data %>% 
  filter(in.phen, in.seed, in.demo) %>%
  # Need size data
  filter(!is.na(No.leaves), !is.na(Leaf.length), No.leaves > 0, Leaf.length > 0) %>%
  # Hmm... okay, separae() and across() doesn't work, so I guess I can merge
  # this in with the data frame above, summarised by mean phen date per plant
  merge(
    y = phen %>% 
      group_by(plantid, Year) %>% 
      summarise(
        mean.phen = mean(phen.julian, na.rm = TRUE), 
        sept.phen = length(unique(phen.julian))
      ) %>%
      ungroup(),
    all.x = TRUE, all.y = FALSE
  )

seed = rbind(
  seed,
  seed %>%
    group_by(Year = factor(Year), plantid) %>%
    mutate(miss.umbel = ifelse(phen.umbels < n(), 0, phen.umbels - n())) %>%
    ungroup() %>%
    distinct(Year, plantid, .keep_all = TRUE) %>%
    uncount(miss.umbel) %>%
    mutate(no.seeds = 0)
) %>%
  mutate(phen.c = mean.phen - round(mean(mean.phen))) %>%
  mutate(size = log(No.leaves * Leaf.length))
