# Visualizations of layouts of plants over time
# Has code for scraping out coordinates, converting them to numeric form, and
# plotting x vs. y coordinates
# So far just has flowering plants in plots
# SN - init 23 jan 2024

library(ggplot2)
library(dplyr)
library(tidyr)

# Read in demo data and do some processing
demo = read.csv('01_data_cleaning/out/demo_imputed_survival.csv') %>%
  # Give me flowering plants only
  filter(!is.na(No.umbels), No.umbels > 0) %>%
  # Split out plantid into component parts (plot column is redundant)
  separate(plantid, into = c('tag', 'plot', 'coord'), sep = '_', remove = FALSE) %>%
  select(-plot) %>%
  # Add in treatment
  merge(read.csv('00_raw_data/plot_treatments.csv') %>% rename(plot = Plot))
  
head(demo)
nrow(demo)

# Give me only the plants where the Y coord is listed alphabetically
demo %>%
  filter(grepl('[A-Z]{1}', coord)) %>% nrow()
# we lose ~20 records here but that's okay

demo = demo %>% filter(grepl('[A-Z]$', coord))

# Okay... now actually separating out x and y coordinate
# Oh actually this is easy lol using gsub

demo = demo %>%
  mutate(
    xcoor = as.numeric(gsub('[A-Z]+$', '', coord)),
    ylett = gsub('^[0-9]{1,2}', '', coord)
  )

unique(demo$xcoor) # cool... kinda surprised by that zero though...
unique(demo$ylett) # looks alright with me

# Convert y coord into numeric coordinate
demo = merge(
  x = demo,
  y = data.frame(ylett = LETTERS, ycoor = 1:length(LETTERS)),
  all.x = TRUE
) %>%
  select(-ylett, coord)

# Cool! Let's try stuff now...

demo %>%
  filter(Plot %in% 1) %>%
  ggplot(aes(x = xcoor, y = ycoor)) +
  geom_point(position = position_jitter(height = 0.25, width = 0.25)) +
  facet_wrap(~ Year)
# wonder if I could do all plots like this...
length(unique(demo$Plot))
# 13 x 8 plot... hmm...

demo %>%
  ggplot(aes(x = xcoor, y = ycoor, colour = trt)) +
  geom_point(position = position_jitter(height = 0.125, width = 0.125)) +
  scale_colour_manual(values = c('black', 'red', 'blue')) +
  facet_grid(rows = vars(Plot), cols = vars(Year))
# less insightful than I'd like
# (should 2020 even be included here...?)

#