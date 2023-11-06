##### Script for playing around with demography data for survival estimates
##### SN - init 6 nov 2023

##############################################
##### Setup/preparation
##############################################

library(ggplot2)
library(dplyr)
library(tidyr)

rm(list = ls())

### Read in demo data

demo = merge(
  # Demo data
  x = read.csv('01_data_cleaning/out/demo_table.csv'),
  # Treatment data
  y = read.csv('00_raw_data/plot_treatments.csv') %>% rename(Plot = plot)
)

head(demo)
tail(demo)
nrow(demo)

### Add in sampling data:

##### Remove plants from coordinates that were not sampled (acc'd to metadata)

demo = demo %>%
  filter(!(Year %in% 2020 & Plot %in% c(1:2, 13:15) & grepl('_\\d{2}\\D$', plantid))) %>%
  filter(!(Year %in% 2020 & Plot %in% 3 & (grepl('_\\d{2}\\D$', plantid) | grepl('\\_[789]\\D', plantid))))

# (I checked... there are five plants with non-empty records that get removed in the above loop)
# (might be an issue in coordinate entry)
# (I think it's a good idea to remove these to be safe)

##############################################
##### Assessing gaps from non-detection or mis-specification
##############################################

# Problem here is that we have plants that can be "dead" (leaf count = 0) but
# re-appear as alive later in the dataset. This could be detection probabilities
# varying (so living plants are missed), dormancy, or the plant could truly be
# dead. Prior assessment of the data suggests that 2/3 or so of the plants that
# were marked as "dead" at some point were observed later. Here I'll try to
# assess that a little more closely.

### How many plants marked as "not detected"
demo %>% filter(!detected) %>% nrow()
# good number
demo %>% 
  group_by(Plot, Year, detected = paste0(ifelse(detected, '', 'not'), 'detected')) %>% 
  summarise(n = n()) %>%
  pivot_wider(names_from = detected, values_from = n)
# looks like a non-zero number of non-detections in each year...

demo %>% filter(!detected & Plot %in% 1)
# ugh...

### Look at deaths by year
demo %>%
  filter(detected) %>%
  group_by(Plot, Year, obs.alive = ifelse(obs.alive, 'living', 'dead.q')) %>%
  summarise(n = n()) %>%
  pivot_wider(names_from = obs.alive, values_from = n)
# oh there are some NAs...

demo %>% filter(is.na(obs.alive)) # oh yeah... what's up with these
# 2023 records: all were alive with flowers eaten off (call these alive
# 2016 records: did not check... but will assume the same thing for this prelim analysis

demo = demo %>%
  mutate(obs.alive = ifelse(is.na(obs.alive), TRUE, obs.alive))

demo %>%
  filter(detected) %>%
  group_by(Plot, Year, obs.alive = ifelse(obs.alive, 'living', 'dead')) %>%
  summarise(n = n()) %>%
  pivot_wider(names_from = obs.alive, values_from = n, values_fill = 0) # %>% View()

# Let's just plot, over time, how many dead plants are observed per plot...

demo %>%
  filter(detected) %>%
  group_by(Plot, trt, Year) %>%
  summarise(n.dead = sum(!obs.alive)) %>%
  ggplot(aes(x = Year, y = n.dead, group = Plot, colour = trt)) + 
  geom_line() +
  scale_colour_manual(values = c('black', 'red', 'blue'))
# okay this is cool, possibly
# spike in 2019, it seems like
# another "spike" in 2022, although this could be due to re-checks

# Dead data frame

demo.dead = demo %>%
  group_by(plantid) %>%
  filter(any(!obs.alive & detected))

head(demo.dead)

# Do a case when?

demo.dead = demo.dead %>%
  # Get rid of observations before first observed year of death...
  filter(Year >= min(Year[!obs.alive & detected])) %>%
  mutate(
    dead.status = case_when(
      n() < 2 ~ 'seen.once',
      any(obs.alive & detected) ~ 'obs.alive',
      sum(!detected) == (n() - 1) ~ 'not.detected',
      .default = 'maybe.actually.dead'
    )
  ) %>%
  arrange(plantid, Year)

# No-detects?
demo.dead %>% distinct(plantid, dead.status) %>% group_by(dead.status) %>% summarise(n = n())
# hmm... good number of cases of each...

demo.dead %>% filter(dead.status %in% 'seen.once')
demo.dead %>% filter(dead.status %in% 'seen.once') %>% group_by(Year) %>% summarise(n = n())
# ugh... some 2022 plants in here that would have been great to confirm in 2023

demo.dead %>% filter(dead.status %in% 'not.detected')
# oh actually I suppose these could all just be fail-to-finds...

# the observed alive plants... hmm... bummer, but it happens
demo.dead %>% filter(dead.status %in% 'obs.alive')
# there are actually some huge gaps in here...
# can we get the distribution of gaps?

demo.dead.gaps = demo.dead %>%
  filter(dead.status %in% 'obs.alive') %>% 
  # should still be grouped by plantid
  # gonna get conservative here and chop out the years after the first time detected alive
  filter(Year <= min(Year[obs.alive & detected])) %>%
  summarise(n.not.obs = sum(!obs.alive))

demo.dead.gaps %>% group_by(n.not.obs) %>% summarise(n = n())
# mostly one year, some two, a handful more than that!

demo.dead.newrecs = demo.dead %>%
  filter(dead.status %in% 'obs.alive') %>%
  filter(Year <= min(Year[obs.alive & detected])) %>%
  mutate(n.not.obs = sum(!obs.alive)) %>%
  arrange(plantid, desc(Year)) %>%
  distinct(plantid, .keep_all = TRUE)

# Distribution of years (might suggest variation in sampling effort)
table(demo.dead.newrecs$Year) 
# most in 2021, surely mostly due to incomplete sampling 2020
# but still a good number in other years

# Gap size vs. year of re-detection
with(demo.dead.newrecs, table(Year, n.not.obs))
# lots of plants reappearing in 2022, including after long gaps

# Size at re-detection
table(demo.dead.newrecs$no.leaves) # mostly one leaf... regression is possible though
hist(demo.dead.newrecs$leaf.leng) # lol nice symmetry!

# Size at re-detection vs. gap size
demo.dead.newrecs %>%
  ggplot(aes(x = no.leaves,y = leaf.leng, fill = n.not.obs)) +
  geom_point(shape = 21, size = 3, position = position_jitter(width = 0.25)) +
  scale_fill_viridis_c()

demo.dead.newrecs %>%
  ggplot(aes(x = n.not.obs,y = leaf.leng, fill = no.leaves)) +
  geom_point(shape = 21, size = 3, position = position_jitter(width = 0.25)) +
  scale_fill_viridis_c()
# very tempted to say plants not observed after a certain date are new recruits...
# (of course - this gap depends on sampling!)

demo.dead.newrecs %>%
  mutate(logsize = log(no.leaves * leaf.leng)) %>%
  ggplot(aes(x = n.not.obs, y = logsize, fill = trt)) +
  geom_point(shape = 21, size = 3, position = position_jitter(width = 0.25))

# Reset this data frame for another plot
demo.dead.newrecs.range = demo.dead %>%
  filter(dead.status %in% 'obs.alive') %>%
  group_by(Plot, trt, plantid) %>%
  summarise(y0 = min(Year), y1 = min(Year[obs.alive & detected]))

head(demo.dead.newrecs.range)  

ggplot(demo.dead.newrecs.range %>% mutate(gapsize = y1 - y0), aes(colour = trt)) +
  geom_point(aes(x = y0, y = plantid), shape = 21) +
  geom_point(aes(x = y1, y = plantid)) +
  geom_segment(aes(x = y0, xend = y1, y = plantid, yend = plantid, colour = trt)) +
  scale_colour_manual(values = c('black', 'red', 'blue'), 'treatment') +
  facet_grid(cols = vars(gapsize)) +
  labs(x = 'Year', y = 'Plant') +
  theme(
    panel.grid = element_blank(), 
    panel.background = element_blank(),
    axis.text.y = element_blank()
    # axis.text.x = element_text(angle = 90)
  )

# A lot of the two-year gaps have 2020 as one of the missing years...
# (we know sampling in 2020 was spotty)

demo.dead.newrecs.range %>%
  filter(y0 < 2020 & y1 > 2020) %>%
  mutate(gap.size = y1 - y0) %>%
  group_by(Plot, trt, gap.size = paste0('gap', gap.size)) %>%
  summarise(n = n()) %>%
  pivot_wider(names_from = gap.size, values_from = n, values_fill = 0)

# looks like a lot of issues in plot 5, 15 in particular
# 14 to a lesser extent
# otherwise problems are sparse/rare...

