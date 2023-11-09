##### Script for playing around with demography data for survival estimates
##### SN - init 6 nov 2023

##############################################
##### Setup/preparation
##############################################

library(ggplot2)
library(dplyr)
library(tidyr)
library(lme4)

rm(list = ls())

### Inverse logit function
ilogit = function(x) (1+exp(-x))^-1

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

# ##### Remove plants from coordinates that were not sampled (acc'd to metadata)
# 
# demo = demo %>%
#   filter(!(Year %in% 2020 & Plot %in% c(1:2, 5, 13:15) & grepl('_\\d{2}\\D$', plantid))) %>%
#   filter(!(Year %in% 2020 & Plot %in% 3 & (grepl('_\\d{2}\\D$', plantid) | grepl('\\_[789]\\D', plantid))))

# (I checked... there are five plants with non-empty records that get removed in the above loop)
# (might be an issue in coordinate entry)
# (I think it's a good idea to remove these to be safe)

# (actually I decided to leave these records in because they aren't actually doing anything...)

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

demo %>%
  filter(Year %in% 2020, Plot %in% c(5, 14, 15)) %>%
  arrange(Plot, plantid)

# How many of these are simply "missed"?

demo.dead.newrecs.range %>%
  filter(y0 < 2020 & y1 > 2020) %>%
  # Of these, how many were marked "dead" in 2019?
  group_by(is.2019 = y0 %in% 2019) %>%
  summarise(n = n()) 
# many of them!

demo.dead.newrecs.range %>%
  filter(y0 < 2020 & y1 > 2020) %>%
  # Of these, how many were marked "dead" in 2019?
  group_by(is.2019 = y0 %in% 2019) %>%
  summarise(n = n())

demo.dead.newrecs.range %>%
  filter(y0 < 2020 & y1 > 2020) %>%
  # Plants last seen in 2019
  filter(y0 %in% 2019) %>%
  # how many of these have coordinates that are >9?
  separate(plantid, into = c('tag', 'plot', 'coord'), sep = '_') %>%
  mutate(past.nine = grepl('\\d{2}', coord)) %>%
  group_by(plot, past.nine) %>%
  summarise(n = n()) %>%
  pivot_wider(names_from = past.nine, values_from = n, values_fill = 0)
# nearly all of the plants in plots 5, 15, 14 are plants in areas not measured
# plots 4, 6, 7, 12 should have been completely sampled

# Let's look at plants that were missed in 2019 *and* 2021
# (I'm tempted to say these ones should be separate plants...)
demo.dead.newrecs.range %>%
  filter(y0 < 2020 & y1 > 2021) %>%
  # Plants last seen in 2019
  filter(y0 %in% 2019) %>%
  # how many of these have coordinates that are >9?
  separate(plantid, into = c('tag', 'plot', 'coord'), sep = '_') %>%
  mutate(past.nine = grepl('\\d{2}', coord)) %>%
  group_by(plot, past.nine) %>%
  summarise(n = n()) %>%
  pivot_wider(names_from = past.nine, values_from = n, values_fill = 0)
# but I guess the question is were these plants seen in 2021? or skipped/missed

demo %>%
  filter(Year %in% 2021) %>%
  filter(plantid %in% (demo.dead.newrecs.range %>%
                         filter(y0 < 2020 & y1 > 2021) %>%
                         # Plants last seen in 2019
                         filter(y0 %in% 2019) %>%
                         pull(plantid)))
# okay it looks like they mostly were checked in 2021... let's say they died]
# (but then wouldn't we also have to assume that plants with no subsequent
# "detected = TRUE" records were never searched?

demo.dead %>% filter(dead.status %in% 'not.detected') # %>% View()
# ugh... suppose it's possible they were entered in mistake once and never re-checked
# but also probably fair to assume they were just not seen
# (maybe check to see if there were other plants entered at some coords... argh maybe later)


##############################################
##### Taking a shot at "fixing" the gaps
##############################################

### Rules:
# - if gap is only one year (i.e., found alive year after being dead), impute alive
# - if gap is two years *and* coordinate was not surveyed in 2020, impute alive
# - if gap is two years and coordinate *was* surveyed in 2020, assign new plant
# - if gap is three years or longer, assign new plant

# Think a case_when could be good for handling this

demo.dead.final = rbind(
  demo.dead %>%
    ungroup() %>%
    filter(!dead.status %in% 'obs.alive') %>%
    distinct(plantid, dead.status),
  demo.dead %>%
    filter(dead.status %in% 'obs.alive') %>%
    summarise(y0 = min(Year), y1 = min(Year[obs.alive])) %>%
    mutate(dead.status = case_when(
      y1 == (y0 + 1) ~ 'gap1',
      ((y1 == 2021) & (y0 == 2019)) &
        # catching plants in plots 1, 2, 5, 13, 14, 15 with coordinate > 9
        (((grepl('\\_[125]\\_', plantid) | grepl('\\_1[345]\\_', plantid)) &
            grepl('\\_\\d{2}\\D$', plantid)) |
           # *or* catching plants in plot 3 with coordinate > 6
           ((grepl('\\_3\\_', plantid)) & 
              (grepl('\\_[789]\\D$', plantid) | grepl('\\_\\d{2}\\D$', plantid)))) ~ 'gap.2020',
      .default = paste0('gap', (y1 - y0))
    )
  ) %>%
    select(-c(y0, y1))
)

demo.dead.final %>% group_by(dead.status) %>% summarise(n = n())

demo.imputed = merge(
  x = demo,
  y = demo.dead.final,
  all.x = TRUE, all.y = TRUE
) %>%
  # any plants that were never seen dead will be NA after this merge - 
  # add a status here
  mutate(dead.status = ifelse(is.na(dead.status), 'not.dead', dead.status))

nrow(demo.imputed)
demo.imputed %>% 
  distinct(plantid, dead.status) %>% 
  group_by(dead.status) %>% 
  summarise(n = n())
# looks okay to me

#### Run the code below to check that the imputing step is working appropriately
#### (I judge that it is)
# demo.imputed %>%
#   # Separate out plants with gaps in records
#   filter(grepl('gap', dead.status)) %>%
#   # For each plant, get the first year observed dead and the first year after that observed alive
#   group_by(plantid) %>%
#   mutate(
#     y0 = min(Year[detected & !obs.alive]), 
#     y1 = min(Year[detected & obs.alive & (Year > y0)])
#   ) %>%
#   # Ungroup, get tag in separate column
#   ungroup() %>%
#   separate(plantid, into = c('tag', 'plot', 'coord'), sep = '_') %>%
#   # Impute "alive" for plants during their gap if gap is one year or is due to sampling in 2020
#   mutate(
#     alive = ifelse(
#       (dead.status %in% c('gap1', 'gap.2020')) & (Year >= y0 & Year < y1),
#       TRUE,
#       obs.alive
#     )
#   ) %>%
#   mutate(
#     tag2 = ifelse(
#       grepl('gap[2345]$', dead.status) & (Year >= y1),
#       paste0(tag, 'b'),
#       tag
#     )
#   ) %>% 
#   group_by(tag, plot, coord) %>%
#   filter(any(tag2 != tag) | any(alive != obs.alive)) %>% 
#   arrange(plot, tag, Year) %>% 
#   print(n = nrow(.))


demo.imputed = rbind(
  # Plants where we don't need to do any fixing/correcting (yet)
  demo.imputed %>%
    filter(!grepl('gap', dead.status)) %>%
    # 'alive' column is imputed based on above
    mutate(alive = obs.alive),
  # Plants where we do need correcting - make these corrections here
  demo.imputed %>%
    # Separate out plants with gaps in records
    filter(grepl('gap', dead.status)) %>%
    # For each plant, get the first year observed dead and the first year after that observed alive
    group_by(plantid) %>%
    mutate(
      y0 = min(Year[detected & !obs.alive]), 
      y1 = min(Year[detected & obs.alive & (Year > y0)])
    ) %>%
    # Ungroup, get tag in separate column
    ungroup() %>%
    separate(plantid, into = c('tag', 'plot', 'coord'), sep = '_') %>%
    # Impute "alive" for plants during their gap if gap is one year or is due to sampling in 2020
    mutate(
      alive = ifelse(
        (dead.status %in% c('gap1', 'gap.2020')) & (Year >= y0 & Year < y1),
        TRUE,
        obs.alive
      )
    ) %>%
    mutate(
      tag = ifelse(
        grepl('gap[2345]$', dead.status) & (Year >= y1),
        paste0(tag, 'b'),
        tag
      )
    ) %>%
    unite(plantid, c(tag, plot, coord), sep = '_') %>%
    select(-c(y0, y1))
)

nrow(demo.imputed)

# IT WORKED ON THE FIRST TRY hell yeah

##############################################
##### Tidy up data frame to get *one* dead record per plant + relevant cols
##############################################

demo.imputed = demo.imputed %>% arrange(plantid, Year)

head(demo.imputed)

### Now that I've gone back and separated old plants from new recruits,
# I could just say "if the last record is alive, it was always alive"
# let's just do this (and do it before removing 2023 plants!)
demo.imputed %>%
  group_by(plantid) %>%
  # annoyingly this only works if I add a new field for alive.last...
  mutate(alive.last = alive[Year == max(Year)]) %>%
  mutate(alive2 = ifelse(alive.last, TRUE, alive)) %>%
  filter(any(alive2 != alive)) # cool -  this works

demo.imputed = demo.imputed %>%
  group_by(plantid) %>%
  # annoyingly this only works if I add a new field for alive.last...
  mutate(alive.last = alive[Year == max(Year)]) %>%
  mutate(alive = ifelse(alive.last, TRUE, alive)) %>%
  # tidy up
  select(-alive.last) %>%
  ungroup()

### First: suppose we should get rid of 2023 records (due to imputing... sad and not ideal!)
# (shit this also would mean no 2023 kernel in general... shit man)
# Also, upon looking at the 2022 data - the 'seen.once' plants in 2022 were, mysteriously, not visited in 2021
# (ahhhhh)
# could fix this in later steps but for now just remove all of these plants...
demo.imputed = demo.imputed %>% 
  filter(Year < 2023) %>%
  filter(!(Year %in% 2022 & dead.status %in% 'seen.once'))

nrow(demo.imputed)
length(unique(demo.imputed$plantid)) # down to 979 plants

### Now, after we trimmed out some records (to minimize effects of censoring)
### let's remove plants that appear only once

demo.imputed = demo.imputed %>%
  group_by(plantid) %>%
  filter(n() > 1) %>%
  ungroup()

nrow(demo.imputed)
length(unique(demo.imputed$plantid)) # down to 870 plants...

### Now: for each plant with dead records, get *only* the first dead record
# (... can we do this with only the "alive" column? do we need "detected"?)

# how many times are plants observed as "dead"?
demo.imputed %>% group_by(plantid) %>% summarise(n.dead = sum(!alive)) %>% group_by(n.dead) %>% summarise(n = n())
# hmm...

# One way to maybe do this is to separate out the alive and dead records
# then take the records that are dead and only get the first (chronologically)
rbind(
  demo.imputed %>% filter(alive),
  demo.imputed %>% filter(!alive) %>% arrange(Year) %>% distinct(plantid, .keep_all = TRUE)
) %>%
  # Check this by going through each plant and getting the maximum gap in years created
  # Ideally all of these would be only 1
  arrange(plantid, Year) %>%
  group_by(plantid) %>%
  summarise(mdiff = max(diff(Year))) %>%
  group_by(mdiff) %>%
  summarise(n = n())
# some badness happening in here...
# ~25% of these have some kind of problem...

# some of them give -Inf which probably means only one record
# (would this record have to be for a dead plant then...?)
# but plenty with a gap greater than that...

rbind(
  demo.imputed %>% filter(alive),
  demo.imputed %>% filter(!alive) %>% arrange(Year) %>% distinct(plantid, .keep_all = TRUE)
) %>%
  group_by(plantid) %>%
  filter(n() < 2)
# yeah these are plants that were only seen dead...

# Solve the above by getting rid of plants that have all records as dead...
demo.imputed = demo.imputed %>%
  group_by(plantid) %>%
  filter(any(alive)) %>%
  ungroup()

# Now look at the rest
rbind(
  demo.imputed %>% filter(alive),
  demo.imputed %>% filter(!alive) %>% arrange(Year) %>% distinct(plantid, .keep_all = TRUE)
) %>%  
  arrange(plantid, Year) %>%
  group_by(plantid) %>%
  # ydiff column will tell us the gaps
  mutate(ydiff = c(NA, diff(Year)))

# ugh... a lot of these are cases where the plant was missed and then found dead after

rbind(
  demo.imputed %>% filter(alive),
  demo.imputed %>% filter(!alive) %>% arrange(Year) %>% distinct(plantid, .keep_all = TRUE)
) %>%  
  arrange(plantid, Year) %>%
  group_by(plantid) %>%
  # ydiff column will tell us the gaps
  mutate(ydiff = c(NA, diff(Year))) %>%
  ungroup() %>%
  filter(ydiff > 1 & !alive) %>%
  nrow()
# 31 of these plants
# come on man...

# Okay. Well, with an uncertain death date, I guess I need to get rid of these.
# All of this stuff is probably introducing bias...
demo.imputed = rbind(
  demo.imputed %>% filter(alive),
  demo.imputed %>% filter(!alive) %>% arrange(Year) %>% distinct(plantid, .keep_all = TRUE)
) %>%  
  arrange(plantid, Year) %>%
  group_by(plantid) %>%
  # ydiff column will tell us the gaps
  mutate(ydiff = c(NA, diff(Year))) %>%
  ungroup() %>%
  # Get rid of columns where there is a gap before mortality
  filter(!(ydiff > 1 & !alive)) %>%
  select(-ydiff)

nrow(demo.imputed)
length(unique(demo.imputed$plantid))
sum(!demo.imputed$alive) # 230 records of mortality...

# Gaps in survival... will want to resolve these
# Options are to merge... ugh

demo.impute.scaffold = demo.imputed %>%
  group_by(plantid) %>%
  summarise(y0 = min(Year), y1 = max(Year)) %>%
  uncount(y1 - y0 + 1) %>%
  group_by(plantid) %>%
  mutate(Year = y0 + (0:(n()-1))) %>%
  select(-c(y0, y1)) %>%
  ungroup()

merge(
  demo.imputed, demo.impute.scaffold,
  all.x = TRUE, all.y = TRUE
  ) %>%
  arrange(plantid, Year) %>% filter(is.na(alive))
# oh... only seven plants affected.
# maybe don't worry about these for now?


##############################################
##### Merge to get previous/next year data in here
##############################################

data.for.model = merge(
  # Demo table with relevant columns
  demo.imputed %>%
    mutate(size = log(no.leaves * leaf.leng)) %>%
    select(plantid, Plot, trt, Year, flowering, no.umbels, size, alive),
  demo.imputed %>%
    select(plantid, Year, alive) %>%
    mutate(Year = Year - 1),
  by = c('plantid', 'Year'), 
  all.x = TRUE, all.y = FALSE,
  suffixes = c('', '.next')
)

nrow(data.for.model)
head(data.for.model)

demo.imputed %>% head(8)
# looks okay

# Now, only relevant rows...
# - remove plants that are dead
# - remove plants that are alive but missing size or flowering
# - remove plants in final obs. eyar (i.e., alive.next is NA)
data.for.model = data.for.model %>%
  filter(alive & !is.na(size) & !is.na(flowering)) %>%
  filter(!is.na(alive.next))

head(data.for.model)
nrow(data.for.model)
length(unique(data.for.model$plantid))
#
# welp... let's try it folks

##############################################
##### Try making some plots...
##############################################

data.for.model %>%
  ggplot(aes(x = size, y = alive.next)) +
  geom_point(
    aes(
      colour = trt, shape = flowering
    ),
    size = 3, alpha = 0.5, position = position_jitter(height = 0.25)
  ) +
  scale_shape_manual(values = c(21, 19)) +
  scale_colour_manual(values = c('black', 'red', 'blue')) +
  facet_wrap(~ Year)

data.for.model %>%
  ggplot(aes(x = size, group = trt, fill = trt)) +
  geom_histogram(position = 'identity', alpha = 0.5, binwidth = 0.25) +
  scale_fill_manual(values = c('black', 'red', 'blue')) +
  facet_wrap(~ Year)
# not seeing big size discrepancies from year to year

##############################################
##### Try fitting some models...
##############################################

### Null model
surv.0 = glmer(
  formula = alive.next ~ (1 | Plot),
  data = data.for.model,
  family = 'binomial'
)

# singularity if including plant-level random effect (ugh)

### Year random effect
surv.y = glmer(
  formula = alive.next ~ (1 | Year) + (1 | Plot),
  data = data.for.model,
  family = 'binomial'
)

AIC(surv.y, surv.0)
# go with the year model
summary(surv.y)
ilogit(2.2986) # baseline 90% survival
ranef(surv.y)
hist(unlist(ranef(surv.y)$Plot))
hist(unlist(ranef(surv.y)$Year))
plot(unlist(ranef(surv.y)$Year))
# interesting...

### Previous size
surv.y.s = glmer(
  formula = alive.next ~ size + (1 | Year) + (1 | Plot),
  data = data.for.model,
  family = 'binomial'
)

summary(surv.y.s)
AIC(surv.y.s, surv.y)
# neat - keep survival in!
# size... what is the unit of size? multiplying size by e?
# anyway
hist(unlist(ranef(surv.y.s)$Plot))
hist(unlist(ranef(surv.y.s)$Year))
plot(unlist(ranef(surv.y.s)$Year))


### Flowering
surv.y.s.f = glmer(
  formula = alive.next ~ size + flowering + (1 | Year) + (1 | Plot),
  data = data.for.model,
  family = 'binomial'
)

AIC(surv.y.s.f, surv.y.s)
# no effect of flowering after accounting for size!
# interesting

### Treatment
surv.y.s.t = glmer(
  formula = alive.next ~ size + trt + (1 | Year) + (1 | Plot),
  data = data.for.model,
  family = 'binomial'
)

AIC(surv.y.s.t, surv.y.s)
# treatment effect no longer present

summary(surv.y.s.t)

### Size effects varying by year?

### Previous size
surv.ys = glmer(
  formula = alive.next ~ trt + (size | Year) + (1 | Plot),
  data = data.for.model,
  family = 'binomial'
)

AIC(surv.ys, surv.y.s, surv.y.s.t)
# no good evidence of size-year interactions

### Size-treatment interaction
surv.y.st = glmer(
  formula = alive.next ~ size * trt + (1 | Year) + (1 | Plot),
  data = data.for.model,
  family = 'binomial'
)

AIC(surv.y.st, surv.y.s.t)
# no interaction

### Treatment effects varying by year
surv.yt.s = glmer(
  formula = alive.next ~ size + (trt | Year) + (1 | Plot),
  data = data.for.model,
  family = 'binomial'
)

AIC(surv.yt.s, surv.y.s.t)
# nope.

##############################################
##### Plot results (from this model at least)
##############################################

data.scaffold = expand.grid(
  trt = unique(data.for.model$trt),
  Year = unique(data.for.model$Year),
  size = (0:60)/10
) %>%
  mutate(
    pred.s = predict(surv.y.s.t, newdata = ., re.form = ~ (1 | Year)) %>% ilogit()
  )

head(data.scaffold)

data.for.model %>%
  ggplot(aes(x = size, y = as.numeric(alive.next))) +
  geom_point(
    aes(colour = trt),
    size = 2, alpha = 0.25, position = position_jitter(height = 0.25)
  ) +
  geom_line(
    data = data.scaffold,
    aes(x = size, y = pred.s, group = trt, colour = trt)
  ) +
  scale_colour_manual(values = c('black', 'red', 'blue')) +
  labs(x = 'Size', y = 'Survival') +
  facet_wrap(~ Year)

