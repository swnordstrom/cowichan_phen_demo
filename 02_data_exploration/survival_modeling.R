### Another script for looking at survival
### Here, the dataset is cleaner (I did more tag reconciliation)
### so there should be fewer record gaps.
### In the script I'll impute survival for short gaps and split plants with
### longer gaps.
### SN - init 21 Nov. 2023

library(ggplot2)
library(dplyr)
library(tidyr)
library(lme4)

rm(list = ls())

demo = merge(
  x = read.csv('01_data_cleaning/out/demo_postcombine.csv'),
  y = read.csv('00_raw_data/plot_treatments.csv') %>% rename(Plot = plot)
) %>%
  # Get rid of unnecessary columns
  select(-c(Xcoor, Ycoor, YrTag, Tag, Plot)) %>%
  arrange(plantid, Year)

head(demo)
nrow(demo)

##### Add "observed" column 

demo = demo %>%
  mutate(
    obs.alive = case_when(
      (No.umbels > 0) & !is.na(No.umbels) ~ TRUE,
      (!No.leaves) & !is.na(No.leaves) ~ FALSE,
      No.leaves > 0 & !is.na(No.leaves) ~ TRUE,
      is.na(No.leaves) ~ FALSE,
      .default = NA
    )
  )

demo %>% group_by(obs.alive) %>% summarise(n = n())

head(demo)

dead.maybe = demo %>%
  group_by(plantid) %>%
  filter(any(!obs.alive)) %>%
  ungroup()

nrow(dead.maybe)

dead.maybe %>%
  group_by(plantid) %>%
  summarise(
    min.dead = min(Year[!obs.alive]),
    gap.size = max(diff(Year))
  ) %>%
  group_by(gap.size) %>%
  summarise(n = n())
# whoa...?
# how can this be true...

dead.maybe %>%
  group_by(plantid) %>%
  filter(max(diff(Year)) > 1) %>%
  print(n = nrow(.))
# ah this is literally just missing from df, not counting NAs
# three cases, two of which (3585 and 3807) have gaps that don't affect
# mortality estimates, but one (3844) does


# Try this: how many plants have any records where there are dead records before
# the *latest* alive record
dead.maybe %>%
  group_by(plantid) %>%
  summarise(bad = any(Year[!obs.alive] < max(Year[obs.alive]))) %>%
  group_by(bad) %>%
  summarise(n = n())
# man... woof! goodness

# Take a look at some of these
dead.maybe %>%
  group_by(plantid) %>%
  filter(any(Year[!obs.alive] < max(Year[obs.alive])))

# Crop out any records before the first observed mortality
dead.maybe %>%
  group_by(plantid) %>%
  filter(any(Year[!obs.alive] < max(Year[obs.alive]))) %>%
  filter(Year >= min(Year[!obs.alive])) %>%
  summarise(gap.size = min(Year[obs.alive]) - min(Year)) %>%
  group_by(gap.size) %>%
  summarise(n = n())
  
dead.maybe %>%
  group_by(plantid) %>%
  filter(any(Year[!obs.alive] < max(Year[obs.alive]))) %>%
  filter(Year >= min(Year[!obs.alive])) %>%
  summarise(y0 = min(Year), y1 = min(Year[obs.alive])) %>%
  mutate(
    gap.len = y1 - y0,
    in.2020 = y1 >= 2020 & y0 <= 2020
  ) %>%
  group_by(gap.len, in.2020) %>%
  summarise(n = n())
# seems pretty easy...

dead.maybe %>%
  group_by(plantid) %>%
  filter(any(Year[!obs.alive] < max(Year[obs.alive]))) %>%
  filter(Year >= min(Year[!obs.alive])) %>%
  filter(min(Year[obs.alive]) - min(Year) > 1, !(min(Year) < 2020 & min(Year[obs.alive] > 2020)))
# looking at a few of these (not every single one), they look like they could be new plants

dead.maybe %>%
  group_by(plantid) %>%
  filter(any(Year[!obs.alive] < max(Year[obs.alive]))) %>%
  filter(Year >= min(Year[!obs.alive])) %>%
  filter(min(Year[obs.alive]) - min(Year) > 1, !(min(Year) < 2020 & min(Year[obs.alive] > 2020))) %>%
  distinct(plantid) %>%
  nrow()
# 31 plants

# Try to do this in one line...

demo2 = rbind(
  # All alive plants
  demo %>% 
    group_by(plantid) %>% 
    filter(all(obs.alive)) %>%
    ungroup() %>%
    mutate(surv = obs.alive),
  # Dead plants, but without any resurrections
  demo %>%
    group_by(plantid) %>%
    filter(any(!obs.alive)) %>%
    filter(!any(Year[obs.alive] > min(Year[!obs.alive]))) %>%
    ungroup() %>%
    mutate(surv = obs.alive),
  # Dead plants, where gap includes 2020 or is only one year
  # IMPUTE survival here
  demo %>%
    group_by(plantid) %>%
    filter(any(!obs.alive)) %>%
    filter(any(Year[obs.alive] > min(Year[!obs.alive]))) %>%
    mutate(
      # First year observed not alive
      y0 = min(Year[!obs.alive]),
      # First year observed alive after death (resurrected)
      y1 = min(Year[obs.alive & Year > y0]),
      # Gap size (number of years between y0 and y1)
      gapsize = y1 - y0,
      # Numeric for if record gap includes 2020
      in.2020 = as.numeric(y1 > 2020 & y0 <= 2020)
    ) %>%
    ungroup() %>%
    # Separate out plantid into tag, plot, coords
    separate(plantid, into = c('Tag', 'Plot', 'Coord'), sep = '_') %>%
    mutate(
      # add 'b' to tag if gap 
      Tag = ifelse((gapsize > (1 + in.2020)) & (Year >= y1), 
                   paste0(Tag, 'b'), 
                   Tag),
      # Create 'alive' column sensu above
      surv = ifelse((gapsize <= (1 + in.2020) & (Year < y1 & Year >= y0)),
                    TRUE,
                    obs.alive)
    ) %>%
    # Delete old columns
    select(-c(y0, y1, gapsize, in.2020)) %>%
    # Add notes
    mutate(
      proc.note = ifelse(grepl('b', Tag), 
                         paste0(proc.note, '; tag modified for gap'),
                         proc.note),
      proc.note = ifelse(surv != obs.alive,
                         paste0(proc.note, '; survival imputed based on subsq. record'),
                         proc.note),
      edited = ifelse(grepl('b', Tag) | surv != obs.alive, TRUE, edited)
    ) %>%
    unite(Tag, Plot, Coord, col = 'plantid', sep = '_') %>%
    select(c(names(demo), surv))
) %>%
  arrange(plantid, Year)

nrow(demo2)
head(demo2)
length(unique(demo2$plantid))

# Remaining plants with gaps?
demo2 %>% group_by(plantid) %>% filter(any(diff(surv) > 0)) %>% distinct(plantid) %>% nrow()
# that's still a lot... guess there are many plants with multiple record gaps

# A while loop would work great for this... but is it worth doing... probably not.
# (I checked - there are no plants with gaps large enough to need tag additions)

### Pass 1
demo2 = rbind(
  # All alive plants
  demo2 %>% 
    group_by(plantid) %>% 
    filter(all(surv)) %>%
    ungroup(),
  # Dead plants, but without any resurrections
  demo2 %>%
    group_by(plantid) %>%
    filter(any(!surv)) %>%
    filter(!any(Year[surv] > min(Year[!surv]))) %>%
    ungroup(),
  # Dead plants, where gap includes 2020 or is only one year
  # IMPUTE survival here
  demo2 %>%
    group_by(plantid) %>%
    filter(any(!surv)) %>%
    filter(any(Year[surv] > min(Year[!surv]))) %>%
    mutate(
      # First year observed not alive
      y0 = min(Year[!surv]),
      # First year observed alive after death (resurrected)
      y1 = min(Year[surv & Year > y0]),
      # Gap size (number of years between y0 and y1)
      gapsize = y1 - y0,
      # Numeric for if record gap includes 2020
      in.2020 = as.numeric(y1 > 2020 & y0 <= 2020)
    ) %>%
    ungroup() %>%
    # Separate out plantid into tag, plot, coords
    mutate(
      # Make notes and add edited flag
      proc.note = ifelse((gapsize <= (1 + in.2020) & (Year < y1 & Year >= y0)),
                         paste0(proc.note, '; survival imputed based on subsq. record'),
                         proc.note),
      edited = ifelse((gapsize <= (1 + in.2020) & (Year < y1 & Year >= y0)), TRUE, edited),
      # Modify 'surv' column to impute survival
      surv = ifelse((gapsize <= (1 + in.2020) & (Year < y1 & Year >= y0)),
                    TRUE,
                    surv)
    ) %>%
    # Delete old columns
    select(-c(y0, y1, gapsize, in.2020)) %>%
    # Add notes
    select(names(demo2))
) %>%
  arrange(plantid, Year)

### Pass 2
demo2 = rbind(
  # All alive plants
  demo2 %>% 
    group_by(plantid) %>% 
    filter(all(surv)) %>%
    ungroup(),
  # Dead plants, but without any resurrections
  demo2 %>%
    group_by(plantid) %>%
    filter(any(!surv)) %>%
    filter(!any(Year[surv] > min(Year[!surv]))) %>%
    ungroup(),
  # Dead plants, where gap includes 2020 or is only one year
  # IMPUTE survival here
  demo2 %>%
    group_by(plantid) %>%
    filter(any(!surv)) %>%
    filter(any(Year[surv] > min(Year[!surv]))) %>%
    mutate(
      # First year observed not alive
      y0 = min(Year[!surv]),
      # First year observed alive after death (resurrected)
      y1 = min(Year[surv & Year > y0]),
      # Gap size (number of years between y0 and y1)
      gapsize = y1 - y0,
      # Numeric for if record gap includes 2020
      in.2020 = as.numeric(y1 > 2020 & y0 <= 2020)
    ) %>%
    ungroup() %>%
    # Separate out plantid into tag, plot, coords
    mutate(
      # Make notes and add edited flag
      proc.note = ifelse((gapsize <= (1 + in.2020) & (Year < y1 & Year >= y0)),
                         paste0(proc.note, '; survival imputed based on subsq. record'),
                         proc.note),
      edited = ifelse((gapsize <= (1 + in.2020) & (Year < y1 & Year >= y0)), TRUE, edited),
      # Modify 'surv' column to impute survival
      surv = ifelse((gapsize <= (1 + in.2020) & (Year < y1 & Year >= y0)),
                    TRUE,
                    surv)
    ) %>%
    # Delete old columns
    select(-c(y0, y1, gapsize, in.2020)) %>%
    # Add notes
    select(names(demo2))
) %>%
  arrange(plantid, Year)

demo2 %>% group_by(plantid) %>% filter(any(diff(surv) > 0)) %>% distinct(plantid) %>% nrow()
# cool

# I want to add Plot as a covariate - scrape that out of the plantid

demo2 = demo2 %>%
  separate(col = plantid, into = c("Tag", "Plot", "Coord"), sep = "_", remove = FALSE) %>%
  select(-c(Tag, Coord))

head(demo2)

demo.mod = merge(
  x = demo2 %>% filter(surv) %>% select(-surv),
  y = demo2 %>% mutate(Year = Year - 1) %>% select(plantid, Year, surv),
  by = c("plantid", "Year")
)

nrow(demo.mod)
head(demo.mod)

demo.mod %>% group_by(plantid) %>% summarise(n.dead = sum(!surv)) %>% distinct(n.dead)
# shouldn't be any 2s...
# (multiple gaps - code currently does not account for this)

# Finally - want to get rid of records without appropriate size info
# OH lol this will just get rid of all of the records I imputed won't it
# lmao
# (oh no actually this is fine - we wanted those for the response in the next year)
# (will also want to get rid of a handful of imputed records with zero leaves)
# Then go ahead and add a size column
demo.mod = demo.mod %>% 
  filter(!(is.na(No.leaves) | is.na(Leaf.length) | !No.leaves | !Leaf.length)) %>%
  mutate(size = log(No.leaves * Leaf.length)) %>%
  # Add a field for flowering too
  mutate(flowering = (No.umbels > 0) & !is.na(No.umbels))

quantile(demo.mod$size)
hist(demo.mod$size)
# possibly a little skewed, but this looks normal-like to me

table(demo.mod$flowering, useNA = 'always')

# Okay, well, anyway

#####
##### Start models
#####

table(demo.mod$surv) # mortality is rare!

### Null model
s_0 = glmer(
  formula = surv ~ 1 + (1 | Plot),
  data = demo.mod,
  family = 'binomial'
)

s_0

### Model with size
s_s = glmer(
  formula = surv ~ size + (1 | Plot),
  data = demo.mod,
  family = 'binomial'
)

s_s
AIC(s_s, s_0)
# Include size

### Model with year and size
s_y.s = glmer(
  formula = surv ~ size + (1 | Year) + (1 | Plot),
  data = demo.mod,
  family = 'binomial'
)

s_y.s
AIC(s_y.s, s_s) # wow... only a marginal improvement from adding year???

### Model with yearly effects of size
s_ys = glmer(
  formula = surv ~ (size | Year) + (1 | Plot),
  data = demo.mod,
  family = 'binomial'
)

s_ys
AIC(s_ys, s_y.s, s_s)
# no evidence for yearly effects of size

### Models with effect of flowering
s_y.f.s = glmer(
  formula = surv ~ size + flowering + (1 | Year) + (1 | Plot),
  data = demo.mod,
  family = 'binomial'
)

s_f.s = glmer(
  formula = surv ~ size + flowering + (1 | Plot),
  data = demo.mod,
  family = 'binomial'
)

s_y.f.s
s_f.s
AIC(s_y.f.s, s_f.s, s_y.s, s_s)
# best supported model includes flowering, not necessarily year

### Treatment effects
s_y.f.s.t = glmer(
  formula = surv ~ size + flowering + trt + (1 | Year) + (1 | Plot),
  data = demo.mod,
  family = 'binomial'
)

s_f.s.t = glmer(
  formula = surv ~ size + flowering + trt + (1 | Plot),
  data = demo.mod,
  family = 'binomial'
)

s_y.f.s.t
s_f.s.t
AIC(s_y.f.s, s_f.s, s_y.f.s.t, s_f.s.t)
# No evidence of treatment effect

### NOTE none of these account for possible overdispersion...

# So best performing model is/are size, flowering, maybe year...

summary(s_y.f.s)

ranef(s_y.f.s)$Plot %>%
  mutate(Plot = row.names(.)) %>%
  merge(demo.mod %>% distinct(Plot, trt)) %>%
  ggplot(aes(x = `(Intercept)`, fill = trt)) +
  geom_histogram(binwidth = 0.25) +
  scale_fill_manual(values = c('black', 'red', 'blue'))
# looks like *super* slight evidence of benefit of irrigation
# (but obviously not strong enough to appear in model...)

ranef(s_f.s)$Plot %>%
  mutate(Plot = row.names(.)) %>%
  merge(demo.mod %>% distinct(Plot, trt)) %>%
  ggplot(aes(x = `(Intercept)`, fill = trt)) +
  geom_histogram(binwidth = 0.25) +
  scale_fill_manual(values = c('black', 'red', 'blue'))
# same here

ranef(s_y.f.s)$Year %>%
  mutate(Year = row.names(.)) %>%
  ggplot(aes(x = Year, y = `(Intercept)`)) +
  geom_point()
# Interesting... seemingly downward negative trend...
# (although model here is not super well supported)

with(demo.mod, table(Year, surv))
# huh - raw data counts really look identical between 2021 and 2022
# interesting

demo.mod %>%
  mutate(surv = as.numeric(surv)) %>%
  ggplot(aes(x = size, y = surv, colour = trt)) +
  geom_point(alpha = 0.5, position = position_jitter(height = 0.25)) +
  scale_colour_manual(values = c('black', 'red', 'blue')) +
  facet_wrap(~ Year)
# Yeah... 2016 numbers are just whacky

# get predictions from models

s_f.s_pred = expand.grid(size = (5:50)/10, flowering = c(TRUE, FALSE)) %>%
  mutate(pred = predict(s_f.s, newdata = ., re.form = NA)) %>%
  mutate(pred = 1 / (1+exp(-pred)))

demo.mod %>%
  mutate(surv = as.numeric(surv)) %>%
  ggplot(aes(x = size, y = surv)) +
  geom_line(
    data = s_f.s_pred,
    aes(x = size, y = pred, group = flowering, linetype = flowering)
  ) +
  geom_point(aes(colour = trt), alpha = 0.5, position = position_jitter(height = 0.25)) +
  scale_colour_manual(values = c('black', 'red', 'blue')) +
  facet_wrap(~ Year)

##### Should look at distributions of size components

demo.mod %>%
  ggplot(aes(x = No.leaves, y = Leaf.length)) +
  geom_point() +
  facet_wrap(~ Year)
# whoa these leaf length sin 2016 are small...

demo.mod %>%
  mutate(Year = factor(Year)) %>%
  ggplot(aes(x = No.leaves, group = Year, fill = Year)) +
  geom_histogram(binwidth = 1, alpha = 0.5, position = 'identity')
# not as helpful as I hoped
# although you can still see 2016 is different

demo.mod %>%
  group_by(Year, No.leaves, trt) %>%
  summarise(n = n()) %>%
  group_by(Year, trt) %>%
  mutate(p = n / sum(n)) %>%
  ggplot(aes(x = No.leaves, y = Year, fill = p)) +
  geom_tile() +
  scale_fill_viridis_c() +
  facet_wrap(~ trt)
# curious...

demo.mod %>%
  mutate(Year = factor(Year)) %>%
  ggplot(aes(x = Leaf.length, group = Year, fill = Year)) +
  geom_histogram(binwidth = 1, alpha = 0.5, position = 'identity')
# definitely looks like a shift over time...
# but this could be ageing!

demo.mod %>%
  mutate(Year = factor(Year)) %>%
  ggplot(aes(x = Leaf.length, group = Year, colour = Year)) +
  geom_density(position = 'identity')
# wowza
# yeah this definitely looks like ageing
# argh...
# (but also - the leaf lengths are much smaller in 2016!)

#####
##### Reduced models
#####

# Try models without 2016 and 2022
# 2016: size measurements are bad
# 2022: final status (2023) is somewhat uncertain

demo.reduced = demo.mod %>% filter(Year %in% 2017:2021)

table(demo.reduced$surv)

### Null model
r_0 = glmer(
  formula = surv ~ (1 | Plot),
  data = demo.reduced,
  family = 'binomial'
)

### Model with size
r_s = glmer(
  formula = surv ~ size + (1 | Plot),
  data = demo.reduced,
  family = 'binomial'
)

AIC(r_s, r_0)

### Model with size + year
r_y.s = glmer(
  formula = surv ~ size + (1 | Year) + (1 | Plot),
  data = demo.reduced,
  family = 'binomial'
)

AIC(r_y.s, r_s)
# No year effects!!!

r_ys = glmer(
  formula = surv ~ (size | Year) + (1 | Plot),
  data = demo.reduced,
  family = 'binomial'
)
# singularity

### Model with treatment

r_s.t = glmer(
  formula = surv ~ size + trt + (1 | Plot),
  data = demo.reduced,
  family = 'binomial'
)

AIC(r_s.t, r_s)
# no treatment effect...

### Model with flowering
r_s.f = r_y.s = glmer(
  formula = surv ~ size + flowering + (1 | Plot),
  data = demo.reduced,
  family = 'binomial'
)

AIC(r_s.f, r_s)
# no effect of flowering either?
# lmao
# woof

ranef(r_s)$Plot %>%
  mutate(Plot = row.names(.)) %>%
  merge(demo.mod %>% distinct(Plot, trt)) %>%
  ggplot(aes(x = `(Intercept)`, fill = trt)) +
  geom_histogram(binwidth = 0.25) +
  scale_fill_manual(values = c('black', 'red', 'blue'))
# leptokurtotic?
# does kind of show mild trend with treatment, once again

summary(r_s)

demo.reduced %>%
  mutate(surv = as.numeric(surv)) %>%
  ggplot(aes(x = size, y = surv)) +
  geom_point(alpha = 0.5, position = position_jitter(width = 0, height = 0.1)) +
  facet_grid(rows = vars(Year), cols = vars(trt))

# Model fits
r_s_pred = data.frame(size = (5:50)/10) %>%
  mutate(pred = predict(r_s, newdata = data.frame(size = (5:50)/10), re.form = NA)) %>%
  mutate(pred = 1 / (1+exp(-pred)))

demo.reduced %>%
  mutate(surv = as.numeric(surv)) %>%
  ggplot(aes(x = size, y = surv)) +
  geom_line(
    data = r_s_pred,
    aes(x = size, y = pred)
  ) +
  geom_point(alpha = 0.5, position = position_jitter(width = 0, height = 0.1)) +
  facet_grid(rows = vars(Year), cols = vars(trt))

range(r_s_pred$pred)
# interesting

### Try to compare size measurements!

r_lenl = glmer(
  formula = surv ~ log(Leaf.length) + (1 | Plot),
  data = demo.reduced,
  family = 'binomial'
)

r_no.l = glmer(
  formula = surv ~ log(No.leaves) + (1 | Plot),
  data = demo.reduced,
  family = 'binomial'
)

AIC(r_s, r_lenl, r_no.l)
# size works best
# okie dokie

### Compare measurements...

merge(r_s_pred, s_f.s_pred, by = 'size', suffixes = c('reduced', 'full')) %>%
  ggplot(aes(x = size)) +
  geom_line(aes(y = predreduced)) +
  geom_line(aes(y = predfull, group = flowering, colour = flowering))
# hmm... okay these differences are kinda sizeable lmao
# gosh

