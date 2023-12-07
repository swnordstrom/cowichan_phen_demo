# Script for initial assessment and compiling/aggregating raw climate data.
# The climate data was compiled by Jenna - currently backed up on Dropbox.

library(ggplot2)
library(dplyr)
library(tidyr)

rm(list = ls())

clim.mo = read.csv('00_raw_data/climate/cgop_weather_monthly_interp.csv')
clim.dy = read.csv('00_raw_data/climate/cgop_weather_daily_interp.csv')
clim.hr = read.csv('00_raw_data/climate/cgop_weather_hourly.csv')

##### Look at hourly data frame

head(clim.hr) # nice!
nrow(clim.hr) # wow!
str(clim.hr)

length(unique(clim.hr$Time))

clim.hr = clim.hr %>%
  separate(Date, into = c('Year', 'Month', 'Day'), sep = '-', remove = FALSE) %>%
  mutate(Date = as.Date(Date)) %>%
  separate(Time, into = c('Hour', 'Minute', 'Second'), sep = ':')
# slow...

head(clim.hr)
unique(clim.hr$Minute) # ah... some NAs... also weird minutes
# yeah... this is a lot of into. maybe not needed.

##### Look at daily data frame

head(clim.dy) # ooh... cool
# would be nice if there was a flag for imputing though...
# ah well.

clim.dy = clim.dy %>%
  separate(Date, into = c('Year', 'Month', 'Day'), sep = '-', remove = FALSE) %>%
  mutate(Date = as.Date(Date)) %>%
  mutate(jDate = as.Date(paste('2020', Month, Day, sep = '-'), format = '%Y-%m-%d'))

head(clim.dy)

nrow(clim.dy)
clim.dy %>% group_by(Year) %>% summarise(n = n()) # wow!
# (wait why are there 365 records for 2023...)
# are they NAs?

clim.dy %>% filter(Year %in% 2023) %>% tail() 
# yes

# Look at variables over time!
clim.dy %>%
  pivot_longer(contains('Temp'), names_to = 'tempStat', values_to = 'value') %>%
  mutate(tempStat = gsub('Temp\\_C', '', tempStat)) %>%
  ggplot(aes(x = jDate, y = value, group = Year)) +
  geom_line(aes(colour = Year)) +
  facet_wrap(~ tempStat, nrow = 3)

# Oh there are some NAs
# where?
clim.dy %>%
  filter(is.na(total_precip_mm) | is.na(minTemp_C) | is.na(maxTemp_C) | is.na(AveTemp_C))
# interesting...
# mostly from 2012 though!
# OH it looks like all of these are 2012 or post-cleaning in 2023!!

clim.dy %>%
  filter(Year %in% 2013:2022) %>%
  filter(is.na(total_precip_mm) | is.na(minTemp_C) | is.na(maxTemp_C) | is.na(AveTemp_C))
# hell yeah Jenna!  

# Deviations from mean

# (first plot the means, removing 2012 maybe just because these seem spottier)
daily.means = clim.dy %>%
  filter(Year > 2012) %>%
  group_by(jDate) %>%
  summarise(
    across(
      c(total_precip_mm, minTemp_C, AveTemp_C, maxTemp_C),
      .fns = list(
        mean = ~ mean(.x, na.rm = TRUE), 
        var = ~ var(.x, na.rm = TRUE), 
        n = ~ sum(!is.na(.x))
      ),
      .names = '{.col}.{.fn}'
    )
  ) %>%
  # Now convert to one row per variable per day  
  # there's probably a more elegant way to do this with one pivot...
  pivot_longer(-jDate, names_to = 'varStat', values_to = 'varval') %>%
  separate(varStat, into = c('varb', 'stat'), sep = '\\.') %>%
  pivot_wider(names_from = stat, values_from = varval)

head(daily.means)

daily.means %>%
  ggplot(aes(x = jDate)) +
  geom_line(aes(y = mean)) +
  geom_ribbon(
    aes(ymin = mean - 2*sqrt(var/n), ymax = mean + 2*sqrt(var/n)),
    alpha = 0.25
  ) +
  facet_wrap(~ varb)
# fairly noisy

# Get z-scores
clim.dy.z = clim.dy %>%
  filter(Year > 2012) %>%
  group_by(jDate) %>%
  mutate(
    across(
      c(total_precip_mm, minTemp_C, AveTemp_C, maxTemp_C),
      .fn = function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
    )
  )

head(clim.dy.z)
# hell yeah!!!

clim.dy.z %>%
  select(-c(Date, Month, Day)) %>%
  pivot_longer(-c(Year, jDate), names_to = 'varb', values_to = 'z') %>%
  ggplot(aes(x = jDate, y = z, group = Year)) +
  geom_line(aes(y = 0), linetype = 2, colour = 'gray77') +
  geom_line(aes(colour = Year)) +
  facet_wrap(~ varb)
# ah... precipitation is not going to work for z-scores...
# everything else looks okay though
# (not sure what to do about precipitation... zeros are inevitable)
# maybe monthly precipitation?
# are there any months with zero precip?
# (also - why are there NAs here?)

# Monthly cumulative precip?

clim.prec.month = clim.dy %>%
  filter(Year %in% 2013:2022) %>%
  group_by(Year = as.numeric(Year), Month = as.numeric(Month)) %>%
  summarise(total_precip = sum(total_precip_mm, na.rm = TRUE))

clim.prec.month %>%
  ggplot(aes(x = Year + (Month-1)/12, y = total_precip)) +
  geom_point(aes(color = total_precip > 0)) +
  scale_colour_manual(values = c('red', 'black')) +
  scale_x_continuous(breaks = 2013:2023) +
  labs(x = 'Year', y = 'Cumulative precipitation') +
  theme(legend.position = 'none')
# ah... there is a one in here!
# (wonder if I can use order statistics to get an expected min of a lognormal...)

clim.prec.month %>%
  mutate(Year = factor(Year)) %>%
  ggplot(aes(x = Month, y = total_precip, group = Year)) +
  geom_line(aes(colour = Year), linewidth = 0.5) +
  scale_x_continuous(breaks = 2*(1:6))

# Precip regimes... looks like Nov-Jan is wet, Feb-Aug is dry, Sept-Oct starts
# getting wet again

# Freeze dates - when is the latest in each year?

clim.dy %>%
  filter(as.numeric(Month) < 9, minTemp_C < 0) %>%
  group_by(Year) %>%
  summarise(last.freeze = max(jDate)) %>%
  ungroup() %>%
  ggplot(aes(x = Year, y = last.freeze)) +
  geom_point(size = 3)
# so often in April, sometimes in March

##### Some possible predictors

### Precip
# - Early season precip (Feb - April)
# - Growing season precip (April - June)
# - Summer precip (June - August)
# - Fall-Winter precip (Sept - Jan)
### Temperatures - use same ranges for aggregation?
# - Early season mean (Feb - April)
# - Early season max (Feb - April)
# - Early season min? (Feb - April)
# - Growing season mean (April - June)
# - Growing season max (April - June)
# - Growing season min (April - June)
# - Summer mean (June - August)
# - Summer max (June - August)
# - Winter mean (Sept - Jan?)
# - Winter max? (Sept - Jan)
# - Last freeze date

### What is a "year"?
# - Growing season of year is April - June of current year (t)
# - so, let year t be June 15 t-1 to June 14 t?
# (should think about this in terms of recruitment... )

# cut.breaks = data.frame(
#   break.day = as.Date(
#     c(paste0('2020-', c('02-01', '04-15', '06-15', '09-01')), '2021-01-01'),
#     format = '%Y-%m-%d'
#   ),
#   label = c('winter1', 'early', 'growing', 'summer', 'winter2')
# )

### Bin data together

# I'll use the cut() function for this
# First - define the breaks for the cut
# these will be the start dates to each one of these time-windows

cut.breaks = as.Date(
    c(paste0('2020-', c('01-01', '02-01', '04-15', '06-15', '09-01')), '2021-01-01'),
    format = '%Y-%m-%d'
  )

# Test out cuts on the jDate column of our dataset
cut(clim.dy$jDate, breaks = cut.breaks, right = FALSE)
# Do any of the
sum(is.na(cut(clim.dy$jDate, breaks = cut.breaks, right = FALSE)))
# looks like this will work...

# Bin dates into time windows (but keep everything else as needed)
clim.summ = clim.dy %>%
  # Convert year, month and day to numeric (for filtering)
  mutate(across(c(Year, Month, Day), as.numeric)) %>%
  # Get only dates in windows with all days present (and having data) in data
  # frame
  filter(
    # all data in 2014-2022
    Year %in% 2014:2022 | 
    # Data in February and later in 2013
    (Year %in% 2013 & Month > 1) | 
    # Data before June 15 in 2023
    (Year %in% 2023 & Month < 6) | (Year %in% 2023 & Month %in% 6 & Day %in% 1:14)
  ) %>%
  # Add the date period (window) using cut()
  mutate(period.start = cut(jDate, breaks = cut.breaks, right = FALSE)) %>%
  # Add some text labels to each window
  # (cut() above will just assign the date - not super legible)
  merge(y = data.frame(
    period.start = cut.breaks,
    # (I don't think winter3 will be assigned to anything...)
    pd.label = c('winter1', 'early', 'grow', 'summer', 'winter2', 'winter3')
    )
  ) %>%
  # Change the `date.pd` column to be in the proper year
  mutate(period.start = as.Date(paste0(Year, gsub('2020', '', period.start)))) %>%
  # Assign summer and winter to the successive year
  mutate(Year = Year + as.logical(pd.label %in% c('summer', 'winter2', 'winter3'))) %>%
  # get rid of numeral after winter[123] to make a single 'winter' label
  # and change the start date for winter to be consistent for all groups
  mutate(pd.label = gsub('\\d', '', pd.label)) %>%
  group_by(Year, pd.label) %>%
  mutate(period.start = min(period.start)) %>%
  ungroup()

# Check period labels and start dates
clim.summ %>% distinct(Year, pd.label, period.start) %>% print(n = nrow(.))

# Get 
clim.summ = merge(
  clim.summ %>%
    group_by(Year, pd.label, period.start) %>%
    summarise(
      # Means of each temperature measure
      meanTemp = mean(AveTemp_C),
      maxxTemp = mean(maxTemp_C),
      minnTemp = mean(minTemp_C),
      # Sum of precipitation
      prec.sum = sum(total_precip_mm)
    ) %>%
    ungroup(),
  clim.summ %>%
    filter(minTemp_C < 0, pd.label %in% c('early', 'grow')) %>%
    group_by(Year) %>%
    summarise(last.freeze = max(jDate)) %>%
    ungroup()
)

head(clim.summ)
nrow(clim.summ)

# Arrange by period start so these are in chrono order
clim.summ = clim.summ %>% arrange(period.start)

# Make a wide version of this data frame for evaluating correlations
clim.wide = clim.summ %>%
  # (the window start dates aren't necessary)
  select(-period.start) %>%
  pivot_longer(c(prec.sum, contains("Temp")), names_to = "stat", values_to = "statVal") %>%
  pivot_wider(names_from = c(pd.label, stat), values_from = statVal)

clim.corr.df = clim.wide[complete.cases(clim.wide),] %>%
  select(-Year) %>%
  mutate(last.freeze = as.numeric(last.freeze)) %>%
  cor() %>%
  as.data.frame() %>%
  mutate(vara = row.names(.)) %>%
  pivot_longer(-vara, names_to = 'varb', values_to = 'corr')

clim.corr.df %>%
  ggplot(aes(x = vara, y = varb)) +
  geom_tile(aes(alpha = abs(corr), fill = corr)) +
  scale_fill_gradient2(low = 'blue', high = 'red', mid = 'white', midpoint = 0) +
  theme(axis.text.x = element_text(angle = 90), legend.position = 'none')
# interesting... temperatures highly correlated - not surprising

clim.corr.df %>%
  filter(vara != varb) %>%
  arrange(desc(abs(corr)))
# damn wtf... super strong correlations between growing season precipitaton and
# mean/max temperatures

clim.corr.df %>%
  filter(vara != varb) %>%
  filter(abs(corr) > 0.5) %>%
  arrange(desc(abs(corr))) %>%
  print(n = nrow(.))

clim.wide %>%
  ggplot(aes(x = grow_meanTemp, y = grow_prec.sum)) +
  geom_label(aes(label = Year))
# whoa...

clim.dy %>%
  filter(as.numeric(Year) %in% 2020:2022, as.numeric(Month) %in% 3:6) %>%
  select(Year, jDate, AveTemp_C, total_precip_mm) %>%
  # arrange(Year, jDate) %>%
  # mutate(total_precip_mm = cumsum(total_precip_mm)) %>%
  pivot_longer(c(AveTemp_C, total_precip_mm), names_to = 'varb', values_to = 'value') %>%
  ggplot(aes(x = jDate, y = value, group = Year, colour = factor(Year))) +
  geom_line() +
  facet_wrap(~ varb)
# oh wow 
# does look like 2022 was colder... and had a lot more rain
# oh lmao there was a huge rain right on april 1
# man okay so it seems like the breaks for periods can be pretty influential

# What is the precip correlation like?
clim.corr.df %>%
  filter(grepl('prec\\.sum', vara)) %>%
  ggplot(aes(x = vara, y = varb)) +
  geom_tile(aes(alpha = abs(corr), fill = corr)) +
  scale_fill_gradient2(low = 'blue', high = 'red', mid = 'white', midpoint = 0) +
  theme(axis.text.x = element_text(angle = 90), legend.position = 'none')
# early season - weak pos. correlation with temperature (makes sense? warmer
#   temp means less snow, more rain?)
# growing season - super strong *negative* correlations with growing season
#   (and positive correlation with last freeze date...)
#   also some negative correlation with winter? autocorrelation?
# summer - strong-ish negative correlation with temperatures
# winter - strong-ish negative correlation with temperatures...

clim.corr.df %>%
  filter(grepl('prec\\.sum', vara)) %>%
  arrange(vara, desc(abs(corr))) %>%
  print(n = nrow(.))
# summer and winter precip are correlated...
# summer correlated with summer max temp and mean temp...
# winter correlated with winter min temperature, but not much else

clim.corr.df %>%
  filter(grepl('grow', vara) & grepl('grow', varb)) %>%
  filter(corr < 1) %>%
  print(n = nrow(.))
# Jesus all of these are super correlated
# guess I should only take one of these variables...
# mean temp, just because it makes sense?

clim.corr.df %>%
  filter(grepl('summer', vara) & grepl('summer', varb)) %>%
  filter(corr < 1) %>%
  arrange(desc(abs(corr))) %>%
  print(n = nrow(.))
# sould do min temp, precip? I hate min temp as a predictor though
# precip and... mean temp? still pretty correlated...

clim.corr.df %>%
  filter(grepl('early', vara) & grepl('early', varb)) %>%
  filter(corr < 1) %>%
  arrange(desc(abs(corr))) %>%
  print(n = nrow(.))
# precip and mean temp ~ 0.6!
# max temp... might as well include in here
# mean and min are super correlated... so maybe
# precip and min, which is ~0.3?

clim.corr.df %>%
  filter(grepl('winter', vara) & grepl('winter', varb)) %>%
  filter(corr < 1) %>%
  arrange(desc(abs(corr))) %>%
  print(n = nrow(.))
# temperatures all pretty correlated,
# precip correlated with max (makes sense)
# not correlated with mean or min?
# mean and precip?

clim.corr.df %>%
  filter(grepl('last', vara)) %>%
  filter(corr < 1) %>%
  arrange(desc(abs(corr))) %>%
  print(n = nrow(.))
# last freeze date is corr with a lot of growing season measures
# correlated with early measures too (makes sense)
# hmm...

clim.corr.df %>%
  filter(
    vara %in% c('grow_mean.temp', 'early_prec.sum', 'early_minnTemp', 
                'winter_meanTemp', 'summer_prec.sum', 'summer_meanTemp',
                'last.freeze'),
    varb %in% c('grow_mean.temp', 'early_prec.sum', 'early_minnTemp', 
                'winter_meanTemp', 'summer_prec.sum', 'summer_meanTemp',
                'last.freeze')
  ) %>%
  filter(vara != varb) %>%
  arrange(desc(abs(corr))) %>%
  print(n = nrow(.))
# these are considerably less correlated, which is good
# summer measures are correlated which is annoying
# last freeze date corelated with winter temperature, which kinda makes sense

clim.wide %>%
  ggplot(aes(x = summer_meanTemp, y = summer_prec.sum)) +
  geom_label(aes(label = Year))
# 2014 an outlier from an otherwise obvious trend
# ugh

clim.wide %>%
  ggplot(aes(x = winter_meanTemp, y = last.freeze)) +
  geom_label(aes(label = Year))
# interesting... much weaker pattern here

# (Old code when I was subsetting variables - note this was on a different
# version of df that was entirely in wide form and was missing period start
# dates)
# # Subset out relevant variables
# clim.subs = clim.summ %>%
#   select(Year, last.freeze, early_prec.sum, early_minnTemp,
#          grow_meanTemp, summer_meanTemp, winter_meanTemp)
# 
# head(clim.subs)
# tail(clim.subs)

# Export it
write.csv(
  clim.summ,
  row.names = FALSE,
  '01_data_cleaning/out/climate_summary.csv'
)

######## Other thoughts, etc.

# wonder if I should also do some de-trending? does it matter???
# z-scoring... would allow me to do polynomial tests...

##### What does a PCA on all of these look like?
# first... make sure that everything is in fact roughly normally distributed...

head(clim.summ)

clim.wide %>%
  mutate(last.freeze = as.numeric(last.freeze)) %>%
  pivot_longer(-Year, names_to = 'var', values_to = 'val') %>%
  filter(!is.na(val)) %>%
  ggplot() +
  geom_density(aes(x = val)) +
  facet_wrap(~ var, scales = 'free_x')
# meh... not in love with this

clim.wide %>%
  mutate(last.freeze = as.numeric(last.freeze)) %>%
  pivot_longer(-Year, names_to = 'var', values_to = 'val') %>%
  filter(!is.na(val)) %>%
  group_by(var) %>%
  mutate(zval = (val - mean(val)) / sd(val)) %>%
  ungroup() %>%
  ggplot() +
  geom_density(aes(x = zval)) +
  facet_wrap(~ var)

# kinda leptokurtotic? but could be much worse.

clim.pca = clim.wide[complete.cases(clim.wide),] %>%
  mutate(last.freeze = as.numeric(last.freeze)) %>%
  select(-Year) %>%
  prcomp(scale. = TRUE)

summary(clim.pca)
clim.pca$rotation %>% round(digits = 3)

clim.pca.loadings = clim.pca$rotation %>%
  as.data.frame() %>%
  mutate(varb = row.names(.)) %>%
  pivot_longer(-varb, names_to = 'pcaxis', values_to = 'loading') 

clim.pca.loadings %>%
  mutate(pcaxis = as.numeric(gsub('PC', '', pcaxis))) %>%
  ggplot(aes(x = varb, y = pcaxis, fill = loading)) +
  geom_tile() +
  scale_fill_gradient2(low = 'blue', high = 'red', mid = 'white', midpoint = 0) +
  scale_y_continuous(breaks = 1:10) +
  theme(axis.text.x = element_text(angle = 90))
