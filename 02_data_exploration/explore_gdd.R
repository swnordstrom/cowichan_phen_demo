# Assessment of cumulative GDD from year to year
# looking at data 2015 - 2023
# (future version may also include temperature)
# SN - init 4 apr 2024

library(ggplot2)
library(dplyr)
library(tidyr)

clim = read.csv('00_raw_data/climate/cgop_weather_daily_interp.csv') %>%
  separate(Date, into = c('Year', 'Month', 'Day'), remove = FALSE, sep = '-') %>%
  mutate(
    Date = as.Date(Date),
    across(c(Year, Month, Day), as.numeric)
  ) %>%
  filter(Year > 2014) %>%
  mutate(jday = as.numeric(Date) - as.numeric(as.Date(paste(Year, '01-01', sep = '-'), format = '%Y-%m-%d')))

jfma = clim %>% filter(Month < 5)

head(jfma)

jfma %>%
  mutate(Year = factor(Year)) %>%
  ggplot(aes(x = jday, y = AveTemp_C, group = Year)) +
  geom_line(aes(colour = Year))

nrow(jfma)

jfma.gdd = merge(
  jfma,
  expand.grid(jday = 0:max(jfma$jday), thresh = 2:10)
) %>%
  arrange(Year, thresh, jday) %>%
  mutate(gdd = ifelse(AveTemp_C > thresh, AveTemp_C - thresh, 0)) %>%
  group_by(Year, thresh) %>%
  mutate(cumul.gdd = cumsum(gdd)) %>%
  ungroup()

jfma.gdd %>%
  ggplot(aes(x = jday, y = cumul.gdd, group = thresh, colour = thresh)) +
  geom_line() +
  facet_wrap(~ Year)

jfma.gdd %>%
  mutate(Year = factor(Year)) %>%
  mutate(newmonth = as.numeric(Day < 2)) %>%
  ggplot() +
  geom_segment(
    data = jfma.gdd %>% filter(Year %in% 2021) %>% filter(Day < 2),
    aes(x = jday, y = 0, xend = jday, yend = max(jfma.gdd$cumul.gdd)),
    linetype = 2, linewidth = 0.1
  ) +
  geom_line(aes(x = jday, y = cumul.gdd, colour = Year)) +
  facet_wrap(~ thresh)

jfma.gdd %>%
  filter(Day < 2) %>%
  mutate(Month = factor(Month)) %>%
  ggplot(aes(x = Year, y = cumul.gdd, colour = Month)) +
  geom_point() +
  geom_line(aes(group = Month)) +
  facet_wrap(~ thresh, scales = 'free_y') +
  scale_x_continuous(breaks = 2014:2023) +
  labs(x = 'Year', y = 'Growing degree days on first day of month') +
  theme(
    axis.text.x = element_text(angle = 90),
    legend.position = 'bottom'
  )

ggsave('02_data_exploration/climate/gdd_compare.png', width = 8, height = 5)

### Look at GDD through the end of June

jfmamj = clim %>% filter(Month < 7)

head(jfmamj)

jfmamj %>%
  mutate(Year = factor(Year)) %>%
  ggplot(aes(x = jday, y = AveTemp_C, group = Year)) +
  geom_line(aes(colour = Year))

nrow(jfma)

jfmamj.gdd = merge(
  jfma,
  expand.grid(jday = 0:max(jfma$jday), thresh = 2:10)
) %>%
  arrange(Year, thresh, jday) %>%
  mutate(gdd = ifelse(AveTemp_C > thresh, AveTemp_C - thresh, 0)) %>%
  group_by(Year, thresh) %>%
  mutate(cumul.gdd = cumsum(gdd)) %>%
  ungroup()

jfmamj.gdd %>%
  ggplot(aes(x = jday, y = cumul.gdd, group = thresh, colour = thresh)) +
  geom_line() +
  facet_wrap(~ Year)

jfmamj.gdd %>%
  mutate(Year = factor(Year)) %>%
  ggplot() +
  geom_segment(
    data = jfmamj.gdd %>% filter(Year %in% 2021) %>% filter(Day < 2),
    aes(x = jday, y = 0, xend = jday, yend = max(jfmamj.gdd$cumul.gdd)),
    linetype = 2, linewidth = 0.1
  ) +
  geom_line(aes(x = jday, y = cumul.gdd, colour = Year)) +
  facet_wrap(~ thresh)

jfmamj.gdd %>%
  mutate(Month = factor(Month)) %>%
  arrange(Year, Month, desc(Day)) %>%
  distinct(Year, Month, thresh, .keep_all = TRUE) %>%
  ggplot(aes(x = Year, y = cumul.gdd, colour = Month)) +
  geom_point() +
  geom_line(aes(group = Month)) +
  facet_wrap(~ thresh, scales = 'free_y') +
  scale_x_continuous(breaks = 2014:2023) +
  labs(x = 'Year', y = 'Growing degree days on first day of month') +
  theme(
    axis.text.x = element_text(angle = 90),
    legend.position = 'bottom'
  )

# From this perspective 2022 was maybe a very cool (cold) year?
# and/or growing season
# 2017 also was below average
# At a lower threshold though it was a warm growing season

jfmamj.gdd %>%
  mutate(Month = factor(Month), Year = factor(Year)) %>%
  arrange(Year, Month, desc(Day)) %>%
  distinct(Year, Month, thresh, .keep_all = TRUE) %>%
  ggplot(aes(x = Month, y = cumul.gdd, colour = Year)) +
  geom_point() +
  geom_line(aes(group = Year)) +
  facet_wrap(~ thresh, scales = 'free_y') +
  labs(x = 'Year', y = 'Growing degree days on first day of month') +
  theme(legend.position = 'bottom')

# Also suggests that 2022 may have been an abnormally cold year
# (or at least growing season - may and june look very cold)

# (what about GDD on day of flowering?)

phen = read.csv('01_data_cleaning/out/phenology_all_cleaned.csv')

phen.gdd = phen %>%
  group_by(year, plot, init.doy) %>%
  summarise(n.init = n()) %>%
  merge(
    x = jfmamj.gdd %>% select(-c(total_precip_mm, minTemp_C, maxTemp_C, AveTemp_C)), 
    y = .,
    by.x = c('Year', 'jday'), by.y = c('year', 'init.doy')
  ) %>%
  merge(read.csv('00_raw_data/plot_treatments.csv')) %>%
  ungroup()

head(phen.gdd)

phen.gdd %>%
  mutate(Year = factor(Year)) %>%
  ggplot(aes(x = cumul.gdd, y = n.init)) +
  geom_point(aes(colour = trt, shape = Year)) +
  # geom_line(aes(colour = trt, linetype = Year)) +
  scale_colour_manual(values = c('black', 'red', 'blue')) +
  facet_wrap(~ thresh, scales = 'free_x')
# Not a useful figure

phen.gdd %>%
  mutate(Year = factor(Year)) %>%
  ggplot(aes(x = jday, y = cumul.gdd, size = n.init, colour = trt)) +
  geom_point() +
  scale_colour_manual(values = c('black', 'red', 'blue')) +
  facet_wrap(~ thresh, scales = 'free_y')

# Can't think of a good way to visualize this...
# Or how/if to use this for analysis...

phen.gdd %>%
  filter(thresh < 8) %>%
  uncount(n.init) %>%
  ggplot(aes(x = cumul.gdd)) +
  geom_density(aes(group = interaction(Year, trt), colour = trt)) +
  scale_colour_manual(values = c('black', 'red', 'blue')) +
  facet_wrap(~ thresh)
# ehhhhh not compelling

jfmamj.gdd %>%
  filter(Year %in% 2021:2023, jday < max(phen.gdd$jday)) %>%
  mutate(Month = factor(Month), Year = factor(Year)) %>%
  ggplot(aes(x = jday, y = cumul.gdd, colour = Year)) +
  geom_point(
    data = phen.gdd %>% mutate(Year = factor(Year)) %>% uncount(n.init),
    alpha = 0.01, shape = 21
  ) +
  geom_line(aes(group = Year)) +
  facet_wrap(~ thresh, scales = 'free_y') +
  theme(legend.position = 'bottom')
# Could maybe be useful...
