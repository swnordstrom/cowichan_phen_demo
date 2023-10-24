### Simple script to combine all of the "multiple-umbel" datasets into one csv
### SN - init 23 Oct 2023

library(ggplot2)
library(dplyr)
library(tidyr)

rm(list = ls())

data.dir = '00_raw_data/lomatium_demography/'

multi.umbel.list = paste0(data.dir, grep('[Uu]mble', dir(data.dir), value = TRUE)) %>%
  lapply(FUN = read.csv)

sapply(multi.umbel.list, names)

shared.cols = Reduce(intersect, lapply(multi.umbel.list, names))

multi.umbel.raw = multi.umbel.list %>%
  lapply(FUN = function(df) df[,shared.cols]) %>%
  do.call(what = rbind)

multi.umbel.raw

str(multi.umbel.raw)

# all appear correct/okay
# although note column... probably missing notes from certain years

# can I get notes in here?
# or should I just collapse all of the columns after 'names'?
# (tidyverse is annoying about doing this... not going to worry about it yet)

# I should add a plant ID though

multi.umbel.proc = multi.umbel.raw %>%
  mutate(plantid = paste0(Tag, "_", Plot, "_", Xcoor, Ycoor))

write.csv(
  multi.umbel.proc,
  file = '01_data_cleaning/out/multi_umbel_cleaned.csv',
  row.names = FALSE
)
