# Script for assessing Lomatium demography data entry
# data entered May 2024 (SN)
# script init June 2024 (SN)

# --- Setup ---

# Packages
library(ggplot2)
library(dplyr)
library(tidyr)

# Clear namespace
rm(list = ls())

# -----------------------------------------------------------------------------------
# Code used to make data entry sheet
# do not rerun
# -----------------------------------------------------------------------------------
# a23 = read.csv('00_raw_data/lomatium_demography/2023_Lomatium_Demography_Data.csv')
# 
# head(a23)
# tail(a23)
# 
# names(a23)
# 
# a23 %>% arrange(Plot, Xcoor, Ycoor, Tag)
# 
# a23 %>% filter(grepl('ag\\spull', notes))
# 
# a24 = a23 %>% 
#   arrange(Plot, Xcoor, Ycoor) %>%
#   filter(is.na(tag.pulled)) %>%
#   select(-c(YrTag, problem.tag, X, X.1, X.2, X.3)) %>%
#   mutate(
#     across(
#       c(No.leaves, Leaf.length, No.umbels, Stalk_Height, umbel.diam, notes, tag.pulled),
#       function(x) NA
#     ),
#     Year = 2024
#   )
# 
# write.csv(
#   a24,
#   row.names = FALSE, na = '',
#   file = '00_raw_data/lomatium_demography/2024_Lomatium_Demography_Data_out.csv'
# )
# -----------------------------------------------------------------------------------

# --- Read in entered data ---

d24 = read.csv('00_raw_data/lomatium_demography/2024_Lomatium_Demography_Data.csv')

head(d24)

# Things to check
# - Umbel diameter counts match umbel number column - DONE
# - Alphanumeric columns - DONE
# - Leaf number range - DONE
# - Missing decimals - DONE
# - Leaf length range - DONE
# - Umbel number range - DONE
# - Tag pulled column has all tag pulls - DONE

# Alnum columns
str(d24)
# Looks good


# - Leaf number range
hist(d24$No.leaves)
table(d24$No.leaves)
d24 %>% filter(No.leaves > 10)
# Looks okay to me


# - Leaf length range
range(d24$Leaf.length, na.rm = TRUE)
# hmm...
hist(d24$Leaf.length)
d24 %>% filter(Leaf.length < 5 | Leaf.length > 30)
d24 %>% filter(is.na(Leaf.length), No.leaves > 0)
# (above line should all return notes for why umbel counts/measures are missing)


# - Umbel number range
table(d24$No.umbels, useNA = 'always')
d24 %>% filter(is.na(No.umbels), No.leaves > 0)
# (above line should all return notes for why umbel counts/measures are missing)


# - Tag pulled column
d24 %>% group_by(tag.pulled) %>% summarise(n = n())
# 30 tags pulled (...), everything else blank
d24 %>% filter(grepl('tag\\spull', notes), is.na(tag.pulled))
# (above line should return empty df)
d24 %>% filter(grepl('[Pp]ull', notes)) %>% arrange(tag.pulled)
# looks fine to me


# - Umbel counts matching
umbel.df = d24 %>%
  filter(!is.na(No.umbels), No.umbels > 1) %>%
  mutate(umbel.diam = gsub('\\s', '', umbel.diam)) %>%
  select(Plot, Tag, Xcoor, Ycoor, No.umbels, umbel.diam, notes) %>%
  # might need to change the 9 in the line below to some other integer
  separate(umbel.diam, into = paste0('u', 1:9), sep = ';', fill = 'right') %>%
  pivot_longer(starts_with('u'), names_to = 'u', values_to = 'diam') %>%
  select(-u) %>%
  filter(!is.na(diam))

# Umbel measurement enduing with a decimal
umbel.df %>% filter(grepl('\\.$', diam))
# (none - good!)

# Looking for places where the No.umbels col doesn't match the number of
# diameters provided
umbel.df %>%
  group_by(Plot, Tag, Xcoor) %>%
  filter(No.umbels != n()) # %>% View()
  
umbel.df %>%
  group_by(diam) %>% summarise(n.cases = n()) %>%
  mutate(diam.numeric = as.numeric(diam)) %>%
  print(n = nrow(.))
# the warning is fine - the character NA converted to NA type

# --- Look at NPs and NPNTs ---

pnp.leaves = d24 %>%
  mutate(
    p.status = case_when(
      grepl('[Nn][Pp]', notes) & !grepl('[Nn][Pp][Nn][Tt]', notes) ~ 'np',
      grepl('[Nn][Pp][Nn][Tt]', notes) ~ 'npnt',
      !grepl('[Nn][Pp]', notes) ~ 'p'
    ),
    leaves = case_when(
      is.na(No.leaves) ~ 'na',
      No.leaves > 0 ~ '>0',
      !No.leaves ~ '0'
    )
  )

pnp.leaves %>%
  group_by(p.status, leaves) %>%
  summarise(n = n())

# - All NP/NPNT have zero leaes (good)
# - Couple plants with zero leaves (hmm...)
# - Two cases where a plant has NA leaves

pnp.leaves %>% filter(p.status %in% 'p', !leaves %in% '>0')

# # First NPs only
# d24 %>%
#   filter(grepl('[Nn][Pp]', notes)) %>%
#   filter(!grepl('[Nn][Pp][Nn][Tt]', notes)) # %>% View()

# --- Some visualizations

# Leaf counts
d24 %>%
  filter(No.leaves > 0, !is.na(No.leaves)) %>%
  ggplot(aes(x = No.leaves)) +
  geom_histogram(binwidth = 0.5) +
  facet_wrap(~ Plot)
# - Plots 5, 6, 7, and 15 seem to have larger plants... 
#   (n.b. didn't find new plants here... may just be sampling effect)
# - Really a low fewer plants in 2
# - This isn't a random sample really so don't take too much from it...

# Leaf lengths
d24 %>%
  filter(!is.na(Leaf.length)) %>%
  ggplot(aes(x = Leaf.length)) + 
  geom_histogram(binwidth = 0.5) +
  facet_wrap(~ Plot)
# Shorter leaves in 1, 13, 14
#   (1 and 13 had lots of new plants, 14 had lots of new plants last year)
# Longer leaves in 7, 9, 15 (snowberry and grass)

# Umbel counts
d24 %>%
  filter(!is.na(No.umbels)) %>%
  ggplot(aes(x = No.umbels)) +
  geom_histogram(binwidth = 1) +
  facet_wrap(~ Plot)
# Of course, lots of new-ish plants in 1, 13, 14 so many zeros

# --- Reexport with tag removed column blanks converted to FALSE ---

d24 %>%
  mutate(tag.pulled = ifelse(is.na(tag.pulled), FALSE, tag.pulled)) # %>%
  # group_by(tag.pulled) %>% summarise(n = n())
  # write.csv()

