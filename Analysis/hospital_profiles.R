library(ppcor)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(lubridate)
library(reshape2)
library(FSA)

staph_isolates = read.csv(here::here("Clean", "staph_isolates.csv")) %>%
  mutate(date = as_date(date))

# #replace S by 1, R by 2. This way, a mean value of 2 between two samples for a given abx implies a change in resistance
# #that approach allows me to ignore NA values in the comparison
# staph_isolates_sr = staph_isolates[,-c(1:5)]
# staph_isolates_sr[staph_isolates_sr == "S"] = 1
# staph_isolates_sr[staph_isolates_sr == "R"] = 3
# 
# staph_isolates_sr = apply(staph_isolates_sr, c(1,2), as.numeric)
# 
# staph_isolates[,-c(1:5)] = staph_isolates_sr
# 
# #assumption: only keep antibiotics with more than 50 results (S or I)
# good_cols = which(colSums(staph_isolates_sr, na.rm = T) > 50)
# good_rows = which(rowSums(staph_isolates_sr, na.rm = T) > 0)

#grouping to count how many profiles occur per year
staph_isolates_profiles = staph_isolates %>%
  mutate(date = floor_date(date, "year")) %>%
  select(-c(1:3)) %>%
  group_by_all() %>%
  count() %>%
  ungroup

#but this treats NA as a difference (ie 1 NA 1 is different from 1 NA NA)
#problem: can't apply a second filtering
#because e.g. NA 3 NA could equally belong to NA 3 1 or NA 3 3, even though these 2 profiles are different
#solution is to just exclude profiles with less than x resistances tested

staph_isolates_profiles$n_tested = apply(staph_isolates_profiles[,-c(1,2,ncol(staph_isolates_profiles))],
                                         1, function(x) sum(!(is.na(x))))

staph_isolates_profiles %>%
  ggplot() +
  geom_histogram(aes(n_tested, fill = SpeciesName), position = "identity", alpha = 0.5)

staph_isolates_profiles %>%
  filter(n_tested >= 13) %>%
  group_by(date, SpeciesName) %>%
  summarise(n_profiles = n()) %>%
  ggplot() +
  geom_line(aes(date, n_profiles, colour = SpeciesName)) +
  theme_bw() +
  theme(legend.position = "bottom")

staph_isolates_profiles %>% 
  filter(n_tested >= 13) %>%
  group_by(date, SpeciesName) %>% 
  filter(n == max(n)) %>% 
  View()

#some further filtering to identify "most common profile" by removing columns always NAs
#then remove columns always equal to S or R
