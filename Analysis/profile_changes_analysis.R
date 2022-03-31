
library(dplyr)
library(tidyverse)
library(ggplot2)
library(lubridate)
library(reshape2)

staph_isolates = read.csv(here::here("Clean", "staph_isolates.csv")) %>%
  mutate(date = as_date(date))

profile_changes = read.csv(here::here("Clean", "profile_changes.csv")) %>%
  mutate(first_date = as_date(first_date)) %>%
  mutate(second_date = as_date(second_date))

profile_changes_accurate = profile_changes %>%
  filter(same_hosp == T & same_sample == F & possible_nos == 0) 

table(profile_changes$possible_nos)

profile_changes %>%
  group_by(project_id, second_date) %>%
  summarise(nos = sum(possible_nos)) %>%
  mutate(nos = (nos>0)) %>%
  select(nos) %>%
  pull %>%
  table

length(unique(profile_changes_accurate$project_id))

##SUMMARY STATS

#total patients
length(unique(staph_isolates$project_id))
#total isolates
nrow(staph_isolates)
#total tests
sum(apply(staph_isolates[,6:59], 2, function(x) length(which(!is.na(x)))))

#total patients with at least one change
length(unique(changing_profiles$project_id))
#total isolates for patients with at least one change
staph_isolates %>% filter(project_id %in% unique(changing_profiles$project_id)) %>% nrow(.)
#total tests for patients with at least one changes
staph_isolates %>% filter(project_id %in% unique(changing_profiles$project_id)) %>% select(6:59) %>% apply(., 2, function(x) length(which(!is.na(x)))) %>% sum
#isolates with at least one change
length(unique(interesting_samples))/2
#changes
nrow(changing_profiles)
#work out prop of tests which are a change
#to prove that it's not just the most tested resistances that change
#tot_samples = apply(staph_isolates[,6:59], 2, function(x) length(which(!is.na(x))))
tot_samples = staph_isolates %>% filter(project_id %in% unique(profile_changes$project_id)) %>% select(6:59) %>% apply(., 2, function(x) length(which(!is.na(x))))
tot_samples[order(tot_samples)]
change_samples = table(profile_changes$antibiotic)
change_samples[order(change_samples)]
round(change_samples/tot_samples[names(change_samples)]*100, 2)
sum(change_samples)/sum(tot_samples)*100

#accurate patients with at least one change
length(unique(profile_changes_accurate$project_id))
#accurate changes
nrow(profile_changes_accurate)
#work out prop of tests which are a change
#to prove that it's not just the most tested resistances that change
#tot_samples = apply(staph_isolates[,6:59], 2, function(x) length(which(!is.na(x))))
tot_samples = staph_isolates %>% filter(project_id %in% unique(profile_changes_accurate$project_id)) %>% select(6:59) %>% apply(., 2, function(x) length(which(!is.na(x))))
tot_samples[order(tot_samples)]
change_samples = table(profile_changes_accurate$antibiotic)
change_samples[order(change_samples)]
round(change_samples/tot_samples[names(change_samples)]*100, 2)
sum(change_samples)/sum(tot_samples)*100

#gains/losses per abx
table(profile_changes_accurate$antibiotic, profile_changes_accurate$change)

table(profile_changes_accurate$source)/table(staph_isolates$SpecimenType)[names(table(profile_changes_accurate$source))]*100

#boxplot of durations between samples with a change
boxplot(int_length(interval(profile_changes_accurate$first_date, profile_changes_accurate$second_date))/60/60/24)


profile_changes_accurate %>%
  group_by(project_id, second_lab_id) %>%
  summarise(n = n()) %>%
  select(n) %>%
  pull %>%
  table