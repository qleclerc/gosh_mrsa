library(ppcor)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(lubridate)
library(reshape2)
options(dplyr.summarise.inform = FALSE)

staph_isolates = read.csv(here::here("Clean", "staph_isolates.csv")) %>%
  mutate(date = as_date(date))

#only keep patients with more than one sample
good_ids = staph_isolates %>%
  count(project_id) %>%
  filter(n > 1) %>%
  select(project_id) %>%
  pull

staph_isolates_profiles = staph_isolates %>%
  filter(project_id %in% good_ids)

interesting_samples = c()

#for each row, if the patient is the same, compare species
for(i in 1:(nrow(staph_isolates_profiles)-1)){
  
  if(staph_isolates_profiles$project_id[i] != staph_isolates_profiles$project_id[i+1]) next

  if(staph_isolates_profiles$SpeciesName[i] != staph_isolates_profiles$SpeciesName[i+1]){
    interesting_samples = c(interesting_samples,
                            staph_isolates_profiles$CultureIsoID[i],
                            staph_isolates_profiles$CultureIsoID[i+1])
  }
}

changing_profiles = staph_isolates_profiles %>%
  filter(CultureIsoID %in% interesting_samples)

profile_changes = data.frame(project_id = "test_id",
                             first_source = "test",
                             second_source = "test",
                             first_date = changing_profiles$date[1],
                             second_date = changing_profiles$date[1],
                             change = 0)

for(i in 1:(nrow(changing_profiles)-1)){
  
  if(changing_profiles$project_id[i] != changing_profiles$project_id[i+1]) next
  # if(changing_profiles$SpecimenType[i] != changing_profiles$SpecimenType[i+1]) next
  
  #rechecking as sometimes multiple samples with same lab id, but no change
  #eg sample 1 mssa, sample 2 mssa and mrsa, don't want to record mssa -> mssa change
  if(changing_profiles$SpeciesName[i] != changing_profiles$SpeciesName[i+1]){
    
    profile_changes = rbind(profile_changes,
                            data.frame(project_id = changing_profiles$project_id[i],
                                       first_source = changing_profiles$SpecimenType[i],
                                       second_source = changing_profiles$SpecimenType[i+1],
                                       first_date = changing_profiles$date[i],
                                       second_date = changing_profiles$date[i+1],
                                       change = changing_profiles$SpeciesName[i]))
  }
}

#remove the 1st test row used to setup the dataframe
profile_changes = profile_changes[-1,]

#this new dataset should contain all the s aureus changes between samples

#need to look whether events are during a single hospitalisation event, or if patients were discharged in the meantime
#if discharged in the meantime, reinfection by someone else would be an explanation
admissions = read.csv(here::here("Data", "combined_patient_ward_stays.csv")) %>%
  mutate(start_datetime = as_date(start_datetime)) %>%
  mutate(end_datetime = as_date(end_datetime)) %>%
  arrange(project_id, start_datetime) %>%
  filter(project_id %in% unique(profile_changes$project_id)) %>%
  filter(!is.na(start_datetime)) %>%
  filter(!is.na(end_datetime))

profile_changes$same_hosp = F
profile_changes$los = 0

for(i in 1:nrow(profile_changes)){
  
  admissions_i = admissions %>%
    filter(project_id == profile_changes$project_id[i])
  
  if(dim(admissions_i)[1] == 0) next
  
  for(j in 1:nrow(admissions_i)){
    
    if(profile_changes$first_date[i] >= admissions_i$start_datetime[j] &&
       profile_changes$second_date[i] <= admissions_i$end_datetime[j]){
      profile_changes$same_hosp[i] = T
      profile_changes$los[i] = admissions_i$end_datetime[j] - admissions_i$start_datetime[j]
    }
  }
  
}


#add column to say if varying profiles come from the same date
profile_changes = profile_changes %>%
  mutate(same_date = (first_date == second_date))

#add column to say if varying profiles come from the same sample source
profile_changes = profile_changes %>%
  mutate(same_source = (first_source == second_source))

#add column to check if could be nosocomial
#decided if there is a matching resistance profile in any patient within 30 days before the change is detected
profile_changes$possible_nos = 0
profile_changes$possible_nos_ward = ""

for(i in 1:(nrow(profile_changes))){
  
  #which hospitals admissions are recorded for patient i
  hosp_i = admissions %>%
    filter(project_id == profile_changes$project_id[i])
  
  if(nrow(hosp_i) == 0) next
  
  #which hospital admissions align with the change detected
  hosp_i = hosp_i %>%
    mutate(valid = (profile_changes$second_date[i] %within% interval(start_datetime, end_datetime))) %>%
    filter(valid == T)
  
  if(nrow(hosp_i) == 0) next
  
  #which other patients were also present in same ward as patient i within the last 30 days
  hosp_j = admissions %>%
    filter(ward_code == hosp_i$ward_code[1]) %>%
    mutate(valid = int_overlaps(interval(start_datetime, end_datetime), interval((profile_changes$second_date[i]-30), profile_changes$second_date[i]))) %>%
    filter(valid == T) %>%
    filter(project_id != profile_changes$project_id[i])
  
  if(nrow(hosp_j) == 0) next
  
  #were these patients positive for mssa or mrsa within 30 days before the change in patient i
  test_samples = staph_isolates_profiles %>%
    filter(project_id %in% unique(hosp_j$project_id)) %>%
    filter(date <= profile_changes$second_date[i] & date > (profile_changes$second_date[i]-30)) %>%
    group_by(project_id, SpeciesName) %>%
    summarise(n = n()) %>%
    ungroup %>%
    select(SpeciesName) %>%
    pull
  
  if(length(test_samples) == 0) next
  
  sample_i = profile_changes$change[i]
  
  profile_changes$possible_nos[i] = profile_changes$possible_nos[i] + sum(test_samples == sample_i)
  profile_changes$possible_nos_ward[i] = hosp_i$ward_code[1]
  
}

write.csv(profile_changes, here::here("Clean", "mrsa_mssa_diversity.csv"), row.names = F)
