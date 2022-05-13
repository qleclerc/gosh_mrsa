library(ppcor)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(lubridate)
library(reshape2)


staph_isolates = read.csv(here::here("Clean", "staph_isolates.csv")) %>%
  mutate(date = as_date(date))

#replace S by 1, R by 2. This way, a mean value of 2 between two samples for a given abx implies a change in resistance
#that approach allows me to ignore NA values in the comparison
staph_isolates_sr = staph_isolates[,-c(1:5)]
staph_isolates_sr[staph_isolates_sr == "S"] = 1
staph_isolates_sr[staph_isolates_sr == "R"] = 3
staph_isolates_sr = apply(staph_isolates_sr, c(1,2), as.numeric)

staph_isolates[,-c(1:5)] = staph_isolates_sr

good_rows = which(rowSums(staph_isolates_sr, na.rm = T) > 0)

staph_isolates_profiles = staph_isolates %>%
  filter(row_number() %in% good_rows)

#only keep patients with more than one sample
good_ids = staph_isolates_profiles %>%
  count(project_id) %>%
  filter(n > 1) %>%
  select(project_id) %>%
  pull

staph_isolates_profiles = staph_isolates_profiles %>%
  filter(project_id %in% good_ids)


interesting_samples_mrsa = c()

#IMPORTANT: filter MRSA/MSSA!! Otherwise can just compare an mrsa with an mssa and miss the next mrsa isolate

staph_isolates_profiles_mrsa = staph_isolates_profiles %>%
  filter(SpeciesName == "Methicillin-Resistant Staphylococcus aureus")

#heavy lifting here
#for each row, if the patient, species and source are the same, compare the resistance profiles
#if there is any 2 in the comparison, that indicates a difference, and is recorded
for(i in 1:(nrow(staph_isolates_profiles_mrsa)-1)){
  
  if(staph_isolates_profiles_mrsa$project_id[i] != staph_isolates_profiles_mrsa$project_id[i+1]) next
  # if(staph_isolates_profiles$SpecimenType[i] != staph_isolates_profiles$SpecimenType[i+1]) next
  
  #extract the two profiles compared
  tt = staph_isolates_profiles_mrsa[c(i:(i+1)),-c(1:5)]
  #work out the means (neat trick: if comparing a value with NA, the mean will stay the same as the non-NA value... so, not either 1 or 3, so won't be flagged)
  test_means = apply(tt, 2, mean)
  #if any value is 2, that indicates a change in resistance profile, record it
  #note in theory there would be duplication if sample 1 is diff from sample 2, and sample 2 is diff from sample 3 (sample 2 would be recorded twice)
  if(any(test_means == 2, na.rm = T)) interesting_samples_mrsa = c(interesting_samples_mrsa,
                                                                   staph_isolates_profiles_mrsa$CultureIsoID[i],
                                                                   staph_isolates_profiles_mrsa$CultureIsoID[i+1])
}

interesting_samples_mssa = c()

staph_isolates_profiles_mssa = staph_isolates_profiles %>%
  filter(SpeciesName == "Methicillin-Susceptible Staphylococcus aureus")

#heavy lifting here
#for each row, if the patient, species and source are the same, compare the resistance profiles
#if there is any 2 in the comparison, that indicates a difference, and is recorded
for(i in 1:(nrow(staph_isolates_profiles_mssa)-1)){
  
  if(staph_isolates_profiles_mssa$project_id[i] != staph_isolates_profiles_mssa$project_id[i+1]) next
  # if(staph_isolates_profiles$SpecimenType[i] != staph_isolates_profiles$SpecimenType[i+1]) next
  
  #extract the two profiles compared
  tt = staph_isolates_profiles_mssa[c(i:(i+1)),-c(1:5)]
  #work out the means (neat trick: if comparing a value with NA, the mean will stay the same as the non-NA value... so, not either 1 or 3, so won't be flagged)
  test_means = apply(tt, 2, mean)
  #if any value is 2, that indicates a change in resistance profile, record it
  #note in theory there would be duplication if sample 1 is diff from sample 2, and sample 2 is diff from sample 3 (sample 2 would be recorded twice)
  if(any(test_means == 2, na.rm = T)) interesting_samples_mssa = c(interesting_samples_mssa,
                                                                   staph_isolates_profiles_mssa$CultureIsoID[i],
                                                                   staph_isolates_profiles_mssa$CultureIsoID[i+1])
}


#this is a list of samples, where: samples are from the same patient, the same source (eg nose), the same species (mrsa or mssa),
# but show a different resistant profile. No limits on dates currently, just have to be subsequent in the data (ie sample is different from
# the one immediately before it chronologically)
changing_profiles_mrsa = staph_isolates_profiles %>%
  filter(CultureIsoID %in% interesting_samples_mrsa)

changing_profiles_mssa = staph_isolates_profiles %>%
  filter(CultureIsoID %in% interesting_samples_mssa)


profile_changes = data.frame(project_id = "test_id",
                             species = "test", 
                             first_source = "test",
                             second_source = "test",
                             first_date = changing_profiles_mssa$date[1],
                             first_lab_id = "test",
                             second_date = changing_profiles_mssa$date[1],
                             second_lab_id = "test",
                             antibiotic = "test",
                             change = 0)

for(i in 1:(nrow(changing_profiles_mrsa)-1)){
  
  if(changing_profiles_mrsa$project_id[i] != changing_profiles_mrsa$project_id[i+1]) next
  # if(changing_profiles$SpecimenType[i] != changing_profiles$SpecimenType[i+1]) next
  
  #extract the two profiles compared
  tt = changing_profiles_mrsa[c(i:(i+1)),-c(1:5)]
  #work out the means (neat trick: if comparing a value with NA, the mean will stay the same as the non-NA value... so, not either 1 or 3, so won't be flagged)
  check = which(apply(tt, 2, mean) == 2)
  #needs a quick check because of the way this dataset was compiled
  #2 subsequent rows may not correspond to 2 different profiles!
  #eg if sample 1 and sample 2 are diff, but sample 2 and 3 are not, and sample 3 and 4 are, the dataset would still list: sample 1 2 3 4
  #yet, comparison between samples 2 and 3 in this loop would crash because they're identical
  if(length(check) == 0) next
  
  tt = tt %>%
    select(all_of(check))
  
  for(j in 1:ncol(tt)){
    profile_changes = rbind(profile_changes,
                            data.frame(project_id = changing_profiles_mrsa$project_id[i],
                                       species = changing_profiles_mrsa$SpeciesName[i],
                                       first_source = changing_profiles_mrsa$SpecimenType[i],
                                       second_source = changing_profiles_mrsa$SpecimenType[i+1],
                                       first_date = changing_profiles_mrsa$date[i],
                                       first_lab_id = changing_profiles_mrsa$CultureIsoID[i],
                                       second_date = changing_profiles_mrsa$date[i+1],
                                       second_lab_id = changing_profiles_mrsa$CultureIsoID[i+1],
                                       antibiotic = colnames(tt)[j],
                                       change = tt[1,j] - tt[2,j]))
  }
  
}

for(i in 1:(nrow(changing_profiles_mssa)-1)){
  
  if(changing_profiles_mssa$project_id[i] != changing_profiles_mssa$project_id[i+1]) next
  # if(changing_profiles$SpecimenType[i] != changing_profiles$SpecimenType[i+1]) next
  
  #extract the two profiles compared
  tt = changing_profiles_mssa[c(i:(i+1)),-c(1:5)]
  #work out the means (neat trick: if comparing a value with NA, the mean will stay the same as the non-NA value... so, not either 1 or 3, so won't be flagged)
  check = which(apply(tt, 2, mean) == 2)
  #needs a quick check because of the way this dataset was compiled
  #2 subsequent rows may not correspond to 2 different profiles!
  #eg if sample 1 and sample 2 are diff, but sample 2 and 3 are not, and sample 3 and 4 are, the dataset would still list: sample 1 2 3 4
  #yet, comparison between samples 2 and 3 in this loop would crash because they're identical
  if(length(check) == 0) next
  
  tt = tt %>%
    select(all_of(check))
  
  for(j in 1:ncol(tt)){
    profile_changes = rbind(profile_changes,
                            data.frame(project_id = changing_profiles_mssa$project_id[i],
                                       species = changing_profiles_mssa$SpeciesName[i],
                                       first_source = changing_profiles_mssa$SpecimenType[i],
                                       second_source = changing_profiles_mssa$SpecimenType[i+1],
                                       first_date = changing_profiles_mssa$date[i],
                                       first_lab_id = changing_profiles_mssa$CultureIsoID[i],
                                       second_date = changing_profiles_mssa$date[i+1],
                                       second_lab_id = changing_profiles_mssa$CultureIsoID[i+1],
                                       antibiotic = colnames(tt)[j],
                                       change = tt[1,j] - tt[2,j]))
  }
  
}

#remove the 1st test row used to setup the dataframe
profile_changes = profile_changes[-1,]

#this new dataset should contain all the resistance changes between isolates

equivalence_table = read.csv(here::here("Clean", "equivalence_table.csv"))[,-1]
colnames(equivalence_table)[1] = "antibiotic"
equivalence_table$antibiotic[equivalence_table$antibiotic == "Fusidic_acid"] = "Fucidin"
equivalence_table$antibiotic[equivalence_table$antibiotic == "Co-Trimoxazole"] = "Cotrimoxazole"

equivalence_table = rbind(equivalence_table,
                          c("antibiotic" = "Gent.Cipro", "class" = "Aminoglycoside"))
equivalence_table = rbind(equivalence_table,
                          c("antibiotic" = "Amik.Fluclox", "class" = "Aminoglycoside"))
equivalence_table = rbind(equivalence_table,
                          c("antibiotic" = "Piperacillin...Tazobactam", "class" = "Penicillin"))
equivalence_table = rbind(equivalence_table,
                          c("antibiotic" = "Pip.Taz.Cipro", "class" = "Penicillin"))
equivalence_table = rbind(equivalence_table,
                          c("antibiotic" = "Penicillin", "class" = "Penicillin"))
equivalence_table = rbind(equivalence_table,
                          c("antibiotic" = "Neomycin", "class" = "Aminoglycoside"))
equivalence_table = rbind(equivalence_table,
                          c("antibiotic" = "Cephradine", "class" = "Cephalosporin"))
equivalence_table = rbind(equivalence_table,
                          c("antibiotic" = "Syncercid", "class" = "Streptogramin"))
equivalence_table = rbind(equivalence_table,
                          c("antibiotic" = "Augmentin", "class" = "Penicillin"))
equivalence_table = rbind(equivalence_table,
                          c("antibiotic" = "Naladixic.Acid", "class" = "Fluoroquinolone"))

#add class of change
profile_changes = profile_changes %>%
  left_join(equivalence_table, "antibiotic") %>%
  filter(!is.na(project_id)) %>% 
  distinct()

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

#hereherehere!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
for(i in 1:nrow(profile_changes)){
  
  admissions_i = admissions %>%
    filter(project_id == profile_changes$project_id[i])
  
  if(dim(admissions_i)[1] == 0) next
  
  for(j in 1:nrow(admissions_i)){
    
    if(profile_changes$first_date[i] >= admissions_i$start_datetime[j] && profile_changes$second_date[i] <= admissions_i$end_datetime[j]) profile_changes$same_hosp[i] = T
    
  }
  
}


#add column to say if varying profiles come from the same single sample (matching times)
profile_changes = profile_changes %>%
  mutate(same_date = (first_date == second_date))

profile_changes = profile_changes %>%
  mutate(same_source = (first_source == second_source))

#add column to check if could be nosocomial
#decided if there is a matching resistance profile in any patient within 1 week before the change is detected

#check if patient was in ward with anyone else with s aureus (null hypothesis)

staph_isolates_profiles[is.na(staph_isolates_profiles)] = 0
profile_changes$possible_nos = 0
for(i in 1:(nrow(profile_changes))){
  
  hosp_i = admissions %>%
    filter(project_id == profile_changes$project_id[i])
  
  if(nrow(hosp_i) == 0) next
  
  hosp_i = hosp_i %>%
    mutate(valid = (profile_changes$second_date[i] %within% interval(start_datetime, end_datetime))) %>%
    filter(valid == T)
  
  if(nrow(hosp_i) == 0) next
  
  #only take patients which were hospitalised in the same ward within the 7 days before the 2nd sample was taken for the patient of interest
  hosp_j = admissions %>%
    filter(ward_code == hosp_i$ward_code[1]) %>%
    mutate(valid = int_overlaps(interval(start_datetime, end_datetime), interval((profile_changes$second_date[i]-30), profile_changes$second_date[i]))) %>%
    filter(valid == T) %>%
    filter(project_id != profile_changes$project_id[i])
  
  if(nrow(hosp_j) == 0) next
  
  #from those patients, only look at any sample taken 30 days before the 2nd sample was taken for the patient of interest
  test_samples = staph_isolates_profiles %>%
    filter(project_id %in% unique(hosp_j$project_id)) %>%
    filter(date <= profile_changes$second_date[i] & date > (profile_changes$second_date[i]-30)) %>%
    select(-c(2:5))
  
  if(nrow(test_samples) == 0) next
  
  sample_i = staph_isolates_profiles %>%
    filter(CultureIsoID == profile_changes$second_lab_id[i]) %>%
    .[1,] %>%
    select(-c(1:5))
  
  already_checked_patient = c()
  
  for(j in 1:nrow(test_samples)){
    if(sum(test_samples[j,-1] - sample_i, na.rm = T) == 0 & !(test_samples[j,-1] %in% already_checked_patient)){
      #extra check to avoid double counting (otherwise if a patient is tested twice and matches target twice, would be double counted)
      already_checked_patient = c(already_checked_patient, test_samples[j,-1])
      profile_changes$possible_nos[i] = profile_changes$possible_nos[i] + 1
    }
  }
  
}


#add column to check if there was antibiotic usage
antibio_data = read.csv(here::here("Clean", "antibio_data.csv")) %>%
  mutate(date = as_date(start_datetime))

profile_changes$antibiotic_use = FALSE
profile_changes$any_antibiotic = FALSE

for(i in 1:(nrow(profile_changes))){
  
  antibio_i = antibio_data %>%
    filter(project_id == profile_changes$project_id[i])
  
  if(nrow(antibio_i) == 0) next
  
  antibio_any = antibio_i %>%
    mutate(valid = (profile_changes$second_date[i] %within% interval(date-30, date))) %>%
    filter(valid == T)
  
  if(nrow(antibio_any) != 0) profile_changes$any_antibiotic[i] = TRUE
  
  antibio_name = profile_changes$class[i]
  
  antibio_i = antibio_i %>%
    filter(class == antibio_name)
  
  if(nrow(antibio_i) == 0) next
  
  antibio_i = antibio_i %>%
    mutate(valid = (profile_changes$second_date[i] %within% interval(date-30, date))) %>%
    filter(valid == T)
  
  if(nrow(antibio_i) == 0) next
  
  profile_changes$antibiotic_use[i] = TRUE
  
}


profile_changes = profile_changes %>%
  mutate(change = replace(change, change == -2, "R")) %>%
  mutate(change = replace(change, change == "2", "S"))

write.csv(profile_changes, here::here("Clean", "res_profiles_diversity.csv"), row.names = F)
