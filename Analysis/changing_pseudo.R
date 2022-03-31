library(ppcor)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(lubridate)
library(reshape2)


#using antibiograms as proxy to identify mrsa
antibiograms = read.csv(here::here("Data" ,"labanalysis_patient_microbiology_positives_sensitivitys_pivot.csv"))

#assumption: species names in that list only are s aureus
pseudo_antibio = antibiograms %>%
  filter(SpeciesName %in% c("Pseudomonas aeruginosa"))

pseudo_antibio[pseudo_antibio == ""] = NA

pseudo_antibio$date = ymd_hms(pseudo_antibio$start_datetime)

good_cols = colnames(pseudo_antibio)[apply(pseudo_antibio, 2, function(x) !(all(is.na(x))))]
good_cols = good_cols[!grepl("MIC", good_cols)]

pseudo_antibio = pseudo_antibio %>%
  select(good_cols) %>%
  select(1, 2, 7, 84, c(13:83)) %>%
  arrange(project_id, date) %>%
  filter(date < "2021-12-01")

#some general plots
#baseline incidence
pseudo_antibio %>%
  mutate(date = floor_date(date, "month")) %>%
  group_by(project_id, date) %>%
  summarise(n = 1) %>%
  ungroup() %>%
  group_by(date) %>%
  summarise(n = sum(n)) %>%
  mutate(prop = n/sum(n)) %>%
  ggplot() +
  geom_line(aes(date, n)) +
  theme_bw()


#replace S by 1, R by 2. This way, a mean value of 2 between two samples for a given abx implies a change in resistance
#that approach allows me to ignore NA values in the comparison
pseudo_antibio_sr = pseudo_antibio[,-c(1:4)]
pseudo_antibio_sr[pseudo_antibio_sr == "S"] = 1
pseudo_antibio_sr[pseudo_antibio_sr == "R"] = 3
pseudo_antibio_sr = apply(pseudo_antibio_sr, c(1,2), as.numeric)

pseudo_antibio[,-c(1:4)] = pseudo_antibio_sr

#assumption: only keep antibiotics with more than 50 results (S or I)
good_cols = which(colSums(pseudo_antibio_sr, na.rm = T) > 50)
good_rows = which(rowSums(pseudo_antibio_sr, na.rm = T) > 0)

pseudo_antibio_profiles = pseudo_antibio %>%
  select(c(1:4,(good_cols+4)))

pseudo_antibio_profiles = pseudo_antibio_profiles %>%
  filter(row_number() %in% good_rows)

#only keep patients with more than one sample
good_ids = pseudo_antibio_profiles %>%
  group_by(project_id) %>%
  summarise(n = n()) %>%
  filter(n > 1) %>%
  select(project_id) %>%
  pull

pseudo_antibio_profiles = pseudo_antibio_profiles %>%
  filter(project_id %in% good_ids)

#assumption: if specimen source is NA, assume it is "*Unspecified" as that value already existed at baseline
pseudo_antibio_profiles$SpecimenSource[is.na(pseudo_antibio_profiles$SpecimenSource)] = "*Unspecified"

interesting_samples = c()

#heavy lifting here
#for each row, if the patient, species and source are the same, compare the resistance profiles
#if there is any 2 in the comparison, that indicates a difference, and is recorded
for(i in 1:(nrow(pseudo_antibio_profiles)-1)){
  
  if(pseudo_antibio_profiles$project_id[i] != pseudo_antibio_profiles$project_id[i+1]) next
  if(pseudo_antibio_profiles$SpecimenSource[i] != pseudo_antibio_profiles$SpecimenSource[i+1]) next
  
  #extract the two profiles compared
  tt = pseudo_antibio_profiles[c(i:(i+1)),-c(1:4)]
  #work out the means (neat trick: if comparing a value with NA, the mean will stay the same as the non-NA value... so, not either 1 or 3, so won't be flagged)
  test_means = apply(tt, 2, mean)
  #if any value is 2, that indicates a change in resistance profile, record it
  #note in theory there would be duplication if sample 1 is diff from sample 2, and sample 2 is diff from sample 3 (sample 2 would be recorded twice)
  if(any(test_means == 2, na.rm = T)) interesting_samples = c(interesting_samples,
                                                              pseudo_antibio_profiles$project_lab_id[i],
                                                              pseudo_antibio_profiles$project_lab_id[i+1])
}


#this is a list of samples, where: samples are from the same patient, the same source (eg nose), the same species (mrsa or mssa),
# but show a different resistant profile. No limits on dates currently, just have to be subsequent in the data (ie sample is different from
# the one immediately before it chronologically)
changing_profiles = pseudo_antibio_profiles %>%
  filter(project_lab_id %in% interesting_samples)




profile_changes = data.frame(project_id = "test_id",
                             source = "test",
                             first_date = changing_profiles$date[1],
                             first_lab_id = "test",
                             second_date = changing_profiles$date[1],
                             second_lab_id = "test",
                             antibiotic = "test",
                             change = 0)

for(i in 1:(nrow(changing_profiles)-1)){
  
  if(changing_profiles$project_id[i] != changing_profiles$project_id[i+1]) next
  if(changing_profiles$SpecimenSource[i] != changing_profiles$SpecimenSource[i+1]) next
  
  #extract the two profiles compared
  tt = changing_profiles[c(i:(i+1)),-c(1:4)]
  #work out the means (neat trick: if comparing a value with NA, the mean will stay the same as the non-NA value... so, not either 1 or 3, so won't be flagged)
  check = which(apply(tt, 2, mean) == 2)
  #needs a quick check because of the way this dataset was compiled
  #2 subsequent rows may not correspond to 2 different profiles!
  #eg if sample 1 and sample 2 are diff, but sample 2 and 3 are not, and sample 3 and 4 are, the dataset would still list: sample 1 2 3 4
  #yet, comparison between samples 2 and 3 in this loop would crash because they're identical
  if(length(check) == 0) next
  
  tt = tt %>%
    select(check)
  
  for(j in 1:ncol(tt)){
    profile_changes = rbind(profile_changes,
                            data.frame(project_id = changing_profiles$project_id[i],
                                       source = changing_profiles$SpecimenSource[i],
                                       first_date = changing_profiles$date[i],
                                       first_lab_id = changing_profiles$project_lab_id[i],
                                       second_date = changing_profiles$date[i+1],
                                       second_lab_id = changing_profiles$project_lab_id[i+1],
                                       antibiotic = colnames(tt)[j],
                                       change = tt[1,j] - tt[2,j]))
  }
  
}

#remove the 1st test row used to setup the dataframe
profile_changes = profile_changes[-1,]

#this new dataset should contain all the resistance changes between isolates

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
  mutate(same_sample = (first_date == second_date))


#add column to check if could be nosocomial
#decided if there is a matching resistance profile in any patient within 1 week before the change is detected

#check if patient was in ward with anyone else with s aureus (null hypothesis)

pseudo_antibio_profiles[is.na(pseudo_antibio_profiles)] = 0
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
    mutate(valid = int_overlaps(interval(start_datetime, end_datetime), interval((profile_changes$second_date[i]-(7*24*60*60)), profile_changes$second_date[i]))) %>%
    filter(valid == T) %>%
    filter(project_id != profile_changes$project_id[i])
  
  if(nrow(hosp_j) == 0) next
  
  #from those patients, only look at any sample taken 7 days before the 2nd sample was taken for the patient of interest
  test_samples = pseudo_antibio_profiles %>%
    filter(project_id %in% unique(hosp_j$project_id)) %>%
    filter(date <= profile_changes$second_date[i] & date > (profile_changes$second_date[i]-(7*24*60*60))) %>%
    select(-c(2:5))
  
  if(nrow(test_samples) == 0) next
  
  sample_i = pseudo_antibio_profiles %>%
    filter(project_lab_id == profile_changes$second_lab_id[i]) %>%
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

profile_changes_accurate = profile_changes %>%
  filter(same_hosp == T & same_sample == F & possible_nos == 0) 

#change 2/-2 to S/R
profile_changes_accurate = profile_changes_accurate %>%
  mutate(change = replace(change, change == -2, "R")) %>%
  mutate(change = replace(change, change == "2", "S"))

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
length(unique(pseudo_antibio$project_id))
#total isolates
nrow(pseudo_antibio)
#total tests
sum(apply(pseudo_antibio[,5:75], 2, function(x) length(which(!is.na(x)))))

#total patients with at least one change
length(unique(changing_profiles$project_id))
#total isolates for patients with at least one change
pseudo_antibio %>% filter(project_id %in% unique(changing_profiles$project_id)) %>% nrow(.)
#total tests for patients with at least one changes
pseudo_antibio %>% filter(project_id %in% unique(changing_profiles$project_id)) %>% select(5:75) %>% apply(., 2, function(x) length(which(!is.na(x)))) %>% sum
#isolates with at least one change
length(unique(interesting_samples))/2
#changes
nrow(changing_profiles)
#work out prop of tests which are a change
#to prove that it's not just the most tested resistances that change
#tot_samples = apply(pseudo_antibio[,6:59], 2, function(x) length(which(!is.na(x))))
tot_samples = pseudo_antibio %>% filter(project_id %in% unique(profile_changes$project_id)) %>% select(5:75) %>% apply(., 2, function(x) length(which(!is.na(x))))
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
#tot_samples = apply(pseudo_antibio[,6:59], 2, function(x) length(which(!is.na(x))))
tot_samples = pseudo_antibio %>% filter(project_id %in% unique(profile_changes_accurate$project_id)) %>% select(5:75) %>% apply(., 2, function(x) length(which(!is.na(x))))
tot_samples[order(tot_samples)]
change_samples = table(profile_changes_accurate$antibiotic)
change_samples[order(change_samples)]
round(change_samples/tot_samples[names(change_samples)]*100, 2)
sum(change_samples)/sum(tot_samples)*100

#gains/losses per abx
table(profile_changes_accurate$antibiotic, profile_changes_accurate$change)

table(profile_changes_accurate$source)

#boxplot of durations between samples with a change
boxplot(int_length(interval(profile_changes_accurate$first_date, profile_changes_accurate$second_date))/60/60/24)
hist(int_length(interval(profile_changes_accurate$first_date, profile_changes_accurate$second_date))/60/60/24, breaks = 40)

profile_changes_accurate %>%
  group_by(project_id, second_lab_id) %>%
  summarise(n = n()) %>%
  select(n) %>%
  pull %>%
  table
