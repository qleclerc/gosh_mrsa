library(ppcor)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(lubridate)
library(forecast)
library(reshape2)


#demographics = xap.read_table("caboodle_patient_demographics")
#demographics dataset: ID, DoB, DoD, death Y/N, ethnicity

#diagnoses = xap.read_table("caboodle_patient_diagnoses")
#diagnostics (not for all patients): ID, time, code, name of diagnostic, further information

#admissions = xap.read_table("caboodle_patient_hospital_admissions")
#admissions (not for all patients): ID, time

#problems = xap.read_table("caboodle_patient_problemlist")
#note that's all the problems they ever had, not reason for admission
#useful problems$diag_name
#inf_probs = c("Staphylococcal infection, unspecified site",
#"Bacterial infection, unspecified",
#"Other bacterial infections of unspecified site",
#"Local infection of skin and subcutaneous tissue, unspecified",
#"Post-traumatic wound infection, not elsewhere classified")
#p_f = problems %>% filter(diag_name %in% inf_probs)

medications = read.csv(here::here("Data", "caboodle_patient_selected_medication_admins_all.csv"))
medications2 = read.csv(here::here("Data", "datalake_patient_selected_medication_admins_all.csv"))
#full details of medications received: time, dose, name, seq (to see order of meds when a single med is given more than once)
#medications_intended = xap.read_table("caboodle_patient_selected_medication_orders_all")
#like above, but the theory (ie the orders for medication, which may/may not have then been given to patients)

medications_short = medications %>%
  mutate(date = as_date(start_datetime)) %>%
  mutate(date = floor_date(date, "day")) %>%
  filter(grepl("CHLORHEX", drug_name) | grepl("MUPIRO", drug_name)) %>%
  mutate(drug_name=replace(drug_name, grepl("CHLORHEX", drug_name), "Chlorhexidine")) %>%
  mutate(drug_name=replace(drug_name, grepl("MUPIRO", drug_name), "Mupirocin")) %>%
  select(project_id, date, drug_name, route_name) %>%
  filter(route_name %in% c("Topical", "Nasal"))


medications2_short = medications2 %>%
  mutate(date = as_date(start_datetime)) %>%
  mutate(date = floor_date(date, "day")) %>%
  filter(grepl("CHLORHEX", drug_name) | grepl("MUPIRO", drug_name)) %>%
  mutate(drug_name=replace(drug_name, grepl("CHLORHEX", drug_name), "Chlorhexidine")) %>%
  mutate(drug_name=replace(drug_name, grepl("MUPIRO", drug_name), "Mupirocin")) %>%
  select(project_id, date, drug_name, route_name) %>%
  filter(route_name %in% c("Topical", "Intranasal", "Nostril (left)", "Nostrils (both)"))

rbind(medications_short, medications2_short) %>%
  mutate(date = floor_date(date, "month")) %>%
  group_by(date, drug_name) %>%
  summarise(n = n()) %>%
  ggplot() +
  geom_line(aes(x = date, y = n, colour = drug_name)) +
  theme_bw()

#meds_simple = xap.read_table("caboodle_patient_selected_medication_comeds_all")
#simplified list of meds, with patient ID, day, and all meds received that day (separated by commas)

#surgery = xap.read_table("caboodle_patient_theatre_list")
#list of surgery operations: time, planned/emergency, procedure name, specialty name

#stays = xap.read_table("caboodle_patient_ward_stays")
#ward stays: entry and exit time, LoS, ward name, icu ward boolean (1 or 0)
#ICU wards are seahorse (Paediatric Intensive Care), flamingo (cardiac critical care), alligator (Cardiac intensive and high dependency care), dolphin (Neonatal Intensive Care)


#admissions_all = xap.read_table("combined_patient_hospital_admissions")
#all hospital admissions for patients since 2000
#admissions_all %>%
#mutate(date = format(start_datetime, format = "%Y")) %>%
#group_by(date) %>%
#summarise(n = n()) %>%
#ungroup %>%
#ggplot() +
#geom_point(aes(x = date, y = n))

#stays_all = xap.read_table("combined_patient_ward_stays")
#all patient ward stays: entry and exit time, ward name, icu ward boolean (1 or 0)
#note ward names likely changed through time...
#and no LoS estimate, need to calculate based on dates

#datalake_diagnoses = xap.read_table("datalake_patient_diagnoses")
#diagnostics (not for all patients): ID, time, code, name of diagnostic, further information

#microbiology = xap.read_table("labanalysis_patient_microbiology_all")



#using antibiograms as proxy to identify mrsa
antibiograms = read.csv(here::here("Data", "labanalysis_patient_microbiology_positives_sensitivitys_pivot.csv"))

staph_antibio = antibiograms %>%
  filter(SpeciesName %in% c("Staphylococcus aureus", "Staphylococcus sp.", "Staphylococcus sp", "Methicillin-Resistant Staphylococcus aureus")) %>%
  mutate(date = ymd_hms(start_datetime))

staph_antibio$SpeciesName[staph_antibio$Flucloxacillin == "R"] = "Methicillin-Resistant Staphylococcus aureus"
staph_antibio$SpeciesName[staph_antibio$Flucloxacillin == "I"] = "Methicillin-Resistant Staphylococcus aureus"
staph_antibio$SpeciesName[staph_antibio$Oxacillin == "R"] = "Methicillin-Resistant Staphylococcus aureus"
staph_antibio$SpeciesName[staph_antibio$Oxacillin == "I"] = "Methicillin-Resistant Staphylococcus aureus"
staph_antibio$SpeciesName[staph_antibio$Cefoxitin == "R"] = "Methicillin-Resistant Staphylococcus aureus"
staph_antibio$SpeciesName[staph_antibio$Cefoxitin == "I"] = "Methicillin-Resistant Staphylococcus aureus"
staph_antibio$SpeciesName[staph_antibio$SpeciesName != "Methicillin-Resistant Staphylococcus aureus"] = "Methicillin-Susceptible Staphylococcus aureus"


staph_antibio_short = staph_antibio %>%
  mutate(date = floor_date(date, "day")) %>%
  group_by(project_id, date, SpeciesName) %>%
  mutate(n = n()) %>%
  select(project_id, date, SpeciesName, n) %>%
  arrange(project_id, date, SpeciesName)

length(unique(staph_antibio_short$project_id))
sum(staph_antibio_short$n)

staph_antibio_short %>%
  mutate(date = floor_date(date, "month")) %>%
  group_by(project_id, date) %>%
  summarise(inf = length(unique(SpeciesName))) %>%
  filter(inf > 1) %>%
  #select(inf) %>%
  #pull() %>%
  #table(.)
  select(project_id) %>%
  pull() %>%
  unique(.) %>%
  length(.)

staph_antibio_distinct = staph_antibio_short %>%
  distinct() %>%
  mutate(n = 1) %>%
  group_by(date, SpeciesName) %>%
  summarise(cs = sum(n))

staph_antibio_distinct %>%
  #mutate(date = year(date)) %>%
  mutate(date = floor_date(date, "month")) %>%
  #filter(date > as_date("2014-01-01")) %>%
  group_by(date, SpeciesName) %>%
  summarise(cs = sum(cs)) %>%
  ggplot(aes(date, cs, colour = SpeciesName)) +
  geom_point() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_smooth(method = "lm")

staph_antibio_distinct %>%
  #mutate(date = year(date)) %>%
  mutate(date = floor_date(date, "month")) %>%
  #filter(date > as_date("2014-01-01")) %>%
  group_by(date, SpeciesName) %>%
  summarise(cs = sum(cs)) %>%
  mutate(freq = cs/sum(cs)) %>%
  ggplot(aes(date, freq, colour = SpeciesName)) +
  geom_point() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


staph_antibio_clean = staph_antibio %>%
  select(good_cols) %>%
  select(c(3, 13:66)) %>%
  melt(.,id.vars = "start_datetime") %>%
  filter(!is.na(value)) %>%
  filter(value %in% c("S", "R")) %>%
  mutate(date = floor_date(start_datetime, "month")) %>%
  arrange(date) %>%
  group_by(date, variable, value) %>%
  summarise(n = n()) %>%
  mutate(prop = n/sum(n)) %>%
  mutate(variable = as.character(variable))


test_variables = unique(staph_antibio_clean$variable)

for(i in c(1,10,19,28,37,44)){
  
  pp = staph_antibio_clean %>%
    filter(value == "R") %>%
    filter(variable %in% test_variables[(i):(i+8)]) %>%
    ggplot() +
    geom_line(aes(date, prop, colour = variable)) +
    theme_bw()
  
  plot(pp)
  
}

staph_antibio_clean %>%
  filter(value == "R") %>%
  filter(variable %in% c("Penicillin", "Tetracycline", "Trimethoprim")) %>%
  ggplot() +
  geom_line(aes(date, prop, colour = variable)) +
  scale_y_continuous(limits = c(0,1)) +
  theme_bw()


abx_classes = data.frame(variable = c("Amik.Fluclox", "Amikacin", "Ampicillin", 
                                      "Augmentin", "Cefotaxime", "Ceftazidime", "Cephradine", "Chloramphenicol", 
                                      "Ciprofloxacin", "Colistin", "Erythromycin", "Flucloxacillin", 
                                      "Fucidin", "Gent.Ceftaz", "Gent.Cipro", "Gent.Pip.Taz", "Gentamicin", 
                                      "Mupirocin", "Naladixic.Acid", "Neomycin", "Nitrofurantoin", 
                                      "Penicillin", "Piperacillin...Tazobactam", "Pip.Taz.Amik", "Pip.Taz.Cipro", 
                                      "Rifampicin", "Teicoplanin", "Tetracycline", "Trimethoprim", 
                                      "Vancomycin", "Sulphonamide", "Oxacillin", "Pristinomycin", "Tobramycin", 
                                      "Mecillinam", "Syncercid", "Ceftriaxone", "Septrin", "Cefuroxime", 
                                      "Linezolid", "Clindamycin", "Timentin", "Imipenem", "Meropenem", 
                                      "Tigecycline", "Daptomycin", "Doxycycline", "Minocycline", "Fosfomycin", 
                                      "Metronidazole", "Cefoxitin", "Co.Trimoxazole..Septrin."),
                         class = c("Mixed", 
                                   "Aminoglycosides", "Penicillins", "Penicillins", "Cephalosporins", 
                                   "Cephalosporins", "Cephalosporins", "Chloramphenicol", "Fluoroquinolones", 
                                   "Polypetides", "Macrolides", "Penicillins", "Fucidin", "Mixed", 
                                   "Mixed", "Mixed", "Aminoglycosides", "Mupirocin", "Fluoroquinolones", 
                                   "Aminoglycosides", "Nitrofurans", "Penicillins", "Mixed", "Mixed", 
                                   "Mixed", "Antimycobacterials", "Glycopeptides", "Tetracyclines", 
                                   "Trimethoprim", "Glycopeptides", "Sulphonamides", "Penicillins", "Pristinomycin", 
                                   "Aminoglycosides", "Penicillins", "Syncercid", "Cephalosporins", 
                                   "Sulfonamides", "Cephalosporins", "Oxazolidinones", "Lincosamides", 
                                   "Mixed", "Carbapenems", "Carbapenems", "Tigecycline", "Lipopetides", 
                                   "Tetracyclines", "Tetracyclines", "Fosfomycin", "Metronidazole", "Cephalosporins", 
                                   "Sulfonamides"))

staph_antibio_class = left_join(staph_antibio_clean, abx_classes, by = "variable")

test_variables = unique(staph_antibio_class$class)

for(i in c(1,5,9,13,17,21,23)){
  
  pp = staph_antibio_class %>%
    group_by(date, class, value) %>%
    summarise(n = sum(n)) %>%
    mutate(prop = n/sum(n)) %>%
    filter(sum(n) > 10) %>% #only keep isolates with more than 10 tested in total
    filter(value == "R") %>%
    filter(class %in% test_variables[i:(i+3)]) %>%
    ggplot() +
    geom_line(aes(date, prop, colour = class)) +
    theme_bw() +
    scale_y_continuous(limits = c(0,1))
  
  plot(pp)
  
}


## testing below and search suggests that if an isolate is resistant to either flucloxacillin, oxacillin or cefoxitin, it is MRSA
#mrsa_ids = micro_staph %>%
#filter(SpeciesGroupName == "Methicillin-Resistant Staphylococcus aureus") %>%
#select(CultureIsoID) %>%
#pull
#let's try with cultureisoid 274645_1, D590197_1, A438456_1
#
#
#resist_mrsa = antibiograms %>%
#filter(CultureIsoID == "A438456_1") %>%
#t(.) %>%
#as.data.frame(.) %>%
#filter(V1 == "R") %>%
#rownames(.)
#
#for(i in mrsa_ids){
#    
#    if(!(i %in% unique(antibiograms$CultureIsoID))) next
#    
#    i_antibio = antibiograms %>%
#        filter(CultureIsoID == i) %>%
#        t(.) %>%
#        as.data.frame(.) %>%
#        filter(V1 == "R") %>%
#        rownames(.)
#
#    if(length(intersect(i_antibio, resist_mrsa)) == 0){
#    
#        cat(which(resist_mrsa == i))
#        break
#        
#    }
#    resist_mrsa = intersect(i_antibio, resist_mrsa)
#    
#}

#looking at changing resistance profiles

good_cols = colnames(staph_antibio)[apply(staph_antibio, 2, function(x) !(all(is.na(x))))]
good_cols = good_cols[!grepl("MIC", good_cols)]

staph_antibio_profiles = staph_antibio %>%
  dplyr::select(good_cols) %>%
  dplyr::select(1, 2, 7, 11, 67, c(13:66)) %>%
  arrange(project_id, date)
#dplyr::select(c("Chloramphenicol", "Ciprofloxacin", "Erythromycin", "Flucloxacillin", "Fucidin", "Mupirocin", "Penicillin", "Tetracycline", "Trimethoprim"))

staph_antibio_sr = staph_antibio_profiles[,-c(1:5)]
staph_antibio_sr[staph_antibio_sr == "S"] = 1
staph_antibio_sr[staph_antibio_sr == "R"] = 3
options(warn = -1)
staph_antibio_sr = apply(staph_antibio_sr, c(1,2), as.numeric)

staph_antibio_profiles[,-c(1:5)] = staph_antibio_sr

good_cols = which(colSums(staph_antibio_sr, na.rm = T) > 50)
good_rows = which(rowSums(staph_antibio_sr, na.rm = T) > 0)

staph_antibio_profiles = staph_antibio_profiles %>%
  dplyr::select(c(1:5,(good_cols+5)))

staph_antibio_profiles = staph_antibio_profiles %>%
  filter(row_number() %in% good_rows)

good_ids = staph_antibio_profiles %>%
  group_by(project_id) %>%
  summarise(n = n()) %>%
  filter(n > 1) %>%
  dplyr::select(project_id) %>%
  pull

staph_antibio_profiles = staph_antibio_profiles %>%
  filter(project_id %in% good_ids)

staph_antibio_profiles$SpecimenSource[is.na(staph_antibio_profiles$SpecimenSource)] = "*Unspecified"

interesting_samples = c()

for(i in 1:(nrow(staph_antibio_profiles)-1)){
  
  if(staph_antibio_profiles$project_id[i] != staph_antibio_profiles$project_id[i+1]) next
  if(staph_antibio_profiles$SpeciesName[i] != staph_antibio_profiles$SpeciesName[i+1]) next
  if(staph_antibio_profiles$SpecimenSource[i] != staph_antibio_profiles$SpecimenSource[i+1]) next
  
  tt = staph_antibio_profiles[c(i:(i+1)),-c(1:5)]
  test_means = apply(tt, 2, mean)
  if(any(test_means == 2, na.rm = T)) interesting_samples = c(interesting_samples,
                                                              staph_antibio_profiles$project_lab_id[i],
                                                              staph_antibio_profiles$project_lab_id[i+1])
}

changing_profiles = staph_antibio_profiles %>%
  filter(project_lab_id %in% interesting_samples)

#cor_profiles = pcor(staph_antibio_profiles)
#cor_profiles_estim = cor_profiles$estimate
#colnames(cor_profiles_estim) = colnames(staph_antibio_profiles)
#rownames(cor_profiles_estim) = colnames(staph_antibio_profiles)

#cor_profiles_estim %>%
#melt(.) %>%
#ggplot() + 
#  geom_tile(aes(x=Var1, y=Var2, fill=value)) +
#  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.1))














