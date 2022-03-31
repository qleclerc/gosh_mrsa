library(dplyr)
library(lubridate)

## 1) Clean S aureus isolates data ####

#using antibiograms as proxy to identify mrsa
antibiograms = read.csv(here::here("Data" ,"labanalysis_patient_microbiology_positives_sensitivitys_pivot.csv"))

#assumption: species names in that list only are s aureus
staph_isolates = antibiograms %>%
  filter(SpeciesName %in% c("Staphylococcus aureus", "Staphylococcus sp.", "Staphylococcus sp", "Methicillin-Resistant Staphylococcus aureus"))

staph_isolates[staph_isolates == ""] = NA

staph_isolates$date = ymd_hms(staph_isolates$start_datetime)

#assumption: resistance to flucloxacillin, oxacillin or cefoxitin is indicative of mrsa
staph_isolates$SpeciesName[staph_isolates$Flucloxacillin == "R"] = "Methicillin-Resistant Staphylococcus aureus"
staph_isolates$SpeciesName[staph_isolates$Flucloxacillin == "I"] = "Methicillin-Resistant Staphylococcus aureus"
staph_isolates$SpeciesName[staph_isolates$Oxacillin == "R"] = "Methicillin-Resistant Staphylococcus aureus"
staph_isolates$SpeciesName[staph_isolates$Oxacillin == "I"] = "Methicillin-Resistant Staphylococcus aureus"
staph_isolates$SpeciesName[staph_isolates$Cefoxitin == "R"] = "Methicillin-Resistant Staphylococcus aureus"
staph_isolates$SpeciesName[staph_isolates$Cefoxitin == "I"] = "Methicillin-Resistant Staphylococcus aureus"
staph_isolates$SpeciesName[staph_isolates$SpeciesName != "Methicillin-Resistant Staphylococcus aureus"] = "Methicillin-Susceptible Staphylococcus aureus"

good_cols = colnames(staph_isolates)[apply(staph_isolates, 2, function(x) !(all(is.na(x))))]
good_cols = good_cols[!grepl("MIC", good_cols)]

staph_isolates = staph_isolates %>%
  select(all_of(good_cols)) %>%
  select(1, 2, 6, 11, 67, c(13:66)) %>%
  arrange(project_id, date) %>%
  filter(date < "2021-12-01") %>%
  filter(date >= "2000-02-01")

tt = staph_isolates[,6:59]
tt[is.na(tt)] = "NA"
tt[tt != "S" & tt != "R"] = NA
staph_isolates[,6:59] = tt

write.csv(staph_isolates, here::here("Clean", "staph_isolates.csv"), row.names = F)

