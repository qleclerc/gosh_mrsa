
#note: this script cleans the raw GOSH data
#this data is the property of GOSH and is not available on this GitHub repository


library(dplyr)
library(lubridate)

## 1) Clean S aureus isolates data ####

#using antibiograms as proxy to identify mrsa
antibiograms = read.csv(here::here("Data" ,"labanalysis_patient_microbiology_positives_sensitivitys_pivot.csv"))

#assumption: species names in that list only are s aureus
staph_isolates = antibiograms %>%
  filter(SpeciesName %in% c("Staphylococcus aureus", "Methicillin-Resistant Staphylococcus aureus"))

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

#remove columns with only NA values
good_cols = colnames(staph_isolates)[apply(staph_isolates, 2, function(x) !(all(is.na(x))))]
#remove MIC columns
good_cols = good_cols[!grepl("MIC", good_cols)]

staph_isolates = staph_isolates %>%
  select(all_of(good_cols)) %>%
  select(1, 6, 10, 11, 67, c(13:66)) %>%
  arrange(project_id, date) %>%
  filter(date < "2021-12-01") %>%
  filter(date >= "2000-02-01")

#remove anything that's not S or R
#(nb this is extremely rare, so acceptable)
tt = staph_isolates[,6:59]
tt[is.na(tt)] = "NA"
tt[tt != "S" & tt != "R"] = NA
staph_isolates[,6:59] = tt

#septrin == cotrimoxazole, so merge and remove cotrimoxazole
for(i in 1:nrow(staph_isolates)){
  if(is.na(staph_isolates$Septrin[i])){
    staph_isolates$Septrin[i] = staph_isolates$Co.Trimoxazole..Septrin.[i]
  }
}

staph_isolates = staph_isolates %>%
  select(colnames(staph_isolates)[apply(staph_isolates, 2, function(x) !(all(is.na(x))))]) %>%
  select(-(Co.Trimoxazole..Septrin.)) %>%
  rename(Cotrimoxazole = Septrin)
  
write.csv(staph_isolates, here::here("Clean", "staph_isolates.csv"), row.names = F)

