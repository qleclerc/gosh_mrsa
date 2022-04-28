
library(dplyr)
library(ggplot2)

data = read.csv(here::here("Data", "caboodle_patient_selected_medication_admins_all.csv")) %>%
  select(project_id, start_datetime, drug_name)
data2 = read.csv(here::here("Data", "datalake_patient_selected_medication_admins_all.csv")) %>%
  select(project_id, start_datetime, drug_name)

antibio_data = rbind(data, data2) %>%
  mutate(start_datetime = as_date(start_datetime))


equivalence_table = data.frame(listed = unique(antibio_data$drug_name),
                               name = unique(antibio_data$drug_name))

for(i in 1:nrow(equivalence_table)){
  
  cat("\nName in data:", equivalence_table$listed[i], "\n")
  
  cat("\nSuggested name data:", equivalence_table$name[i], "\n")
  
  change_name = readline("Change suggested name? ")
  
  if(change_name != ""){
    
    equivalence_table$name[i] = change_name
    
  }
  
}

write.csv(equivalence_table, here::here("Clean", "equivalence_table.csv"), row.names = F)

colnames(antibio_data)[3] = "listed"

antibio_data = right_join(antibio_data, equivalence_table, by = "listed")

write.csv(antibio_data, here::here("Clean", "antibio_data.csv"), row.names = F)
