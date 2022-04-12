
data = read.csv(here::here("Data", "caboodle_patient_selected_medication_admins_all.csv")) %>%
  select(project_id, start_datetime, drug_name)
data2 = read.csv(here::here("Data", "datalake_patient_selected_medication_admins_all.csv")) %>%
  select(project_id, start_datetime, drug_name)

data = rbind(data, data2)

tt = data[grepl("EICOPLAN", data$drug_name),]
nrow(tt)/nrow(data)*100

length(unique(tt$project_id))/length(unique(data$project_id))*100

tt %>%
  mutate(date = as_date(start_datetime)) %>%
  mutate(date = floor_date(date, "month")) %>%
  group_by(date) %>%
  summarise(n = n()) %>%
  ggplot() +
  geom_line(aes(date, n))
