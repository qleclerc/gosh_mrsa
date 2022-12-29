
library(tidyverse)
library(ggplot2)
library(lubridate)
library(reshape2)

admissions = read.csv(here::here("Data", "combined_patient_ward_stays.csv")) %>%
  mutate(start_datetime = as_date(start_datetime)) %>%
  mutate(end_datetime = as_date(end_datetime)) %>%
  arrange(project_id, start_datetime) %>%
  filter(!is.na(start_datetime)) %>%
  filter(!is.na(end_datetime))

i = 1
to_remove = c()

while(i < nrow(admissions)){
  
  i_step = 1
  
  repeat{
    
    if(admissions$project_id[i] == admissions$project_id[i+i_step] &&
       admissions$end_datetime[i] == admissions$start_datetime[i+i_step]){
      admissions$end_datetime[i] = admissions$end_datetime[i+i_step]
      to_remove = c(to_remove, i+i_step)
      i_step = i_step+1
    }
    else{
      i = i+i_step
      break
    }
  }
}

combined_admissions = admissions[-to_remove,c(1:3)]

write.csv(combined_admissions, here::here("Clean", "combined_admissions.csv"), row.names = F)
