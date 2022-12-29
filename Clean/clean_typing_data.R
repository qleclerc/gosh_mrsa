
library(dplyr)
library(lubridate)
library(ggplot2)
library(RColorBrewer)

data1 = read.csv(here::here("Data", "omni_patient_selected_lab_components_all.csv"))
data2 = read.csv(here::here("Data", "caboodle_patient_selected_lab_components_all.csv"))

setdiff(colnames(data1), colnames(data2))

data2 = data2 %>% 
  select(-c("SpecimenType", "Method", "PathologyType", "SourceType",
            "Section", "SubSection", "LCD", "labcomponent_header",
            "component_type", "component_subtype", "component_datatype",
            "LCRF", "Unit", "low_range", "high_range", "flag", "abnormal",
            "encounter_key", "reference_range"))

data1 = data1 %>% 
  select(-c("specimen_code", "site_code", "site_name",
            "specimen_comment_text", "working_diagnosis",
            "test_set_code", "test_code", "result_code",
            "unit_code", "unit_name", "lower_range", "upper_range", "reference_range")) %>%
  rename(SpecimenSource = specimen_name,
         collected_datetime = collected,
         received_datetime = received,
         verified_datetime = authorised,
         component_basename = test_name,
         TestName = test_set_name,
         Value = result_text,
         NumericValue = result_numeric) %>%
  mutate(LabTestKey = NA,
         component_name = NA,
         ResultStatus = NA) %>%
  select(colnames(data2))

data = rbind(data2, data1) %>%
  arrange(start_datetime) %>%
  filter(TestName %in% c("STAU_OLD",
                         "STAU_OLD2",
                         "Staph aureus type/toxin testing",
                         "Staph. aureus Type/Toxin - 1",
                         "Staph. aureus Type/Toxin - 2",
                         "Staph. aureus Type/Toxin - 3",
                         "TYPE")) %>%
  filter(!(component_name %in% c("Date Final Report Received",
                                 "Date Interim Report Received",
                                 "Date final report received",
                                 "Reference Number","Result From",
                                 "Organism Sent","Isolate Number"))) %>%
  mutate(Value = stringr::str_to_lower(Value)) %>%
  mutate(across(matches("datetime"), as_date))

length(unique(data$project_id))
length(unique(data$Value))

data$Value[grepl("not detected", data$Value)] = "negative"
data$Value[grepl("absent", data$Value)] = "negative"
data$Value[grepl("neg", data$Value)] = "negative"
data$Value[data$Value == "-"] = "negative" #assumption

data$Value[grepl("gene detected", data$Value)] = "positive"
data$Value[grepl("gene present", data$Value)] = "positive"
data$Value[data$Value == "+"] = "positive" #assumption

sort(table(data$Value), decreasing = T)[1:100]

sort(table(data$component_basename), decreasing = T)



## only keep st data ###
#22, 30, 239 (asian one)

data_cc = data %>%
  filter(component_basename %in% c("BRES2", "MLST", "MLST Clonal Complex",
                                   "MLST Profile", "MLST ST", "MLSTAP", "MLSTCC",
                                   "STP3", "SPATYPE1","SPATYPE2","SPATYPE3",
                                   "SPARS1","SPARS2","SPARS3")) %>%
  mutate(Value = gsub(" ", "", Value)) %>%
  mutate(Value = gsub("-", "", Value)) %>%
  mutate(Value = gsub("st", "", Value)) %>%
  mutate(Value = gsub("spa", "", Value)) %>%
  mutate(Value = gsub("type", "", Value)) %>%
  mutate(Value = gsub(":", "", Value)) %>%
  mutate(Value = gsub("\\*", "", Value)) %>%
  mutate(Value = gsub("cc", "", Value)) %>%
  mutate(Value = gsub("\\.", ",", Value)) %>%
  filter(Value != "diinctrain") %>%
  filter(Value != "nontypable") %>%
  filter(Value != "widespreadsporadic") %>%
  filter(Value != "skinassociatedpattern") %>%
  filter(Value != "lytic")


data_cc$Value[data_cc$Value == "1111111"] = "1,1,1,1,1,1,1"
data_cc$Value[data_cc$Value == "1,1,1,1,1,1,1,"] = "1,1,1,1,1,1,1"
data_cc$Value[data_cc$Value == "2222632"] = "2,2,2,2,6,3,2"

spatypes = openxlsx::read.xlsx(here::here("Clean", "spatypes.xlsx")) %>%
  mutate(mlst = gsub("ST-","",mlst)) %>%
  mutate(mlst = gsub(",.*", "", mlst)) %>%
  mutate(rep = gsub("-", "", rep))
  
spatypes$mlst[spatypes$spa=="t304"] = 8
spatypes$mlst[spatypes$spa=="t688"] = 5
spatypes$mlst[spatypes$spa=="t267"] = 97
spatypes$mlst[spatypes$spa=="t548"] = 5
spatypes$mlst[spatypes$spa=="t657"] = 772
spatypes$mlst[spatypes$spa=="t786"] = 88
spatypes$mlst[spatypes$spa=="t852"] = 22
spatypes$mlst[spatypes$spa=="t020"] = 22
spatypes$mlst[spatypes$spa=="t437"] = 59
spatypes$mlst[spatypes$spa=="t690"] = 88



mlstdata = read.delim(here::here("Clean","bigsdb.txt"))

mlstdata$Value = apply(mlstdata[,-c(1,9)], 1, function(x) paste0(x, collapse = ","))

mlstdata = mlstdata %>%
  select(ST, clonal_complex, Value) %>%
  mutate(ST = as.character(ST)) %>%
  mutate(clonal_complex = gsub("CC", "", clonal_complex))

mlstdata$clonal_complex[mlstdata$ST=="80"] = "80"
mlstdata$clonal_complex[mlstdata$ST=="88"] = "88"
mlstdata$clonal_complex[mlstdata$ST=="59"] = "59"

data_cc$ST = "None"
data_cc$CC = "None"

for(i in 1:nrow(data_cc)){
  
  #STs or CCs
  #note: assuming that the value is an ST, so might bias towards representative CCs
  if(data_cc$component_basename[i] %in% c("MLST", "MLST Clonal Complex",
                                          "MLST ST", "MLSTCC")){
    mlstdata_index = mlstdata %>%
      filter(ST == data_cc$Value[i])
    if(nrow(mlstdata_index) > 0){
      data_cc$ST[i] = as.character(data_cc$Value[i])
      data_cc$CC[i] = as.character(mlstdata_index$clonal_complex)
    }
  }
  
  
  #mlst
  if(data_cc$component_basename[i] %in% c("MLST Profile", "MLSTAP")){
    mlstdata_index = mlstdata %>%
      filter(Value == data_cc$Value[i])
    if(nrow(mlstdata_index) > 0){
      data_cc$ST[i] = as.character(mlstdata_index$ST)
      data_cc$CC[i] = as.character(mlstdata_index$clonal_complex)
    }
  }
  
  #spa raw seq
  if(data_cc$component_basename[i] %in% c("SPARS1", "SPARS2", "SPARS3")){
    spatypes_index = spatypes %>%
      filter(rep == data_cc$Value[i])
    if(nrow(spatypes_index) > 0){
      data_cc$Value[i] = spatypes_index$spa
    }
  }
  
  #spa types
  #will catch those from loop above
  if(grepl("t[0-9]", data_cc$Value[i])){
    spatypes_index = spatypes %>%
      filter(spa == data_cc$Value[i])
    if(nrow(spatypes_index) > 0){
      data_cc$ST[i] = as.numeric(spatypes_index$mlst)
      
      mlstdata_index = mlstdata %>%
        filter(ST == data_cc$ST[i])
      
      if(nrow(mlstdata_index) > 0){
        data_cc$CC[i] = as.character(mlstdata_index$clonal_complex)
        data_cc$ST[i] = as.character(data_cc$ST[i])
      }
    }
  }
  
  
  
  #random h000a things
  if(grepl("h", data_cc$Value[i])){
    
    data_cc$Value[i] = gsub("\\D", "", data_cc$Value[i])
    
    mlstdata_index = mlstdata %>%
      filter(ST == data_cc$Value[i])
    if(nrow(mlstdata_index) > 0){
      data_cc$ST[i] = as.character(data_cc$Value[i])
      data_cc$CC[i] = as.character(mlstdata_index$clonal_complex)
    }
  }
  
  
  #mrsa15
  if(grepl("a15", data_cc$Value[i])){
    data_cc$ST[i] = "22"
    data_cc$CC[i] = "22"
  }
  
  
  #mrsa16
  if(grepl("a16", data_cc$Value[i])){
    data_cc$ST[i] = "36"
    data_cc$CC[i] = "30"
  }
  
  
}

data_cc$CC[is.na(data_cc$CC)] = "None"
data_cc$CC[data_cc$CC == ""] = "None"

data_cc$ST[is.na(data_cc$ST)] = "None"

#adding STs that are not CCs
data_cc$CC[data_cc$ST != "None" & data_cc$CC == "None"] = paste0("ST",
                                                              data_cc$ST[data_cc$ST != "None" & 
                                                                           data_cc$CC == "None"])

100 - (sum(data_cc$ST == "None")/nrow(data_cc)*100)
data_cc %>% filter(ST == "None") %>% select(Value) %>% table %>% sort(.,decreasing = T)

write.csv(data_cc, here::here("Clean", "typing_data.csv"), row.names = F)

# data_cc_sum = data_cc %>%
#   mutate(CC = replace(CC, grepl("ST", CC), "Other")) %>%
#   filter(CC != "0") %>%
#   select(project_id, start_datetime, CC) %>%
#   distinct() %>%
#   mutate(start_datetime = floor_date(start_datetime, "year")) %>%
#   group_by(start_datetime, CC) %>%
#   summarise(n = n())
# 
# ggplot() +
#   geom_col(data = data_cc_sum, aes(start_datetime, n, fill = CC), position = "stack") +
#   theme_bw() +
#   labs(y = "Number of typed isolates", x = "Year") +
#   geom_text(data = data_cc_sum %>%
#               group_by(start_datetime) %>%
#               summarise(total = sum(n)),
#             aes(start_datetime, total, label = total), vjust = -0.5) +
#   geom_text(aes(x = as.Date("2005-01-01"), y = 200,
#                 label = "Numbers above bars correspond to total\nnumber of typed isolates for that year")) +
#   scale_fill_brewer(palette = "Paired") +
#   theme(axis.text = element_text(size = 12),
#         axis.title = element_text(size = 12),
#         legend.text = element_text(size = 12),
#         legend.title = element_text(size = 12))
# 
# ggsave("CC_fig.png")
# 
# data_other_sum = data_cc %>%
#   filter(grepl("ST", CC)) %>%
#   select(project_id, start_datetime, CC) %>%
#   distinct() %>%
#   mutate(start_datetime = floor_date(start_datetime, "year")) %>%
#   group_by(start_datetime, CC) %>%
#   summarise(n = n())
# 
# ggplot() +
#   geom_col(data = data_other_sum, aes(start_datetime, n, fill = CC), position = "fill") +
#   theme_bw() +
#   labs(y = "Proportion of typed isolates", x = "Year", fill = "Other STs") +
#   geom_text(data = data_other_sum %>%
#               group_by(start_datetime) %>%
#               summarise(total = sum(n)),
#             aes(start_datetime, 1, label = total), vjust = -0.5)
# 
