

library(dplyr)
library(ggplot2)
library(lubridate)

res_profiles_diversity = read.csv(here::here("Clean", "res_profiles_diversity.csv")) %>%
  mutate(first_date = as_date(first_date),
         second_date = as_date(second_date))

staph_isolates = read.csv(here::here("Clean", "staph_isolates.csv")) %>%
  mutate(date = as_date(date))

mrsa_mssa_diversity = read.csv(here::here("Clean", "mrsa_mssa_diversity.csv")) %>%
  mutate(first_date = as_date(first_date),
         second_date = as_date(second_date))

mrsanoabx = table(mrsa_mssa_diversity %>%
                    filter(change == "Methicillin-Resistant Staphylococcus aureus") %>%
                    filter(same_hosp == T) %>%
                    mutate(delay = as.numeric(second_date-first_date)) %>%
                    filter(delay > 2) %>%
                    select(any_antibiotic) %>% pull)[1]
mrsaabx = table(mrsa_mssa_diversity %>%
                  filter(change == "Methicillin-Resistant Staphylococcus aureus") %>%
                  filter(same_hosp == T) %>%
                  mutate(delay = as.numeric(second_date-first_date)) %>%
                  filter(delay > 2) %>%
                  select(any_antibiotic) %>% pull)[2]

mssanoabx = table(mrsa_mssa_diversity %>%
                    filter(change == "Methicillin-Susceptible Staphylococcus aureus") %>%
                    filter(same_hosp == T) %>%
                    mutate(delay = as.numeric(second_date-first_date)) %>%
                    filter(delay > 2) %>%
                    select(any_antibiotic) %>% pull)[1]
mssaabx = table(mrsa_mssa_diversity %>%
                  filter(change == "Methicillin-Susceptible Staphylococcus aureus") %>%
                  filter(same_hosp == T) %>%
                  mutate(delay = as.numeric(second_date-first_date)) %>%
                  filter(delay > 2) %>%
                  select(any_antibiotic) %>% pull)[2]

chisq.test(matrix(c(mrsanoabx,mrsaabx,
                    mssanoabx,mssaabx),
                  byrow = T, nrow = 2))

mrsaabx/(mrsaabx+mrsanoabx)
mssaabx/(mssaabx+mssanoabx)


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

antibio_data = read.csv(here::here("Clean", "antibio_data.csv")) %>%
  mutate(date = as_date(start_datetime)) %>%
  filter(class %in% c("Aminoglycoside", "Carbapenem", "Cephalosporin",
                      "Fluoroquinolone", "Fosfomycin", "Fusidic_acid",
                      "Glycopeptide", "Lincosamide", "Macrolide", "Metronidazole",
                      "Nitrofuran", "Penicillin", "Polypeptide", "Sulfonamide",
                      "Tetracycline", "Rifamycin"))

admissions = read.csv(here::here("Data", "combined_patient_ward_stays.csv")) %>%
  mutate(start_datetime = as_date(start_datetime)) %>%
  mutate(end_datetime = as_date(end_datetime)) %>%
  arrange(project_id, start_datetime) %>%
  filter(project_id %in% unique(res_profiles_diversity$project_id)) %>%
  filter(!is.na(start_datetime)) %>%
  filter(!is.na(end_datetime))

sort(table(antibio_data$class))

antibio_data %>%
  group_by(class) %>%
  summarise(n = length(unique(project_id))) %>%
  arrange(-n) %>%
  mutate(n = n/22206) %>%
  mutate(class = replace(class, class == "Fusidic_acid", "Fucidin")) %>%
  mutate(class = factor(class, levels = class)) %>%
  ggplot() +
  geom_col(aes(class, n)) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.title = element_text(size = 12)) +
  labs(x = "", y = "Proportion of patients exposed to antibiotic")

ggsave(here::here("Figures", "suppfig7.png"))


#patients with change in a hosp, who had any abx
abxchange = res_profiles_diversity %>%
  filter(same_hosp == T) %>%
  filter(any_antibiotic == T) %>%
  nrow
#patients with change in a hosp, who had no abx 
noabxchange = res_profiles_diversity %>%
  filter(same_hosp == T) %>%
  filter(any_antibiotic == F) %>%
  nrow

#look for patients with no change who had abx anytime between isolate 1 date and isolate 1 date + 60 days
#because 60 days is cut off to see vast majority of changes (Fig 7d)
# nochange = staph_isolates %>%
#   filter(!project_id %in% unique(res_profiles_diversity$project_id)) %>%
#   select(project_id, date)
# nochange$any_abx = FALSE
# 
# for(i in 1:nrow(nochange)){
#   
#   if(nrow(antibio_data %>%
#           filter(project_id == nochange$project_id[i]) %>%
#           filter(date %within% interval(nochange$date[i], nochange$date[i]+60))) > 0) nochange$any_abx[i] = T
#   
# }
# 
# write.csv(nochange, here::here("Clean", "antibiotic_post_isolate.csv"), row.names = F)  

nochange = read.csv(here::here("Clean", "antibiotic_post_isolate.csv"))

#patients with no change, who had any abx
abxnochange = nochange %>%
  filter(any_abx == T) %>%
  nrow
#patients with no change, who had no abx
noabxnochange = nochange %>%
  filter(any_abx == F) %>%
  nrow

#isolates with/without sus test data x isolates from nose,throat,skin/not
chisq.test(matrix(c(abxchange,abxnochange,
                    noabxchange,noabxnochange),
                  byrow = T, nrow = 2))

abxchange/(abxchange+abxnochange)
noabxchange/(noabxchange+noabxnochange)


antibio_data %>%
  group_by(project_id) %>%
  summarise(n = length(unique(class))) %>%
  select(n) %>%
  pull %>% table/20000*100

antibio_data %>%
  mutate(date = floor_date(date, "month")) %>%
  group_by(date, class) %>%
  summarise(n = n()) %>%
  ggplot() +
  geom_line(aes(date, n, colour = class)) +
  theme_bw() +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y")

antibio_data %>%
  filter(class == "Fluoroquinolone") %>%
  mutate(date = floor_date(date, "month")) %>%
  group_by(date) %>%
  summarise(n = n()) %>%
  ggplot() +
  geom_line(aes(date, n)) +
  theme_bw() +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y")
