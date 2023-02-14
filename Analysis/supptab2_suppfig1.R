library(dplyr)
library(lubridate)
library(openxlsx)
library(ggplot2)
library(ggtext)

staph_isolates = read.csv(here::here("Clean", "staph_isolates.csv")) %>% 
  mutate(date = floor_date(as_date(date), "year"))

demographics = read.csv(here::here("Data", "caboodle_patient_demographics.csv")) %>%
  mutate(birth_date = floor_date(as_date(birth_date), "year"))

admissions = read.csv(here::here("Clean", "combined_admissions.csv")) %>%
  mutate(start_datetime = as_date(start_datetime)) %>%
  mutate(end_datetime = as_date(end_datetime)) %>%
  mutate(los = end_datetime-start_datetime) %>%
  mutate(start_datetime = floor_date(start_datetime, "year")) %>%
  filter(start_datetime >= as_date("2000-01-01")) %>%
  filter(los > 1)

# descriptive stats
clean_tab = data.frame(Year = c(2000:2021), Unique_p = 0, Female = 0, Age_median = 0,
                       Age_IQR = "0-0", Ethnicity = 0, Top1sample = "a (0%)",
                       Top2sample = "a (0%)", Top3sample = "a (0%)", Adm = 0, LoS_median = 0,
                       LoS_IQR = "0-0")

dates = as_date(seq.Date(as.Date("2000-01-01"), by = "year", length.out = 23))

for(i in 1:(length(dates)-1)){
  
  iso_i = filter(staph_isolates, date >= dates[i] & date < dates[i+1])
  demo_i = filter(demographics, project_id %in% iso_i$project_id)
  adm_i = filter(admissions, start_datetime >= dates[i] & start_datetime < dates[i+1])
  
  clean_tab[i,2] = length(unique(iso_i$project_id))
  clean_tab[i,3] = round(sum(demo_i$sex=="F")/nrow(demo_i)*100, 2)
  clean_tab[i,4] = median(interval(demo_i$birth_date,dates[i])/years(1))
  clean_tab[i,5] = paste0(quantile(interval(demo_i$birth_date,dates[i])/years(1), 0.25),
                          "-",quantile(interval(demo_i$birth_date,dates[i])/years(1), 0.75))
  clean_tab[i,6] = round((nrow(demo_i)-
                            sum(demo_i$ethnicity_name=="White British")-
                            sum(demo_i$ethnicity_name=="")-
                            sum(demo_i$ethnicity_name=="Prefer Not To Say"))/
                           (nrow(demo_i)-
                              sum(demo_i$ethnicity_name=="")-
                              sum(demo_i$ethnicity_name=="Prefer Not To Say"))*100,2)
  clean_tab[i,7] = paste0(names(sort(table(iso_i$SpecimenType), decreasing = T))[1], " ",
                          round(sort(table(iso_i$SpecimenType), decreasing = T)[1]/nrow(iso_i)*100,2))
  clean_tab[i,8] = paste0(names(sort(table(iso_i$SpecimenType), decreasing = T))[2], " ",
                          round(sort(table(iso_i$SpecimenType), decreasing = T)[2]/nrow(iso_i)*100,2))
  clean_tab[i,9] = paste0(names(sort(table(iso_i$SpecimenType), decreasing = T))[3], " ",
                          round(sort(table(iso_i$SpecimenType), decreasing = T)[3]/nrow(iso_i)*100,2))
  clean_tab[i,10] = nrow(adm_i)
  clean_tab[i,11] = median(adm_i$los)
  clean_tab[i,12] = paste0(quantile(adm_i$los, 0.25),
                           "-",quantile(adm_i$los, 0.75))
  
}

write.xlsx(clean_tab, here::here("Clean", "supptab2.xlsx"))

#todo: remove sample source, clean formatting of median IQR columns

cpalette = c(RColorBrewer::brewer.pal(12, "Paired"), "grey", "black")

staph_isolates %>%
  filter(SpecimenType %in% names(sort(table(SpecimenType),decreasing=T))[1:15]) %>% 
  mutate(SpecimenType = replace(SpecimenType, SpecimenType=="Skin Swab", "Skin swab")) %>%
  group_by(date, SpecimenType) %>%
  summarise(n=n()) %>%
  group_by(date) %>%
  mutate(n=n/sum(n)) %>%
  ggplot() +
  geom_col(aes(date, n, fill = SpecimenType), position = "stack") +
  scale_fill_manual(values = cpalette) +
  theme_bw() +
  labs(x = "Year", y = "Proportion of isolates", fill = "Swab type:") +
  theme(axis.text = element_text(size=12),
        axis.title.y = element_markdown(size=12),
        axis.title.x = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12)) +
  scale_x_date(breaks = as.Date(c("2000-01-01", "2002-01-01", "2004-01-01",
                                  "2006-01-01", "2008-01-01", "2010-01-01",
                                  "2012-01-01", "2014-01-01", "2016-01-01",
                                  "2018-01-01", "2020-01-01")),
               date_labels = "%Y",
               limits = as.Date(c("1999-06-01", "2021-07-01")))

ggsave(here::here("Figures", "suppfig1.png"))

