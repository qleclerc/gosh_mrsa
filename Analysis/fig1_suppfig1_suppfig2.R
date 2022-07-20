#time-series analysis of s aureus cases
#tldr: no seasonality

library(tidyr)
library(dplyr)
library(lubridate)
library(ggplot2)
library(ggtext)
library(scales)
library(cowplot)
library(forecast)
library(zoo)

mrsa_col = "#EFC000FF"
mssa_col = "#0073C2FF"

staph_isolates = read.csv(here::here("Clean", "staph_isolates.csv")) %>% 
  mutate(date = as_date(date))

#some general plots
#baseline incidence
p1 = staph_isolates %>%
  mutate(date = floor_date(date, "month")) %>%
  group_by(date) %>%
  summarise(n = n()) %>%
  ggplot() +
  geom_line(aes(date, n), size = 0.8) +
  theme_bw() +
  labs(y = "Number of S. aureus isolates", x = "Time (months)") +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12)) +
  scale_x_date(limits = as.Date(c("2000-02-01", "2021-11-01")),
               date_breaks = "2 years", date_labels = "%Y") +
  geom_vline(xintercept = as.Date("2020-03-26"), linetype = 2, colour = "green4", size = 1)



#prop MRSA/MSSA
p2 = staph_isolates %>%
  mutate(date = floor_date(date, "month")) %>%
  group_by(date, SpeciesName) %>%
  summarise(n = n()) %>%
  mutate(prop = n/sum(n)) %>%
  ggplot() +
  geom_line(aes(date, prop, colour = SpeciesName), size = 0.8) +
  theme_bw() +
  labs(y = "Proportion of MSSA and MRSA isolates", x = "Time (months)", colour = "") +
  scale_x_date(limits = as.Date(c("2000-02-01", "2021-11-01")),
               date_breaks = "2 years", date_labels = "%Y") +
  geom_vline(xintercept = as.Date("2020-03-26"), linetype = 2, colour = "green4", size = 1) +
  scale_colour_manual(values = c(mrsa_col, mssa_col),
                      breaks = c("Methicillin-Resistant Staphylococcus aureus",
                                 "Methicillin-Susceptible Staphylococcus aureus"),
                      labels = c("Methicillin-Resistant *Staphylococcus aureus*",
                                 "Methicillin-Susceptible *Staphylococcus aureus*")) +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.text = element_markdown(size=12))

plot_grid(plot_grid(p1 + theme(legend.position = "none"),
                    NULL,
                    p2 + theme(legend.position = "none"),
                    nrow = 3, labels = c("a)", "", "b)"), label_size = 12,
                    rel_heights = c(0.7, 0.05, 1), vjust = c(1,0)),
          get_legend(p2 + theme(legend.position = "bottom")),
          nrow = 2, rel_heights = c(1,0.05))

ggsave(here::here("Figures", "fig1.png"), heigh = 8, width = 10)


# isolates per patient

staph_isolates %>%
  mutate(date = floor_date(date, "year")) %>%
  group_by(date, SpeciesName, project_id) %>%
  summarise(n = n()) %>%
  ggplot() +
  geom_boxplot(aes(date, n, group = interaction(date, SpeciesName), colour = SpeciesName),
               outlier.shape = NA, size = 0.8) +
  theme_bw() +
  scale_y_continuous(limits = c(0,8), breaks = seq(0,8,2)) +
  labs(x = "Time (years)", y = "Number of isolates per patient", colour = "") +
  theme(legend.position = "bottom",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 12)) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  scale_colour_manual(values = c(mrsa_col, mssa_col))

ggsave(here::here("Figures", "suppfig1.png"))


# ethnicity
data = read.csv(here::here("Data", "caboodle_patient_demographics.csv"))

admissions = read.csv(here::here("Data", "combined_patient_ward_stays.csv")) %>%
  mutate(start_datetime = as_date(start_datetime)) %>%
  mutate(end_datetime = as_date(end_datetime)) %>%
  arrange(project_id, start_datetime) %>%
  filter(!is.na(start_datetime)) %>%
  filter(!is.na(end_datetime))

admissions = right_join(admissions, data, "project_id")

table(data$ethnicity_name)

p1 = admissions %>%
  mutate(start_datetime = floor_date(start_datetime, "month")) %>%
  filter(start_datetime > as.Date("2000-01-01")) %>%
  filter(start_datetime < as.Date("2021-12-01")) %>%
  filter(ethnicity_name != "" & ethnicity_name != "Prefer Not To Say") %>%
  mutate(ethnicity_name = replace(ethnicity_name, ethnicity_name != "White British", "Other Groups")) %>%
  group_by(start_datetime) %>%
  count(ethnicity_name) %>%
  mutate(n = n/sum(n)) %>%
  ggplot() +
  geom_line(aes(start_datetime, n, colour = ethnicity_name), size = 0.8) +
  geom_smooth(data = . %>%
                filter(start_datetime > "2014-01-01" & start_datetime < "2020-03-26"),
              aes(start_datetime, n, colour = ethnicity_name), method = "lm") +
  theme_bw() +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  labs(x = "Time (months)", y = "Proportion of patients", colour = "Ethnicity:") +
  theme(legend.position = "bottom",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12)) +
  geom_vline(xintercept = as.Date("2020-03-26"), linetype = 2, colour = "green4", size = 1) +
  geom_vline(xintercept = as.Date("2014-01-01"), linetype = 3, colour = "black", size = 1)

p2 = staph_isolates %>%
  mutate(date = floor_date(date, "month")) %>%
  group_by(date, SpeciesName) %>%
  summarise(n = n()) %>%
  mutate(prop = n/sum(n)) %>%
  filter(SpeciesName == "Methicillin-Resistant Staphylococcus aureus") %>%
  ggplot() +
  geom_line(aes(date, prop), size = 0.8, colour = mrsa_col) +
  geom_smooth(data = . %>%
                filter(date > "2014-01-01" & date < "2020-03-26"),
              aes(date, prop), method = "lm", colour = mrsa_col) +
  theme_bw() +
  labs(y = "Proportion of MRSA isolates", x = "Time (months)", colour = "") +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12)) +
  scale_x_date(limits = as.Date(c("2000-02-01", "2021-11-01")),
               date_breaks = "2 years", date_labels = "%Y") +
  geom_vline(xintercept = as.Date("2020-03-26"), linetype = 2, colour = "green4", size = 1) +
  geom_vline(xintercept = as.Date("2014-01-01"), linetype = 3, colour = "black", size = 1)


plot_grid(p2, p1, ncol = 1, labels = c("a)", "b)"), hjust = 0)

ggsave(here::here("Figures", "suppfig2.png"))



# descriptive stats
cat(length(unique(staph_isolates$project_id)), "patients with a positive S aureus isolate between",
    as.character(min(staph_isolates$date)), "and", as.character(max(staph_isolates$date)))
cat(staph_isolates %>%
      filter(SpeciesName == "Methicillin-Resistant Staphylococcus aureus") %>%
      select(project_id) %>% pull %>% unique %>% length,
    "patients with a positive MRSA isolate")
cat(staph_isolates %>%
      filter(SpeciesName == "Methicillin-Susceptible Staphylococcus aureus") %>%
      select(project_id) %>% pull %>% unique %>% length,
    "patients with a positive MSSA isolate")
(19777+3506-22206)
19777 - 1077
3506 - 1077

(19777+3506-22206)/22206
(19777 - 1077)/22206
(3506 - 1077)/22206

staph_isolates %>%
  count(project_id) %>%
  filter(n > 1) %>%
  select(n) %>% pull %>% sum

staph_isolates %>%
  count(project_id) %>%
  filter(n > 1) %>%
  nrow()/22206*100

staph_isolates %>%
  filter(SpeciesName == "Methicillin-Resistant Staphylococcus aureus") %>%
  count(project_id) %>%
  filter(n > 1) %>%
  select(n) %>% pull %>% sum

staph_isolates %>%
  filter(SpeciesName == "Methicillin-Resistant Staphylococcus aureus") %>%
  count(project_id) %>%
  filter(n > 1) %>%
  nrow()

staph_isolates %>%
  filter(SpeciesName == "Methicillin-Susceptible Staphylococcus aureus") %>%
  count(project_id) %>%
  filter(n > 1) %>%
  select(n) %>% pull %>% sum

staph_isolates %>%
  filter(SpeciesName == "Methicillin-Susceptible Staphylococcus aureus") %>%
  count(project_id) %>%
  filter(n > 1) %>%
  nrow()


table(staph_isolates$SpeciesName)

#lm
lm_dat = staph_isolates %>%
  mutate(date = floor_date(date, "year")) %>%
  group_by(date, SpeciesName, project_id) %>%
  summarise(n = n()) %>%
  select(date, SpeciesName, n)

lm_res = lm(n ~ date + SpeciesName, data = lm_dat)
summary(lm_res)
