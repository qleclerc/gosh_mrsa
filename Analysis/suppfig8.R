

library(dplyr)
library(ggplot2)
library(lubridate)
library(cowplot)

admissions = read.csv(here::here("Clean", "combined_admissions.csv")) %>%
  mutate(start_datetime = as_date(start_datetime)) %>%
  mutate(end_datetime = as_date(end_datetime)) %>%
  arrange(project_id, start_datetime) %>%
  filter(!is.na(start_datetime)) %>%
  filter(!is.na(end_datetime)) %>%
  mutate(delay = end_datetime - start_datetime) %>%
  filter(delay > 2)

median(admissions$delay)
quantile(admissions$delay)

pa = ggplot(admissions) +
  geom_histogram(aes(delay, y = after_stat(density)), binwidth = 1, colour = "white") +
  geom_vline(aes(xintercept = median(delay)), linetype = "dashed", linewidth = 1) +
  scale_x_continuous(breaks = seq(0,120,20)) + 
  coord_cartesian(xlim = c(0,80)) +
  theme_bw() +
  labs(x = "Hospital length of stay (days)", y = "Proportion") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))

table(admissions$delay>80)/nrow(admissions)

mrsa_mssa_diversity = read.csv(here::here("Clean", "mrsa_mssa_diversity.csv")) %>%
  mutate(first_date = as_date(first_date),
         second_date = as_date(second_date))

mrsa_changes = mrsa_mssa_diversity %>%
  filter(same_date == F & change == "Methicillin-Resistant Staphylococcus aureus") %>%
  filter(same_hosp == T) %>%
  mutate(delay = as.numeric(second_date-first_date)) %>%
  filter(delay > 2) %>% 
  select(project_id, los) %>%
  distinct()

median(mrsa_changes$los)

pb = ggplot(mrsa_changes) +
  geom_histogram(aes(los, y = 7*after_stat(density)), colour= "white", binwidth = 7) +
  geom_vline(aes(xintercept = median(los)), linetype = "dashed", linewidth = 1) +
  scale_x_continuous(breaks = seq(0,700,50)) + 
  coord_cartesian(xlim = c(0,700)) +
  theme_bw() +
  labs(x = "Hospital length of stay (days)", y = "Proportion") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))

table(mrsa_changes$los>700)/nrow(mrsa_changes)

mssa_changes = mrsa_mssa_diversity %>%
  filter(same_date == F & change == "Methicillin-Susceptible Staphylococcus aureus") %>%
  filter(same_hosp == T) %>%
  mutate(delay = as.numeric(second_date-first_date)) %>%
  filter(delay > 2) %>% 
  select(project_id, los) %>%
  distinct()

pc = ggplot(mssa_changes) +
  geom_histogram(aes(los, y = 7*after_stat(density)), colour= "white", binwidth = 7) +
  geom_vline(aes(xintercept = median(los)), linetype = "dashed", linewidth = 1) +
  scale_x_continuous(breaks = seq(0,700,50)) + 
  coord_cartesian(xlim = c(0,700)) +
  theme_bw() +
  labs(x = "Hospital length of stay (days)", y = "Proportion") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))

median(mssa_changes$los)

table(mssa_changes$los>700)/nrow(mssa_changes)

plot_grid(NULL,pa,NULL,pb,NULL,pc,
          ncol = 1,
          rel_heights = c(0.1,1,0.05,1,0.05,1),
          labels = c("","a) All patients", "",
                     "b) Patients with changes from MSSA to MRSA", "",
                     "c) Patients with changes from MRSA to MSSA"),
          hjust = 0, vjust=0)

ggsave(here::here("Figures", "suppfig8.png"), height = 6)
