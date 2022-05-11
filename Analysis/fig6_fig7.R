
library(tidyverse)
library(ggplot2)
library(scales)
library(lubridate)
library(cowplot)

mrsa_mssa_diversity = read.csv(here::here("Clean", "mrsa_mssa_diversity.csv")) %>%
  mutate(first_date = as_date(first_date),
         second_date = as_date(second_date))

res_profiles_diversity = read.csv(here::here("Clean", "res_profiles_diversity.csv")) %>%
  mutate(first_date = as_date(first_date),
         second_date = as_date(second_date))

isolates_diversity = read.csv(here::here("Clean", "isolates_diversity.csv")) %>%
  mutate(date = as_date(date))

staph_isolates = read.csv(here::here("Clean", "staph_isolates.csv")) %>%
  mutate(date = as_date(date))


# DIVERSITY ####
#how many patients with both mrsa and mssa isolates ever detected?
length(unique(mrsa_mssa_diversity$project_id))
#how many patients had mrsa and mssa isolates detected on the same day?
mrsa_mssa_diversity %>%
  filter(same_date == T) %>%
  select(project_id) %>% pull %>% unique %>% length
#sample source
mrsa_mssa_diversity %>%
  filter(same_date == T) %>%
  select(same_source) %>% pull %>% table %>% prop.table


#how many patients with resistance profile diversity?
length(unique(res_profiles_diversity$project_id))
#how many patients with resistance profile diversity on the same day?
res_profiles_diversity %>%
  filter(same_date == T) %>%
  select(project_id) %>% pull %>% unique %>% length
#sample source
res_profiles_diversity %>%
  filter(same_date == T) %>%
  select(project_id, second_lab_id, same_source) %>%
  distinct() %>%
  select(same_source) %>% pull %>% table %>% prop.table
#diversity by MRSA/MSSA isolates
res_profiles_diversity %>%
  filter(same_date == T) %>%
  select(species) %>% pull %>% table


#fig6
pa = right_join(mrsa_mssa_diversity %>%
                  filter(same_date == T) %>%
                  mutate(date = floor_date(second_date, "year")) %>%
                  select(project_id, date) %>%
                  distinct() %>%
                  count(date),
                staph_isolates %>%
                  mutate(date = floor_date(date, "year")) %>%
                  select(project_id, date) %>%
                  distinct() %>%
                  count(date), by = "date") %>%
  mutate(prop = n.x/n.y*100)  %>%
  ggplot() +
  geom_line(aes(date, prop)) +
  theme_bw() +
  labs(x = "Time (years)", y = "Patients with MRSA-MSSA\ndiversity (%)") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12)) +
  scale_y_continuous(limits = c(0,4), breaks = seq(0,4,1)) +
  scale_x_date(breaks = as.Date(c("2000-01-01", "2002-01-01", "2004-01-01",
                                  "2006-01-01", "2008-01-01", "2010-01-01",
                                  "2012-01-01", "2014-01-01", "2016-01-01",
                                  "2018-01-01", "2020-01-01")), date_labels = "%Y",
               limits = as.Date(c("2000-01-01", "2021-01-01")))

pb = right_join(res_profiles_diversity %>%
                  filter(same_date == T) %>%
                  mutate(date = floor_date(second_date, "year")) %>%
                  select(project_id, date, species) %>%
                  distinct() %>%
                  group_by(date, species) %>%
                  summarise(n = n()),
                staph_isolates %>%
                  mutate(date = floor_date(date, "year")) %>%
                  select(project_id, date, SpeciesName) %>%
                  rename(species = SpeciesName) %>%
                  distinct() %>%
                  group_by(date, species) %>%
                  summarise(n = n()), by = c("date", "species")) %>%
  mutate(prop = n.x/n.y*100)  %>%
  ggplot() +
  geom_line(aes(date, prop, colour = species)) +
  theme_bw() +
  labs(x = "Time (years)", y = "Patients with resistance profile\ndiversity (%)",
       colour = "") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.position = "bottom") +
  scale_y_continuous(limits = c(0,12), breaks = seq(0,12,2)) +
  scale_x_date(breaks = as.Date(c("2000-01-01", "2002-01-01", "2004-01-01",
                                  "2006-01-01", "2008-01-01", "2010-01-01",
                                  "2012-01-01", "2014-01-01", "2016-01-01",
                                  "2018-01-01", "2020-01-01")), date_labels = "%Y",
               limits = as.Date(c("2000-01-01", "2021-01-01")))

#number of unique profiles identified when diversity is identified
pc = isolates_diversity %>%
  filter(n_profiles > 1) %>%
  group_by(species, n_profiles) %>%
  summarise(n = n()) %>%
  mutate(n = n/sum(n)*100) %>%
  ggplot() +
  geom_col(aes(as.factor(n_profiles), n, fill = species)) +
  facet_wrap(~species) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 12),
        legend.position = "none") +
  labs(x = "Number of unique resistance profiles identified within the same patient",
       y = "Patients with diversity (%)") +
  scale_y_continuous(limits = c(0,100), breaks = seq(0,100,20))

#number of differences identified between isolates with different profiles
pd = res_profiles_diversity %>%
  filter(same_date == T) %>%
  group_by(project_id, species, first_date) %>%
  summarise(n_diffs = length(unique(antibiotic))) %>%
  group_by(species, n_diffs) %>%
  summarise(n = n()) %>%
  mutate(n = n/sum(n)*100) %>%
  ggplot() +
  geom_col(aes(as.factor(n_diffs), n, fill = species)) +
  facet_wrap(~species) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 12),
        legend.position = "none") +
  labs(x = "Number of differing resistances between isolates within the same patient",
       y = "Patients with diversity (%)") +
  scale_y_continuous(limits = c(0,100), breaks = seq(0,100,20))


plot_grid(plot_grid(pa, NULL, pb,
                    ncol = 1,
                    rel_heights = c(0.8,0.01,1),
                    labels = c("a)", "", "b)")),
          NULL,
          pc,
          NULL,
          pd,
          ncol = 1,
          rel_heights = c(1,0.01,0.4,0.01,0.4),
          labels = c("", "", "c)", "", "d)"))

ggsave(here::here("Figures", "fig6.png"), height = 14, width = 11)




# EVOLUTION ####
#changes from mssa to mrsa
mrsa_changes = mrsa_mssa_diversity %>%
  filter(same_date == F & change == "Methicillin-Resistant Staphylococcus aureus")

nrow(mrsa_changes)
length(unique(mrsa_changes$project_id))

#changes from mssa to mrsa < 365 days
mrsa_changes = mrsa_mssa_diversity %>%
  filter(same_date == F & change == "Methicillin-Resistant Staphylococcus aureus") %>%
  mutate(delay = second_date-first_date) %>%
  mutate(delay = as.numeric(delay, units = "days")) %>%
  filter(delay <= 365)

nrow(mrsa_changes)
length(unique(mrsa_changes$project_id))

#changes from mssa to mrsa < 365 days & > 1 day
mrsa_changes = mrsa_mssa_diversity %>%
  filter(same_date == F & change == "Methicillin-Resistant Staphylococcus aureus") %>%
  mutate(delay = second_date-first_date) %>%
  mutate(delay = as.numeric(delay, units = "days")) %>%
  filter(delay <= 365 & delay > 1)

nrow(mrsa_changes)
length(unique(mrsa_changes$project_id))


pa = ggplot(mrsa_changes) +
  geom_histogram(aes(delay, y = 7*..density..), binwidth = 7, colour = "white") +
  geom_vline(xintercept = median(mrsa_changes$delay), linetype = "dashed") +
  theme_bw() +
  scale_x_continuous(breaks = seq(0, 365, 30)) +
  labs(x = "Delay between MSSA and MRSA isolates (days)", y = "Proportion") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))

#changes from mssa to mrsa, potential nosocomial
mrsa_changes = mrsa_mssa_diversity %>%
  filter(same_date == F & change == "Methicillin-Resistant Staphylococcus aureus") %>%
  mutate(delay = second_date-first_date) %>%
  mutate(delay = as.numeric(delay, units = "days")) %>%
  filter(delay <= 365 & delay > 1) %>%
  filter(same_source == T & same_hosp == T) %>%
  filter(possible_nos > 0)
nrow(mrsa_changes)
length(unique(mrsa_changes$project_id))

pb = mrsa_changes %>%
  mutate(second_date = floor_date(second_date, "year")) %>%
  group_by(second_date) %>%
  summarise(n = n()) %>%
  ggplot() +
  geom_bar(aes(second_date, n), stat="identity") +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12)) +
  labs(x = "Detection date for potential MRSA nosocomial transmission (year)",
       y = "Count")

#changes from mssa to mrsa, mysterious events
mrsa_changes = mrsa_mssa_diversity %>%
  filter(same_date == F & change == "Methicillin-Resistant Staphylococcus aureus") %>%
  mutate(delay = second_date-first_date) %>%
  mutate(delay = as.numeric(delay, units = "days")) %>%
  filter(delay <= 365 & delay > 1) %>%
  filter(same_source == T & same_hosp == T) %>%
  filter(possible_nos == 0)
nrow(mrsa_changes)
length(unique(mrsa_changes$project_id))



#changes in resistance profiles
res_profiles_changes = res_profiles_diversity %>%
  filter(same_date == F)
length(unique(res_profiles_changes$project_id)) #number of patients
length(unique(res_profiles_changes$project_id))/length(unique(staph_isolates$project_id)) #proportion

res_profiles_changes = res_profiles_diversity %>%
  filter(same_date == F & same_hosp == T) %>%
  mutate(delay = second_date-first_date) %>%
  mutate(delay = as.numeric(delay, units = "days")) %>%
  filter(delay <= 365 & delay > 1)
length(unique(res_profiles_changes$project_id)) #number of patients
length(unique(res_profiles_changes$project_id))/length(unique(staph_isolates$project_id)) #proportion

pc = ggplot(res_profiles_changes) +
  geom_histogram(aes(delay, y = 7*..density.., group = species, fill = species), binwidth = 7,
                 colour = "white", position = "identity", alpha = 0.5) +
  geom_vline(xintercept = median(res_profiles_changes$delay), linetype = "dashed") +
  theme_bw() +
  scale_x_continuous(breaks = seq(0, 365, 30)) +
  labs(x = "Delay between isolates with differing resistance profile (days)", y = "Proportion",
       fill = "") +
  theme(legend.position = "bottom",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 12)) +
  guides(fill = guide_legend(override.aes = list(alpha = 1)))

pd = right_join(res_profiles_changes %>%
                  mutate(date = floor_date(second_date, "year")) %>%
                  select(project_id, date, species) %>%
                  distinct() %>%
                  group_by(date, species) %>%
                  summarise(n = n()),
                staph_isolates %>%
                  mutate(date = floor_date(date, "year")) %>%
                  select(project_id, date, SpeciesName) %>%
                  rename(species = SpeciesName) %>%
                  distinct() %>%
                  group_by(date, species) %>%
                  summarise(n = n()),
                by = c("date", "species")) %>%
  mutate(prop = n.x/n.y*100)  %>%
  ggplot() +
  geom_line(aes(date, prop, colour = species)) +
  theme_bw() +
  labs(x = "Time (years)", y = "Patients with changes in isolates resistance profile (%)",
       colour = "") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.position = "bottom")


plot_grid(plot_grid(pa, NULL, pb,
                    nrow = 1, rel_widths = c(1,0.05,1),
                    labels = c("a)", "", "b)")),
          NULL,
          plot_grid(pc + theme(legend.position = "none"),
                    NULL,
                    pd + theme(legend.position = "none"),
                    nrow = 1, rel_widths = c(1,0.05,1),
                    labels = c("c)", "", "d)")),
          get_legend(pd),
          nrow = 4,
          rel_heights = c(1,0.05,1,0.1))

ggsave(here::here("Figures", "fig7.png"))


#total tests
sum(apply(staph_isolates[,6:56], 2, function(x) length(which(!is.na(x)))))
#total isolates for patients with at least one change
staph_isolates %>% filter(project_id %in% unique(res_profiles_changes$project_id)) %>% nrow(.)
#total tests for patients with at least one changes
staph_isolates %>% filter(project_id %in% unique(res_profiles_changes$project_id)) %>% select(6:56) %>% apply(., 2, function(x) length(which(!is.na(x)))) %>% sum
#work out prop of tests which are a change
#to prove that it's not just the most tested resistances that change
#tot_samples = apply(staph_isolates[,6:59], 2, function(x) length(which(!is.na(x))))
tot_samples = staph_isolates %>% filter(project_id %in% unique(res_profiles_changes$project_id)) %>% select(6:56) %>% apply(., 2, function(x) length(which(!is.na(x))))
tot_samples[order(tot_samples)]
change_samples = table(res_profiles_changes$antibiotic)
change_samples[order(change_samples)]
round(change_samples/tot_samples[names(change_samples)]*100, 2)
sum(change_samples)/sum(tot_samples)*100

#gains/losses per abx
table(res_profiles_changes$antibiotic, res_profiles_changes$change) %>% prop.table(margin = 1) %>% round(2)
