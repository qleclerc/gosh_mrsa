
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



#fig5
pa = right_join(mrsa_mssa_diversity %>%
                  filter(same_date == T) %>%
                  mutate(date = floor_date(first_date, "month")) %>%
                  select(project_id, date) %>%
                  distinct() %>%
                  count(date),
                staph_isolates %>%
                  mutate(date = floor_date(date, "month")) %>%
                  select(project_id, date) %>%
                  distinct() %>%
                  count(date), by = "date") %>%
  mutate(prop = n.x/n.y*100)  %>%
  ggplot() +
  geom_line(aes(date, prop)) +
  theme_bw() +
  labs(x = "Time (months)", y = "Percentage of patients with MRSA-MSSA diversity (%)") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))

pb = right_join(res_profiles_diversity %>%
                  filter(same_date == T) %>%
                  mutate(date = floor_date(first_date, "month")) %>%
                  select(project_id, date, species) %>%
                  distinct() %>%
                  group_by(date, species) %>%
                  summarise(n = n()),
                staph_isolates %>%
                  mutate(date = floor_date(date, "month")) %>%
                  select(project_id, date, SpeciesName) %>%
                  rename(species = SpeciesName) %>%
                  distinct() %>%
                  group_by(date, species) %>%
                  summarise(n = n()), by = c("date", "species")) %>%
  mutate(prop = n.x/n.y*100)  %>%
  ggplot() +
  geom_line(aes(date, prop, colour = species)) +
  theme_bw() +
  labs(x = "Time (months)", y = "Percentage of patients with resistance profile diversity (%)",
       colour = "") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.position = "bottom")

plot_grid(pa, NULL, pb,
          ncol = 1,
          rel_heights = c(1,0.05,1),
          labels = c("a)", "", "b)"))

ggsave(here::here("Figures", "fig5.png"), height = 11, width = 14)




# EVOLUTION ####
#wards where change was detected
mrsa_mssa_diversity %>%
  filter(same_date == F & change == "Methicillin-Resistant Staphylococcus aureus") %>%
  select(possible_nos_ward) %>% pull %>% table() %>% sort()

mrsa_mssa_diversity %>%
  mutate(second_date = floor_date(second_date, "year")) %>%
  group_by(second_date) %>%
  summarise(n = n()) %>%
  ggplot() +
  geom_bar(aes(second_date, n), stat="identity") +
  theme_bw()

table(mrsa_mssa_diversity$possible_nos)

mrsa_mssa_diversity %>%
  filter(possible_nos > 0) %>%
  mutate(second_date = floor_date(second_date, "year")) %>%
  group_by(second_date) %>%
  summarise(n = n()) %>%
  ggplot() +
  geom_bar(aes(second_date, n), stat="identity") +
  theme_bw()


res_profiles_diversity %>%
  filter(same_date == T) %>%
  group_by(project_id, species, second_date) %>%
  summarise(n = n()) %>%
  mutate(date = floor_date(second_date, "year")) %>%
  group_by(date, species) %>%
  summarise(mean_n = median(n),
            sd_n = sd(n)) %>%
  ggplot() +
  geom_line(aes(x = date, y = mean_n, colour = species)) +
  theme_bw() +
  theme(legend.position = "bottom")

#total patients
length(unique(staph_isolates$project_id))
#total isolates
nrow(staph_isolates)
#total tests
sum(apply(staph_isolates[,6:56], 2, function(x) length(which(!is.na(x)))))

#total patients with at least one change
length(unique(res_profiles_diversity$project_id))
#total isolates for patients with at least one change
staph_isolates %>% filter(project_id %in% unique(res_profiles_diversity$project_id)) %>% nrow(.)
#total tests for patients with at least one changes
staph_isolates %>% filter(project_id %in% unique(res_profiles_diversity$project_id)) %>% select(6:56) %>% apply(., 2, function(x) length(which(!is.na(x)))) %>% sum
#work out prop of tests which are a change
#to prove that it's not just the most tested resistances that change
#tot_samples = apply(staph_isolates[,6:59], 2, function(x) length(which(!is.na(x))))
tot_samples = staph_isolates %>% filter(project_id %in% unique(res_profiles_diversity$project_id)) %>% select(6:56) %>% apply(., 2, function(x) length(which(!is.na(x))))
tot_samples[order(tot_samples)]
change_samples = table(res_profiles_diversity$antibiotic)
change_samples[order(change_samples)]
round(change_samples/tot_samples[names(change_samples)]*100, 2)
sum(change_samples)/sum(tot_samples)*100

#accurate patients with at least one change
length(unique(res_profiles_diversity$project_id))
#accurate changes
nrow(res_profiles_diversity)
#work out prop of tests which are a change
#to prove that it's not just the most tested resistances that change
#tot_samples = apply(staph_isolates[,6:59], 2, function(x) length(which(!is.na(x))))
tot_samples = staph_isolates %>% filter(project_id %in% unique(res_profiles_diversity$project_id)) %>% select(6:59) %>% apply(., 2, function(x) length(which(!is.na(x))))
tot_samples[order(tot_samples)]
change_samples = table(res_profiles_diversity$antibiotic)
change_samples[order(change_samples)]
round(change_samples/tot_samples[names(change_samples)]*100, 2)
sum(change_samples)/sum(tot_samples)*100

#gains/losses per abx
table(res_profiles_diversity$antibiotic, res_profiles_diversity$change)

table(res_profiles_diversity$source)/table(staph_isolates$SpecimenType)[names(table(res_profiles_diversity$first_source))]*100

#boxplot of durations between samples with a change
boxplot(int_length(interval(res_profiles_diversity$first_date, res_profiles_diversity$second_date))/60/60/24)


res_profiles_diversity %>%
  group_by(project_id, second_lab_id) %>%
  summarise(n = n()) %>%
  select(n) %>%
  pull %>%
  table
