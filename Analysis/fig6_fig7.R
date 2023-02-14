
library(tidyverse)
library(ggplot2)
library(ggtext)
library(scales)
library(lubridate)
library(cowplot)

mrsa_col = "#EFC000FF"
mssa_col = "#0073C2FF"

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
mrsa_mssa_diversity %>%
  filter(same_date == T) %>%
  nrow
#sample source
mrsa_mssa_diversity %>%
  filter(same_date == T) %>%
  select(same_source) %>% pull %>% table %>% prop.table


#patients with different isolates detected
res_profiles_diversity %>%
  filter(same_date == T) %>%
  select(project_id) %>% pull %>% unique %>% length
res_profiles_diversity %>%
  filter(species == "Methicillin-Resistant Staphylococcus aureus") %>%
  select(project_id) %>% pull %>% unique %>% length
res_profiles_diversity %>%
  filter(species == "Methicillin-Susceptible Staphylococcus aureus") %>%
  select(project_id) %>% pull %>% unique %>% length
#how many patients with resistance profile diversity on the same day?
res_profiles_diversity %>%
  filter(species == "Methicillin-Resistant Staphylococcus aureus") %>%
  filter(same_date == T) %>%
  select(project_id) %>% pull %>% unique %>% length
res_profiles_diversity %>%
  filter(species == "Methicillin-Resistant Staphylococcus aureus") %>%
  filter(same_date == T) %>%
  group_by(project_id, first_lab_id) %>%
  summarise(n = n()) %>%
  nrow
res_profiles_diversity %>%
  filter(species == "Methicillin-Susceptible Staphylococcus aureus") %>%
  filter(same_date == T) %>%
  select(project_id) %>% pull %>% unique %>% length
res_profiles_diversity %>%
  filter(species == "Methicillin-Susceptible Staphylococcus aureus") %>%
  filter(same_date == T) %>%
  group_by(project_id, first_lab_id) %>%
  summarise(n = n()) %>%
  nrow

#sample source
res_profiles_diversity %>%
  filter(same_date == T) %>%
  select(project_id, second_lab_id, same_source) %>%
  distinct() %>%
  select(same_source) %>% pull %>% table %>% prop.table


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
  mutate(prop = n.x/n.y)  %>%
  ggplot() +
  geom_line(aes(date, prop), size = 0.8) +
  theme_bw() +
  labs(x = "Time (years)", y = "Proportion of patients\nwith MRSA-MSSA diversity") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12)) +
  scale_y_continuous(limits = c(0,0.04), breaks = seq(0,0.04,0.01)) +
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
  mutate(prop = n.x/n.y)  %>%
  ggplot() +
  geom_line(aes(date, prop, colour = species), size = 0.8) +
  theme_bw() +
  labs(x = "Time (years)", y = "Proportion of patients with\nphenotypic resistance diversity",
       colour = "") +
  scale_y_continuous(limits = c(0,0.14), breaks = seq(0,0.14,0.02)) +
  scale_x_date(breaks = as.Date(c("2000-01-01", "2002-01-01", "2004-01-01",
                                  "2006-01-01", "2008-01-01", "2010-01-01",
                                  "2012-01-01", "2014-01-01", "2016-01-01",
                                  "2018-01-01", "2020-01-01")), date_labels = "%Y",
               limits = as.Date(c("2000-01-01", "2021-01-01"))) +
  scale_colour_manual(values = c(mrsa_col, mssa_col),
                      breaks = c("Methicillin-Resistant Staphylococcus aureus",
                                 "Methicillin-Susceptible Staphylococcus aureus"),
                      labels = c("Methicillin-Resistant *Staphylococcus aureus*",
                                 "Methicillin-Susceptible *Staphylococcus aureus*")) +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.text = element_markdown(size=12),
        legend.position = "bottom")

#number of unique profiles identified when diversity is identified
pc = isolates_diversity %>%
  filter(n_profiles > 1) %>%
  group_by(species, n_profiles) %>%
  summarise(n = n()) %>%
  mutate(n = n/sum(n)) %>%
  ggplot() +
  geom_col(aes(as.factor(n_profiles), n, fill = species), position = position_dodge(0.9)) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 12),
        legend.position = "none") +
  labs(x = "Unique antibiograms recorded within the same patient",
       y = "Proportion") +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.2)) +
  scale_fill_manual(values = c(mrsa_col, mssa_col))

#number of differences identified between isolates with different profiles
pd = res_profiles_diversity %>%
  filter(same_date == T) %>%
  group_by(project_id, species, first_date) %>%
  summarise(n_diffs = length(unique(antibiotic))) %>%
  group_by(species, n_diffs) %>%
  summarise(n = n()) %>%
  mutate(n = n/sum(n)) %>%
  ggplot() +
  geom_col(aes(as.factor(n_diffs), n, fill = species), position = position_dodge(0.9)) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 12),
        legend.position = "none") +
  labs(x = "Differing resistances between antibiograms for the same patient",
       y = "Proportion") +
  scale_y_continuous(limits = c(0,1,0.2), breaks = seq(0,1,0.2)) +
  scale_fill_manual(values = c(mrsa_col, mssa_col))


#which changes for each abx
#R vs S changes
res_m = res_profiles_diversity %>%
  filter(same_date == T) %>%
  group_by(species, antibiotic) %>%
  summarise(n = n()) %>%
  group_by(species) %>%
  mutate(prop = n/sum(n)) %>%
  tibble() %>%
  complete(species, antibiotic)
res_m[is.na(res_m)] = 0

top10mrsa = res_m %>%
  filter(species == "Methicillin-Resistant Staphylococcus aureus") %>%
  group_by(antibiotic) %>%
  summarise(n = sum(n)) %>%
  arrange(-n) %>%
  select(antibiotic) %>% pull %>% .[1:10]

top10mssa = res_m %>%
  filter(species == "Methicillin-Susceptible Staphylococcus aureus") %>%
  group_by(antibiotic) %>%
  summarise(n = sum(n)) %>%
  arrange(-n) %>%
  select(antibiotic) %>% pull %>% .[1:10]

pe = ggplot(res_m %>%
              filter(species == "Methicillin-Resistant Staphylococcus aureus") %>%
              filter(antibiotic %in% top10mrsa) %>%
              mutate(antibiotic = factor(antibiotic, levels = top10mrsa))) +
  geom_col(aes(antibiotic, prop),
           colour = "white", fill = mrsa_col) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.title = element_text(size = 12),
        legend.position = "none") +
  labs(x = "", y = "Proportion of\ndiffering resistances")

pf = ggplot(res_m %>%
              filter(species == "Methicillin-Susceptible Staphylococcus aureus") %>%
              filter(antibiotic %in% top10mssa) %>%
              mutate(antibiotic = factor(antibiotic, levels = top10mssa))) +
  geom_col(aes(antibiotic, prop),
           colour = "white", fill = mssa_col) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.title = element_text(size = 12),
        legend.position = "none") +
  labs(x = "", y = "Proportion of\ndiffering resistances")


plot_grid(plot_grid(pa, NULL, pb + theme(legend.position = "none"),
                    ncol = 1,
                    rel_heights = c(0.8,0.01,0.8),
                    labels = c("a)", "", "b)"), hjust = 0),
          NULL,
          plot_grid(pc, NULL, pd, nrow = 1, labels = c("c)", "", "d)"),
                    rel_widths = c(1,0.05,1), hjust = 0),
          NULL,
          plot_grid(pe, NULL, pf, nrow = 1, labels = c("e)", "", "f)"),
                    rel_widths = c(1,0.05,1), hjust = 0),
          get_legend(pb),
          ncol = 1,
          rel_heights = c(0.7,0.01,0.4,0.01,0.4,0.1))

ggsave(here::here("Figures", "fig6.png"), height = 14, width = 11)




# EVOLUTION ####
#work out the complete delay distribution for all isolates (not only different ones)
#only keep patients with more than one sample
good_ids = staph_isolates %>%
  count(project_id) %>%
  filter(n > 1) %>%
  select(project_id) %>%
  pull

staph_isolates_profiles = staph_isolates %>%
  filter(project_id %in% good_ids)

all_delays = data.frame()

#for each row, if the patient is the same, compare species
for(i in 1:(nrow(staph_isolates_profiles)-1)){
  
  if(staph_isolates_profiles$project_id[i] != staph_isolates_profiles$project_id[i+1]) next
  
  all_delays = rbind(all_delays,
                     data.frame(species1 = staph_isolates_profiles$SpeciesName[i],
                                species2 = staph_isolates_profiles$SpeciesName[i+1],
                                delay = staph_isolates_profiles$date[i+1] - staph_isolates_profiles$date[i]))
}

all_delays$same_species = (all_delays$species1 == all_delays$species2)
all_delays$delay = as.numeric(all_delays$delay)

sum(all_delays$delay > 130)/sum(all_delays$delay > 2)

pb = ggplot(all_delays %>% filter(delay > 2)) +
  geom_histogram(aes(delay, y = 7*..density..), binwidth = 7, colour = "white") +
  geom_vline(xintercept = median(all_delays$delay), linetype = "dashed", size = 1) +
  scale_x_continuous(breaks = seq(0,120,20)) + 
  coord_cartesian(xlim = c(0,130)) +
  theme_bw() +
  labs(x = "Delay between any subsequent isolates (days)", y = "Proportion") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))

#all changes
mrsa_mssa_diversity %>%
  filter(same_date == F) %>%
  select(project_id) %>% pull %>% unique %>% length
mrsa_mssa_diversity %>%
  filter(same_date == F) %>%
  nrow
#changes from mrsa to mssa
mrsa_mssa_diversity %>%
  filter(same_date == F & change == "Methicillin-Susceptible Staphylococcus aureus") %>%
  nrow
mrsa_mssa_diversity %>%
  filter(same_date == F & change == "Methicillin-Susceptible Staphylococcus aureus") %>%
  select(project_id) %>% pull %>% unique %>% length

#changes from mssa to mrsa
mrsa_mssa_diversity %>%
  filter(same_date == F & change == "Methicillin-Resistant Staphylococcus aureus") %>%
  nrow
mrsa_mssa_diversity %>%
  filter(same_date == F & change == "Methicillin-Resistant Staphylococcus aureus") %>%
  select(project_id) %>% pull %>% unique %>% length

#changes from mssa to mrsa outside hosp
mrsa_mssa_diversity %>%
  filter(same_date == F & change == "Methicillin-Resistant Staphylococcus aureus") %>%
  filter(same_hosp == F) %>%
  nrow

#interesting changes from mssa to mrsa inside hosp
mrsa_changes = mrsa_mssa_diversity %>%
  filter(same_date == F & change == "Methicillin-Resistant Staphylococcus aureus") %>%
  filter(same_hosp == T) %>%
  mutate(delay = as.numeric(second_date-first_date))

table(mrsa_changes$delay <= 2)

mrsa_changes = mrsa_changes %>%
  filter(delay > 2)


nrow(mrsa_changes)
length(unique(mrsa_changes$project_id))

median(mrsa_changes$delay)
table(mrsa_changes$delay)

mrsa_changes %>% filter(any_antibiotic == T) %>% nrow()
mrsa_changes %>% filter(any_antibiotic == T) %>% select(project_id) %>% pull %>% unique %>% length

median(mrsa_changes %>% select(project_id, los) %>% distinct() %>% select(los) %>% pull)
quantile(mrsa_changes %>% select(project_id, los) %>% distinct() %>% select(los) %>% pull)

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


pa = ggplot(mrsa_changes) +
  geom_histogram(aes(delay, y = 7*..density..), binwidth = 7, colour = "white") +
  geom_vline(xintercept = median(mrsa_changes$delay), linetype = "dashed", size = 1) +
  scale_x_continuous(breaks = seq(0,120,20)) + 
  coord_cartesian(xlim = c(0,130)) +
  scale_y_continuous(breaks = seq(0,0.6,0.1), limits = c(0,0.6)) +
  theme_bw() +
  labs(x = "Delay between MSSA and MRSA isolates (days)", y = "Proportion") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))

#changes from mssa to mrsa, potential nosocomial
mrsa_changes %>%
  filter(possible_nos == 0) %>%
  select(project_id) %>%
  pull %>% unique %>% length
mrsa_changes = mrsa_changes %>%
  filter(possible_nos > 0)

nrow(mrsa_changes)
length(unique(mrsa_changes$project_id))

pc = mrsa_changes %>%
  mutate(second_date = floor_date(second_date, "year")) %>%
  group_by(second_date) %>%
  summarise(n = n()) %>%
  ggplot() +
  geom_bar(aes(second_date, n), stat="identity") +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12)) +
  labs(x = "Time (years)",
       y = "Number of potential\nnosocomial-acquired MRSA")


#changes from mssa to mrsa, mysterious events
mrsa_changes_mys = mrsa_mssa_diversity %>%
  filter(same_date == F & change == "Methicillin-Resistant Staphylococcus aureus") %>%
  filter(same_hosp == T) %>%
  mutate(delay = as.numeric(second_date-first_date)) %>%
  filter(delay > 2) %>%
  filter(possible_nos == 0) %>%
  filter(any_antibiotic == F)

nrow(mrsa_changes_mys)
length(unique(mrsa_changes_mys$project_id))



#changes in resistance profiles
#how many patients with resistance profile diversity on different day?
res_profiles_diversity %>%
  filter(species == "Methicillin-Resistant Staphylococcus aureus") %>%
  filter(same_date == F) %>%
  select(project_id) %>% pull %>% unique %>% length
res_profiles_diversity %>%
  filter(species == "Methicillin-Resistant Staphylococcus aureus") %>%
  filter(same_date == F) %>%
  group_by(project_id, first_lab_id) %>%
  summarise(n = n()) %>%
  nrow
res_profiles_diversity %>%
  filter(species == "Methicillin-Susceptible Staphylococcus aureus") %>%
  filter(same_date == F) %>%
  select(project_id) %>% pull %>% unique %>% length
res_profiles_diversity %>%
  filter(species == "Methicillin-Susceptible Staphylococcus aureus") %>%
  filter(same_date == F) %>%
  group_by(project_id, first_lab_id) %>%
  summarise(n = n()) %>%
  nrow


res_profiles_changes = res_profiles_diversity %>%
  filter(same_date == F & same_hosp == T) %>%
  mutate(delay = as.numeric(second_date-first_date)) %>%
  filter(delay > 2)

length(unique(res_profiles_changes$project_id))
res_profiles_changes %>%
  group_by(project_id, first_lab_id) %>%
  summarise(n = n()) %>%
  nrow

#mystery changes
res_profiles_changes %>%
  filter(possible_nos == 0) %>%
  filter(antibiotic_use == F) %>%
  group_by(project_id, first_lab_id) %>%
  summarise(n = n()) %>%
  nrow
res_profiles_changes %>%
  filter(possible_nos > 0) %>%
  group_by(project_id, first_lab_id) %>%
  summarise(n = n()) %>%
  nrow
res_profiles_changes %>%
  filter(antibiotic_use == T) %>%
  group_by(project_id, first_lab_id) %>%
  summarise(n = n()) %>%
  nrow



res_profiles_changes %>%
  filter(species == "Methicillin-Resistant Staphylococcus aureus") %>%
  select(project_id) %>% pull %>% unique %>% length
res_profiles_changes %>%
  filter(species == "Methicillin-Resistant Staphylococcus aureus") %>%
  group_by(project_id, first_lab_id) %>%
  summarise(n = n()) %>%
  nrow
res_profiles_changes %>%
  filter(species == "Methicillin-Susceptible Staphylococcus aureus") %>%
  select(project_id) %>% pull %>% unique %>% length
res_profiles_changes %>%
  filter(species == "Methicillin-Susceptible Staphylococcus aureus") %>%
  group_by(project_id, first_lab_id) %>%
  summarise(n = n()) %>%
  nrow

res_profiles_changes %>%
  filter(species == "Methicillin-Resistant Staphylococcus aureus") %>%
  filter(possible_nos > 0) %>%
  select(project_id) %>% pull %>% unique %>% length
res_profiles_changes %>%
  filter(species == "Methicillin-Resistant Staphylococcus aureus") %>%
  filter(possible_nos > 0) %>%
  group_by(project_id, first_lab_id) %>%
  summarise(n = n()) %>%
  nrow
res_profiles_changes %>%
  filter(species == "Methicillin-Susceptible Staphylococcus aureus") %>%
  filter(possible_nos > 0) %>%
  select(project_id) %>% pull %>% unique %>% length
res_profiles_changes %>%
  filter(species == "Methicillin-Susceptible Staphylococcus aureus") %>%
  filter(possible_nos > 0) %>%
  group_by(project_id, first_lab_id) %>%
  summarise(n = n()) %>%
  nrow

res_profiles_changes %>%
  filter(species == "Methicillin-Resistant Staphylococcus aureus") %>%
  filter(antibiotic_use == T) %>%
  select(project_id) %>% pull %>% unique %>% length
res_profiles_changes %>%
  filter(species == "Methicillin-Resistant Staphylococcus aureus") %>%
  filter(antibiotic_use == T) %>%
  group_by(project_id, first_lab_id) %>%
  summarise(n = n()) %>%
  nrow
res_profiles_changes %>%
  filter(species == "Methicillin-Susceptible Staphylococcus aureus") %>%
  filter(antibiotic_use == T) %>%
  select(project_id) %>% pull %>% unique %>% length
res_profiles_changes %>%
  filter(species == "Methicillin-Susceptible Staphylococcus aureus") %>%
  filter(antibiotic_use == T) %>%
  group_by(project_id, first_lab_id) %>%
  summarise(n = n()) %>%
  nrow

#mystery changes
res_profiles_changes %>%
  filter(antibiotic_use == F,
         same_hosp == T)

nrow(res_profiles_changes)
length(unique(res_profiles_changes$project_id)) #number of patients
length(unique(res_profiles_changes$project_id))/length(unique(staph_isolates$project_id)) #proportion

changes_delays = res_profiles_changes %>%
  group_by(project_id, first_lab_id, species) %>%
  summarise(delay = mean(delay)) %>%
  select(delay) %>% pull

median(changes_delays)
sum(changes_delays > 130)

pd = res_profiles_changes %>%
  group_by(project_id, first_lab_id, species) %>%
  summarise(delay = mean(delay)) %>%
  ggplot() +
  geom_histogram(aes(delay, y = 7*..density.., group = species, fill = species), binwidth = 7,
                 colour = "white", position = "identity", alpha = 0.5) +
  geom_vline(xintercept = median(res_profiles_changes$delay), linetype = "dashed", size = 1) +
  scale_x_continuous(breaks = seq(0,120,20)) +
  coord_cartesian(xlim = c(0,130)) +
  scale_y_continuous(breaks = seq(0,0.6,0.1), limits = c(0,0.6)) +
  theme_bw() +
  labs(x = "Delay before change in detected phenotypic resistance (days)", y = "Proportion",
       fill = "") +
  theme(legend.position = "bottom",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 12)) +
  guides(fill = guide_legend(override.aes = list(alpha = 1))) +
  scale_fill_manual(values = c(mrsa_col, mssa_col))

pe = right_join(res_profiles_changes %>%
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
  mutate(prop = n.x/n.y)  %>%
  ggplot() +
  geom_line(aes(date, prop, colour = species)) +
  theme_bw() +
  labs(x = "Time (years)", y = "Proportion of patients with changes\nin phenotypic resistance",
       colour = "") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.position = "bottom") +
  scale_colour_manual(values = c(mrsa_col, mssa_col))


#R vs S changes
res_m = res_profiles_changes %>%
  group_by(species, antibiotic, change, antibiotic_use) %>%
  summarise(n = n()) %>%
  group_by(species) %>%
  mutate(prop = n/sum(n)) %>%
  tibble() %>%
  complete(species, antibiotic, change, antibiotic_use)
res_m[is.na(res_m)] = 0

top10mrsa = res_m %>%
  filter(species == "Methicillin-Resistant Staphylococcus aureus") %>%
  group_by(antibiotic) %>%
  summarise(n = sum(n)) %>%
  arrange(-n) %>%
  select(antibiotic) %>% pull %>% .[1:10]

top10mssa = res_m %>%
  filter(species == "Methicillin-Susceptible Staphylococcus aureus") %>%
  group_by(antibiotic) %>%
  summarise(n = sum(n)) %>%
  arrange(-n) %>%
  select(antibiotic) %>% pull %>% .[1:10]

pf = ggplot(res_m %>%
              filter(species == "Methicillin-Resistant Staphylococcus aureus") %>%
              filter(antibiotic %in% top10mrsa) %>%
              mutate(antibiotic = factor(antibiotic, levels = top10mrsa))) +
  geom_col(data = . %>%
             filter(antibiotic_use == F),
           aes(antibiotic, prop, fill = change), colour = "white", position = position_dodge(), alpha = 0.5) +
  geom_col(data = . %>%
             filter(antibiotic_use == T),
           aes(antibiotic, prop, fill = change), colour = "white", position = position_dodge()) +
  scale_y_continuous(breaks = seq(0,0.08,0.02), limits = c(0,0.095)) +
  geom_text(data = . %>%
              filter(antibiotic_use == F),
            aes(antibiotic, prop, fill = change, label = change),
            position = position_dodge(0.9), vjust = -0.3) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.title = element_text(size = 12),
        legend.position = "none") +
  labs(x = "", y = "Proportion of changes in MRSA") +
  scale_fill_manual(values = c(mrsa_col, mrsa_col))

pg = ggplot(res_m %>%
              filter(species == "Methicillin-Susceptible Staphylococcus aureus") %>%
              filter(antibiotic %in% top10mssa) %>%
              mutate(antibiotic = factor(antibiotic, levels = top10mssa))) +
  geom_col(data = . %>%
             filter(antibiotic_use == F),
           aes(antibiotic, prop, fill = change), colour = "white", position = position_dodge(), alpha = 0.5) +
  geom_col(data = . %>%
             filter(antibiotic_use == T),
           aes(antibiotic, prop, fill = change), colour = "white", position = position_dodge()) +
  scale_y_continuous(breaks = seq(0,0.08,0.02), limits = c(0,0.095)) +
  geom_text(data = . %>%
              filter(antibiotic_use == F),
            aes(antibiotic, prop, fill = change, label = change),
            position = position_dodge(0.9), vjust = -0.3) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.title = element_text(size = 12),
        legend.position = "none") +
  labs(x = "", y = "Proportion of changes in MSSA") +
  scale_fill_manual(values = c(mssa_col, mssa_col))


plot_grid(plot_grid(plot_grid(pa, NULL, pb, NULL, pd + theme(legend.position = "none"),
                              ncol = 1, rel_heights = c(1,0.05,1,0.05,1),
                              labels = c("a)", "", "b)", "","d)"), align = "v",
                              hjust = 0),
                    NULL,
                    plot_grid(pc,
                              NULL,
                              pe + theme(legend.position = "none"),
                              ncol = 1, rel_heights = c(1,0.05,1),
                              labels = c("c)", "", "e)"), hjust = 0, vjust = c(1.5, 0, -0.5)),
                    ncol = 3,
                    rel_widths = c(1,0.05,1)), 
          NULL,
          pf,
          pg,
          ncol = 1,
          rel_heights = c(1,0.01,0.6,0.6),
          labels = c("", "", "f)", "g)"), vjust = 0, hjust = 0)

ggsave(here::here("Figures", "fig7.png"), height = 14, width = 10)


