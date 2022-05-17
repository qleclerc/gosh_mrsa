
library(dplyr)
library(tidyr)
library(lubridate)
library(reshape2)
library(ggplot2)
library(scales)
library(cowplot)

mrsa_col = "#EFC000FF"
mssa_col = "#0073C2FF"

staph_isolates = read.csv(here::here("Clean", "staph_isolates.csv"), stringsAsFactors = F) %>%
  mutate(date = as_date(date))

#the line to remove all isolates with NAs only (ie nothing tested at all)
#staph_isolates = staph_isolates[!apply(staph_isolates[,6:59], 1, function(x) all(is.na(x))),]

mrsa_isolates_m = staph_isolates %>%
  filter(SpeciesName == "Methicillin-Resistant Staphylococcus aureus") %>%
  select(-c(1:4)) %>%
  melt(id.vars = "date") %>%
  filter(!is.na(value)) %>%
  mutate(date = floor_date(date, "month")) %>%
  group_by(date, variable) %>%
  summarise(n = n()) %>%
  filter(variable %in% c("Chloramphenicol", "Mupirocin", "Tetracycline",
                         "Teicoplanin", "Vancomycin", "Syncercid", "Cotrimoxazole",
                         "Linezolid", "Clindamycin", "Fosfomycin", "Cefoxitin"))

mssa_isolates_m = staph_isolates %>%
  filter(SpeciesName == "Methicillin-Susceptible Staphylococcus aureus") %>%
  select(-c(1:4)) %>%
  melt(id.vars = "date") %>%
  filter(!is.na(value)) %>%
  mutate(date = floor_date(date, "month")) %>%
  group_by(date, variable) %>%
  summarise(n = n()) %>%
  filter(variable %in% c("Chloramphenicol", "Mupirocin", "Tetracycline",
                         "Teicoplanin", "Vancomycin", "Syncercid", "Cotrimoxazole",
                         "Linezolid", "Clindamycin", "Fosfomycin", "Cefoxitin"))

mrsa_per_day = staph_isolates %>%
  filter(SpeciesName == "Methicillin-Resistant Staphylococcus aureus") %>%
  mutate(date = floor_date(date, "month")) %>%
  group_by(date) %>%
  summarise(total = n())

mrsa_isolates_m = left_join(mrsa_isolates_m, mrsa_per_day, by = "date")
mrsa_isolates_m$species = "Methicillin-Resistant Staphylococcus aureus"

mssa_per_day = staph_isolates %>%
  filter(SpeciesName == "Methicillin-Susceptible Staphylococcus aureus") %>%
  mutate(date = floor_date(date, "month")) %>%
  group_by(date) %>%
  summarise(total = n())

mssa_isolates_m = left_join(mssa_isolates_m, mssa_per_day, by = "date")
mssa_isolates_m$species = "Methicillin-Susceptible Staphylococcus aureus"

all_isolates = rbind(mrsa_isolates_m, mssa_isolates_m) %>%
  tibble() %>%
  complete(date, variable, species)

p1 = ggplot(all_isolates %>%
              filter(variable %in% c("Mupirocin"))) +
  geom_line(aes(date, n/total, colour = species)) +
  scale_y_continuous(limits = c(0,1)) +
  facet_wrap(~variable, nrow = 1) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 12)) +
  labs(x = "Time (months)", y = "Proportion of isolates tested") +
  scale_x_date(limits = as.Date(c("2000-02-01", "2021-11-01"))) +
  scale_colour_manual(values = c(mrsa_col, mssa_col))

p2 = ggplot(all_isolates %>% 
              filter(variable %in% c("Fosfomycin", "Cotrimoxazole")) %>%
              mutate(variable = factor(variable, levels = c("Cotrimoxazole", "Fosfomycin")))) +
  geom_line(aes(date, n/total, colour = species)) +
  scale_y_continuous(limits = c(0,1)) +
  facet_wrap(~variable, nrow = 1) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 12)) +
  labs(x = "Time (months)", y = "Proportion of isolates tested", colour = "") +
  scale_x_date(limits = as.Date(c("2000-02-01", "2021-11-01"))) +
  scale_colour_manual(values = c(mrsa_col, mssa_col))

p3 = ggplot(all_isolates %>% 
              filter(variable %in% c("Clindamycin", "Cefoxitin"))) +
  geom_line(aes(date, n/total, colour = species)) +
  scale_y_continuous(limits = c(0,1)) +
  facet_wrap(~variable, nrow = 1) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 12)) +
  labs(x = "Time (months)", y = "Proportion of isolates tested") +
  scale_x_date(limits = as.Date(c("2000-02-01", "2021-11-01"))) +
  scale_colour_manual(values = c(mrsa_col, mssa_col))

p4 = ggplot(all_isolates %>% 
              filter(variable %in% c("Chloramphenicol", "Tetracycline"))) +
  geom_line(aes(date, n/total, colour = species)) +
  scale_y_continuous(limits = c(0,1)) +
  facet_wrap(~variable, ncol = 1) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 12)) +
  labs(x = "Time (months)", y = "Proportion of isolates tested") +
  scale_x_date(limits = as.Date(c("2000-02-01", "2021-11-01"))) +
  scale_colour_manual(values = c(mrsa_col, mssa_col))

p5 = ggplot(all_isolates %>% 
              filter(variable %in% c("Linezolid", "Syncercid", "Teicoplanin", "Vancomycin"))) +
  geom_line(aes(date, n/total, colour = species)) +
  scale_y_continuous(limits = c(0,1)) +
  facet_wrap(~variable, ncol = 1) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 12)) +
  labs(x = "Time (months)", y = "Proportion of isolates tested") +
  scale_x_date(limits = as.Date(c("2000-02-01", "2021-11-01"))) +
  scale_colour_manual(values = c(mrsa_col, mssa_col))

plot_grid(plot_grid(plot_grid(plot_grid(p1 + theme(legend.position = "none"),
                                        NULL,
                                        p4 + theme(legend.position = "none"),
                                        ncol = 1, rel_heights = c(0.5, 0.05, 1),
                                        labels = c("a)", "", "d)"), hjust = 0, label_size = 12),
                              NULL,
                              plot_grid(p2 + theme(legend.position = "none"),
                                        NULL,
                                        p3 + theme(legend.position = "none"),
                                        nrow = 3, rel_heights = c(1,0.05,1),
                                        labels = c("b)", "", "c)"), hjust = 0, label_size = 12),
                              ncol = 3, rel_widths = c(0.7,0.05,1)),
                    NULL,
                    p5 + theme(legend.position = "none"),
                    rel_widths = c(1,0.05,0.4), nrow = 1, labels = c("", "", "e)")),
          get_legend(p1+theme(legend.position = "bottom")),
          ncol = 1, rel_heights = c(1,0.05))

ggsave(here::here("Figures","fig5.png"), height = 10, width = 14)





mrsa_isolates_m = staph_isolates %>%
  filter(SpeciesName == "Methicillin-Resistant Staphylococcus aureus") %>%
  select(-c(1:4)) %>%
  melt(id.vars = "date") %>%
  filter(!is.na(value)) %>%
  mutate(date = floor_date(date, "month")) %>%
  group_by(date, variable) %>%
  summarise(n = n())

mssa_isolates_m = staph_isolates %>%
  filter(SpeciesName == "Methicillin-Susceptible Staphylococcus aureus") %>%
  select(-c(1:4)) %>%
  melt(id.vars = "date") %>%
  filter(!is.na(value)) %>%
  mutate(date = floor_date(date, "month")) %>%
  group_by(date, variable) %>%
  summarise(n = n())

mrsa_per_day = staph_isolates %>%
  filter(SpeciesName == "Methicillin-Resistant Staphylococcus aureus") %>%
  mutate(date = floor_date(date, "month")) %>%
  group_by(date) %>%
  summarise(total = n())

mrsa_isolates_m = left_join(mrsa_isolates_m, mrsa_per_day, by = "date")
mrsa_isolates_m$species = "Methicillin-Resistant Staphylococcus aureus"

mssa_per_day = staph_isolates %>%
  filter(SpeciesName == "Methicillin-Susceptible Staphylococcus aureus") %>%
  mutate(date = floor_date(date, "month")) %>%
  group_by(date) %>%
  summarise(total = n())

mssa_isolates_m = left_join(mssa_isolates_m, mssa_per_day, by = "date")
mssa_isolates_m$species = "Methicillin-Susceptible Staphylococcus aureus"

all_isolates = rbind(mrsa_isolates_m, mssa_isolates_m) %>%
  tibble() %>%
  complete(date, variable, species)

sp1 = ggplot(all_isolates) +
  geom_line(aes(date, n/total, colour = species)) +
  scale_y_continuous(limits = c(0,1)) +
  facet_wrap(~variable) +
  theme_bw() +
  labs(x = "Time (months)", y = "Proportion of isolates tested", colour = "") +
  scale_x_date(limits = as.Date(c("2000-02-01", "2021-11-01"))) +
  scale_colour_manual(values = c(mrsa_col, mssa_col)) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(here::here("Figures", "suppfig5.png"), sp1, height = 12, width = 12)
